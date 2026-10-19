#!/usr/bin/env python3
"""timsTOF (diaPASEF) counterpart of fragment_intensity_eval.py's `extract` (docs/testing.md, "Fragment-intensity
accuracy"): observed fragment intensities read from the raw Bruker .d at each identified precursor's apex, written
in the same per-ion Parquet layout, so that `fragment_intensity_eval.py score` scores every predictor on identical
precursors and ion sets.

Needs Bruker's timsdata library (libtimsdata.so / timsdata.dll, from Bruker's TDF SDK; DIA-NN ships a copy in its
installation directory), passed with --timsdata. It is not part of this repository and is loaded only at run time.
numpy + pyarrow otherwise.

For every precursor of a DIA-NN report at Q.Value <= --qvalue (optionally only those also in --also-in reports):
  - the diaPASEF isolation window that holds its m/z and whose scan (1/K0) range holds its apex 1/K0 (report IM);
  - the frame of that window group closest to the apex RT and --apex-halfwidth frames of the same group either side;
  - in each frame, the scans within +-(--im-halfwidth) 1/K0 of the apex IM, clipped to the window's scan range;
  - raw TOF peaks of those scans summed (Bruker timsdata library: tims_read_scans_v2, per-frame index -> m/z);
  - every theoretical b and y ion (no neutral loss, fragment z 1 .. min(--max-fragment-charge, precursor z),
    ordinals 1 .. n-1, inside the frame's m/z range): intensity = sum of the raw peaks within --ppm (profile-like TOF
    data, so a sum, not the single most intense centroid of the Orbitrap extractor). Unobserved ions are kept as 0.
  - flags as in fragment_intensity_eval.py: ambiguous (another ion of the precursor within 2 x ppm) and shared (the
    ion lies within ppm of a b/y ion, z <= 3, of another identified precursor of the same window whose RT.Start..
    RT.Stop overlaps the frames used and whose apex 1/K0 is within 2 x im-halfwidth).

usage: tims_extract.py --report report.parquet (DIA-NN 2.x, with IM) --d run.d --timsdata libtimsdata.so --out run.ions.parquet
         [--also-in other_report.parquet ...] [--ppm 15] [--im-halfwidth 0.02] [--apex-halfwidth 1]
         [--max-fragment-charge 3] [--qvalue 0.01] [--workers 16]
"""
import argparse, ctypes, os, sqlite3, sys
from multiprocessing import Pool
import numpy as np, pyarrow as pa, pyarrow.parquet as pq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fragment_intensity_eval as fie

A = None   # args, set in main / worker init
_dll = _h = None
_k0 = {}

def open_tims(lib, d):
    try:
        dll = ctypes.cdll.LoadLibrary(lib)
    except OSError as e:
        raise SystemExit(f"cannot load Bruker's timsdata library {lib!r} ({e}); pass --timsdata with the path to "
                         "libtimsdata.so / timsdata.dll (Bruker TDF SDK, or the copy in DIA-NN's directory)")
    dll.tims_open.argtypes = [ctypes.c_char_p, ctypes.c_uint32]; dll.tims_open.restype = ctypes.c_uint64
    dll.tims_read_scans_v2.argtypes = [ctypes.c_uint64, ctypes.c_int64, ctypes.c_uint32, ctypes.c_uint32, ctypes.c_void_p, ctypes.c_uint32]
    dll.tims_read_scans_v2.restype = ctypes.c_uint32
    for f in ("tims_index_to_mz", "tims_scannum_to_oneoverk0"):
        getattr(dll, f).argtypes = [ctypes.c_uint64, ctypes.c_int64, ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.c_uint32]
        getattr(dll, f).restype = ctypes.c_uint32
    dll.tims_set_num_threads.argtypes = [ctypes.c_uint32]; dll.tims_set_num_threads(1)   # parallel over processes instead
    h = dll.tims_open(d.encode(), 0)
    if h == 0: raise SystemExit(f"tims_open failed for {d}")
    return dll, h

def convert(fn, frame, x):
    x = np.ascontiguousarray(x, np.float64); out = np.empty_like(x)
    fn(_h, frame, x.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), out.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), len(x))
    return out

def read_scans(frame, s0, s1):
    """(tof indices, intensities) of scans s0..s1-1 of one frame, concatenated"""
    n = 1 << 16
    while True:
        buf = np.empty(n, np.uint32)
        need = _dll.tims_read_scans_v2(_h, frame, s0, s1, buf.ctypes.data, 4 * n)
        if need == 0: raise RuntimeError(f"tims_read_scans_v2 failed (frame {frame})")
        if need <= 4 * n: break
        n = need // 4 + 1
    counts = buf[:s1 - s0].astype(np.int64); pos = s1 - s0; idx, inten = [], []
    for c in counts:
        idx.append(buf[pos:pos + c]); inten.append(buf[pos + c:pos + 2 * c]); pos += 2 * c
    return (np.concatenate(idx) if idx else np.zeros(0, np.uint32)), (np.concatenate(inten) if inten else np.zeros(0, np.uint32))

def k0_table(frame, nscans):
    if frame not in _k0: _k0[frame] = convert(_dll.tims_scannum_to_oneoverk0, frame, np.arange(nscans, dtype=np.float64))
    return _k0[frame]

def init(args):
    global A, _dll, _h
    A = args; _dll, _h = open_tims(args.timsdata, args.d)

def work(chunk):
    """chunk: list of precursor dicts; returns column dict of ion rows"""
    out = {k: [] for k in ("pid", "modseq", "precursor_charge", "ion_type", "ordinal", "fragment_charge", "mz", "intensity",
                           "apex_intensity", "observed_mz", "ambiguous", "shared", "n_scans", "apex_rt", "window_lo", "window_hi",
                           "apex_im", "im_lo", "im_hi")}
    for p in chunk:
        mzs_all, ints_all, apex_mz, apex_int = [], [], None, None
        for f in p["frames"]:
            k0 = k0_table(f, A.nscans)
            s0 = int(np.searchsorted(-k0, -(p["im"] + A.im_halfwidth))); s1 = int(np.searchsorted(-k0, -(p["im"] - A.im_halfwidth), side="right"))
            s0, s1 = max(s0, p["scan_begin"]), min(s1, p["scan_end"])
            if s1 <= s0: continue
            ti, it = read_scans(f, s0, s1)
            if not len(ti): continue
            u, inv = np.unique(ti, return_inverse=True); s = np.bincount(inv, weights=it.astype(np.float64))
            mz = convert(_dll.tims_index_to_mz, f, u.astype(np.float64))
            mzs_all.append(mz); ints_all.append(s)
            if f == p["apex_frame"]: apex_mz, apex_int = mz, s
        if mzs_all:
            mz = np.concatenate(mzs_all); it = np.concatenate(ints_all); o = np.argsort(mz); mz, it = mz[o], it[o]
        else:
            mz, it = np.zeros(0), np.zeros(0)
        cs = np.concatenate([[0.0], np.cumsum(it)])
        theo = fie.theoretical_ions(p["rm"], min(A.max_fragment_charge, p["z"]), A.mz_lo, A.mz_hi)
        tm = np.array([t[3] for t in theo]); order = np.argsort(tm); amb = np.zeros(len(theo), bool); stm = tm[order]
        close = np.diff(stm) <= stm[:-1] * 2 * A.ppm * 1e-6; amb[order[:-1][close]] = True; amb[order[1:][close]] = True
        tol = tm * A.ppm * 1e-6
        a = np.searchsorted(mz, tm - tol); b = np.searchsorted(mz, tm + tol, side="right"); tot = cs[b] - cs[a]
        if apex_mz is not None:
            ca = np.concatenate([[0.0], np.cumsum(apex_int)])
            ap = ca[np.searchsorted(apex_mz, tm + tol, side="right")] - ca[np.searchsorted(apex_mz, tm - tol)]
        else: ap = np.zeros(len(tm))
        omz = p["other_mz"]
        for i, (t, o_, z, m) in enumerate(theo):
            shared = False
            if tot[i] > 0 and len(omz):
                j = np.searchsorted(omz, m - tol[i]); shared = bool(j < len(omz) and omz[j] <= m + tol[i])
            if tot[i] > 0:
                seg = slice(a[i], b[i]); obs_mz = float(np.sum(mz[seg] * it[seg]) / tot[i])
            else: obs_mz = float("nan")
            for k, v in (("pid", p["pid"]), ("modseq", p["seq"]), ("precursor_charge", p["z"]), ("ion_type", t), ("ordinal", o_),
                         ("fragment_charge", z), ("mz", m), ("intensity", float(tot[i])), ("apex_intensity", float(ap[i])),
                         ("observed_mz", obs_mz), ("ambiguous", bool(amb[i])), ("shared", shared), ("n_scans", len(p["frames"])),
                         ("apex_rt", p["apex_rt"]), ("window_lo", p["wlo"]), ("window_hi", p["whi"]), ("apex_im", p["im"]),
                         ("im_lo", p["im"] - A.im_halfwidth), ("im_hi", p["im"] + A.im_halfwidth)):
                out[k].append(v)
    return out

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--report", required=True); ap.add_argument("--d", required=True); ap.add_argument("--timsdata", required=True)
    ap.add_argument("--out", required=True); ap.add_argument("--also-in", action="append", default=[])
    ap.add_argument("--ppm", type=float, default=15.0); ap.add_argument("--im-halfwidth", type=float, default=0.02)
    ap.add_argument("--apex-halfwidth", type=int, default=1); ap.add_argument("--max-fragment-charge", type=int, default=3)
    ap.add_argument("--qvalue", type=float, default=0.01); ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--limit", type=int, default=0, help="test: only the first N placed precursors")
    args = ap.parse_args(argv)
    if not os.path.isfile(os.path.join(args.d, "analysis.tdf")):
        raise SystemExit(f"{args.d}: not a Bruker .d directory (no analysis.tdf)")
    db = sqlite3.connect(f"file:{args.d}/analysis.tdf?mode=ro", uri=True)
    g = dict(db.execute("select Key, Value from GlobalMetadata"))
    args.mz_lo, args.mz_hi = float(g["MzAcqRangeLower"]), float(g["MzAcqRangeUpper"])
    fr = np.array(db.execute("select f.Id, f.Time, i.WindowGroup from Frames f join DiaFrameMsMsInfo i on f.Id = i.Frame order by f.Id").fetchall())
    args.nscans = int(db.execute("select max(NumScans) from Frames").fetchone()[0])
    win = db.execute("select WindowGroup, ScanNumBegin, ScanNumEnd, IsolationMz, IsolationWidth from DiaFrameMsMsWindows").fetchall()
    init(args)
    # 1/K0 of every scan, from the first diaPASEF frame (used only to place precursors in windows)
    k0 = k0_table(int(fr[0, 0]), args.nscans)
    W = [dict(g=int(w[0]), s0=int(w[1]), s1=int(w[2]), lo=w[3] - w[4] / 2, hi=w[3] + w[4] / 2, k_hi=float(k0[int(w[1])]), k_lo=float(k0[int(w[2]) - 1])) for w in win]
    gframes = {gg: (fr[fr[:, 2] == gg, 0].astype(int), fr[fr[:, 2] == gg, 1] / 60.0) for gg in sorted({w["g"] for w in W})}
    rows = fie.read_report(args.report, args.qvalue)
    t = pq.read_table(args.report, columns=["Precursor.Id", "IM"]).to_pydict()
    im = dict(zip(t["Precursor.Id"], t["IM"]))
    keep = None
    for other in args.also_in:
        names = pq.read_schema(other).names
        o = pq.read_table(other, columns=["Precursor.Id", "Q.Value"] + (["Decoy"] if "Decoy" in names else [])).to_pydict()
        dec = o.get("Decoy", [0] * len(o["Precursor.Id"]))
        s = {pid for pid, q, d in zip(o["Precursor.Id"], o["Q.Value"], dec) if q <= args.qvalue and not d}
        keep = s if keep is None else keep & s
    prec, skipped = [], 0
    for r in rows:
        if keep is not None and r["Precursor.Id"] not in keep: continue
        try: rm = fie.residue_masses(r["Modified.Sequence"])
        except KeyError: skipped += 1; continue
        z = int(r["Precursor.Charge"]); pmz = (rm.sum() + fie.H2O + z * fie.PROTON) / z; pim = float(im[r["Precursor.Id"]])
        cand = [w for w in W if w["lo"] <= pmz <= w["hi"] and w["k_lo"] <= pim <= w["k_hi"]]
        if not cand: skipped += 1; continue
        w = min(cand, key=lambda w: abs((w["lo"] + w["hi"]) / 2 - pmz) + abs((w["k_lo"] + w["k_hi"]) / 2 - pim))
        ids, rts = gframes[w["g"]]; k = int(np.argmin(np.abs(rts - r["RT"])))
        sel = ids[max(0, k - args.apex_halfwidth):k + args.apex_halfwidth + 1]
        prec.append(dict(pid=fie.norm_id(r["Precursor.Id"]), seq=r["Modified.Sequence"], z=z, rm=rm, pmz=pmz, im=pim, win=W.index(w),
                         scan_begin=w["s0"], scan_end=w["s1"], frames=[int(x) for x in sel], apex_frame=int(ids[k]), apex_rt=float(rts[k]),
                         rt0=r["RT.Start"], rt1=r["RT.Stop"], t0=float(rts[max(0, k - args.apex_halfwidth)]), t1=float(rts[min(len(rts) - 1, k + args.apex_halfwidth)]),
                         wlo=w["lo"], whi=w["hi"]))
    if args.limit: prec = prec[:args.limit]
    print(f"report {len(rows)} precursors at q <= {args.qvalue}; placed {len(prec)}; skipped {skipped}; filtered by --also-in: {keep is not None}", flush=True)
    by_win = {}
    for j, p in enumerate(prec): by_win.setdefault(p["win"], []).append(j)
    theo3 = [np.array([t[3] for t in fie.theoretical_ions(p["rm"], min(3, p["z"]))]) for p in prec]
    for j, p in enumerate(prec):
        others = [q for q in by_win[p["win"]] if q != j and prec[q]["rt0"] <= p["t1"] and prec[q]["rt1"] >= p["t0"]
                  and abs(prec[q]["im"] - p["im"]) <= 2 * args.im_halfwidth]
        p["other_mz"] = np.sort(np.concatenate([theo3[q] for q in others])) if others else np.zeros(0)
    chunks = [prec[i::args.workers * 8] for i in range(args.workers * 8)]
    with Pool(args.workers, initializer=init, initargs=(args,)) as pool:
        parts = pool.map(work, chunks)
    cols = {k: sum((x[k] for x in parts), []) for k in parts[0]}
    tab = pa.table(cols).replace_schema_metadata({"ppm": str(args.ppm), "im_halfwidth": str(args.im_halfwidth), "apex_halfwidth": str(args.apex_halfwidth),
                                                  "max_fragment_charge": str(args.max_fragment_charge), "qvalue": str(args.qvalue), "report": args.report, "d": args.d})
    pq.write_table(tab, args.out)
    obs = np.array(cols["intensity"]) > 0
    print(f"wrote {args.out}: {len(prec)} precursors, {len(obs)} theoretical ions, {int(obs.sum())} observed ({obs.mean():.1%})", flush=True)

if __name__ == "__main__":
    main()
