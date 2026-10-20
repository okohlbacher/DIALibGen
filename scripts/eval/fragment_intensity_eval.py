#!/usr/bin/env python3
"""Predicted against observed fragment intensities, observed in the raw MS2 spectra.

The project's standard fragment-intensity benchmark (docs/testing.md, "Fragment-intensity accuracy"). It replaces
scoring against a search engine's fragment export, which carries only the ions that engine chose to quantify
(DIA-NN 1.7.12: six singly charged y ions per precursor) and so cannot see b ions or higher fragment charges.

Two steps, both deterministic, numpy + pyarrow only:

  extract   For every precursor of a DIA-NN report at Q.Value <= --qvalue: find the MS2 scans whose isolation
            window holds its m/z, take the scan closest to its apex RT and the --apex-halfwidth scans either
            side in the same window, and read every theoretical b and y ion (no neutral loss; fragment charge
            1 .. min(--max-fragment-charge, precursor charge); ordinals 1 .. n-1; inside the scan range) as the
            most intense centroid within --ppm, summed over those scans. Unobserved ions are kept with
            intensity 0. Two flags per ion:
              ambiguous  another theoretical ion of the same precursor lies within 2 x ppm (excluded always)
              shared     the matched peak is also within ppm of a theoretical b/y ion (z <= 3) of another
                         identified precursor co-eluting (RT.Start..RT.Stop covers a used scan) in the same window
            Output: one Parquet row per theoretical ion.

  score     Scores one or more libraries (and optionally a replicate run's extraction, as if it were a
            predictor) against an extraction, on identical precursors and identical ion sets:
              ion sets      union  -- the b/y no-loss ions any of the given libraries stores
                            theo2  -- every theoretical b/y ion, z <= min(2, precursor z)
                            theo3  -- every theoretical b/y ion, z <= min(3, precursor z)
              conventions   zero     -- unobserved ions enter as observed 0
                            excluded -- only observed ions (> 0) enter
              variants      all, unshared (ions flagged shared removed)
            An ion a library does not store is predicted 0 under every convention; a predictor with no intensity
            on any scored ion gets cosine 0 (a miss), never a dropped precursor. A precursor is scored when
            it is in every predictor and has >= 3 observed ions in the set (so every predictor, convention and
            set within a variant is scored on the same precursors). Cosine and spectral angle
            (1 - 2 acos(cos) / pi).
            Outputs: <prefix>.summary.tsv, <prefix>.breakdown.tsv (ion type x fragment charge x precursor
            charge: observable and observed ions, share of observed and of predicted intensity, cosine within
            the group), <prefix>.per_precursor.parquet.

Library layouts read: DIALibGen Parquet (compact, list columns), DIA-NN 2.x Parquet (one row per fragment),
DIA-NN / DIALibGen TSV. Precursor ids are compared after mapping modification names to UniMod ids.

usage:
  fragment_intensity_eval.py extract --report report.tsv|parquet --mzml run.mzML --out run.ions.parquet
         [--ppm 10] [--apex-halfwidth 1] [--max-fragment-charge 3] [--qvalue 0.01]
  fragment_intensity_eval.py score --ions run.ions.parquet --library NAME=PATH [--library NAME=PATH ...]
         [--replicate NAME=run2.ions.parquet] --out-prefix PREFIX
"""
import argparse, base64, json, math, re, sys, zlib
import numpy as np
import pyarrow as pa, pyarrow.csv as pcsv, pyarrow.compute as pc, pyarrow.parquet as pq

PROTON, H2O = 1.007276467, 18.010564684
AA = dict(A=71.037113805, R=156.101111050, N=114.042927470, D=115.026943065, C=103.009184505, E=129.042593135,
          Q=128.058577540, G=57.021463735, H=137.058911875, I=113.084064015, L=113.084064015, K=128.094963050,
          M=131.040484645, F=147.068413945, P=97.052763875, S=87.032028435, T=101.047678505, W=186.079312980,
          Y=163.063328575, V=99.068413945, U=150.953633405)
UNIMOD = {"1": 42.010564684, "4": 57.021463735, "7": 0.984015595, "21": 79.966330927, "35": 15.994914620,
          "121": 114.042927470}
NAME2UNIMOD = {"Acetyl": "1", "Carbamidomethyl": "4", "Deamidated": "7", "Phospho": "21", "Oxidation": "35",
               "GG": "121"}
GROUPS = [(t, z) for t in "by" for z in (1, 2, 3)]


def norm_id(s):
    """Precursor id with modification names mapped to UniMod ids, so DIALibGen and DIA-NN ids compare."""
    return re.sub(r"\(([A-Za-z]+)\)", lambda m: f"(UniMod:{NAME2UNIMOD[m.group(1)]})"
                  if m.group(1) in NAME2UNIMOD else m.group(0), s)


def residue_masses(modseq):
    """Residue masses with modifications; an N-terminal modification is added to the first residue."""
    masses, pending, i = [], 0.0, 0
    s = modseq.strip("_")
    while i < len(s):
        c = s[i]
        if c in "([":
            j = s.index(")" if c == "(" else "]", i)
            name = s[i + 1:j]
            key = name.split(":", 1)[1] if name.startswith("UniMod:") else NAME2UNIMOD.get(name)
            if key not in UNIMOD:
                raise KeyError(name)
            if masses:
                masses[-1] += UNIMOD[key]
            else:
                pending += UNIMOD[key]
            i = j + 1
            continue
        if c.isalpha():
            masses.append(AA[c] + pending); pending = 0.0
        i += 1
    return np.array(masses)


def theoretical_ions(rm, zmax, mz_lo=-1.0, mz_hi=1e9):
    """[(type, ordinal, z, mz)] for b and y, ordinals 1..n-1, z 1..zmax, inside [mz_lo, mz_hi]."""
    n = len(rm)
    pre = np.cumsum(rm)[:-1]
    suf = np.cumsum(rm[::-1])[:-1]
    out = []
    for z in range(1, zmax + 1):
        for o in range(1, n):
            b = (pre[o - 1] + z * PROTON) / z
            y = (suf[o - 1] + H2O + z * PROTON) / z
            if mz_lo <= b <= mz_hi:
                out.append(("b", o, z, b))
            if mz_lo <= y <= mz_hi:
                out.append(("y", o, z, y))
    return out


# ---- mzML ----------------------------------------------------------------------------------------------------------
_ACC = re.compile(r'accession="(MS:\d+)"')
_VAL = re.compile(r'\bvalue="([^"]*)"')
_UNIT = re.compile(r'unitName="([^"]*)"')


def read_ms2(path):
    """All MS2 spectra: dict of arrays (rt in minutes, isolation window) and per-spectrum peak arrays.

    A line-oriented reader for indexed or plain mzML (one element per line, as ProteoWizard, OpenMS and
    ThermoRawFileParser write it); supports 32/64-bit floats, zlib or no compression. Numpress is refused."""
    meta = {"index": [], "rt": [], "lo": [], "hi": [], "target": [], "scan_lo": [], "scan_hi": []}
    peaks = []
    cur = None
    arr_kind = arr_bits = arr_zlib = None
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            s = line.lstrip()
            if s.startswith("<spectrum "):
                m = re.search(r'\bindex="(\d+)"', s)
                cur = {"index": int(m.group(1)), "ms": None, "rt": None, "target": None, "lo": 0.0, "hi": 0.0,
                       "scan_lo": None, "scan_hi": None, "mz": None, "int": None}
                continue
            if cur is None:
                continue
            if s.startswith("<cvParam"):
                m = _ACC.search(s)
                if not m:
                    continue
                acc = m.group(1)
                v, u = _VAL.search(s), _UNIT.search(s)
                val, unit = (v.group(1) if v else None), (u.group(1) if u else None)
                if acc == "MS:1000511":
                    cur["ms"] = int(val)
                    if cur["ms"] != 2:            # skip MS1 quickly
                        cur = None
                elif acc == "MS:1000016":
                    cur["rt"] = float(val) / (60.0 if unit in ("second", None) else 1.0)
                elif acc == "MS:1000827":
                    cur["target"] = float(val)
                elif acc == "MS:1000828":
                    cur["lo"] = float(val)
                elif acc == "MS:1000829":
                    cur["hi"] = float(val)
                elif acc == "MS:1000501":
                    cur["scan_lo"] = float(val)
                elif acc == "MS:1000500":
                    cur["scan_hi"] = float(val)
                elif acc == "MS:1000514":
                    arr_kind = "mz"
                elif acc == "MS:1000515":
                    arr_kind = "int"
                elif acc == "MS:1000523":
                    arr_bits = 64
                elif acc == "MS:1000521":
                    arr_bits = 32
                elif acc == "MS:1000574":
                    arr_zlib = True
                elif acc == "MS:1000576":
                    arr_zlib = False
                elif acc in ("MS:1002312", "MS:1002313", "MS:1002314"):
                    raise SystemExit("numpress-compressed mzML is not supported; convert without numpress")
                continue
            if s.startswith("<binaryDataArray "):
                arr_kind = arr_bits = arr_zlib = None
                continue
            if s.startswith("<binary>"):
                raw = base64.b64decode(s[8:s.index("</binary>")])
                if arr_zlib and raw:          # an empty array (encodedLength="0", SCIEX mzML) has no zlib stream
                    raw = zlib.decompress(raw)
                v = np.frombuffer(raw, dtype=np.float64 if arr_bits == 64 else np.float32)
                if arr_kind in ("mz", "int"):
                    cur[arr_kind] = v
                continue
            if s.startswith("</spectrum>"):
                if cur["ms"] == 2 and cur["target"] is not None:
                    for k in ("index", "rt", "target", "scan_lo", "scan_hi"):
                        meta[k].append(cur[k] if cur[k] is not None else np.nan)
                    meta["lo"].append(cur["target"] - cur["lo"]); meta["hi"].append(cur["target"] + cur["hi"])
                    mz = cur["mz"] if cur["mz"] is not None else np.zeros(0)
                    it = cur["int"] if cur["int"] is not None else np.zeros(0)
                    order = np.argsort(mz, kind="stable")
                    peaks.append((np.asarray(mz, np.float64)[order], np.asarray(it, np.float64)[order]))
                cur = None
    return {k: np.asarray(v) for k, v in meta.items()}, peaks


def best_peak(mz_arr, int_arr, target, ppm):
    """(intensity, mz) of the most intense centroid within ppm of target, or (0, nan)."""
    tol = target * ppm * 1e-6
    a, b = np.searchsorted(mz_arr, target - tol), np.searchsorted(mz_arr, target + tol, side="right")
    if b <= a:
        return 0.0, float("nan")
    k = a + int(np.argmax(int_arr[a:b]))
    return float(int_arr[k]), float(mz_arr[k])


# ---- report -------------------------------------------------------------------------------------------------------
def read_report(path, qvalue):
    cols = ["Precursor.Id", "Modified.Sequence", "Precursor.Charge", "Q.Value", "RT", "RT.Start", "RT.Stop",
            "Precursor.Mz"]
    if path.endswith(".parquet"):
        t = pq.read_table(path, columns=cols)
    else:
        t = pcsv.read_csv(path, parse_options=pcsv.ParseOptions(delimiter="\t"),
                          convert_options=pcsv.ConvertOptions(include_columns=cols))
    t = t.filter(pc.less_equal(t["Q.Value"], qvalue))
    seen, rows = set(), []
    for r in t.to_pylist():
        if r["Precursor.Id"] in seen:
            continue
        seen.add(r["Precursor.Id"]); rows.append(r)
    return rows


def extract(args):
    rows = read_report(args.report, args.qvalue)
    print(f"report: {len(rows)} precursors at Q.Value <= {args.qvalue}", flush=True)
    meta, peaks = read_ms2(args.mzml)
    nscan = len(meta["rt"])
    print(f"mzML: {nscan} MS2 spectra", flush=True)
    # isolation windows: scans grouped by window, each group in RT order
    win_key = np.round(meta["target"], 3)
    windows = {}
    for i, w in enumerate(win_key):
        windows.setdefault(w, []).append(i)
    wins = sorted(windows)
    win_lo = np.array([np.median(meta["lo"][windows[w]]) for w in wins])
    win_hi = np.array([np.median(meta["hi"][windows[w]]) for w in wins])
    win_scans = [np.array(sorted(windows[w], key=lambda i: meta["rt"][i])) for w in wins]
    win_rt = [meta["rt"][s] for s in win_scans]
    print(f"isolation windows: {len(wins)} (median width {np.median(win_hi - win_lo):.3f} Th)", flush=True)

    # precursors: masses, window, apex scans
    prec, skipped = [], 0
    for r in rows:
        try:
            rm = residue_masses(r["Modified.Sequence"])
        except KeyError:
            skipped += 1; continue
        z = int(r["Precursor.Charge"])
        pmz = (rm.sum() + H2O + z * PROTON) / z
        inside = np.nonzero((win_lo <= pmz) & (pmz <= win_hi))[0]
        if len(inside) == 0:
            skipped += 1; continue
        w = int(inside[np.argmin(np.abs((win_lo[inside] + win_hi[inside]) / 2 - pmz))])
        k = int(np.argmin(np.abs(win_rt[w] - r["RT"])))
        lo, hi = max(0, k - args.apex_halfwidth), min(len(win_scans[w]), k + args.apex_halfwidth + 1)
        prec.append(dict(pid=norm_id(r["Precursor.Id"]), seq=r["Modified.Sequence"], z=z, rm=rm, pmz=pmz, win=w,
                         apex=int(win_scans[w][k]), scans=win_scans[w][lo:hi], rt=r["RT"], rt0=r["RT.Start"],
                         rt1=r["RT.Stop"], report_mz=r["Precursor.Mz"]))
    print(f"precursors placed: {len(prec)}, skipped (unknown modification or no window): {skipped}", flush=True)
    if prec:
        drt = np.array([abs(meta["rt"][p["apex"]] - p["rt"]) for p in prec]) * 60.0
        print(f"apex scan vs report RT: median |dRT| {np.median(drt):.2f} s, max {drt.max():.2f} s", flush=True)
    dmz = np.array([abs(p["pmz"] - p["report_mz"]) for p in prec if p["report_mz"] == p["report_mz"]])
    if len(dmz):
        print(f"computed vs report precursor m/z: max |diff| {dmz.max():.4f}", flush=True)

    # co-eluting identified precursors per window, for the shared flag
    by_win = {}
    for j, p in enumerate(prec):
        by_win.setdefault(p["win"], []).append(j)
    theo_all = [theoretical_ions(p["rm"], min(3, p["z"])) for p in prec]

    out = {k: [] for k in ("pid", "modseq", "precursor_charge", "ion_type", "ordinal", "fragment_charge", "mz",
                           "intensity", "apex_intensity", "observed_mz", "ambiguous", "shared", "n_scans",
                           "apex_rt", "window_lo", "window_hi")}
    for j, p in enumerate(prec):
        scans = p["scans"]
        rts = meta["rt"][scans]
        scan_lo = np.nanmax([meta["scan_lo"][scans].min(), 0.0])
        scan_hi = np.nanmin([meta["scan_hi"][scans].max(), 1e9])
        theo = [t for t in theoretical_ions(p["rm"], min(args.max_fragment_charge, p["z"]), scan_lo, scan_hi)]
        mzs = np.array([t[3] for t in theo])
        order = np.argsort(mzs)
        amb = np.zeros(len(theo), bool)
        smz = mzs[order]
        close = np.diff(smz) <= smz[:-1] * 2 * args.ppm * 1e-6
        amb[order[:-1][close]] = True; amb[order[1:][close]] = True
        others = [q for q in by_win.get(p["win"], []) if q != j and prec[q]["rt0"] <= rts.max()
                  and prec[q]["rt1"] >= rts.min()]
        omz = np.sort(np.array([t[3] for q in others for t in theo_all[q]])) if others else np.zeros(0)
        for i, (t, o, z, mz) in enumerate(theo):
            tot, apex_int, best_mz, best_int = 0.0, 0.0, float("nan"), -1.0
            for s in scans:
                v, m = best_peak(peaks[s][0], peaks[s][1], mz, args.ppm)
                tot += v
                if s == p["apex"]:
                    apex_int = v
                if v > best_int and v > 0:
                    best_int, best_mz = v, m
            shared = False
            if tot > 0 and len(omz):
                tol = best_mz * args.ppm * 1e-6
                a = np.searchsorted(omz, best_mz - tol)
                shared = a < len(omz) and omz[a] <= best_mz + tol
            for k, v in (("pid", p["pid"]), ("modseq", p["seq"]), ("precursor_charge", p["z"]), ("ion_type", t),
                         ("ordinal", o), ("fragment_charge", z), ("mz", mz), ("intensity", tot),
                         ("apex_intensity", apex_int), ("observed_mz", best_mz), ("ambiguous", bool(amb[i])),
                         ("shared", bool(shared)), ("n_scans", len(scans)), ("apex_rt", float(meta["rt"][p["apex"]])),
                         ("window_lo", float(win_lo[p["win"]])), ("window_hi", float(win_hi[p["win"]]))):
                out[k].append(v)
    t = pa.table(out)
    t = t.replace_schema_metadata({"ppm": str(args.ppm), "apex_halfwidth": str(args.apex_halfwidth),
                                   "max_fragment_charge": str(args.max_fragment_charge), "qvalue": str(args.qvalue),
                                   "report": args.report, "mzml": args.mzml})
    pq.write_table(t, args.out)
    obs = np.array(out["intensity"]) > 0
    print(f"wrote {args.out}: {len(prec)} precursors, {len(obs)} theoretical ions, {obs.sum()} observed", flush=True)


# ---- libraries ----------------------------------------------------------------------------------------------------
TYPE_ENUM = {2: "b", 5: "y"}          # DIALibGen FragmentType: Unknown=0, A, B, C, X, Y, Z


def read_library(path, wanted):
    """{precursor id: {(type, ordinal, z): intensity}} for b/y no-loss ions of the wanted (normalised) ids."""
    lib = {}
    if path.endswith(".parquet"):
        names = pq.read_schema(path).names
        cols = ["Precursor.Id", "Relative.Intensity", "Fragment.Type", "Fragment.Charge", "Fragment.Series.Number",
                "Fragment.Loss.Type"] + (["Decoy"] if "Decoy" in names else [])
        t = pq.read_table(path, columns=cols)
        # normalise each distinct id once, then filter in Arrow: a library has millions of rows
        raw = [u for u in pc.unique(t["Precursor.Id"]).to_pylist() if norm_id(u) in wanted]
        t = t.filter(pc.is_in(t["Precursor.Id"], value_set=pa.array(raw, t["Precursor.Id"].type)))
        compact = pa.types.is_list(t.schema.field("Fragment.Type").type)
        for r in t.to_pylist():
            if r.get("Decoy"):
                continue
            pid = norm_id(r["Precursor.Id"])
            if compact:
                d = lib.setdefault(pid, {})
                for it, ty, ch, se, lo in zip(r["Relative.Intensity"], r["Fragment.Type"], r["Fragment.Charge"],
                                              r["Fragment.Series.Number"], r["Fragment.Loss.Type"]):
                    if lo == 0 and ty in TYPE_ENUM:
                        d[(TYPE_ENUM[ty], int(se), int(ch))] = float(it)
            elif r["Fragment.Loss.Type"] in ("noloss", "", None) and r["Fragment.Type"] in ("b", "y"):
                lib.setdefault(pid, {})[(r["Fragment.Type"], int(r["Fragment.Series.Number"]),
                                         int(r["Fragment.Charge"]))] = float(r["Relative.Intensity"])
    else:
        opts = pcsv.ConvertOptions(include_columns=["Precursor.Id", "LibraryIntensity", "FragmentType",
                                                    "FragmentCharge", "FragmentSeriesNumber", "FragmentLossType",
                                                    "Decoy"], include_missing_columns=True)
        reader = pcsv.open_csv(path, parse_options=pcsv.ParseOptions(delimiter="\t"), convert_options=opts,
                               read_options=pcsv.ReadOptions(block_size=1 << 26))
        for batch in reader:
            for pid, it, ty, ch, se, lo, dec in zip(*(batch.column(i).to_pylist() for i in range(7))):
                if dec or ty not in ("b", "y") or lo not in ("noloss", "", None):
                    continue
                pid = norm_id(pid)
                if pid in wanted:
                    lib.setdefault(pid, {})[(ty, int(se), int(ch))] = float(it)
    return lib


def read_ions(path):
    t = pq.read_table(path).to_pydict()
    ions = {}
    for i in range(len(t["pid"])):
        p = ions.setdefault(t["pid"][i], {"z": t["precursor_charge"][i], "ions": {}})
        p["ions"][(t["ion_type"][i], t["ordinal"][i], t["fragment_charge"][i])] = (
            t["intensity"][i], t["ambiguous"][i], t["shared"][i])
    return ions


# ---- scoring ------------------------------------------------------------------------------------------------------
def cosine(a, b):
    """Cosine of observed a and predicted b. A prediction with no intensity on observed ions scores 0 (a complete
    miss, so a library that stores none of them is penalised, not dropped); nothing observed is undefined (nan)."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    if na == 0:
        return float("nan")
    return float(a @ b / (na * nb)) if nb > 0 else 0.0


def angle(c):
    return 1 - 2 * math.acos(min(1.0, max(-1.0, c))) / math.pi if c == c else float("nan")


def ion_set(name, p, union):
    z = p["z"]
    keys = [k for k, (v, amb, sh) in p["ions"].items() if not amb]
    if name == "union":
        return [k for k in keys if k in union]
    zmax = min(2 if name == "theo2" else 3, z)
    return [k for k in keys if k[2] <= zmax]


def score_precursor(p, predictors, keys, convention, unshared):
    """{predictor: cosine} on one precursor, or None when < 3 observed ions remain."""
    if unshared:
        keys = [k for k in keys if not p["ions"][k][2]]
    if convention == "excluded":
        keys = [k for k in keys if p["ions"][k][0] > 0]
    if sum(p["ions"][k][0] > 0 for k in keys) < 3:
        return None
    o = [p["ions"][k][0] for k in keys]
    return {name: cosine(o, [pred.get(k, 0.0) for k in keys]) for name, pred in predictors.items()}, len(keys)


def summarise(x):
    x = np.asarray([v for v in x if v == v], float)
    if not len(x):
        return dict(n=0, median=float("nan"), q25=float("nan"), q75=float("nan"), ge09=float("nan"),
                    sa_median=float("nan"))
    return dict(n=int(len(x)), median=float(np.median(x)), q25=float(np.quantile(x, .25)),
                q75=float(np.quantile(x, .75)), ge09=float((x >= 0.9).mean()),
                sa_median=float(np.median([angle(c) for c in x])))


def write_tsv(path, rows):
    if not rows:
        return
    cols = list(rows[0])
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(f"{r[c]:.6g}" if isinstance(r[c], float) else str(r[c]) for c in cols) + "\n")


def score(args, extra_predictors=None):
    ions = read_ions(args.ions)
    wanted = set(ions)
    preds = {}
    for spec in args.library or []:
        name, path = spec.split("=", 1)
        preds[name] = read_library(path, wanted)
        print(f"library {name}: {len(preds[name])} of {len(wanted)} precursors", flush=True)
    for spec in args.replicate or []:
        name, path = spec.split("=", 1)
        rep = read_ions(path)
        preds[name] = {pid: {k: v[0] for k, v in p["ions"].items()} for pid, p in rep.items()}
        print(f"replicate {name}: {len(preds[name])} precursors", flush=True)
    if extra_predictors:
        preds.update(extra_predictors)
    libs = [n for n in preds if n not in {s.split("=", 1)[0] for s in (args.replicate or [])}]
    common = [pid for pid in ions if all(pid in preds[n] for n in preds)]
    print(f"precursors in every predictor: {len(common)}", flush=True)
    names = list(preds)
    summary, per = [], {"pid": common, "precursor_charge": [ions[pid]["z"] for pid in common]}
    for variant in ("all", "unshared"):
        for set_name in ("union", "theo2", "theo3"):
            for conv in ("zero", "excluded"):
                res = {n: [] for n in names}; nions = []
                for pid in common:
                    p = ions[pid]
                    union = set().union(*(preds[n][pid].keys() for n in libs)) if libs else set()
                    r = score_precursor(p, {n: preds[n][pid] for n in names}, ion_set(set_name, p, union), conv,
                                        variant == "unshared")
                    for n in names:
                        res[n].append(r[0][n] if r else float("nan"))
                    nions.append(r[1] if r else 0)
                tag = f"{variant}.{set_name}.{conv}"
                for n in names:
                    per[f"{n}.{tag}"] = res[n]
                    s = summarise(res[n])
                    ref = np.array(res[names[0]]); cur = np.array(res[n]); ok = ~np.isnan(ref) & ~np.isnan(cur)
                    summary.append(dict(variant=variant, ion_set=set_name, convention=conv, predictor=n, **s,
                                        median_ions=float(np.median([k for k in nions if k])) if any(nions) else 0.0,
                                        paired_median_diff_vs_first=float(np.median(cur[ok] - ref[ok])) if ok.any() else float("nan"),
                                        share_higher_than_first=float((cur[ok] > ref[ok]).mean()) if ok.any() else float("nan")))
    write_tsv(f"{args.out_prefix}.summary.tsv", summary)
    pq.write_table(pa.table(per), f"{args.out_prefix}.per_precursor.parquet")

    # breakdown: ion type x fragment charge x precursor charge (theo3 ions, ambiguous excluded)
    rows = []
    for zp in sorted({ions[pid]["z"] for pid in common}):
        group_pids = [pid for pid in common if ions[pid]["z"] == zp]
        tot_obs = {pid: sum(v for k, (v, a, s) in ions[pid]["ions"].items() if not a) for pid in group_pids}
        tot_pred = {n: {pid: sum(preds[n][pid].values()) for pid in group_pids} for n in names}
        for t, zf in GROUPS:
            if zf > zp:
                continue
            n_theo, n_obs, sh_obs, gcos, gnone = [], [], [], {n: [] for n in names}, {n: [] for n in names}
            sh_pred = {n: [] for n in names}
            for pid in group_pids:
                keys = [k for k, (v, a, s) in ions[pid]["ions"].items() if not a and k[0] == t and k[2] == zf]
                obs = [ions[pid]["ions"][k][0] for k in keys]
                n_theo.append(len(keys)); n_obs.append(sum(v > 0 for v in obs))
                if tot_obs[pid] > 0:
                    sh_obs.append(sum(obs) / tot_obs[pid])
                for n in names:
                    if tot_pred[n][pid] > 0:
                        sh_pred[n].append(sum(v for k, v in preds[n][pid].items() if k[0] == t and k[2] == zf)
                                          / tot_pred[n][pid])
                ok = [k for k, v in zip(keys, obs) if v > 0]
                if len(ok) >= 3:
                    for n in names:
                        pv = [preds[n][pid].get(k, 0.0) for k in ok]
                        gcos[n].append(cosine([ions[pid]["ions"][k][0] for k in ok], pv))
                        gnone[n].append(not any(v > 0 for v in pv))
            row = dict(precursor_charge=zp, ion_type=t, fragment_charge=zf, precursors=len(group_pids),
                       mean_theoretical_ions=float(np.mean(n_theo)), mean_observed_ions=float(np.mean(n_obs)),
                       share_observed_intensity=float(np.mean(sh_obs)) if sh_obs else float("nan"))
            for n in names:
                row[f"share_predicted_{n}"] = float(np.mean(sh_pred[n])) if sh_pred[n] else float("nan")
            row["n_within"] = len(gcos[names[0]]) if names else 0
            for n in names:
                # cosine within the group (precursors with >= 3 observed ions of the group); a predictor with no
                # intensity on any of them scores 0 there, and how often that happens is reported beside it
                row[f"cos_within_{n}"] = float(np.median(gcos[n])) if gcos[n] else float("nan")
                row[f"none_predicted_{n}"] = float(np.mean(gnone[n])) if gnone[n] else float("nan")
            rows.append(row)
    write_tsv(f"{args.out_prefix}.breakdown.tsv", rows)
    for r in summary:
        if r["variant"] == "all":
            print(f"{r['ion_set']:6s} {r['convention']:8s} {r['predictor']:14s} n {r['n']:6d} median {r['median']:.4f} "
                  f"(IQR {r['q25']:.3f}-{r['q75']:.3f}) SA {r['sa_median']:.4f} ions {r['median_ions']:.0f}", flush=True)
    return summary


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    e = sub.add_parser("extract")
    e.add_argument("--report", required=True); e.add_argument("--mzml", required=True)
    e.add_argument("--out", required=True)
    e.add_argument("--ppm", type=float, default=10.0)
    e.add_argument("--apex-halfwidth", type=int, default=1)
    e.add_argument("--max-fragment-charge", type=int, default=3)
    e.add_argument("--qvalue", type=float, default=0.01)
    s = sub.add_parser("score")
    s.add_argument("--ions", required=True)
    s.add_argument("--library", action="append", help="NAME=PATH; the first is the paired reference")
    s.add_argument("--replicate", action="append", help="NAME=ions.parquet of another run of the same sample")
    s.add_argument("--out-prefix", required=True)
    a = ap.parse_args(argv)
    extract(a) if a.cmd == "extract" else score(a)


if __name__ == "__main__":
    main()
