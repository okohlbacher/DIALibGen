#!/usr/bin/env python
"""The raw-spectrum fragment-intensity benchmark (scripts/eval/fragment_intensity_eval.py) on a synthetic run.

Everything the benchmark's numbers rest on, on data whose answer is known:
  - apex-scan choice and the sum over the neighbouring scans of the same isolation window
  - the ppm tolerance (a centroid 15 ppm off is not an observation at 10 ppm)
  - unobserved ions kept as 0, both scoring conventions, ambiguous ions excluded
  - the shared flag: b ions of PEPC*TIDEK are also b ions of the co-eluting PEPC*TIDEKR in the same window
  - DIALibGen compact Parquet, DIA-NN long Parquet and TSV libraries, with modification names normalised
  - a predictor with no intensity on the scored ions scores 0, in the summary and in the per-group breakdown
    (a library without b ions, or without charge-3 fragments, must not turn a group into nan)
  - the timsTOF extractor imports without Bruker's library and refuses a missing one with a clear message

Exits 77 (ctest SKIP) without numpy or pyarrow.
"""
import base64, math, os, sys, tempfile

try:
    import numpy as np
    import pyarrow as pa, pyarrow.parquet as pq
except ImportError:
    sys.exit(77)

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts", "eval"))
import fragment_intensity_eval as fie

A, B = "PEPC(UniMod:4)TIDEK", "PEPC(UniMod:4)TIDEKR"
failures = []


def check(cond, what):
    if not cond:
        failures.append(what)
        print("FAIL:", what)


def ions(seq, zmax):
    return {(t, o, z): mz for t, o, z, mz in fie.theoretical_ions(fie.residue_masses(seq), zmax)}


def pmz(seq, z):
    return (fie.residue_masses(seq).sum() + fie.H2O + z * fie.PROTON) / z


ia, ib = ions(A, 2), ions(B, 2)
# Observed peaks of A: four y ions and two b ions (b ions shared with B), one y ion 15 ppm off.
obs_a = {("y", 3, 1): 100.0, ("y", 4, 1): 400.0, ("y", 5, 1): 250.0, ("y", 6, 1): 50.0,
         ("b", 3, 1): 80.0, ("b", 4, 1): 120.0, ("b", 5, 1): 60.0}
off_ion = ("y", 7, 1)
obs_b = {("y", 2, 1): 30.0, ("y", 3, 1): 60.0, ("y", 4, 1): 90.0}
scale = {0: 0.5, 1: 1.0, 2: 0.25}            # scan 1 is the apex


def spectrum(i):
    peaks = {}
    for k, v in obs_a.items():
        peaks[ia[k]] = v * scale[i]
    peaks[ia[off_ion] * (1 + 15e-6)] = 999.0
    for k, v in obs_b.items():
        peaks[ib[k]] = v * scale[i]
    peaks[1234.5678] = 77.0                  # unrelated
    mz = np.array(sorted(peaks)); it = np.array([peaks[m] for m in mz], np.float32)
    return mz, it


def b64(a):
    return base64.b64encode(a.tobytes()).decode()


def mzml(path):
    L = ['<?xml version="1.0" encoding="utf-8"?>', "<mzML>", "<run>", '<spectrumList count="4">']
    rts = [60.0, 66.0, 72.0]                  # seconds: 1.0, 1.1, 1.2 min
    L += ['<spectrum id="scan=1" index="0" defaultArrayLength="1">',
          '<cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="1" />', "</spectrum>"]
    for i, rt in enumerate(rts):
        mz, it = spectrum(i)
        L += [f'<spectrum id="scan={i + 2}" index="{i + 1}" defaultArrayLength="{len(mz)}">',
              '<cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="2" />',
              f'<cvParam cvRef="MS" accession="MS:1000016" name="scan start time" value="{rt}" unitAccession="UO:0000010" unitName="second" unitCvRef="UO" />',
              '<cvParam cvRef="MS" accession="MS:1000501" name="scan window lower limit" value="150" />',
              '<cvParam cvRef="MS" accession="MS:1000500" name="scan window upper limit" value="2000" />',
              '<cvParam cvRef="MS" accession="MS:1000827" name="isolation window target m/z" value="550" />',
              '<cvParam cvRef="MS" accession="MS:1000828" name="isolation window lower offset" value="100" />',
              '<cvParam cvRef="MS" accession="MS:1000829" name="isolation window upper offset" value="100" />',
              f'<binaryDataArray encodedLength="0">',
              '<cvParam cvRef="MS" accession="MS:1000514" name="m/z array" />',
              '<cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" />',
              '<cvParam cvRef="MS" accession="MS:1000576" name="no compression" />',
              f"<binary>{b64(mz)}</binary>", "</binaryDataArray>",
              f'<binaryDataArray encodedLength="0">',
              '<cvParam cvRef="MS" accession="MS:1000515" name="intensity array" />',
              '<cvParam cvRef="MS" accession="MS:1000521" name="32-bit float" />',
              '<cvParam cvRef="MS" accession="MS:1000576" name="no compression" />',
              f"<binary>{b64(it)}</binary>", "</binaryDataArray>", "</spectrum>"]
    L += ["</spectrumList>", "</run>", "</mzML>"]
    open(path, "w").write("\n".join(L) + "\n")


with tempfile.TemporaryDirectory() as d:
    mzml(f"{d}/run.mzML")
    check(450 < pmz(A, 2) < 650 and 450 < pmz(B, 2) < 650, "both precursors inside the 450-650 window")
    with open(f"{d}/report.tsv", "w") as fh:
        fh.write("Precursor.Id\tModified.Sequence\tPrecursor.Charge\tQ.Value\tRT\tRT.Start\tRT.Stop\tPrecursor.Mz\n")
        fh.write(f"{A}2\t{A}\t2\t0.001\t1.1\t1.0\t1.2\t{pmz(A, 2):.5f}\n")
        fh.write(f"{B}2\t{B}\t2\t0.002\t1.1\t1.0\t1.2\t{pmz(B, 2):.5f}\n")
        fh.write(f"ELVISK2\tELVISK\t2\t0.5\t1.1\t1.0\t1.2\t300\n")            # above the q-value cut
    fie.main(["extract", "--report", f"{d}/report.tsv", "--mzml", f"{d}/run.mzML", "--out", f"{d}/ions.parquet"])
    t = pq.read_table(f"{d}/ions.parquet").to_pylist()
    rows = {(r["pid"], r["ion_type"], r["ordinal"], r["fragment_charge"]): r for r in t}
    check({r["pid"] for r in t} == {f"{A}2", f"{B}2"}, "only precursors at Q.Value <= 0.01")
    ra = rows[(f"{A}2", "y", 4, 1)]
    check(abs(ra["intensity"] - 400.0 * 1.75) < 1e-3, f"apex +-1 scans summed: {ra['intensity']}")
    check(abs(ra["apex_intensity"] - 400.0) < 1e-3, "apex scan is the one nearest the report RT")
    check(rows[(f"{A}2", "y", 7, 1)]["intensity"] == 0.0, "a centroid 15 ppm away is not observed at 10 ppm")
    check(rows[(f"{A}2", "y", 2, 1)]["intensity"] == 0.0, "unobserved ions kept with intensity 0")
    check(rows[(f"{A}2", "b", 3, 1)]["shared"] and rows[(f"{A}2", "b", 4, 1)]["shared"],
          "b ions shared with the co-eluting PEPC*TIDEKR")
    check(not any(rows[(f"{A}2", "y", o, 1)]["shared"] for o in (3, 4, 5, 6)), "y ions of A are not shared")
    check(all(r["fragment_charge"] <= 2 for r in t), "fragment charge capped at the precursor charge")

    # libraries: "exact" stores A's observed ions with the observed ratios; "flat" stores them all equal
    exact = {k: v for k, v in obs_a.items()}
    flat = {k: 1.0 for k in obs_a}
    pa_ids = [f"PEPC(Carbamidomethyl)TIDEK2", f"PEPC(Carbamidomethyl)TIDEKR2"]
    enum = {"b": 2, "y": 5}

    def compact(lib_a, lib_b):
        cols = {"Precursor.Id": pa_ids, "Relative.Intensity": [], "Fragment.Type": [], "Fragment.Charge": [],
                "Fragment.Series.Number": [], "Fragment.Loss.Type": []}
        for lib in (lib_a, lib_b):
            ks = sorted(lib)
            cols["Relative.Intensity"].append([float(lib[k]) for k in ks])
            cols["Fragment.Type"].append([enum[k[0]] for k in ks])
            cols["Fragment.Series.Number"].append([k[1] for k in ks])
            cols["Fragment.Charge"].append([k[2] for k in ks])
            cols["Fragment.Loss.Type"].append([0] * len(ks))
        return pa.table(cols)

    pq.write_table(compact(exact, obs_b), f"{d}/exact.parquet")
    long = {"Precursor.Id": [], "Decoy": [], "Relative.Intensity": [], "Fragment.Type": [], "Fragment.Charge": [],
            "Fragment.Series.Number": [], "Fragment.Loss.Type": []}
    for pid, lib in ((f"{A}2", flat), (f"{B}2", {k: 1.0 for k in obs_b})):
        for k, v in lib.items():
            for c, x in zip(long, (pid, 0, v, k[0], k[2], k[1], "noloss")):
                long[c].append(x)
    pq.write_table(pa.table(long), f"{d}/flat.parquet")
    with open(f"{d}/flat.tsv", "w") as fh:
        fh.write("Precursor.Id\tLibraryIntensity\tFragmentType\tFragmentCharge\tFragmentSeriesNumber\tFragmentLossType\tDecoy\n")
        for k, v in flat.items():
            fh.write(f"{A}2\t{v}\t{k[0]}\t{k[2]}\t{k[1]}\tnoloss\t0\n")
        for k in obs_b:
            fh.write(f"{B}2\t1.0\t{k[0]}\t{k[2]}\t{k[1]}\tnoloss\t0\n")

    # "yonly" stores no b ion; "miss" stores, for A, only an ion that is not observed
    pq.write_table(compact({k: v for k, v in obs_a.items() if k[0] == "y"}, obs_b), f"{d}/yonly.parquet")
    pq.write_table(compact({("y", 2, 1): 1.0}, obs_b), f"{d}/miss.parquet")
    fie.main(["score", "--ions", f"{d}/ions.parquet", "--library", f"exact={d}/exact.parquet",
              "--library", f"flat={d}/flat.parquet", "--library", f"flattsv={d}/flat.tsv",
              "--library", f"yonly={d}/yonly.parquet", "--library", f"miss={d}/miss.parquet",
              "--out-prefix", f"{d}/s"])
    per = pq.read_table(f"{d}/s.per_precursor.parquet").to_pydict()
    ia_row = per["pid"].index(f"{A}2")
    get = lambda col: per[col][ia_row]
    # every observed ion of A is stored by both libraries, so with observed-only ions exact is perfect ...
    check(abs(get("exact.all.union.excluded") - 1.0) < 1e-9, f"exact library scores 1: {get('exact.all.union.excluded')}")
    check(abs(get("exact.all.theo2.zero") - 1.0) < 1e-9, "unobserved ions are 0 in the observation and in exact")
    o = np.array([v * 1.75 for v in obs_a.values()])
    want = float(o.sum() / (np.linalg.norm(o) * math.sqrt(len(o))))
    check(abs(get("flat.all.union.excluded") - want) < 1e-6, f"flat cosine {get('flat.all.union.excluded')} vs {want}")
    check(abs(get("flat.all.union.excluded") - get("flattsv.all.union.excluded")) < 1e-9, "TSV reads as Parquet")
    oy = np.array([obs_a[("y", i, 1)] for i in (3, 4, 5, 6)]) * 1.75
    want_u = float(oy.sum() / (np.linalg.norm(oy) * 2.0))
    check(abs(get("flat.unshared.union.excluded") - want_u) < 1e-6, "unshared drops the shared b ions")
    check(get("miss.all.union.excluded") == 0.0, f"no intensity on any observed ion scores 0: {get('miss.all.union.excluded')}")
    summary = open(f"{d}/s.summary.tsv").read().splitlines()
    check(len(summary) == 1 + 2 * 3 * 2 * 5, f"summary rows {len(summary)}")
    lines = open(f"{d}/s.breakdown.tsv").read().splitlines()
    check(lines[0].startswith("precursor_charge\tion_type\tfragment_charge"), "breakdown header")
    head = lines[0].split("\t")
    rows_b = [dict(zip(head, ln.split("\t"))) for ln in lines[1:]]
    b1 = [r for r in rows_b if r["precursor_charge"] == "2" and r["ion_type"] == "b" and r["fragment_charge"] == "1"]
    check(len(b1) == 1, "one b 1+ row for 2+ precursors")
    if b1:
        r = b1[0]
        # both precursors show b3-b5 (the same peaks); "exact" stores them for A (cosine 1) but not for B (0)
        check(int(r["n_within"]) == 2, f"b 1+ group scored for both precursors: {r['n_within']}")
        check(float(r["cos_within_yonly"]) == 0.0, f"a library without b ions scores 0 on the b group: {r['cos_within_yonly']}")
        check(float(r["none_predicted_yonly"]) == 1.0, "and is counted as predicting none of them")
        check(abs(float(r["cos_within_exact"]) - 0.5) < 1e-9, f"exact on the b group, median of 1 and 0: {r['cos_within_exact']}")
        check(abs(float(r["none_predicted_exact"]) - 0.5) < 1e-9, "exact predicts none of B's b ions")
        check(not any(v == "nan" for k, v in r.items() if k.startswith("cos_within_")), f"no nan in the b 1+ group: {r}")

# the timsTOF extractor: importable without Bruker's library, and a clear refusal without it
import tims_extract
try:
    tims_extract.open_tims(os.path.join(tempfile.gettempdir(), "no-such-libtimsdata.so"), "x.d")
    check(False, "a missing timsdata library is refused")
except SystemExit as e:
    check("timsdata library" in str(e) and "--timsdata" in str(e), f"refusal names the library and the option: {e}")
try:
    tims_extract.main(["--report", "r.parquet", "--d", tempfile.gettempdir(), "--timsdata", "x", "--out", "o.parquet"])
    check(False, "a directory without analysis.tdf is refused")
except SystemExit as e:
    check("analysis.tdf" in str(e), f"refusal names analysis.tdf: {e}")

if failures:
    print(f"{len(failures)} failure(s)")
    sys.exit(1)
print("fragment_intensity_eval: all checks passed")
