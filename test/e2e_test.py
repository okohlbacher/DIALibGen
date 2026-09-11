#!/usr/bin/env python3
"""End-to-end: a FASTA in, libraries out, checked for what they must contain.

Everything else in the suite checks a piece. This runs the tool the way a user
does and then reads the result back, because the failures that matter here are
not crashes -- they are a library that is written successfully and is wrong:
NaN retention times, an ion-mobility column that is not 1/K0, three "different"
decoy methods that are the same method, or a recipe that does not name the
models that produced it.

Usage: e2e_test.py <binary> <proteins.fasta> <model-dir> [workdir]
Exits 77 (ctest SKIP) when pyarrow or the models are unavailable.
"""
import json
import os
import re
import subprocess
import sys
import tempfile

try:
    import numpy as np
    import pyarrow.parquet as pq
except ImportError as e:
    print(f"SKIP: {e}", file=sys.stderr)
    sys.exit(77)

MODELS = ("peptdeep_rt_dynamic.onnx", "peptdeep_ms2_dynamic.onnx", "peptdeep_ccs_dynamic.onnx")


def main(binary, fasta, model_dir, workdir):
    missing = [m for m in MODELS if not os.path.exists(os.path.join(model_dir, m))]
    if missing:
        print(f"SKIP: no models in {model_dir}: {', '.join(missing)}", file=sys.stderr)
        return 77

    env = dict(os.environ, DIALIBGEN_MODEL_DIR=model_dir)
    os.chdir(workdir)
    failures = []

    def check(ok, msg):
        print(("  ok   " if ok else "  FAIL ") + msg)
        if not ok:
            failures.append(msg)

    def run(out, config=None):
        args = [binary, "-in", fasta, "-out", out]
        if config is not None:
            name = out + ".config.json"
            with open(name, "w") as fh:
                json.dump(config, fh)
            args += ["-config", name]
        p = subprocess.run(args, capture_output=True, text=True, env=env)
        if p.returncode != 0:
            raise SystemExit(f"FAIL {out}: exit {p.returncode}\n{p.stderr}")
        # Both streams: TOPPBase's writeLogInfo_ goes to stdout while the
        # predictor's device lines go to stderr, and the digest line read below
        # is on the former.
        return p.stdout + p.stderr

    log = run("lib.parquet")
    run("lib.tsv")
    run("irt.parquet", {"irt_rescale": True})
    run("rev.parquet", {"decoys": "reverse"})
    run("mut.parquet", {"decoys": "mutate"})

    table = pq.read_table("lib.parquet")
    meta = {k.decode(): v.decode() for k, v in (table.schema.metadata or {}).items()}
    cfg = json.loads(meta["odia.config_json"])
    d = table.to_pydict()

    # The digest line and the file must agree. They are produced by different
    # code paths, and a writer that drops rows would otherwise be invisible.
    m = re.search(r"digest:[^\n]*?(\d+) precursors", log)
    if not m:
        raise SystemExit("FAIL: no digest line in the tool's output:\n" + log)
    precursors = int(m.group(1))
    check(table.num_rows == precursors,
          f"the file holds the {precursors} precursors the digest reported")

    # The recipe must name the models ACTUALLY used, including when they were
    # found via DIALIBGEN_MODEL_DIR rather than written in the config -- a
    # library that does not say what produced it is not reproducible from itself.
    check(all(cfg[k].endswith(m) for k, m in
              (("rt_model", MODELS[0]), ("ms2_model", MODELS[1]), ("ccs_model", MODELS[2]))),
          "the embedded recipe names the models that produced the library")
    check(cfg["enzyme"] == "Trypsin/P" and cfg["fixed_modifications"] == ["Carbamidomethyl (C)"],
          "the embedded recipe carries the effective defaults")

    rt = np.array(d["RT"], dtype=float)
    check(np.isfinite(rt).all(), "no NaN retention times")
    check(0.0 <= rt.min() and rt.max() <= 1.0,
          f"RT is the raw 0..1 model domain ({rt.min():.3f}..{rt.max():.3f})")

    ccs = np.array(d["CCS"], dtype=float)
    check(np.isfinite(ccs).all() and 100 < ccs.min(),
          f"CCS is physically plausible ({ccs.min():.0f}..{ccs.max():.0f} A^2)")
    # 1/K0, not CCS again: the two are proportional but not equal, and emitting
    # one under the other's name costs a diaPASEF consumer the whole dimension.
    im = np.array(d["IM"], dtype=float)
    check(np.isfinite(im).all() and 0.4 < im.min() and im.max() < 2.0,
          f"IM is 1/K0, not a copy of CCS ({im.min():.3f}..{im.max():.3f})")

    pmz = np.array(d["Precursor.Mz"], dtype=float)
    lo, hi = cfg["precursor_mz"]
    check(lo <= pmz.min() and pmz.max() <= hi,
          f"precursor m/z inside the configured window ({pmz.min():.1f}..{pmz.max():.1f})")

    inten = [np.asarray(x, dtype=float) for x in d["Relative.Intensity"]]
    check(all(np.isfinite(x).all() for x in inten), "no NaN fragment intensities")
    # NOT "every spectrum peaks at exactly 1.0", on either side. Intensities are
    # relative to the PREDICTED base peak, which is often outside the fragment
    # m/z window or below the top-N cut and so absent from the file (observed
    # per-spectrum maxima run from 0.46 to 1.0); and the stored float32 maximum
    # overshoots 1.0 by ~2e-5, which is numerical, not a fragment brighter than
    # the base peak. The bound is here to catch an unnormalised spectrum -- an
    # order of magnitude out -- not to police the last digit.
    peak = max(x.max() for x in inten)
    check(peak <= 1.001, f"fragment intensities are relative to the base peak (max {peak:.6f})")
    check(min(len(x) for x in inten) >= cfg["fragments"][0],
          f"every precursor keeps at least min_fragments ({min(len(x) for x in inten)})")

    irt = pq.read_table("irt.parquet")
    irt_rt = np.array(irt.to_pydict()["RT"], dtype=float)
    check(irt_rt.max() > 1.5, f"irt_rescale changes the RT domain ({irt_rt.min():.1f}..{irt_rt.max():.1f})")
    irt_meta = {k.decode(): v.decode() for k, v in (irt.schema.metadata or {}).items()}
    # Identical content in two different RT domains must not share a cache key.
    check(meta.get("odia.params") != irt_meta.get("odia.params"),
          "a raw and an iRT library do not share a cache fingerprint")

    def decoy_fragments(path):
        t = pq.read_table(path).to_pydict()
        return [np.asarray(p, dtype=float)
                for p, dec in zip(t["Product.Mz"], t["Decoy"]) if dec]

    rev, mut = decoy_fragments("rev.parquet"), decoy_fragments("mut.parquet")
    check(len(rev) == precursors and len(mut) == precursors,
          "each decoy method produced one decoy per precursor")
    # Compared on FRAGMENT m/z, not on the sequence string: a decoy row carries
    # its target's sequence, and with recompute_decoy_mz off it inherits the
    # target's precursor m/z too, so those are identical by design. The methods
    # were once distinguishable in neither -- all three mapped to Mutate.
    shared = sum(1 for a, b in zip(rev, mut) if len(a) == len(b) and np.allclose(a, b))
    check(shared == 0, f"reverse and mutate are different methods ({shared} shared spectra)")

    with open("lib.tsv") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        rows = sum(1 for _ in fh)
    check(rows > 0 and "Precursor.Mz" in header and "Product.Mz" in header,
          f"the TSV is the DIA-NN dialect with {rows} rows")

    print()
    if failures:
        print("e2e: FAILED")
        for f in failures:
            print("  -", f)
        return 1
    print("e2e: all checks passed")
    return 0


if __name__ == "__main__":
    if len(sys.argv) < 4:
        sys.exit(__doc__)
    work = sys.argv[4] if len(sys.argv) > 4 else tempfile.mkdtemp(prefix="dialibgen-e2e.")
    sys.exit(main(os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2]),
                  os.path.abspath(sys.argv[3]), work))
