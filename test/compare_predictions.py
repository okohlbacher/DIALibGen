#!/usr/bin/env python
"""Compare ODIA's predicted iRT against the independent Python reference.

Tensor agreement is not prediction agreement: dtype, layout, stride or ordering
assumptions can differ while the encoded values match. This compares what the
models actually return.

The sequence list deliberately mixes lengths, so a broken length-grouping or a
lost reordering shows up as wrong values rather than as an error.

Usage: compare_predictions.py <odia_predict_rt> <model.onnx>
"""
import math
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import peptdeep_reference as ref

# Modified and terminal cases, plus every length a tryptic digest can produce
# (min_peptide_length 7 to max_peptide_length 30). Four lengths were covered
# before, and a defect returning the right value for those four and 0.0 for
# every other length passed the whole suite.
_FILLER = "AGVLIFMPWSTYHKRQENDCAGVLIFMPWST"
SEQUENCES = [
    "ELVISLIVESK",
    "PEPTIDEK",
    "PEPTIDER",
    "PEPC(Carbamidomethyl)TIDEK",
    ".(Acetyl)PEPTIDEK",
    "PEPTIDER.(Amidated)",
    "AVVPASLSGQDVGSFAYLTIK",
    "PEPTIDEM(Oxidation)K",
] + [_FILLER[:n] for n in range(7, 31)]

# Set from measurement, not taste. Within one build, repeated calls are
# bit-identical and batch composition moves the result by 1 ULP (~1.5e-8). Across
# the two builds compared here -- the C++ links OpenMS's libonnxruntime.so, the
# Python its own copy -- kernels differ in the last bits and the spread is up to
# 2.7e-7 on values near 0.88, about 3 ULP of float32.
#
# 1e-6 sits above that and well below the systematic bias this test exists to
# catch: at the previous 1e-5 a constant +9e-6 offset on every prediction passed
# as "identical to six decimal places".
TOLERANCE = 1e-6


def main(tool, model):
    env = dict(os.environ, ODIA_ORT_THREADS="1")
    proc = subprocess.run([tool, model, *SEQUENCES], capture_output=True, text=True,
                          env=env)
    if proc.returncode != 0:
        print(f"  FAIL ODIA predictor failed: {proc.stderr.strip()}")
        return 1
    if "provider:" not in proc.stderr:
        print("  FAIL predictor did not report which execution provider it used")
        return 1

    cpp = {}
    for line in proc.stdout.strip().split("\n"):
        seq, value = line.rsplit("\t", 1)
        cpp[seq] = float(value)

    os.environ["ODIA_ORT_THREADS"] = "1"
    python = dict(zip(SEQUENCES, ref.predict_rt(model, SEQUENCES)))

    failures = 0
    for seq in SEQUENCES:
        if seq not in cpp:
            print(f"  FAIL {seq}: no prediction from ODIA")
            failures += 1
            continue
        # An explicit finiteness check, because abs(nan - x) > tol is False:
        # multiplying every prediction by NaN previously passed this test with
        # all eight lines reported "ok".
        if not math.isfinite(cpp[seq]):
            print(f"  FAIL {seq}: ODIA returned {cpp[seq]}")
            failures += 1
            continue
        if not math.isfinite(python[seq]):
            print(f"  FAIL {seq}: reference returned {python[seq]}")
            failures += 1
            continue
        d = abs(cpp[seq] - python[seq])
        if d > TOLERANCE:
            print(f"  FAIL {seq}: ODIA {cpp[seq]:.6f} vs reference {python[seq]:.6f} "
                  f"(delta {d:.2e})")
            failures += 1
        elif len(SEQUENCES) < 12:
            print(f"  ok   {seq}: {cpp[seq]:.6f}")

    print(f"prediction comparison: {failures} failures over {len(SEQUENCES)} sequences")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1], sys.argv[2]))
