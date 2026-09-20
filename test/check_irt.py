#!/usr/bin/env python
"""ODIA's iRT calibration must reproduce OpenMS's own fixture.

The fixture (src/tests/class_tests/openms/data/peptdeep_irt_peptides_predicted.csv)
carries, for the 11 Biognosys standards, both the ONNX model's raw prediction
and the calibrated iRT that AlphaPeptDeep's add_irt_column_to_precursor_df
produces from it. That makes it an independent statement about two separate
things -- our inference and our line -- from an authority neither was written
from.

The line itself is deliberately NOT pinned as a constant anywhere in ODIA: it
depends on the model checkpoint and is refitted at run time. This test pins the
line for THIS checkpoint, so a silent model swap is visible.

usage: check_irt.py <odia_irt_calibration> <rt_model.onnx> <standards.tsv> <fixture.csv>
"""
import math
import subprocess
import sys

TOOL, MODEL, STANDARDS, FIXTURE = sys.argv[1:5]
failures = []


def fail(msg):
    failures.append(msg)
    print(f"  FAIL {msg}")


proc = subprocess.run([TOOL, MODEL, STANDARDS], capture_output=True, text=True)
if proc.returncode != 0:
    print(f"  FAIL tool failed: {proc.stderr.strip()}")
    sys.exit(1)

got = {}
per_peptide = {}
for line in proc.stdout.splitlines():
    parts = line.split("\t")
    if parts[0] == "P":
        per_peptide[parts[1]] = (float(parts[2]), float(parts[3]))
    else:
        got[parts[0]] = float(parts[1])

want = {}
with open(FIXTURE) as f:
    header = f.readline().rstrip("\n").split(",")
    at = {n: i for i, n in enumerate(header)}
    for line in f:
        row = line.rstrip("\n").split(",")
        want[row[at["sequence"]]] = {
            "irt": float(row[at["irt"]]),
            "rt_pred_onnx": float(row[at["rt_pred_onnx"]]),
            "irt_pred_onnx": float(row[at["irt_pred_onnx"]]),
        }

if len(per_peptide) != len(want):
    fail(f"{len(per_peptide)} standards predicted, fixture has {len(want)}")

# The model's raw output. This is inference, and it must agree to float32.
for seq, (raw, _) in per_peptide.items():
    if seq not in want:
        fail(f"{seq} is not in the fixture")
        continue
    if not math.isfinite(raw) or abs(raw - want[seq]["rt_pred_onnx"]) > 2e-6:
        fail(f"{seq}: raw prediction {raw:.9f}, fixture says "
             f"{want[seq]['rt_pred_onnx']:.9f}")

# The calibrated value. This is the line, and it is the thing under test.
# 5e-4 iRT on a 125-unit scale is 4 ppm of the range -- far tighter than the
# 8.76 iRT the worst standard genuinely misses, so this cannot pass by being
# loose.
for seq, (_, calibrated) in per_peptide.items():
    if seq not in want:
        continue
    if not math.isfinite(calibrated) or abs(calibrated - want[seq]["irt_pred_onnx"]) > 5e-4:
        fail(f"{seq}: calibrated {calibrated:.6f}, fixture says "
             f"{want[seq]['irt_pred_onnx']:.6f}")

# The fit for this checkpoint. Not a constant ODIA uses -- it refits every run
# -- but pinned here so a model swap is loud rather than silent.
for name, expected, tol in (("slope", 152.235611, 1e-3),
                            ("intercept", -39.232154, 1e-3),
                            ("peptides", 11, 0),
                            ("max_abs_error", 8.756966, 1e-3)):
    if name not in got:
        fail(f"tool did not report {name}")
    elif not math.isfinite(got[name]) or abs(got[name] - expected) > tol:
        fail(f"{name} is {got[name]}, expected {expected}")

# The calibration must be monotone increasing, or it would reorder the library.
if got.get("slope", 0) <= 0:
    fail(f"slope {got.get('slope')} is not positive; the mapping would invert "
         f"the retention-time order")

print(f"iRT calibration: {len(failures)} failures over {len(per_peptide)} standards")
sys.exit(1 if failures else 0)
