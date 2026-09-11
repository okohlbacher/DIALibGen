#!/usr/bin/env python
"""Compare ODIA's predicted collision cross-sections against the reference.

CCS is the third model kind, and the one the constructor used to reject: it
takes three inputs, matching neither the RT model's two nor the MS2 model's
five, so a guard written for the first two refused it outright.

Batches are mixed-length and mixed-charge on purpose. Charge is the only thing
distinguishing two rows of the same peptide, so a version that used the first
row's charge for the whole batch would agree with a single-peptide fixture.

Usage: compare_ccs.py <odia_predict_ccs> <model.onnx>
"""
import json
import math
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import peptdeep_reference as ref

# The same peptide at several charges, because CCS depends strongly on charge
# and nothing else in this list would notice if it were ignored.
BATCHES = [
    [("ELVISLIVESK", 2)],
    [("ELVISLIVESK", 2), ("ELVISLIVESK", 3), ("ELVISLIVESK", 4)],
    [("PEPTIDEK", 2), ("AVVPASLSGQDVGSFAYLTIK", 3), ("AAAK", 1),
     ("PEPC(Carbamidomethyl)TIDEK", 2), ("PEPTIDEM(Oxidation)K", 3),
     (".(Acetyl)PEPTIDEK", 2)],
]
TOLERANCE = 1e-4

# Literal expected values, in square angstroms.
#
# Agreement with the oracle cannot catch a constant that is wrong in both, and
# the charge scale is exactly that kind of constant: applying 1.0 instead of 0.1
# in the C++ and in the reference at once passes every comparison here. Only a
# literal pin catches it, which is why check_invariants.py pins residue masses.
#
# These are also a sanity floor: a singly-charged tryptic peptide is a few
# hundred square angstroms, so a value in the tens or the thousands is wrong
# whatever the comparison says.
PINNED = [
    ("ELVISLIVESK", 2, 387.3969),
    ("ELVISLIVESK", 3, 477.7039),
    ("PEPTIDEK", 2, 325.3780),
    ("AVVPASLSGQDVGSFAYLTIK", 3, 594.8837),
    ("PEPC(Carbamidomethyl)TIDEK", 2, 356.4488),
]
# Loose enough for the last-bit spread between builds and thread counts, far
# tighter than any scale error.
PINNED_TOLERANCE = 0.05


def main(tool, model):
    failures = 0
    os.environ["ODIA_ORT_THREADS"] = "1"
    seen_spread = 0.0

    for pairs in BATCHES:
        proc = subprocess.run(
            [tool, model] + [f"{s}:{z}" for s, z in pairs],
            capture_output=True, text=True,
            env=dict(os.environ, ODIA_ORT_THREADS="1"))
        if proc.returncode != 0:
            print(f"  FAIL batch of {len(pairs)}: {proc.stderr.strip()}")
            failures += 1
            continue
        got = ref.load_json(proc, "ccs")
        want = ref.predict_ccs(model, [s for s, _ in pairs], [z for _, z in pairs])
        if len(got) != len(want):
            print(f"  FAIL batch of {len(pairs)}: {len(got)} values, expected {len(want)}")
            failures += 1
            continue
        worst = 0.0
        bad = False
        for k, (g, w) in enumerate(zip(got, want)):
            if not math.isfinite(g) or not math.isfinite(w):
                print(f"  FAIL {pairs[k]}: non-finite CCS {g} vs {w}")
                failures += 1
                bad = True
                continue
            worst = max(worst, abs(g - w))
        if bad:
            continue
        if worst > TOLERANCE:
            print(f"  FAIL batch of {len(pairs)}: worst {worst:.2e}")
            failures += 1
        else:
            print(f"  ok   batch of {len(pairs)}: worst {worst:.2e}, "
                  f"CCS {min(got):.1f}-{max(got):.1f}")
        seen_spread = max(seen_spread, max(got) - min(got))

    # A predictor that ignored charge, or returned a constant, would agree with
    # itself everywhere. The values must actually vary.
    if seen_spread < 10.0:
        print(f"  FAIL every predicted CCS is within {seen_spread:.2f} A^2 of every "
              f"other; the comparison cannot distinguish a constant")
        failures += 1

    # The pinned values, which are the only check that survives an error made
    # identically in both implementations.
    pinned = subprocess.run(
        [tool, model] + [f"{s}:{z}" for s, z, _ in PINNED],
        capture_output=True, text=True, env=dict(os.environ, ODIA_ORT_THREADS="1"))
    if pinned.returncode != 0:
        print(f"  FAIL pinned values: {pinned.stderr.strip()}")
        failures += 1
    else:
        bad_pins = 0
        for (seq, charge, want), got in zip(PINNED, ref.load_json(pinned, "pinned ccs")):
            if not math.isfinite(got) or abs(got - want) > PINNED_TOLERANCE:
                print(f"  FAIL {seq} z={charge}: CCS {got} A^2, pinned at {want}")
                bad_pins += 1
        failures += bad_pins
        if not bad_pins:
            print(f"  ok   {len(PINNED)} pinned CCS values")

    # Feeding the wrong model kind must say so rather than fail inside Run.
    for other, kind in ((model.replace("_ccs_", "_rt_"), "retention-time"),
                        (model.replace("_ccs_", "_ms2_"), "MS2")):
        if not os.path.exists(other):
            continue
        proc = subprocess.run([tool, other, "PEPTIDEK:2"],
                              capture_output=True, text=True)
        if proc.returncode == 0:
            print(f"  FAIL the {kind} model was accepted by predictCCS")
            failures += 1
        elif "CCS model takes 3" not in proc.stderr:
            print(f"  FAIL the {kind} model was refused without saying why: "
                  f"{proc.stderr.strip()}")
            failures += 1

    # Charge is what a single-peptide fixture cannot test.
    charge_only = subprocess.run(
        [tool, model, "ELVISLIVESK:2", "ELVISLIVESK:3"],
        capture_output=True, text=True, env=dict(os.environ, ODIA_ORT_THREADS="1"))
    if charge_only.returncode == 0:
        a, b = ref.load_json(charge_only, "charge-only ccs")
        if abs(a - b) < 1.0:
            print(f"  FAIL the same peptide at charge 2 and 3 gave {a:.3f} and "
                  f"{b:.3f}; charge is not reaching the model")
            failures += 1

    print(f"CCS comparison: {failures} failures over {len(BATCHES)} batches")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1], sys.argv[2]))
