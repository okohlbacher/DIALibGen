#!/usr/bin/env python
"""Guard for the PeptDeep reference encoder.

The reference is the oracle ODIA's C++ encoder is validated against, so if it
drifts, everything checked against it is confidently wrong.

The previous version of this file was near-vacuous: of 13 deliberate defects
injected into the encoder, 12 passed. Its "accumulation" assertion tested an
unmodified peptide, nine of eleven composition entries were never touched,
`aa_indices` was checked at one position, the '?' bucket had no coverage, and
`predict_rt` was not exercised at all. Every assertion here exists because a
specific mutation survived without it.

Expected values are written as literals -- absolute indices and pinned model
outputs -- so a permuted element list or a rewritten helper cannot satisfy them
by construction.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import peptdeep_reference as ref

failures = []


def check(cond, msg):
    if not cond:
        failures.append(msg)


# --- element list: absolute positions, not self-referential lookups ----------
check(len(ref.MOD_ELEMENTS) == 109, f"element list is {len(ref.MOD_ELEMENTS)}, must be 109")
for index, symbol in ((0, "C"), (1, "H"), (2, "N"), (3, "O"), (4, "P"), (5, "S"),
                      (38, "Dy"), (104, "2H"), (105, "13C"),
                      (106, "15N"), (107, "18O"), (108, "?")):
    check(ref.MOD_ELEMENTS[index] == symbol,
          f"element {index} must be {symbol!r}, got {ref.MOD_ELEMENTS[index]!r}")
check(len(set(ref.MOD_ELEMENTS)) == 109,
      "element symbols must be unique or index lookup is ambiguous")

# --- every composition entry, against literal (index, count) pairs -----------
# Indices 0-5 are C H N O P S, fixed by the yaml's own comment.
EXPECTED = {
    "Acetyl":          {0: 2, 1: 2, 3: 1},
    "Carbamidomethyl": {0: 2, 1: 3, 2: 1, 3: 1},
    "Oxidation":       {3: 1},
    "GG":              {0: 4, 1: 6, 2: 2, 3: 2},
    "Deamidated":      {1: -1, 2: -1, 3: 1},
    "Phospho":         {1: 1, 3: 3, 4: 1},
    "Amidated":        {1: 1, 2: 1, 3: -1},
    "Dehydrated":      {1: -2, 3: -1},
    "1":               {0: 2, 1: 2, 3: 1},          # UniMod aliases must agree
    "4":               {0: 2, 1: 3, 2: 1, 3: 1},
    "35":              {3: 1},
    "21":              {1: 1, 3: 3, 4: 1},
    "121":             {0: 4, 1: 6, 2: 2, 3: 2},
    "2":               {1: 1, 2: 1, 3: -1},
    "23":              {1: -2, 3: -1},
    "7":               {1: -1, 2: -1, 3: 1},
}
for name, expected in EXPECTED.items():
    if name not in ref.COMPOSITION:
        failures.append(f"{name} missing from COMPOSITION")
        continue
    v = ref.mod_vector(name)
    got = {int(i): int(v[i]) for i in np.nonzero(v)[0]}
    check(got == expected, f"{name} vector is {got}, expected {expected}")
check(set(ref.COMPOSITION) == set(EXPECTED),
      f"COMPOSITION and the expected table disagree on: "
      f"{sorted(set(ref.COMPOSITION) ^ set(EXPECTED))}")

# --- the '?' bucket: unknown elements ACCUMULATE, known ones are assigned -----
ref.COMPOSITION["_probe_one"] = {"Xx": 3}
ref.COMPOSITION["_probe_two"] = {"Xx": 3, "Yy": 4}
check(ref.mod_vector("_probe_one")[108] == 3,
      "an unknown element belongs in slot 108")
check(ref.mod_vector("_probe_two")[108] == 7,
      "two unknown elements must accumulate, not overwrite")
del ref.COMPOSITION["_probe_one"], ref.COMPOSITION["_probe_two"]

# --- aa_indices: every position, not just one --------------------------------
aa, mod_x = ref.encode("PEPC(Carbamidomethyl)TIDEK")
check(list(aa) == [0, 16, 5, 16, 3, 20, 9, 4, 5, 11, 0], f"aa_indices wrong: {list(aa)}")
check(aa.dtype == np.int64, f"aa_indices must be int64, got {aa.dtype}")
check(mod_x.dtype == np.float32, f"mod_x must be float32, got {mod_x.dtype}")
check(mod_x.shape == (11, 109), f"a 9-mer encodes to (11, 109), got {mod_x.shape}")

# --- the three site classes must land on three different rows ----------------
for seq, row in (("(UniMod:1)PEPTIDEK", 0), (".(Acetyl)PEPTIDEK", 0),
                 ("PEPTIDEK(GG)", 8), ("PEPTIDER.(Amidated)", 9)):
    _, m = ref.encode(seq)
    rows = [int(r) for r in np.nonzero(m.any(axis=1))[0]]
    check(rows == [row], f"{seq}: modification belongs on row {row}, got {rows}")

# --- accumulation at a site, on a MODIFIED peptide ---------------------------
_, m = ref.encode("PEPC(UniMod:4)(UniMod:4)TIDEK")
check(float(m[4].sum()) == 14.0,
      f"two modifications on one site must accumulate to 14 atoms, got {m[4].sum()}")

# --- names never contribute residues; "UniMod" contains an M -----------------
residues, _ = ref.parse("(UniMod:1)ADAWEEIRR")
check("".join(residues) == "ADAWEEIRR", f"got residues {''.join(residues)}")

# --- input the encoder must refuse rather than silently mis-encode -----------
for bad in ("", "peptidek"):
    try:
        ref.encode(bad)
        failures.append(f"encode({bad!r}) should have been rejected")
    except ValueError:
        pass

# --- the model itself: pinned outputs, dtypes and length grouping ------------
# Falls back to the OpenMS install's model directory. $ODIA_OPENMS, then
# $CONDA_PREFIX, then the model is simply absent and the block below is skipped.
_prefix = os.environ.get("ODIA_OPENMS") or os.environ.get("CONDA_PREFIX") or ""
MODEL = os.environ.get(
    "ODIA_RT_MODEL",
    os.path.join(_prefix, "share", "OpenMS", "models",
                 "peptdeep_rt_dynamic.onnx") if _prefix else "")
if os.path.exists(MODEL):
    # Obtained by running the shipped model directly, and reproduced
    # independently during review.
    got = ref.predict_rt(MODEL, ["ELVISLIVESK", "PEPTIDEK", "PEPTIDER"])
    for seq, want, actual in zip(("ELVISLIVESK", "PEPTIDEK", "PEPTIDER"),
                                 (0.819559, 0.124775, 0.175533), got):
        check(abs(actual - want) < 1e-4, f"RT({seq}) = {actual:.6f}, expected {want}")

    # Mixed lengths in one call must agree with one call each, which only holds
    # if each length is batched separately -- padding is not inert.
    alone = ref.predict_rt(MODEL, ["ELVISLIVESK"])[0]
    check(abs(got[0] - alone) < 1e-6,
          f"length grouping broken: {got[0]} in a mixed batch vs {alone} alone")

    # Duplicates must not collapse, and order must be preserved.
    dup = ref.predict_rt(MODEL, ["PEPTIDER", "ELVISLIVESK", "PEPTIDER"])
    check(len(dup) == 3 and dup[0] == dup[2] and abs(dup[1] - alone) < 1e-6,
          f"duplicate or ordering handling is wrong: {dup}")
else:
    print(f"  note: {MODEL} absent, model assertions skipped")

for f in failures:
    print(f"  FAIL {f}")
print(f"peptdeep encoding: {len(failures)} failures")
sys.exit(1 if failures else 0)
