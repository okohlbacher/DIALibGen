#!/usr/bin/env python
"""Compare ODIA's C++ PeptDeep encoder against the independent Python reference.

The two were written from doc/04-peptdeep-encoding.md separately -- the C++ from
OpenMS's AASequence, the Python from the raw string -- so agreement is evidence
that the spec was implemented rather than that one copied the other.

Usage: compare_encoders.py <odia_encode_dump> [sequence...]
"""
import json
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import peptdeep_reference as ref

DEFAULT = [
    "PEPTIDEK",                          # unmodified
    "PEPC(Carbamidomethyl)TIDEK",        # residue modification
    ".(Acetyl)PEPTIDEK",                 # N-terminal, row 0
    "PEPTIDER.(Amidated)",               # C-terminal, row nAA+1, negative count
    "PEPTIDEM(Oxidation)K",              # a different residue modification
    ".(Acetyl)PEPC(Carbamidomethyl)TIDEM(Oxidation)K",   # several at once
    "AVVPASLSGQDVGSFAYLTIK",             # long, unmodified
    "C(Carbamidomethyl)C(Carbamidomethyl)PEPTIDEK",      # adjacent modifications
]


def compare_one(cpp, seq):
    """Compare one peptide's encoding. Returns a failure string, or None."""
    aa, mod_x = ref.encode(seq)
    if cpp["aa_indices"] != aa.tolist():
        return (f"{seq}: aa_indices differ\n"
                f"        C++    {cpp['aa_indices']}\n        python {aa.tolist()}")
    py_rows = {str(r): {str(c): float(mod_x[r][c]) for c in np.nonzero(mod_x[r])[0]}
               for r in range(mod_x.shape[0]) if mod_x[r].any()}
    cpp_rows = {r: {c: float(v) for c, v in cols.items()}
                for r, cols in cpp["mod_x"].items()}
    if cpp_rows != py_rows:
        return (f"{seq}: mod_x differs\n"
                f"        C++    {cpp_rows}\n        python {py_rows}")
    return None


def compare_batch(dump, sequences):
    """Encode several equal-length peptides as ONE batch and check every row.

    A one-peptide batch cannot detect a row-indexing defect: deleting the row
    term from the encoder's mod_x index passed every test, while collapsing
    every peptide's modifications onto row 0 and shifting a predicted iRT by
    0.096 -- twice the error that motivated getting terminal placement right.
    """
    proc = subprocess.run([dump, *sequences], capture_output=True, text=True)
    if proc.returncode != 0:
        return [f"batch of {len(sequences)}: C++ failed: {proc.stderr.strip()}"]
    batch = json.loads(proc.stdout)
    if batch["rows"] != len(sequences):
        return [f"batch reported {batch['rows']} rows, expected {len(sequences)}"]
    out = []
    for row, (cpp, seq) in enumerate(zip(batch["peptides"], sequences)):
        bad = compare_one(cpp, seq)
        if bad:
            out.append(f"row {row}: {bad}")
    return out


def main(dump, sequences):
    failures = 0
    for seq in sequences:
        proc = subprocess.run([dump, seq], capture_output=True, text=True)
        if proc.returncode != 0:
            print(f"  FAIL {seq}: C++ encoder failed: {proc.stderr.strip()}")
            failures += 1
            continue
        cpp = json.loads(proc.stdout)["peptides"][0]

        aa, mod_x = ref.encode(seq)
        if cpp["aa_indices"] != aa.tolist():
            print(f"  FAIL {seq}: aa_indices differ\n"
                  f"        C++    {cpp['aa_indices']}\n        python {aa.tolist()}")
            failures += 1
            continue

        py_rows = {str(r): {str(c): float(mod_x[r][c])
                            for c in np.nonzero(mod_x[r])[0]}
                   for r in range(mod_x.shape[0]) if mod_x[r].any()}
        cpp_rows = {r: {c: float(v) for c, v in cols.items()}
                    for r, cols in cpp["mod_x"].items()}
        if cpp_rows != py_rows:
            print(f"  FAIL {seq}: mod_x differs\n"
                  f"        C++    {cpp_rows}\n        python {py_rows}")
            failures += 1
            continue
        print(f"  ok   {seq}")

    # Multi-peptide batches, so every row is checked and not only row 0.
    # Modifications are placed on different rows deliberately: a defect that
    # collapses rows would otherwise be invisible when all rows look alike.
    BATCHES = [
        ["PEPTIDEKKKR", "PEPC(Carbamidomethyl)TIDEKKR", "PEPTIDEM(Oxidation)KKR"],
        [".(Acetyl)PEPTIDEKKKR", "PEPTIDEKKK R".replace(" ", ""), "PEPTIDEKKKR.(Amidated)"],
        ["PEPTIDEKKKR"],
    ]
    for batch in BATCHES:
        bad = compare_batch(dump, batch)
        for message in bad:
            print(f"  FAIL batch: {message}")
        failures += len(bad)
        if not bad:
            print(f"  ok   batch of {len(batch)}")

    print(f"encoder comparison: {failures} failures")
    return 1 if failures else 0


if __name__ == "__main__":
    seqs = sys.argv[2:] or DEFAULT
    sys.exit(main(sys.argv[1], seqs))
