#!/usr/bin/env python
"""Convert an ODIA library's CCS column into the 1/K0 a timsTOF reports.

This lives outside ODIA on purpose. CCS is a property of the ion; 1/K0 is what
a particular instrument measures for it, and the two are related through the
drift gas and that instrument's calibration. ODIA therefore emits square
angstroms and stops, and the conversion happens here -- on the consumer side,
which knows the instrument.

The relation is Mason-Schamp in its reduced form,

    1/K0 = CCS * sqrt(mu) / (z * C),   mu = m_ion * m_gas / (m_ion + m_gas)

with m_ion = mz * z and m_gas the drift gas mass. C folds together the physical
constants and the gas conditions.

**C is fitted, not remembered.** The constant and the drift gas mass are values
this project has no authoritative local copy of, and a plausible-looking number
in the wrong units is worse than no number at all. Given a reference library
that carries predicted 1/K0 for the same precursors -- DIA-NN's does -- C is
determined by the data, and the scatter around the fit says whether the
functional form is right. A tight fit is evidence; a loose one means this
should not be used.

Usage:
  ccs_to_mobility.py fit    <odia.tsv> <reference-with-IM>
  ccs_to_mobility.py apply  <odia.tsv> <out.tsv> --coefficient C [--gas 28.0]
"""
import argparse
import csv
import math
import statistics
import sys

DRIFT_GAS_MASS = 28.0    # N2, the gas every timsTOF here runs.


def reduced_mass(mz, charge, gas):
    m_ion = mz * charge
    return m_ion * gas / (m_ion + gas)


def mobility(ccs, mz, charge, coefficient, gas=DRIFT_GAS_MASS):
    return ccs * math.sqrt(reduced_mass(mz, charge, gas)) / (charge * coefficient)


def coefficient_from(ccs, mz, charge, im, gas=DRIFT_GAS_MASS):
    """The C that would make this one precursor's CCS give this 1/K0."""
    return ccs * math.sqrt(reduced_mass(mz, charge, gas)) / (charge * im)


def read_precursors(path, fields):
    """{(sequence, charge): {field: value}} over the target precursors."""
    import os
    out = {}
    if path.endswith(".parquet"):
        import pyarrow.parquet as pq
        pf = pq.ParquetFile(path)
        names = [f.name for f in pf.schema_arrow]
        want = {k: names.index(v) for k, v in fields.items() if v in names}
        decoy = names.index("Decoy") if "Decoy" in names else None
        seq_i, chg_i = want["sequence"], want["charge"]
        for batch in pf.iter_batches(batch_size=200000):
            cols = {i: batch.column(i).to_pylist() for i in set(want.values()) |
                    ({decoy} if decoy is not None else set())}
            for r in range(batch.num_rows):
                if decoy is not None and cols[decoy][r]:
                    continue
                key = (cols[seq_i][r], int(cols[chg_i][r]))
                if key not in out:
                    out[key] = {k: cols[i][r] for k, i in want.items()}
        return out

    with open(path, newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
        want = {k: header.index(v) for k, v in fields.items() if v in header}
        decoy = header.index("Decoy") if "Decoy" in header else None
        for row in reader:
            if decoy is not None and row[decoy] not in ("0", ""):
                continue
            key = (row[want["sequence"]], int(float(row[want["charge"]])))
            if key not in out:
                out[key] = {k: row[i] for k, i in want.items()}
    return out


def strip_mods(s):
    keep, depth = [], 0
    for c in s:
        if c in "([{":
            depth += 1
        elif c in ")]}":
            depth = max(0, depth - 1)
        elif depth == 0 and c.isalpha():
            keep.append(c)
    return "".join(keep)


def fit(odia_path, reference_path, gas):
    ours = read_precursors(odia_path, {"sequence": "Modified.Sequence",
                                       "charge": "Precursor.Charge",
                                       "mz": "Precursor.Mz", "ccs": "CCS"})
    theirs = read_precursors(reference_path, {"sequence": "Modified.Sequence",
                                              "charge": "Precursor.Charge",
                                              "mz": "Precursor.Mz", "im": "IM"})
    ours = {(strip_mods(k[0]), k[1]): v for k, v in ours.items()}
    theirs = {(strip_mods(k[0]), k[1]): v for k, v in theirs.items()}

    per_charge = {}
    coefficients = []
    for key in set(ours) & set(theirs):
        a, b = ours[key], theirs[key]
        try:
            ccs = float(a["ccs"])
            mz = float(a["mz"])
            im = float(b["im"])
        except (ValueError, KeyError, TypeError):
            continue
        if not (ccs > 0 and im > 0):
            continue
        c = coefficient_from(ccs, mz, key[1], im, gas)
        coefficients.append(c)
        per_charge.setdefault(key[1], []).append(c)

    if len(coefficients) < 100:
        print(f"only {len(coefficients)} usable pairs; not enough to fit")
        return 1

    median = statistics.median(coefficients)
    spread = statistics.median([abs(c - median) / median for c in coefficients])
    print(f"pairs                {len(coefficients):,}")
    print(f"coefficient (median) {median:.4f}")
    print(f"relative spread      {100 * spread:.2f}% (median |c-C|/C)")
    print("per charge:")
    for z in sorted(per_charge):
        v = per_charge[z]
        m = statistics.median(v)
        s = statistics.median([abs(x - m) / m for x in v])
        print(f"  z={z}  n={len(v):>8,}  C={m:.4f}  spread {100 * s:.2f}%")

    # The whole point of the fit is that its tightness is the evidence. A
    # coefficient that varies with charge means the functional form is wrong,
    # not that the constant needs a per-charge table.
    charges = [statistics.median(v) for z, v in per_charge.items() if len(v) > 100]
    if charges:
        drift = (max(charges) - min(charges)) / statistics.median(charges)
        print(f"\ncoefficient drift across charge states: {100 * drift:.2f}%")
        if drift > 0.05:
            print("  WARNING: C should not depend on charge. The Mason-Schamp form "
                  "above is not describing this reference, so do not use it.")
    print(f"\napply with: --coefficient {median:.4f}")
    return 0


def apply(odia_path, out_path, coefficient, gas):
    """Rewrite the library with IM filled in from CCS."""
    converted = skipped = 0
    with open(odia_path, newline="") as fin, open(out_path, "w", newline="") as fout:
        reader = csv.reader(fin, delimiter="\t")
        writer = csv.writer(fout, delimiter="\t", lineterminator="\n")
        header = next(reader)
        for name in ("IM", "CCS", "Precursor.Mz", "Precursor.Charge"):
            if name not in header:
                raise SystemExit(f"library has no {name} column")
        i_im, i_ccs = header.index("IM"), header.index("CCS")
        i_mz, i_z = header.index("Precursor.Mz"), header.index("Precursor.Charge")
        writer.writerow(header)
        for row in reader:
            try:
                ccs = float(row[i_ccs])
                mz = float(row[i_mz])
                z = int(float(row[i_z]))
                if ccs > 0 and z > 0:
                    row[i_im] = f"{mobility(ccs, mz, z, coefficient, gas):.6f}"
                    converted += 1
                else:
                    skipped += 1
            except (ValueError, IndexError):
                skipped += 1
            writer.writerow(row)
    print(f"converted {converted:,} rows, left {skipped:,} without ion mobility")
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["fit", "apply"])
    ap.add_argument("library")
    ap.add_argument("second")
    ap.add_argument("--coefficient", type=float, default=None)
    ap.add_argument("--gas", type=float, default=DRIFT_GAS_MASS)
    args = ap.parse_args()

    if args.mode == "fit":
        return fit(args.library, args.second, args.gas)
    if args.coefficient is None:
        raise SystemExit("apply needs --coefficient; get one from `fit`")
    return apply(args.library, args.second, args.coefficient, args.gas)


if __name__ == "__main__":
    sys.exit(main())
