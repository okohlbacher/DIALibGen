#!/usr/bin/env python
"""Assert the free-cysteine RT correction (doc/30) is exactly what it claims.

Three libraries from one FASTA:
  on   -- CAM-free, correction on
  off  -- CAM-free, correction off
  cam  -- carbamidomethylated, correction on

on - off must be the bucketed offset for every peptide and zero for peptides
without cysteine; cam must be untouched, because the correction counts
UNMODIFIED cysteines and an alkylated one is not one. That self-gating is the
whole reason the flag can default to on, so it is the part worth guarding.

Usage: check_free_cysteine.py <on.parquet> <off.parquet> <cam_on.parquet> <cam_off.parquet>
"""
import sys

import pyarrow.parquet as pq

EXPECTED = {0: 0.0, 1: 0.0343, 2: 0.0645}   # 3+ shares the 3 bucket
THIRD_AND_ABOVE = 0.0580
TOL = 1e-6


def rt_by_sequence(path):
    d = pq.read_table(path, columns=["Modified.Sequence", "RT"]).to_pydict()
    out = {}
    for seq, rt in zip(d["Modified.Sequence"], d["RT"]):
        out.setdefault(seq, rt)
    return out


def expected(n):
    return EXPECTED[n] if n in EXPECTED else THIRD_AND_ABOVE


def main(on, off, cam_on, cam_off):
    a, b = rt_by_sequence(on), rt_by_sequence(off)
    shared = [s for s in a if s in b]
    if not shared:
        raise SystemExit("no shared peptides: the two libraries are unrelated")

    bad = 0
    seen = set()
    for seq in shared:
        n = seq.count("C")          # CAM-free, so every C is a free one
        want = expected(n)
        got = a[seq] - b[seq]
        seen.add(min(n, 3))
        if abs(got - want) > TOL:
            bad += 1
            if bad <= 5:
                print(f"  FAIL {seq}: {n} Cys, expected {want:+.4f}, got {got:+.6f}")
    if bad:
        raise SystemExit(f"{bad}/{len(shared)} peptides carry the wrong offset")
    if not {0, 1} <= seen:
        raise SystemExit(f"fixture exercises only {sorted(seen)} cysteine counts; "
                         "it must contain both cysteine-free and cysteine-bearing peptides")
    print(f"  offsets exact on {len(shared)} peptides, buckets {sorted(seen)}")

    c, d = rt_by_sequence(cam_on), rt_by_sequence(cam_off)
    shared = [s for s in c if s in d]
    with_c = sum(1 for s in shared if "C" in s)
    worst = max((abs(c[s] - d[s]) for s in shared), default=0.0)
    if worst > TOL:
        raise SystemExit(f"alkylated library moved by {worst:.8f}: the correction is "
                         "counting modified cysteines, so it is not self-gating")
    if not with_c:
        raise SystemExit("alkylated fixture has no cysteine peptides, so it proves nothing")
    print(f"  alkylated library unchanged across {with_c} cysteine peptides")


if __name__ == "__main__":
    if len(sys.argv) != 5:
        raise SystemExit(__doc__)
    main(*sys.argv[1:])
