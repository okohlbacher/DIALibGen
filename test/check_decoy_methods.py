#!/usr/bin/env python
"""Assert the decoy methods are actually different from each other.

This exists because they were not. DIALibraryGenerator hand-rolled its own
name->method chain instead of calling parseDecoyMethod, and every name it did
not recognise fell through to Mutate -- so "reverse" and "shuffle" silently
produced mutation decoys and three libraries came out byte-identical. Nothing
failed, nothing warned; the only symptom was 100% pairwise fragment identity.

Usage: check_decoy_methods.py <target.parquet> <m1>=<f1> <m2>=<f2> ...
"""
import sys

import pyarrow.parquet as pq

# A decoy that reproduces most of its target's fragments is not a null. The bar
# is deliberately loose -- some coincidence is expected, since a permuted
# peptide keeps its terminal ions when the permuted region does not reach them.
MAX_TARGET_IDENTITY = 0.25
# Two DIFFERENT methods agreeing this closely means one of them is not running.
MAX_PAIRWISE_IDENTITY = 0.50


def fragments(path, decoy):
    t = pq.read_table(path, columns=["Precursor.Id", "Product.Mz", "Decoy"]).to_pydict()
    return {p: tuple(round(x, 4) for x in mz)
            for p, mz, d in zip(t["Precursor.Id"], t["Product.Mz"], t["Decoy"])
            if bool(d) == decoy}


def identity(a, b, keys, strip):
    same = total = 0
    for k in keys:
        x, y = a[k], b[k[:-6] if strip else k]
        n = min(len(x), len(y))
        same += sum(1 for i in range(n) if x[i] == y[i])
        total += n
    return same / total if total else 0.0


def main(target_path, pairs):
    targets = fragments(target_path, decoy=False)
    got = {m: fragments(p, decoy=True) for m, p in pairs}
    first = next(iter(got.values()))
    keys = [k for k in first if k[:-6] in targets and all(k in g for g in got.values())]
    if len(keys) < 50:
        raise SystemExit(f"only {len(keys)} comparable decoys; the fixture is too small "
                         "for this test to mean anything")
    bad = []
    for m, g in got.items():
        f = identity(g, targets, keys, strip=True)
        print(f"  {m:16s} {100*f:6.2f}% of fragments equal the target's")
        if f > MAX_TARGET_IDENTITY:
            bad.append(f"{m} reproduces {100*f:.1f}% of its target's fragments")
    ms = list(got)
    for i in range(len(ms)):
        for j in range(i + 1, len(ms)):
            f = identity(got[ms[i]], got[ms[j]], keys, strip=False)
            print(f"  {ms[i]:16s} vs {ms[j]:16s} {100*f:6.2f}% identical")
            if f > MAX_PAIRWISE_IDENTITY:
                bad.append(f"{ms[i]} and {ms[j]} agree on {100*f:.1f}% of fragments "
                           "-- one of them is probably not dispatched")
    if bad:
        raise SystemExit("decoy methods are not distinct:\n  " + "\n  ".join(bad))
    print(f"  {len(got)} methods, all distinct, on {len(keys)} decoys")


if __name__ == "__main__":
    if len(sys.argv) < 4:
        raise SystemExit(__doc__)
    main(sys.argv[1], [a.split("=", 1) for a in sys.argv[2:]])
