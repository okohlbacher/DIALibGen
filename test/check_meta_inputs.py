#!/usr/bin/env python
"""The MS2 meta-input constants must match the pinned table, in both languages.

A comparison between ODIA and the Python oracle cannot catch a constant that is
wrong in both -- swapping the charge and NCE scales in each at once passes every
case in compare_ms2.py, and so does relabelling timsTOF as Lumos. Only a literal
pin catches those, which is the same reason data/peptdeep_mod_elements.txt
exists.

The pin is not a third copy of the implementation: it carries the provenance of
each value in its header, and it is the file the header comment cites.
"""
import os
import re
import sys

root = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import peptdeep_reference as ref

pinned = {}
instruments = {}
with open(f"{root}/data/peptdeep_meta_inputs.txt") as f:
    for line in f:
        line = line.rstrip("\n")
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if parts[0] == "instrument":
            instruments[parts[1]] = int(parts[2])
        else:
            pinned[parts[0]] = float(parts[1])

header = open(f"{root}/include/odia/PeptDeepEncoder.h").read()
encoder = open(f"{root}/src/predict/PeptDeepEncoder.cpp").read()

failures = []


def check(what, got, want):
    if got != want:
        failures.append(f"{what} is {got!r}, pinned as {want!r}")


def cpp_constant(name):
    m = re.search(rf"{name}\s*=\s*([0-9.]+)f?;", header)
    return float(m.group(1)) if m else None


check("C++ CHARGE_SCALE", cpp_constant("CHARGE_SCALE"), pinned["charge_scale"])
check("C++ NCE_SCALE", cpp_constant("NCE_SCALE"), pinned["nce_scale"])
check("reference CHARGE_SCALE", ref.CHARGE_SCALE, pinned["charge_scale"])
check("reference NCE_SCALE", ref.NCE_SCALE, pinned["nce_scale"])

# The C++ map is a literal in the source; read it rather than trusting a
# separately maintained list.
cpp_instruments = dict(
    (name, int(index))
    for name, index in re.findall(r'\{"([A-Za-z]+)",\s*(\d+)\}', encoder))
for name, index in instruments.items():
    check(f"C++ instrument {name}", cpp_instruments.get(name.upper(),
                                                        cpp_instruments.get(name)), index)
    check(f"reference instrument {name}", ref.INSTRUMENTS.get(name.upper()), index)

check("C++ unknown instrument", float(re.search(r"return it == known.end\(\) \? (\d+)", encoder).group(1)),
      pinned["unknown_instrument"])
check("reference unknown instrument", float(ref.UNKNOWN_INSTRUMENT),
      pinned["unknown_instrument"])

# An instrument the pin does not name must not be silently mapped to a real
# one; that is the failure mode the "unknown" slot exists to prevent.
for stray in ("Astral", "Orbitrap", "", "qe "):
    if ref.INSTRUMENTS.get(stray.upper()) is not None:
        failures.append(f"reference maps unlisted instrument {stray!r}")

for f in failures:
    print(f"  FAIL {f}")
print(f"meta inputs: {len(failures)} failures "
      f"({len(instruments)} instruments, {len(pinned)} scalars)")
sys.exit(1 if failures else 0)
