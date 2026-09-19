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
aliases = {}
with open(f"{root}/data/peptdeep_meta_inputs.txt") as f:
    for line in f:
        line = line.rstrip("\n")
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if parts[0] == "instrument":
            instruments[parts[1]] = int(parts[2])
        elif parts[0] == "alias":
            aliases[parts[1]] = parts[2]
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

# The INDEX map must know only the five names the model was trained with. An
# alias reaches a trained slot through canonicalInstrument, a layer above this
# one, and does so loudly -- it is logged and recorded in the provenance. What
# must never happen is the index map growing a name of its own.
for stray in ("Astral", "Orbitrap", "", "qe "):
    if ref.INSTRUMENTS.get(stray.upper()) is not None:
        failures.append(f"reference maps unlisted instrument {stray!r}")

# The alias table, pinned against upstream's instrument_group. An alias that
# stops resolving falls through to the UNTRAINED slot, and no library built
# from it looks wrong afterwards.
cpp_aliases = dict(re.findall(r'\{"([A-Za-z0-9+]+)",\s*"([A-Za-z]+)"\}', encoder))
def fold(name):
    return "".join(c for c in name.upper() if c not in "-_ ")
for name, canonical in aliases.items():
    check(f"C++ alias {name}", cpp_aliases.get(fold(name)), canonical)
    if canonical not in instruments:
        failures.append(f"alias {name!r} resolves to {canonical!r}, which is not a pinned instrument")
# Every canonical name must map to itself, or a config naming it exactly would
# be refused by the very check meant to protect it.
for name in instruments:
    check(f"C++ alias {name} (identity)", cpp_aliases.get(fold(name)), name)

for f in failures:
    print(f"  FAIL {f}")
print(f"meta inputs: {len(failures)} failures "
      f"({len(instruments)} instruments, {len(aliases)} aliases, {len(pinned)} scalars)")
sys.exit(1 if failures else 0)
