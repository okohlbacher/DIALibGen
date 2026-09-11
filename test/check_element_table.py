#!/usr/bin/env python
"""The C++ header and the Python data file must not drift.

They are the same table for two languages. Nothing compiled or tested the
header before, so it could diverge from the list the oracle reads with nothing
failing -- and the index IS the feature position.
"""
import os
import re
import sys

root = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

data = [ln.strip() for ln in open(f"{root}/data/peptdeep_mod_elements.txt")
        if ln.strip() and not ln.startswith("#")]
header = open(f"{root}/include/odia/PeptDeepElements.h").read()
body = header[header.index("PEPTDEEP_MOD_ELEMENTS{"):header.index("};")]
from_header = re.findall(r'"([^"]+)"', body)

# The whole table pinned by content, not just its length and endpoints.
# Only 13 of 109 indices were pinned before, so swapping Fe with Zn -- both of
# which appear in real UniMod compositions -- passed every test on BOTH sides,
# since the C++ and the Python read the same file. A shared error is invisible
# to a comparison between them; only a literal pin catches it.
EXPECTED_SHA256 = "085570abe03e183cb34a472ce99cb9fe62de8191503dd3a53a4ae572dabe5475"

failures = []
if len(data) != 109:
    failures.append(f"data file has {len(data)} elements, must be 109")
if from_header != data:
    for i, (a, b) in enumerate(zip(from_header, data)):
        if a != b:
            failures.append(f"header and data diverge at index {i}: {a!r} vs {b!r}")
            break
    else:
        failures.append(f"header has {len(from_header)} elements, data has {len(data)}")
import hashlib
digest = hashlib.sha256("\n".join(data).encode()).hexdigest()
if digest != EXPECTED_SHA256:
    failures.append(f"element table content changed: sha256 {digest}, "
                    f"expected {EXPECTED_SHA256}. If this is a deliberate "
                    f"upstream update, regenerate and update the pin.")
if len(set(data)) != len(data):
    failures.append("duplicate symbols would make index lookup ambiguous")
# 'No' must survive as a string; YAML 1.1 reads a bare No as false.
if "No" not in data:
    failures.append("element 'No' is missing -- a YAML round-trip may have made it False")

for f in failures:
    print(f"  FAIL {f}")
print(f"element table: {len(failures)} failures ({len(data)} elements)")
sys.exit(1 if failures else 0)
