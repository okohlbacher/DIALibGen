# Backlog

Work that is worth doing and is not scheduled. Anything here should be
actionable: what it is, and what it actually touches.

## The macOS CLI takes 5½ minutes on its first run

Measured on an M-series Mac against the published v0.9.0 tarball: **335 s on the
first invocation, 1 s on every one after.** The unsigned v0.2.3 tarball started in
0.4 s. It is not a hang — a `sample` of the process shows dyld parked in
`Loader::mapSegments`, validating a code signature, while syspolicyd logs one XPC
round trip roughly every second. 145 Developer-ID-signed dylibs, evaluated one at
a time, verdict cached afterwards.

So signing made the first run two orders of magnitude slower, and most people will
conclude it has hung and kill it. This is now the worst thing about the macOS
build, and it reaches Homebrew too: the CLI cask installs that same tarball.

**A `.tar.gz` cannot carry a stapled ticket**, so each dylib is checked with Apple
individually. The fix is to ship the CLI in a container that can: a signed,
notarized, **stapled** `.dmg` or `.pkg`, whose ticket covers the whole payload and
is consulted once, offline. FASTag already does exactly this — see its
"Package the CLI as a signed, notarized disk image" step, which is the same problem
solved.

Options, cheapest first:
- Ship `DIALibGen-macos-<arch>.dmg` alongside the tarball and point the
  CLI cask at it. The tarball stays for people who want a plain archive.
- Or a `.pkg`, which installs to a fixed prefix and is what `brew install --cask`
  handles most cleanly.
- Measure before and after. The claim to verify is first-run wall time from a
  fresh download, not that the ticket exists.

Until then, say so in the README: the first run is slow, it is not stuck.

## Establish the licence on the AlphaPeptDeep weights

The project does not ship the three ONNX models and says the terms are "not
established". Researched 2026-09-12; the position is narrower than "unclear":

**The weights are unlicensed, not restricted.** AlphaPeptDeep's code is Apache-2.0
(verified byte-for-byte), but that licence has never been applied to the weights by
anyone at MannLabs. `MANIFEST.in` and the Dockerfile both show the Apache metadata
attaches to a distribution that provably excludes the models — `peptdeep` downloads
them at first run from a GitHub release whose body never uses the word "licence",
and `pretrained_models_v3.zip` contains no LICENSE, NOTICE or COPYING (both ZIP
central directories enumerated). Apache-2.0 §1 defines the licensed Work as
material carrying a notice "included in or attached to" it; these carry none.

Equally, **nothing restricts them**: no non-commercial clause, no academic-use
clause, no citation-as-condition, anywhere in the repo, licence, release, docs or
paper. Training-data terms do not propagate — the ProteomeXchange sets are CC0
except PXD004732, whose EBI terms explicitly add no restrictions of their own.

The ONNX files are a **second** derivation and not MannLabs artefacts at all:
OpenMS produced them with `tools/scripts/export_peptdeep_models_to_onnx.py`, which
carries no copyright line and no SPDX header, unlike every OpenMS C++ file. The PR
that introduced them discusses provenance nowhere. So there is nothing to inherit
from OpenMS either.

**Searching for a file is exhausted — the remaining move is to ask.**
- Wen-Feng Zeng (@jalew188), first author, created the `pre-trained-models`
  release; or `opensource@alphapept.com`, the contact in MannLabs' own CLA. One
  GitHub issue asking "under what licence are the release weights published?" is
  the whole fix.
- Timo Sachsenberg (OpenMS, **Tübingen — same institution**) maintains the PeptDeep
  files and reviewed the PR. He can say whether OpenMS obtained permission before
  uploading to archive.openms.de.

Until one of them answers, `scripts/fetch-models.sh` downloads them on demand,
which is the same thing a source build of OpenMS does.
