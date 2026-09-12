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
- Ship `DIALibraryGenerator-macos-<arch>.dmg` alongside the tarball and point the
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

## Rename the tool and the repository to DIALibGen

`DIALibraryGenerator` is 21 characters to type at a prompt, and the code already
says `DIALIBGEN` wherever it had to be short — `DIALIBGEN_MODEL_DIR`,
`DIALIBGEN_DATA_DIR`, `DIALIBGEN_VERSION`. The rename makes the tool agree with
itself.

197 occurrences across 49 files, but most are textual. The ones that carry
identity, and the order they have to be done in, are what make this more than a
find-and-replace:

**Breaks on rename, plan for it**

- **Release asset names.** `DIALibraryGenerator-macos-arm64.tar.gz` and the
  twelve others are named in `ci.yml`, `windows.yml`, `release-complete`'s
  assertions, both Homebrew casks and the README. A cask pins a digest to a URL,
  so old releases keep old names for ever and only new ones move.
- **Both cask tokens**, and the tap repository itself
  (`homebrew-dialibrarygenerator`). A changed token is a different cask to
  Homebrew: existing users get no upgrade, they get an orphan. Casks support
  `old_tokens` for exactly this — use it rather than asking people to uninstall
  and reinstall.
- **The desktop app's identity.** `tauri.conf.json` carries
  `productName: DIALibraryGenerator` and `identifier: de.openms.dialibrarygenerator`.
  The identifier is what macOS and the MSI use to decide whether an install is an
  upgrade or a second copy; changing it strands the installed app rather than
  replacing it. Decide deliberately whether to change it at all — the identifier
  need not follow the name.
- **`share/DIALibraryGenerator/`**, which `dataDirCandidates()` searches. An
  installed tree from a previous version is not found by a renamed binary, which
  is fine, but the two must move together.
- **The embedded recipe and the cache fingerprint.** Both name the tool. A
  rename should invalidate old caches rather than silently match them; check
  whether `fingerprintParams` needs another version bump.

**Straightforward**

- `src/DIALibraryGenerator.cpp` and `cmake/DIALibraryGeneratorConfig.cmake.in`
  (filenames), `project()` and the target name, the exported CMake package name
  and therefore every `find_package(DIALibraryGenerator)` downstream.
- `gui/package.json`, `gui/package-lock.json`, `gui/src-tauri/Cargo.toml` and
  `Cargo.lock` (crate `dialibrarygenerator-gui`) — all five files the version
  guard already checks, so it will catch a half-done rename the same way it
  catches a half-done bump.
- The doc URL in the tool description, which points at
  `TOPP_DIALibraryGenerator.html` upstream.
- `README.md`, `CHANGELOG.md` (historical entries stay as they were), the tap's
  README, and the `software-signing` skill's per-repo table.

**Sequence.** Rename the repo last: GitHub redirects the old path, so the tap and
any clone keep working while the rest lands. Cut a release under the new names
before pointing the casks at them, or the updater will compute digests for URLs
that do not exist yet.
