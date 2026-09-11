# Changelog

All notable changes to this project are documented here.
This project follows [Semantic Versioning](https://semver.org/).

## [0.1.0] — unreleased

First standalone release. Extracted from
[OpenDIAlyzer](https://github.com/okohlbacher/ODIA) at commit `70a21a2`, where
this code was developed; the history there is not carried over.

### Added
- Standalone CMake project building `odia_library` and the `DIALibraryGenerator`
  tool, installable and consumable via `find_package(DIALibraryGenerator)`.
- `nlohmann/json` and Boost declared as dependencies. Both were previously
  satisfied only by a single shared conda prefix, which hid the requirement.
- `THIRD-PARTY-NOTICES.md`; Apache-2.0 attribution corrected on
  `include/odia/PeptDeepElements.h`, which is derived from AlphaPeptDeep
  constants and was stamped BSD-3-Clause.
- `example/proteins.fasta` and `example/default.json`.
- Tests that need no model or fixture: config round-trip, unknown-key rejection,
  `-write_config`.

### Changed
- Test and script fallbacks no longer hardcode cluster-absolute paths; they
  resolve through `$ODIA_OPENMS` / `$CONDA_PREFIX`.

### Known limitations
Carried over unchanged from extraction — see README, "Known limitations".
No behaviour of the generator was modified in the split.
