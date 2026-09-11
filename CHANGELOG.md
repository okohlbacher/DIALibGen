# Changelog

All notable changes to this project are documented here.
This project follows [Semantic Versioning](https://semver.org/).

## [0.2.2] — 2026-09-11

### Fixed
- **The embedded recipe did not name models supplied via `DIALIBGEN_MODEL_DIR`.**
  It was built before the model search ran, so a library produced that way
  recorded `"rt_model": ""` and could not state what produced it. The cache
  fingerprint was unaffected — it hashes model contents — so nothing but reading
  the recipe back would have shown it.
- **Out-of-range config numbers were accepted.** nlohmann converts a negative
  JSON number to `std::size_t` without complaint, so `"missed_cleavages": -1`
  became 18446744073709551615 and the digest ran on it. A precursor charge of 0
  divided by zero when the m/z was formed, and a `min_relative_intensity`
  outside 0..1 kept everything or nothing. Every numeric key is now
  range-checked, and the tool is tested against its own `-write_config` output
  so a default can never fall outside the ranges it enforces.

### Added
- `test/e2e_test.py`: runs the tool the way a user does and reads the result
  back — RT domain, CCS and 1/K0 ranges, fragment normalisation, decoy methods
  being genuinely different, the embedded recipe, the DIA-NN TSV. Registered
  with ctest and skipped without the models.

### Changed
- The README's build recipe is now exactly the package set CI installs, rather
  than a second list maintained by hand.

## [0.2.1] — 2026-09-11

### Added
- **Homebrew tap**:
  [okohlbacher/homebrew-dialibrarygenerator](https://github.com/okohlbacher/homebrew-dialibrarygenerator),
  with a cask for the desktop app, a cask for the CLI and a formula. The CLI
  cask uses `command_wrapper` rather than a `binary` stanza: OpenMS resolves its
  data directory with `_NSGetExecutablePath`, which does not follow symlinks, so
  a symlink into `$(brew --prefix)/bin` would send the data lookup **and** the
  `@executable_path/../lib` dylib closure to `/opt/homebrew`, where neither is.
  The releases are unsigned, which on current macOS means a Homebrew-installed
  copy is killed by Gatekeeper on first run; the tap's README documents that and
  the tarball route that is unaffected.

### Fixed
- **The desktop app declared macOS 11.0 as its minimum and meant 14.0.** That
  number is enforced by macOS, so the app launched on a system where the CLI
  bundled inside it — built for 13.3, since libc++ shipped the floating-point
  `std::to_chars` there — could not run at all. Tauri's default of 11.0 was
  never true for this app. `brew audit` is what noticed, by comparing the
  plist against the cask's `depends_on macos:`.

## [0.2.0] — 2026-09-11

The release that makes this a standalone tool rather than a tool that happens to
build. Everything in *Fixed* below was a defect an OpenMS-internal tool would
never have, and that this one had only because it links OpenMS from outside.

### Added
- **Desktop app** (`gui/`): Tauri 2 + React, for macOS, Windows and Linux. The
  form is generated from the tool's own `-write_config` output, so it cannot
  offer a parameter the CLI lacks, cannot default one differently, and surfaces
  a newly added parameter with no code change. Named presets, last-used state,
  a model-directory picker, streamed log and cancellation.
- **CI** on Linux x64/arm64 and macOS x64/arm64, plus GUI checks (TypeScript,
  vitest, `cargo test`, `cargo clippy -D warnings`) on all three OSes.
- `test/standalone_test.sh`: the bare-environment gates — runs under `env -i`,
  reports its own version, no Qt network noise, and refuses bad input. Run in CI
  against the **installed** tree, not the build tree.
- `test/write_ctd_test.sh`: `-write_ctd` produces a descriptor that agrees with
  `-write_ini` parameter for parameter, and the registration shim does not leak
  into `--help`.
- `cmake/version-guard.cmake`: every file carrying the version is checked
  against `project(VERSION)` at configure time.

- **Homebrew tap** for macOS:
  [okohlbacher/homebrew-dialibrarygenerator](https://github.com/okohlbacher/homebrew-dialibrarygenerator),
  with a cask for the desktop app and one for the CLI. Neither is signed, so a
  Mac with Gatekeeper enforcing refuses them on first run; the tap's README
  covers that and the release archive, which is unaffected. The CLI cask uses
  `command_wrapper` rather than a `binary` stanza: a symlink into
  `$(brew --prefix)/bin` would make the tool resolve its own `share/` relative
  to `/opt/homebrew/bin`, where it is not.
- `CMAKE_OSX_DEPLOYMENT_TARGET` (13.0). Without it the macOS floor was whatever
  SDK the build machine happened to have, and it moved when a CI runner image
  was retired — which reaches a user as an app that installs, launches and then
  cannot run.
- Release assets have FIXED names (`DIALibraryGenerator-gui-macos-arm64.dmg`,
  …). Tauri spells them with its own version and architecture, which would make
  every README link and every cask URL depend on Tauri's naming.

### Fixed
- **The tool reported OpenMS's version.** `version_` and `verboseVersion_` are
  now set from `project(VERSION)`; `--helphelp` carries both numbers, which is
  what a bug report needs.
- **OpenMS's update check ran on every invocation**, printing
  `QIODevice::read (QNetworkReplyHttpImpl): device not open` to stderr whenever
  it failed — an error that reads like this tool's. Disabled unless the user
  sets `OPENMS_DISABLE_UPDATE_CHECK` themselves.
- **`-threads` was registered twice and the tool's copy was dead.** TOPPBase
  registers it after `registerOptionsAndFlags_`, with a default of 1, and wins
  at parse time — so every run that did not pass `-threads` used one core and
  one inference session. The default is now applied to argv before TOPPBase
  parses it: 0 = all available cores.
- **`-write_ctd`, `-write_cwl` and `-write_json` aborted** with "Requested tool
  'DIALibraryGenerator' does not exist!" — `ToolHandler`'s registry is a
  hard-coded list of official TOPP tools. The tool now registers itself as an
  internal tool for the duration of a descriptor run; OpenMS still writes every
  byte of the descriptor. `-write_cwl`/`-write_json` refuse in this tool's own
  words on an OpenMS without `ENABLE_TDL`, rather than letting OpenMS truncate
  the target and then throw.
- **The data directory was a compile-time absolute path into the source tree**,
  so an installed, packaged or relocated binary resolved its iRT standards out
  of a build directory that may not exist. Resolution is now relative to the
  executable, with `DIALIBGEN_DATA_DIR` as an override; `data/` installs to
  `share/DIALibraryGenerator` (it previously went to a lower-case directory
  nothing looked in).
- **An unknown `decoys` value was silently treated as `mutate`** and then
  written into the provenance under the name that was typed. Refused, with the
  known methods listed.
- **`-out` accepted any extension**, writing a DIA-NN TSV for anything that was
  not `.parquet` — so `library.parqet` produced a TSV under a Parquet name.
- **`schema_version` was accepted at any value and never read.** A config
  written for a schema this build does not know was interpreted with this
  build's meaning for every key. Refused.
- **The fingerprint under-determined the library.** `min_relative_intensity` was
  formatted at three decimals, so `1e-4` and `0.0` shared a key; and no RT-domain
  token was recorded, so a raw-RT and an iRT library of the same content
  fingerprinted identically. Both closed; the key is now `v3`.
- **A missing model died inside the ONNX session constructor** with
  `Load model from  failed`. The tool now searches `DIALIBGEN_MODEL_DIR`, its
  own `share/`, and `share/OpenMS/models`, and names the file it could not find
  along with every directory it looked in.
- **`-write_config` required an `-in` and an `-out`** that it returned before
  reading. Both are optional now and enforced where they are actually needed.

### Changed
- Apache Arrow floor lowered from 23 to **19**, which is what the code actually
  uses (`parquet::arrow::OpenFile`'s `Result` overload). 23 excluded bioconda's
  libarrow 21 — the version OpenMS itself pins — so a consumer linking both
  could not satisfy it.
- `LibraryGenerator::fingerprintParams` takes an `irt_rescale` argument.
- Copyright headers say "the DIALibraryGenerator authors", not "the ODIA
  authors".

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
