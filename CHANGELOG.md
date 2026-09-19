# Changelog

All notable changes to this project are documented here.
This project follows [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Changed — breaking
- **An unknown `instrument` is now refused instead of silently predicting for no
  instrument at all.** The MS2 model one-hot encodes the *index* of the name, and
  a name it did not know indexed `max_instrument_num - 1` — a slot still holding
  its random initialisation. Every slot of the shipped checkpoint carries
  non-zero weight, so nothing downstream could notice: the library predicted,
  wrote and read back looking entirely normal. `"Astral"` hit that slot.

  Anyone who was passing a name that "worked" this way was getting predictions
  from an untrained slot; the error now names every accepted spelling.

### Added
- **Instrument aliases**, taken from upstream's own `instrument_group`: `Astral`,
  `Fusion`, `Eclipse`, `OrbitrapTribrid` → `Lumos`; `QE+`, `QEHF`, `QEHFX`,
  `Exploris` → `QE`; `timsTOF Pro/SCP/HT/Ultra/flex` → `timsTOF`; `TripleTOF`,
  `ZenoTOF` → `SciexTOF`. Case, spaces, `-` and `_` are ignored. The resolution
  is logged and the canonical name is what lands in the provenance.
- **Per-instrument default `nce`**, applied only when the config names none, and
  resolved before `-write_config` so the dumped config is the one that was used:
  `Lumos` (and so `Astral`) 25, everything else upstream's 30. `nce_source` in
  the recipe records whether the operator or the tool chose the number.

  `timsTOF` stays at upstream's 30 although we measured 40 as better on K562
  diaPASEF (spectral angle 0.9041 ± 0.0004 against 0.8939 ± 0.0011 over three
  replicates, +695 precursors end to end): the curve falls about four times more
  steeply above its peak than below, so a default at the measured maximum puts
  every cooler collision-energy ramp on the steep side. The README says to set
  `nce: 40` for a method like ours. Most of the end-to-end win is the label, not
  the NCE — `QE`/30 → `timsTOF`/30 is +1,662 precursors and already passes
  DIA-NN on protein groups.
- `SciexTOF` and `ThermoTOF` warn that they carry no trained weights in the
  shipped checkpoint. Read out of it, only `QE` and `timsTOF` have weights
  outside the meta layer's initialisation bound; `Lumos` is the no-correction
  baseline the other two are deltas from.
- `test/tools/odia_instrument_dump.cpp` and two tests: the whole
  name → canonical → slot → NCE table is pinned, and an unknown instrument in a
  config must fail. `test/check_meta_inputs.py` now validates the alias table
  against the C++ source the same way it already validated the index map.

## [0.10.1] — 2026-09-17

### Added
- **The PeptDeep models ship with the tool.** Every release now carries
  `peptdeep_{rt,ms2,ccs}_dynamic.onnx` in `share/DIALibGen/models`, which is
  already the second entry in the tool's own search order — so an installed copy
  predicts immediately, with nothing set and nothing to download. Their
  provenance and the paper to cite are in
  [THIRD-PARTY-NOTICES.md](THIRD-PARTY-NOTICES.md).

  This also ends a papercut that had no good answer before: a Homebrew upgrade
  replaces the Caskroom directory, so models installed into it were deleted on
  every upgrade and had to be fetched again. Bundled models arrive with each
  version instead.

  `dialibgen-fetch-models` stays, for a build from source, a refresh, or
  `--check`. `DIALIBGEN_MODEL_DIR` and the config's `rt_model`/`ms2_model`/
  `ccs_model` still override the bundled copies.

### Changed
- **The release gates no longer hand the tool its models.** The end-to-end run
  against the installed tree, and the bare-environment run against the staged
  bundle, both stopped setting `DIALIBGEN_MODEL_DIR` — so what they prove is
  that a user's install finds its own models. Both workflows also assert the
  three files are in the install tree, on every platform.

## [0.10.0] — 2026-09-12

### Changed
- **Renamed to DIALibGen**, everywhere: the tool, the CMake project and its
  exported package, the installed binary, `share/DIALibGen`, the desktop app,
  the release assets and both Homebrew casks. The code already said `DIALIBGEN`
  wherever it had to be short — `DIALIBGEN_MODEL_DIR`, `DIALIBGEN_DATA_DIR`,
  `DIALIBGEN_VERSION` — so the tool now agrees with itself, and the name is
  9 characters to type rather than 21.

  What this breaks, and what it does not:
  - **Release assets are named `DIALibGen-*` from this version.** Earlier
    releases keep their own names for ever; a cask pins a digest to a URL, so
    nothing that already exists moves.
  - **Both cask tokens changed**, and carry `old_tokens` so an existing install
    upgrades rather than being orphaned.
  - **The bundle identifier `de.openms.dialibrarygenerator` is unchanged**, on
    purpose. Nobody sees it, and it is what macOS and the MSI use to decide
    whether an install is an upgrade or a second copy — changing it would strand
    every installed copy for no gain.
  - **`find_package(DIALibraryGenerator)` becomes `find_package(DIALibGen)`**,
    and the exported targets file is renamed with it.
  - `DIALIBGEN_*` environment variables are untouched; they were already right.

## [0.9.1] — 2026-09-12

### Added
- **`dialibgen-fetch-models` ships with the tool.** It downloads the three
  AlphaPeptDeep ONNX models from OpenMS's archive, checks each against a pinned
  SHA256, and installs them where the tool already looks — so nothing has to be
  set in the environment afterwards. `--check` verifies an existing set, `--dir`
  and `--prefix` place them elsewhere, and re-running it leaves correct files
  alone.
  It is installed next to the binary because the models are not redistributed
  with this project: a release tarball or a Homebrew cask is the whole of what
  most people have, and pointing them at a script in the git repository is not an
  answer. 0.9.0 shipped the script in the repository only.

### Documented
- **The first run of the macOS CLI takes about five minutes.** Measured at 335 s
  against the published 0.9.0 tarball, 1 s on every run after. macOS validates
  each of the 145 bundled libraries with Apple individually, and a `.tar.gz`
  cannot carry a stapled ticket that would answer for all of them at once. It is
  not stuck. The `.dmg` is stapled and unaffected. See BACKLOG.md for the fix.

## [0.9.0] — 2026-09-12

### Added
- **The macOS builds are signed and notarized.** With
  `Developer ID Application: Oliver Kohlbacher (9WF4NVY9MY)`, the hardened
  runtime and a secure timestamp on all 146 Mach-O files in the bundle. The
  `.dmg` and the `.app` inside it each carry a stapled ticket and so need no
  network; the CLI tarball cannot carry one, because no archive format can, and
  is checked online once on first run. This is what the Homebrew casks needed:
  an unsigned payload installed by a package manager carries a provenance record
  Gatekeeper refuses, clearing the quarantine attribute does not remove it, and
  Homebrew 6 dropped `--no-quarantine`.
  Signing is optional in CI — without the certificate secret the artefacts are
  ad-hoc signed exactly as before, so forks and pull requests still build.

### Fixed
- **The Windows CLI archive and `.msi` could not generate a library.** Neither
  carried `share/OpenMS`, so the tool exited with "Cannot find shared data!
  OpenMS cannot function without it!" on the first mass calculation. `--help`
  worked, which is why the release gate never caught it: it asked only for
  `--help`, and that reads none of OpenMS's data. Both gates now build a real
  library from the bare bundle instead, on every platform.
- **The Windows bundle carried 271 DLLs and 267 MB** — a blind sweep of every
  `*.dll` in conda's `Library/bin`, Qt's and contrib's, including all of MKL,
  gRPC, ICU and a Qt built for debugging. Replaced with the dependency closure
  `dumpbin` walks, plus the sets that are `LoadLibrary`'d rather than imported
  and so appear in no walk. 66 DLLs; the archive is 108 MB, the `.msi` 114 MB
  and the setup `.exe` 84 MB, down from 267/275/168.
- **A failed build leg could publish an incomplete release.** `fail-fast` is off,
  so the surviving legs uploaded and the release went public missing whatever
  the failed leg owed it — while `release-complete`, the job that exists to
  catch exactly that, was SKIPPED by the same failure.
- **A `workflow_dispatch` selected against a tag uploaded to the real release.**
  The ref was tested before the event, so the run whose entire purpose is to
  touch no release took the release branch. Both workflows test the event first.
- **A rejected notarization named only a submission id.** The reasons live behind
  `notarytool log`, a second round trip nobody makes on a step that takes tens
  of minutes to reach again; it is now fetched on failure.

## [0.2.3] — 2026-09-12

### Added
- **Windows x64 builds.** The CLI archive and the `.msi` installer now ship
  alongside Linux and macOS, and the full suite — 23 tests including the
  end-to-end library build against the real models — passes there.

### Fixed
- **The exported CMake target carried build-machine include paths.**
  `odia_library` attached ONNX Runtime's and Boost's include directories as
  bare `PUBLIC` paths, so `find_package(DIALibraryGenerator)` handed a consumer
  absolute paths from whatever machine built it. CMake only *errors* when such
  a path lies inside the source tree, which is why this surfaced on Windows and
  was silent everywhere else.
- **Model paths reached ONNX Runtime in the wrong character type.**
  `Ort::Session` takes `ORTCHAR_T*` — `wchar_t*` on Windows — so
  `std::string::c_str()` could not compile there. Now `std::filesystem::path`,
  which is the right type on every platform and carries a non-ASCII path
  correctly.
- The installed `DIALibraryGeneratorConfig.cmake` also looks for conda-forge's
  `onnxruntime_conda` spelling, without which `find_package` fails on Windows.

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
