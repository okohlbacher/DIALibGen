# Changelog

All notable changes to this project are documented here.
This project follows [Semantic Versioning](https://semver.org/).

## [0.11.2] — 2026-10-11

### Breaking

- **TSV column names.** The `.tsv` library is now DIA-NN's TSV format under
  the column names that both DIA-NN 1.7 and 2.x read: `ModifiedPeptide`,
  `PrecursorCharge`, `PrecursorMz`, `iRT`, `IonMobility`, `ProductMz`,
  `LibraryIntensity`, `FragmentType`, `FragmentCharge`, `FragmentSeriesNumber`,
  `FragmentLossType`, `Decoy`, instead of the dotted names (`Modified.Sequence`,
  `Precursor.Mz`, ...). DIA-NN 1.7 could not read the dotted names at all. New
  columns `UniprotID` (accessions), `ProteinName` (entry names) and
  `Proteotypic` give DIA-NN its protein groups; with `Protein.Group` alone both
  versions reported "0 protein groups" unless searched with `--fasta
  --reannotate`. Modifications are written as UniMod accessions
  (`C(UniMod:4)`, `(UniMod:1)PEPTIDE`) instead of OpenMS names.
  `Precursor.Id`, `Protein.Group` and `CCS` stay. **Scripts that read
  DIALibGen TSVs by the dotted names must switch to the new ones**; DIALibGen
  itself reads both. The Parquet layout is unchanged.
- **Exit code.** Generation now refuses, with exit code 6 (ILLEGAL_PARAMETERS)
  and no output, a requested modification it would silently not apply: a
  protein-terminal specificity (`Acetyl (Protein N-term)`, fixed or variable;
  not supported yet), or a variable modification that leaves no modified
  precursor in the whole library. Until 0.11.1 such runs exited 0 with an
  unmodified library.

### Added

- `-mode convert` rewrites a library in the format of `-out` without
  changing its content, e.g. Parquet to DIA-NN TSV.
- `peptdeep_ms2_pruned.onnx`: the MS2 model without its unused
  modification-neutral-loss branch, shipped beside the three stock models and
  run by default. The b/y channels generation reads are bit-identical to the
  full model's, so the library content is the same; MS2 needs about 15 % less
  CPU. The recipe names and hashes the model that ran (so a library equals
  0.11.1's in content, not byte for byte). `"ms2_model":
  "peptdeep_ms2_dynamic.onnx"` selects the full model. Made by
  `scripts/prune-ms2-modloss.py`; `dialibgen-fetch-models` verifies it.
- Opt-in performance tracing: `ODIA_TRACE=<dir>` writes phase, call, batch and
  ONNX Runtime events with /proc samples, `ODIA_ORT_PROFILE` ONNX Runtime's
  per-operator profile ([tracing](docs/tracing.md)). Off, it costs nothing
  measurable.
- `-machine:require_gpu`: fail instead of running a model on the CPU when its
  CUDA session cannot be created.
- `-machine:tf32 true|false` (default `false`): TensorFloat-32 on NVIDIA GPUs,
  for prediction and CUDA training. Off, a GPU library matches the CPU library
  to float rounding; on is faster but changes ~1 % of fragment sets.
- `-machine:gpu_sessions` (default 4) CUDA sessions per GPU and model, and
  `-machine:gpu_devices` to spread them over several GPUs.

### Changed

- **Generation on the CPU is 6.5 times faster at 96 threads** (whole human-sized
  proteome, 5.8 M precursors, 2 x 96-core EPYC: 1,025 -> 158 s),
  1.7 times at 16 threads (1,018 -> 595 s) and 1.5 times at 8 (1,729 -> 1,156 s),
  with the same library content:
  - one ONNX Runtime session per thread instead of at most 16, sharing one
    copy of the weights, fed 64-row batches costliest first, so that every
    session stays busy to the end of each call;
  - the pruned MS2 graph (above);
  - RT and CCS predicted first and MS2 last, with the library written while
    MS2 runs; nothing appears under the output name before it is complete;
  - fragment m/z without per-ordinal sequence copies, the encoder without
    ModificationsDB's lock, parallel RT parsing, at most 16 digest threads;
  - above 32 threads on a multi-socket machine the run's memory is
    interleaved over the NUMA nodes (-16 % at 96 threads; `ODIA_NUMA_INTERLEAVE`
    overrides).
  The output does not depend on `-threads`. Peak memory at 96 threads is
  9.5 GiB (0.11.1: 14.0), at 16 threads 4.7 GiB (13.9).
- `-threads 0` means one thread per **physical** core of the effective CPU
  budget, not per hardware thread: SMT siblings share a core's arithmetic
  units, and sessions on them ran slower. On Linux that is the cores of the
  affinity mask (cgroup cpusets included), capped by a cgroup CPU quota
  (`docker --cpus`, Kubernetes limits); on Windows the cores of all processor
  groups; on macOS `hw.physicalcpu`. The log says how the count was arrived at.
- Tune's re-prediction (`-tune_predict_sessions 0`, the default) uses
  `-threads` sessions; until 0.11.1 it used up to 16 whatever `-threads` said.
  Set `-threads` (or `-tune_predict_sessions`) to the CPUs it may use.
- When a thread cannot be started (a process or cgroup thread limit), the
  work runs on the threads that exist instead of the process aborting.
- The TSV writer holds at most ~0.7 GB of rendered text, whatever the thread
  count (about 4 GB at 96 threads before), and formats each precursor's
  repeated fields once per precursor; the bytes are unchanged.
- ONNX Runtime thread pools never spin while waiting (where more than one
  intra-op thread is used, spinning cost 29-43 % more CPU); Arrow's thread
  pools follow `-threads` instead of the hardware-thread count. The desktop
  app starts the CLI with `OPENBLAS_NUM_THREADS=1` (OpenBLAS otherwise starts
  a spinning thread per hardware thread, ~14 CPU-s per run on a large node);
  command-line users can set it themselves.
- Residue-less peptide-terminal variable modifications (`Acetyl (N-term)`,
  `Amidated (C-term)`) are applied at the default `max_variable_modifications`
  1; they were dropped silently. The log counts the precursors per variable
  modification. Generation cache key v6.
- The `SciexTOF` warning no longer says its predictions are close to Lumos:
  at the same NCE their median spectral cosine to Lumos is 0.914 (ThermoTOF:
  0.993, which is now its own warning).
- The `-out` help text says what each format is: `.tsv` is DIA-NN's TSV
  format, `.parquet` DIALibGen's own layout, which DIA-NN does not read.
- A model requested on CUDA that falls back to the CPU on a machine with a GPU
  (typically: cuDNN 9 for CUDA 12 is missing, so RT and CCS cannot use the GPU)
  now prints a `WARNING` with the reason and the fix, and Parquet output records
  it as `odia.warnings`. Previously the only sign was `[RT] provider CPU`.
- GPU prediction no longer uses TF32 by default (ONNX Runtime's CUDA provider
  enabled it). GPU libraries now equal CPU libraries to float rounding.
- Digestion, precursor enumeration, MS2 fragment selection (overlapped with
  the next block's prediction) and the Parquet and TSV writers run on
  `-threads` host threads. One H100 builds a whole-proteome library in
  91 s instead of 154 s (0.11.1 with cuDNN, TF32 on), two H100s in 65 s.
  Refine and tune write their output on `-threads` threads as well.

### Fixed

- A CUDA build of ONNX Runtime on a machine without an NVIDIA driver no longer
  reports a GPU fallback: no `WARNING`, no `odia.warnings` in the library.

### Documentation

- Usage: speed, threads and memory; modifications (what is refused, DIA-NN's
  `--mod` for UniMod accessions it does not know, such as `UniMod:2`); the MS2
  model. Tracing guide. GPU section (cuDNN 9 for CUDA 12, fallback warnings,
  TF32, sessions); a tested CUDA LibTorch build recipe. On one H100, tuning
  the RT head and re-predicting a 3.5 M-precursor library takes about a
  minute instead of 12.5 min on 64 CPU threads, with the same held-out
  accuracy.

## [0.11.1] — 2026-10-10

First published 0.11 release. 0.11.0 was validated but not published; its
changes below are part of this release.

### Changed

- Compact Parquet libraries store `Precursor.Id`, `Modified.Sequence` and
  `Protein.Group` as plain strings instead of one global Arrow dictionary.
  Parquet repeated that dictionary as the dictionary page of every row group,
  so readers decoded every string of the library once per row group. Parquet
  now dictionary-encodes each row group until its dictionary page reaches
  1 MiB and writes the rest plain. On a 5.81 M-precursor library the file
  shrinks from 1,038 to 765 MB with identical content; DIALibGen loads it in
  9.9 s instead of 21.0 s and writes it in 13.0 s instead of 24.6 s. pyarrow
  and pandas now see `string` instead of `dictionary<string>` (pandas: `str`
  or `object` instead of `category`); DIALibGen 0.10.1 and 0.11.0 read old
  and new files alike. The flat writer keeps its dictionary builders.
- The Linux AppImages embed AppImage type2-runtime revision 8f39b89, which
  creates the directories it extracts into (`--appimage-extract`,
  extract-and-run) with mode 0700 instead of 0755.
  The build now packs a vendored copy of that pinned runtime instead of
  upstream's mutable `continuous` download.

## [0.11.0] — 2026-09-20

See the [validation record](docs/testing.md#release-verification) for tested
contracts, platform checks and coverage limits.

### Added

- One TOPP executable for generation, observed-value refinement and RT/CCS
  model tuning: `-mode generate|refine|tune`.
- Native `-generation:<key>` options for every generation setting, including
  TOPP INI/CTD export. Explicit CLI/INI settings override JSON configuration.
- Refinement and training from DIALibRefine, with report quality gates,
  provenance sidecars and optional saved ONNX models. Tuned ONNX models embed
  their provenance.
- Whole-library tune mode preserves precursor keys and predicts unseen entries
  without observed-value write-in. Portable training uses CPU; CUDA remains
  available with suitable source-build dependencies.
- Usage/migration documentation, generated parameter reference and a qualified
  summary of the completed 0.10.1 cross-run benchmark.
- Instrument aliases and per-instrument NCE defaults. Unknown instruments are
  refused; the effective recipe records the selected NCE and its source.
- Runtime dependency inventories, upstream license texts and platform-specific
  corresponding-source archives, including desktop dependencies.
- Desktop generation, refinement and fine-tuning modes, with RT/CCS head
  selection, native training settings, model export and mode-aware saved settings.

### Changed

- `-threads` now follows the TOPP default of 1. Set 0 for automatic inference
  parallelism. Training uses its separate `-machine:threads` setting.
- The desktop app explicitly selects the requested mode and uses the CLI's
  defaults and training schema. Fractional NCE and m/z settings remain editable;
  thread counts must be whole numbers. Existing output paths are refused before
  starting a child process.
- User documentation replaces obsolete pre-release notes and same-run tuning
  recommendations. Superseded Python tuning helpers are not part of the tool.
- Portable archives contain the executable, runtime libraries, data, models
  and notices. Headers, static libraries and CMake development files remain
  available through source installation rather than the portable archives.

### Fixed

- Training provenance now records true SHA-256 model hashes. Earlier refiner
  fields with that name held a 64-bit fingerprint; library/report cache
  fingerprints are now explicitly labelled FNV-1a64.
- Met excision now retains N-terminal peptides at every permitted missed-cleavage
  count, including cases where the fully cleaved peptide is below the minimum
  length. Cache fingerprint v5 invalidates libraries containing the old omission.
- Variable modifications no longer duplicate unmodified precursors. Duplicate
  charges, unsupported fragment charges and non-finite model output are refused;
  empty assays are removed consistently across output formats.
- Decoys retain terminal modifications, reject self-identical sequences and use
  a shuffle pinned across standard libraries. Calibration data and decoy mass
  policy now participate in cache keys.
- Atomic library and model writes preserve existing files on failure. Parquet
  loads use bounded record batches; strict TSV/Parquet validation covers nulls,
  truncated rows, fragment annotations, numeric domains and missing mobility.
- Refinement validates the report before training, retains the full training
  recipe in library provenance and counts unique join keys. Terminal-modification
  aliases, min-max RT scaling and fragment-quality counters are corrected.
- Keeping unidentified precursors while replacing RT now requires whole-library
  RT re-prediction in reference-run minutes. This prevents mixed RT units;
  keeping the original RT values remains available with `-no_write_rt`.
- Refinement resolves additional named modifications through OpenMS with residue
  and terminal specificity. Unresolved library tokens produce a warning and a
  provenance count. Compact empirical Parquet references are refused explicitly.
- Fragment-intensity replacement requires matching target/decoy fragment sets,
  including when fragment restriction is disabled.
- Library operations check the 32-bit transition capacity before allocation or
  append, refusing oversized inputs instead of wrapping stored offsets.
- Training validates its public API controls, measures elapsed wall time and
  rejects cross-charge protein-group conflicts before assigning held-out cohorts.
- Windows builds prioritize the Arrow headers belonging to the linked runtime,
  avoiding a Parquet reader ABI mismatch with OpenMS contrib headers.
- Desktop health probes, launch/cancel races, automatic NCE, settings recovery,
  temporary files and IPC dispatch now have regression coverage.
- Update Vite and Vitest to patched versions, removing the reported development
  server and test-server dependency vulnerabilities.

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
  not stuck. The `.dmg` is stapled and unaffected. See the
  [subsequent startup measurements](docs/testing.md#macos-startup).

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
