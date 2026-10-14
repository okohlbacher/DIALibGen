# Changelog

All notable changes to this project are documented here.
This project follows [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Changed

- `docs/report/DIALibGen-report-0.11.4.md` and `.pdf` replace the 0.11.3 report: the laptop
  timings of 0.11.4 (CPU, CoreML with a cold and a warm cache, MPS and CPU tuning), the
  Apple-GPU options of the release, and the unchanged quality, identification and entrapment
  sections; also attached to the 0.11.4 release.

## [0.11.4] — 2026-10-13

### Added

- Apple GPU ([usage](docs/usage.md#apple-gpu)). `-machine:device coreml`
  runs the MS2 model on the GPU through ONNX Runtime's CoreML execution
  provider: one compiled static-shape session per residue count, batches
  padded to 2048 rows, the models' scalar `PRelu` rewritten into `LeakyRelu`
  (the provider's MLProgram builder cannot compile it), compiled models
  cached under `~/Library/Caches/DIALibGen/coreml`. Measured on an M4 Max:
  MS2 15-18 times faster than one CPU thread; on a whole human library the
  CPU's precursors and order, identical RT/CCS, intensities within 4.9e-6,
  8 of 3.49 M fragment sets different, min cosine 0.987. RT and
  CCS stay on the CPU path, loudly: CoreML has no LSTM operator, the
  predictor reads ONNX Runtime's node placement and falls back when an
  arithmetic operator is left on the CPU (`WARNING`, `odia.warnings`;
  `-machine:require_gpu` refuses). `-machine:coreml_units gpu|all|ane|cpu`
  picks the compute units; `-machine:gpu_sessions` is the number of
  concurrent batches there.
- `-machine:device mps`: fine-tuning on the Apple GPU (LibTorch's Metal
  Performance Shaders back end; no operator falls back with LibTorch 2.10),
  with `-tune_predict_gpu` re-predicting through CoreML. Not bit-identical to
  CPU training; the provenance says so.
- `-machine:device` and `-machine:coreml_units` are accepted in `generate`.
- `scripts/coreml_unroll_lstm.py`: the LSTM-free RT and CCS graphs for
  CoreML, as an experiment that is not shipped. The unrolled graphs match
  the CPU's LSTM kernel (RT 6e-7, CCS 4.3e-4 A^2 on CoreML) but run at the
  speed of one CPU thread (16.2 k against 18.9 k rows/s at 16 tokens, 9.5 k
  against 10.2 k at 32), and the CPU path has one session per core, so RT
  and CCS stay on the CPU, overlapped with the GPU's MS2 phase instead.
- `-mode prepare` (macOS): compiles the bundled MS2 model's CoreML shapes
  for `-generation:peptide_length` into the cache once, for an installation
  or before a first run. The cache is capped (`ODIA_COREML_CACHE_MAX`, 4G by
  default): above it the least recently used shapes of other models go when
  a run's predictor ends, never the running model's own. Every CoreML MS2
  phase ends with one line of shapes built, rows run and padding rows.
- `-machine:cores performance` (macOS, every mode): `-threads 0` counts the
  performance cores only (`hw.perflevel0.physicalcpu`, 12 of the M4 Max's
  16). Opt-in, not the default: on a whole human library 12 threads needed
  18 % less CPU time than 16 for the same output, but the wall time went
  either way with the background load (-3 % at full scale, +16 % on a small
  input), so it saves energy and fan noise, not time.
- CTest `coreml_matches_cpu` (skips without the CoreML provider).
- [Building on macOS](docs/building.md#macos-apple-silicon): the conda-forge
  packages carry both back ends; environment and compiler notes.
- `docs/report/DIALibGen-report-0.11.3.md` and `.pdf`: the report on the tool, its usage and
  its performance against DIA-NN (runtime, RT and fragment-intensity accuracy, identifications),
  also attached to the 0.11.3 release.

### Changed

- `generate` on a GPU (CUDA, CoreML): the RT phase (with the iRT
  calibration) and the CCS phase run on a thread of their own while the GPU
  predicts MS2, which does not depend on them; the streaming writer waits
  for both. On a whole human library on an M4 Max the two CPU phases took
  77 s of 16 idle cores before 236 s of GPU time. The GPU phase needs cores
  for its host side, so on CoreML the overlapped phases run on half of
  `-threads` (`ODIA_OVERLAP_RESERVE`): on the Astral library, untimed and
  under a changing load, the phases in order took 238 s, overlapped on all
  16 cores 232 s (the GPU phase stretched from 153 to 229 s), on 12 cores
  208 s and on 8 cores 189 s. The RT and CCS calls are the same calls as
  before whatever the session count, so the library is the same: the Astral
  CoreML library and the one-H100 whole-proteome library are byte-identical
  with and without the overlap (`ODIA_OVERLAP=0` keeps the phases in order).
  The MS2 predictor is now opened before them (trace phase `ms2_sessions`),
  on every provider.
- CoreML: before a model is compiled in the run, a child process
  (`DIALibGen --probe-coreml`) compiles its probe shape and runs one batch,
  once per process and model, so that a Core ML that aborts the process (a
  macOS virtual machine without a usable GPU) leaves the run on the CPU
  with a `WARNING` instead of taking it down; a shape Core ML refuses inside
  the run falls back to a CPU session for its batches, loudly.
- CoreML: the shapes the library needs (its residue counts, each with the
  batch sizes a chunk of that many precursors can take; `-mode prepare`
  builds every shape of the window and keeps none in memory) are compiled,
  or loaded from the cache, in the background from the predictor's
  construction on (`ODIA_COREML_BUILDERS`, 8), in the order the chunks ask
  for them, so they
  overlap the RT and CCS phases instead of being built one after the other
  by the first worker to need each while the others waited inside the MS2
  phase (66 s of a cold-cache whole-proteome run, 16 s warm). Measured:
  Core ML compiles one model at a time whatever the builder count (42
  shapes cold: 105 s with 4, 8 or 16 builders), cached loads 23 s with 8
  builders against 34 s one by one. A worker that needs a shape still being
  built waits for that one alone.
- CoreML: from 22 residues on, batches are padded to 512, 1024 or 2048 rows
  instead of 2048 alone (`ODIA_COREML_BUCKETS`). The rarer long lengths
  padded 15-38 % of their rows, and a padded row of 25-30 residues costs
  three times a short one; on a 5.8 M-precursor library the padding falls
  from 11.5 % to 8.8 % of the rows run, for 18 more compiled shapes (about
  65 MB each, built ahead as above).
- `-machine:device` defaults to `auto`: on Apple silicon fine-tuning runs on
  the GPU (MPS) when LibTorch has the back end and one small tensor summed
  on it comes back right, elsewhere on the CPU as before (CUDA is still
  asked for by name). A LibTorch that reports MPS but cannot use it (a macOS
  virtual machine) trains on the CPU with a log line and the reason under
  `device_note` in the provenance; `-machine:device mps` by name still fails
  loudly. The Astral RT tune on an M4 Max
  took 85 s on MPS against 847 s on 16 CPU threads, with the same held-out
  metric; MPS training is not bit-identical to CPU training, and the
  provenance records the resolved `device` and the requested
  `device_requested`. `-machine:device cpu` trains on the CPU. Prediction on
  the Apple GPU stays opt-in (`-machine:device coreml`): CoreML loses on
  small inputs, where loading the compiled shapes dominates (a 400 k-precursor
  library: 32.5 s against 28.8 s on the CPU), and its MS2 output agrees with
  the CPU's to float rounding rather than exactly, so a default `generate`
  still writes the CPU library on every platform.
- `-machine:threads` defaults to `0`, the platform's fixed count: 16 on Linux
  and Windows as before, **8 on macOS**. On an M4 Max every training
  operator forks LibTorch's OpenMP team and the join waits for the four
  efficiency cores' share: 12.4 s per epoch at 16 threads against 8.1 s at 8
  (4 ties with 8, 12 takes 10.2 s), with 2,299 s of system time in an 847-s
  CPU tune (about 100 s at 8). The tuned weights depend on the count, which
  the provenance records as before; `-machine:threads 16` on a Mac reproduces
  a Linux or Windows model.
- ONNX Runtime's log goes through one process-wide logger; its two per-session
  CoreML notices (how many nodes the provider took, "some nodes were not
  assigned") are replaced by the predictor's own provider line and warning.
- The third-party notices name what of Boost the mzPeak reader compiles in
  (Boost.JSON and Boost.Compat, header-only) and which Boost headers the
  release builds use (1.90 on Linux and macOS, 1.92 on Windows);
  `licenses/README.md` no longer cites a Boost version the builds do not use.

## [0.11.3] — 2026-10-13

### Added (experimental)

- `-run run.mzML` in `refine` and `tune`: built-in identification of the run
  (fragment-evidence candidates, OpenSWATH extraction with stock OpenMS, a
  semi-supervised discriminant, target-decoy FDR), with ion mobility on
  diaPASEF runs; its report (`-out_ids`) then serves as `-ids`. Settings are
  the `-search:` options. Not validated yet: the entrapment gate has not
  passed and the 16 GB / 60 min on 8 cores resources gate is not met; `-ids`
  remains the validated route
  ([design](docs/design/built-in-identification.md)).
  Faster extraction, calibration and run loading by default (nested
  OpenSWATH window loop with a team of up to 16, calibration split by seed,
  hashes and cache removal in the background, a capped decode team): the
  full timsTOF `-run` at 16 threads takes 3,317 s instead of 4,574 s
  (-27.5 %) with the same report
  ([measured](docs/design/built-in-identification.md#-run-efficiency-defaults-0113)).
- `-run` reads mzPeak runs (`.mzpeak` archives or their unpacked directory)
  besides mzML, through the vendored mzpeak-openms reader into the same
  OpenSWATH maps; diaPASEF frames stored whole are split by their windows'
  1/K0 bands. Experimental, as `-run`.

### Added

- `-out_format` (all modes) selects the layout of `-out`; `auto`, the default,
  keeps choosing by extension, and the default `.parquet` and `.tsv` layouts
  are unchanged. The opt-in layouts carry their version
  in their name:
  - `diann2_parquet` (`.parquet`): DIA-NN 2.0's own Parquet library layout,
    the 27 columns, order and types of a library DIA-NN 2.0 writes, including
    `Flags` and `Protein.Ids`. DIA-NN 2.0 refuses DIALibGen's own Parquet
    ("missing the 'Flags' column") and loads this one completely: on a
    3,000-protein library the same protein isoforms, protein groups and
    precursors as from the TSV, and a search of a DIA slice gives an
    identical report. The recipe travels in its metadata, which DIA-NN
    ignores. DIALibGen reads it back.
  - `openswath_pqp` (`.pqp`, also chosen by `auto`) and `openswath_tsv`
    (`.tsv`): OpenSWATH's PQP and assay TSV, written by OpenMS's
    `TransitionPQPFile` and `TransitionTSVFile` from a `TargetedExperiment`.
    `TargetedFileConverter` reads both back with every precursor and
    transition, and `OpenSwathDecoyGenerator` and `OpenSwathWorkflow` run on
    the PQP. Replaces the manual column mapping. The writers hold the library
    as a `TargetedExperiment`: about 0.8 kB (PQP) and 2 kB (TSV) per
    transition.
- Desktop app: an **output layout** choice offers the five layouts (DIALibGen
  Parquet, DIA-NN TSV, DIA-NN 2 Parquet, OpenSWATH PQP, OpenSWATH TSV); it
  sets the extension and passes `-out_format` where the extension does not
  select the layout. The default is unchanged.
- Protein-terminal modifications (`Acetyl (Protein N-term)`,
  `Amidated (Protein C-term)`, fixed or variable) apply to protein-terminal
  peptides, including the N-terminus after initiator-Met excision; the digest
  records where each peptide sits in its proteins. A fixed one keeps a peptide
  that also occurs inside a protein unmodified as well. The protein N-terminal
  acetylated forms equal DIA-NN 2.0's (`--var-mod UniMod:1,42.010565,*n
  --met-excision`) on a 40-protein FASTA. They were refused in 0.11.2.

### Changed

- README restructured: an overview figure, the modes, examples, an instrument
  table and the input/output formats; a [backlog](BACKLOG.md).
- The pruned MS2 model (`peptdeep_ms2_pruned.onnx`) is removed; the full
  model, `peptdeep_ms2_dynamic.onnx`, is the default again. Libraries name the
  full model in their provenance; their content is unchanged. MS2 takes about
  15-19 % more CPU than in 0.11.2.
- **Generation: the host work runs in parallel**, with the same library
  content as 0.11.2 and 0.11.1. Whole human-sized proteome (5.8 M precursors),
  2 x 96-core EPYC 9654, traced, one run each on the same node:

  | version | 16 threads | 96 threads |
  |---|---|---|
  | 0.11.2 | 593.6 s | 157.9 s |
  | 0.11.3 | 682.8 s | 146.2 s |

  0.11.3 is 15 % slower than 0.11.2 at 16 threads and 7.5 % faster at 96:
  0.11.2 ran the pruned MS2 model, which 0.11.3 removed, and without it MS2
  needs about 15-19 % more CPU. On a GPU (16 CPU threads, same node for both
  versions):

  | version | 1 H100 | 2 H100 |
  |---|---|---|
  | 0.11.2 | 88.9 s | 65.2 s |
  | 0.11.3 | 82.4 s | 47.3 s |

  The changes:
  - parallel work runs on a persistent worker pool instead of threads started
    for every step; worker k is bound to NUMA node k mod n (Linux;
    `ODIA_POOL_PIN=0` turns it off). Memory is no longer interleaved over the
    NUMA nodes above 32 threads, which with the pool cost time
    (`ODIA_NUMA_INTERLEAVE=1` still does it);
  - modified sequences are parsed without OpenMS's process-wide ResidueDB
    lock (`ODIA_PARSE_CHECK=1` compares every parse with OpenMS's);
  - the digest merges its peptide map in parallel shards, and enumeration
    appends precursors and fragments in parallel;
  - MS2 and CCS: parsing, fragment selection and appending of the
    neighbouring blocks overlap each prediction;
  - CUDA: the probe session is kept as the first session, and the other
    sessions are built per device concurrently.
- **Libraries load in parallel**: the TSV is parsed on `-threads` threads and
  Parquet columns are decoded in parallel (`-mode convert`, `refine` and
  `tune`). Numbers are parsed with `std::from_chars` (or `strtod_l` in the C
  locale) instead of a stream per number. The loaded library is the same at
  every thread count. Loading a TSV library at 16 threads on 2 x 96-core EPYC
  9654 takes 5.7 s instead of 23.7 s. A report number that underflows to zero (below
  4.9e-324) is now rejected by refine on Linux too.
- C-terminal modifications are written as DIA-NN writes them,
  `PEPTIDEK(UniMod:2)` instead of `PEPTIDEK.(UniMod:2)`, in the TSV and the
  DIA-NN Parquet. DIA-NN 2.0 reads both and reports `PEPTIDEK(UniMod:2)`, so its
  report now matches the library text. OpenMS parses both.
- A requested protein-terminal modification that modifies no precursor is
  refused (exit code 6), as a variable one already was.
- **`-machine:threads` (CPU training threads) defaults to 16** instead of 4.
  The tuned weights depend on the thread count (it fixes the order of
  LibTorch's reductions), so at default settings the tuned models and the
  re-predicted RT/CCS differ from 0.11.2's; `-machine:threads 4` reproduces
  0.11.2's training. The default is a fixed number, independent of `-threads`
  and of the machine, and `.tune.json` and the tuned ONNX metadata record it.
- **CPU training on Linux is 1.5 to 3 times faster per epoch** (RT head, Astral run,
  2 x 96-core EPYC: 16.7 -> 11.3 s at 4 threads, 9.1 -> 4.4 s at 16,
  10.8 -> 3.5 s at 64) with bit-identical tuned weights at the same count: glibc kept handing the
  training tensors back to the system and faulting them in again at every
  step (30-37 % of training CPU was system time). It now keeps freed memory
  while a head trains (about 2 GB more resident during training) and returns
  it afterwards.
- Tune records its training thread count (`threads`, `interop_threads`) in
  `.tune.json` and the tuned ONNX metadata.
  LibTorch's inter-op pool is one thread (tune never uses it). Tracing
  (`ODIA_TRACE`) shows a phase per tuned head. At a given thread count none of
  this changes the tuned weights.

### Fixed

- Refine joins C-terminally modified precursors of a DIA-NN report
  to the library. DIA-NN reports `PEPTIDEK(UniMod:2)`, generation stores
  `PEPTIDEK.(Amidated)`, and until 0.11.2 the two keys differed, so such
  precursors were silently left unmatched. A modification OpenMS knows only
  at the C-terminus now has one key in either spelling; one that can also
  sit on the last residue (`Methyl` on K) is joined to the library's spelling
  when the library holds only one of the two forms.

### Documentation

- Usage: the output layouts, the DIA-NN 2.0 Parquet columns, OpenSWATH export
  and decoys, protein-terminal modifications and the terminal notation of each
  output.

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
