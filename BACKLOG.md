# Backlog

Things that are known, not done, and not forgotten. Each carries the reason it
is here rather than done.

## One-step mode: generate, tune and refine in one invocation

**The owner's next feature candidate** (asked to add it, 2026-10-19).

Adapting a library to a run takes at least two invocations: `generate`, then
`tune`, `refine`, or `refine -tune`, which already tunes before it writes
observed values. Most users want one thing: *a library for this FASTA that
fits this run*.

Proposed command:

```bash
DIALibGen -mode auto -in proteins.fasta -ids first_pass_report.parquet \
  -out adapted.tsv -generation:instrument timsTOF
```

What it would do, in order, with one provenance record covering all stages:

1. Generate the library from the FASTA (or reuse an existing one given with
   `-library`, skipping this step).
2. Tune the RT and CCS models on the report and re-predict the **whole**
   library, so that unidentified precursors are adapted too.
3. Optionally write observed values for the identified precursors while keeping
   the full library. After step 2 the whole library is in the run's RT units,
   which is the only combination in which that write-in is unit-safe.

Already in 0.12.1, so not part of this item: instrument defaults from YAML
files (`-generation:instrument` sets NCE, m/z windows, charges and ion
mobility), `-tune_heads auto` (CCS only when the report has 1/K0), mode
inference without `-mode` (a FASTA generates, a library converts), `-threads`
defaulting to the physical cores, and a 23-option `--help`. Still missing:

| Proposed | Maps to | Default |
|---|---|---|
| `-mode auto`, or inferring it from a FASTA `-in` with `-ids` (today refused) | `generate`, then `tune` or `refine -tune` | |
| `-observed none\|rt\|rt+im` | refinement write-in flags with `-no_filter` | `none` |
| `-effort fast\|default\|thorough` | `train:epochs`, `stop:*`, `cohort:*` presets | `default` |

`-tune_out_models` keeps the models as in `tune`; without it they are temporary.

Why it is here rather than done:

- The report has to come from a first-pass search (or `-run`, experimental),
  and the report's precursors must match the library's digest and
  modification settings. The mode needs a compatibility check that fails
  early with a useful message instead of training on a handful of matches.
- The default for `-observed` needs evidence. The 0.10.1 benchmark found no
  reliable benefit from observed-value write-in, so the proposal defaults to
  `none`.
- `-effort` presets must be derived from measurements, not chosen by feel.
- Acceptance test: one `-mode auto` command produces the same precursors and
  transitions as `generate` followed by `tune` (for `-observed none`) or by
  `refine -tune -no_filter` (for `-observed rt` and `rt+im`): TSV
  byte-identical, Parquet identical apart from provenance metadata. A separate
  `tune` followed by `refine -no_filter` is not the reference; refine refuses
  observed RT with the filter off unless the RT head is tuned in the same call.

## Done

- Fragment intensities scored on b and y ions at every fragment charge, from the
  raw spectra (owner rule, 2026-10-19; `scripts/eval/fragment_intensity_eval.py`,
  docs/testing.md): the report's intensity section and summary row use it (Astral:
  DIALibGen 0.941, DIA-NN 2.0 0.839, replicate 0.974), and its timsTOF section.

- Outputs that DIA-NN reads completely: the TSV carries DIA-NN's protein
  columns since 0.11.2; `-out_format diann2_parquet` writes DIA-NN 2.0's
  Parquet layout (0.11.3).
- Terminal modifications: peptide-terminal ones apply at any
  `max_variable_modifications` (0.11.2), protein-terminal ones apply to
  protein-terminal peptides (0.11.3); a requested modification that modifies
  nothing is refused.
- The SciexTOF warning (0.11.2).
- OpenSWATH export: `-out library.pqp` and `-out_format openswath_tsv`,
  written by OpenMS (0.11.3).
- Persistent worker pool (0.11.3; item 1 below).
- GPU timing on the H100 (item 4 below): `generate` of the whole proteome in
  77.7 s on one H100 and 48.3 s on two (0.12.1, report section 3); the Linux
  packages ship ONNX Runtime's CUDA build but not CUDA or cuDNN (0.12.1).
- Four-arm performance comparison (laptop CPU, laptop with the Apple GPU,
  server CPU, server with one and two H100): the 0.12.1 report
  (`docs/report/`, section 3). Energy was not measured.

## Accepted trade-offs (owner decisions)

- `-run` on Windows: the identification report's `Evidence` differs between
  `-threads 1` and `4` in about the ninth digit (MSVC, Intel OpenMP; Linux and
  macOS write the same bytes); q-values and the refined library are the same,
  so IDs are unaffected. Accepted by the owner 2026-10-19. Evidence:
  `identify_e2e` checks Evidence on Windows within 1e-6 relative (commit 14e74ba).
- Linux packages are about 760-800 MB because they ship ONNX Runtime's CUDA
  build (0.12.1). Accepted by the owner 2026-10-19.
- CoreML (Apple GPU) libraries equal the CPU's to float rounding (intensities
  within 1e-4, RT and CCS identical), not byte for byte. Accepted by the owner
  2026-10-19.

## `generate` needs 1.5 x DIA-NN's CPU time at low thread counts

Measured again for 0.13.0 (report section 3): the same whole-proteome FASTA
and digest (5,808,438 precursors at charges 2-4; DIALibGen keeps 12 fragments
per precursor), both tools one after the other on one quiet 2 x 96-core node
(kim, 20 October 2026), DIALibGen 0.13.0 from the release package:

| threads | DIALibGen 0.13.0 | DIA-NN 2.0 predictor | DIALibGen CPU / DIA-NN CPU |
|---|---|---|---|
| 96 | 151 s | 180 s | 12,924 / 11,448 s |
| 16 | 690 s | 464 s | 10,937 / 6,619 s |

The "several times slower" of 0.10.1 (18 min 26 s against 2 min 43 s) is
gone: at 96 threads DIALibGen is now the faster tool. At 16 threads it still
takes 1.5 x DIA-NN's wall time and 1.65 x its CPU time.
The two do not do identical work (different models, fragment selection and
output format), so this is an end-to-end comparison.

Why it is here rather than done: the remaining gap is CPU cost per
prediction (item 3 of the follow-ups below: profile ONNX Runtime per
operation for the MS2 transformer and the RT/CCS LSTMs); anything that
changes outputs (reduced precision) needs an owner decision. On a GPU the
question does not arise (49 s on two H100).

## Benchmark the ONNX predictors on a GPU (H100)

Status 0.12.1: the whole-proteome timing, tune on CUDA and the GPU-versus-CPU
differences are measured (Done above, report sections 3 and 4). Still open:
per-model throughput settings (batch size, sessions, CUDA streams, fp16) and
the prediction step of `-run` on the GPU.

Every timing so far is CPU-only. On an H100 node, `generate` with the MS2 model
on CUDA took 6:10-6:27 wall for the whole proteome, about 2.7 x the 96-thread
CPU run (a side measurement during a review, two runs sharing one GPU, not a
benchmark). GPU output is deterministic but not byte-identical to CPU output.

To do: measure RT, MS2 and CCS throughput on one H100 (batch size, number of
sessions and CUDA streams, fp32 vs fp16/TF32, device transfer and host-side
encoding/post-processing as the new bottleneck), whole-proteome `generate`
and the prediction steps of `tune` and `-run`, against the tuned CPU path and
DIA-NN on the same node; quantify how far GPU output differs from CPU output.

Why it is here rather than done: the CPU analysis comes first, and its
fixes (session limit, batch sizes) change what the GPU path is compared
against.

Measured 2026-10-10 on one H100 (spock): whole proteome 154 s with RT, MS2
and CCS on CUDA (needs cuDNN 9 for CUDA 12; without it RT and CCS silently
fall back to the CPU: 310-395 s), 126 s on two H100s; ~94 s of it is serial
host work. Decision (owner, 2026-10-10, revised): TF32 is an option and
OFF by default, so GPU libraries match CPU libraries; ON (opt-in) saves
~26 s per proteome on one H100 but changes ~1 % of fragment sets and CCS
by up to ~1.9 A^2.

## Performance follow-ups after 0.11.2 (from the generate performance analysis, 2026-10-11)

0.11.2 carries the scheduling changes, shared read-only data, parallel phases,
the pruned MS2 model, streamed writing during MS2, and the small items
(OpenBLAS threads, error flag, Arrow pool). Still open, in order:

1. **Persistent worker pool** (done in 0.11.3) instead of a `std::thread` per call: node-level
   affinity, half the workers per socket, each worker first-touching its own
   arena. This removes about 5,600 thread creations and 97 glibc arenas (4.3-5.5 GB
   retained) per whole-proteome run, and the NUMA-interleave workaround.
   Estimated at least the interleave gain at 96 sessions (-10 % wall, -20 % CPU).
   Acceptance: no slower than 0.11.2 at 16 threads, and at or below the
   interleaved time at 96, with identical output.
2. **The same speed-ups for `-run` and tune** (feature/builtin-identification):
   IonLadder in `PredictedAssays::rank`, the session cap and chunking in the
   search-decoy prediction, and a larger prefilter block than 8,192 pairs. The
   prediction step is about 1,464 s of the timsTOF `-run`; an estimated -500 s.
   Acceptance: `ids.parquet` identical on a reduced cell and on one full run.
   **Done in 0.11.3** (perf/run-speed, with B2 decoy tokens and IonLadder, and B3:
   the searched pairs keep the prefilter's assays): the full timsTOF run
   1:26:29 -> 1:16:16, the prefilter 1,427 -> 1,028 s, mostly 0.11.2's
   predictor; report table identical (design document, "The prediction port
   and the decoys"). The OpenSWATH cap stays at 8 (B1 re-measured: 16 is
   slower on both runs, a passive wait policy changes nothing). Next for
   `-run`: the extraction straggler window (one window ends each timsTOF call
   alone), the calibration's OpenSWATH call (270 s on 1.8 cores), and the serial
   hashing and loading (about 100 s).
3. **CPU cost per prediction**, the remaining gap to DIA-NN at 8-16 threads
   (0.11.2 at 8 threads: 9,100 vs 6,386 CPU-s; 0.13.0 at 16 threads: 10,937 vs
   6,619 CPU-s, 1.65 x). First profile ONNX Runtime per operation
   (MS2 transformer layers, RT/CCS LSTMs). Anything that changes outputs
   (reduced precision; int8 was rejected) needs an owner decision.
4. **GPU re-benchmark** (done in 0.12.1: 77.7 s per whole proteome on one
   H100). cuDNN 9 is not bundled; the Linux packages ship ONNX Runtime's CUDA
   build, and CUDA 12 and cuDNN 9 come from the system.
5. **Reuse predictions across runs** (`-run` predicts the same library for every run).
6. **Hardware cache counters** (`perf`) once IBMI allows them, to plot cache
   misses per phase. The 2026-10-11 analysis inferred cache behaviour from
   placement and batch-size experiments instead.
