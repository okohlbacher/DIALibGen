# Backlog

Things that are known, not done, and not forgotten. Each carries the reason it
is here rather than done.

## One-step mode: generate, tune and refine in one invocation

Adapting a library to a run takes at least two invocations: `generate`, then
`tune`, `refine`, or `refine -tune`, which already tunes before it writes
observed values. The user has to know which of about ninety options matter. Most users want one thing: *a library for this
FASTA that fits this run*.

Proposed command:

```bash
DIALibGen -mode auto -in proteins.fasta -ids first_pass_report.parquet \
  -out adapted.tsv -instrument timsTOF
```

What it would do, in order, with one provenance record covering all stages:

1. Generate the library from the FASTA (or reuse an existing one given with
   `-library`, skipping this step).
2. Tune the RT and CCS models on the report and re-predict the **whole**
   library, so that unidentified precursors are adapted too.
3. Optionally write observed values for the identified precursors while keeping
   the full library. After step 2 the whole library is in the run's RT units,
   which is the only combination in which that write-in is unit-safe.

Simplified parameters. A small top-level set maps onto the existing namespaced
options; everything else keeps its default and stays reachable through the
namespaced options, INI or JSON.

| Simple option | Maps to | Default |
|---|---|---|
| `-instrument NAME` | `-generation:instrument` and its NCE default | `QE` |
| `-nce N` | `-generation:nce` | instrument default |
| `-charges 2 3 4` | `-generation:precursor_charges` | as `generate` |
| `-missed_cleavages N` | `-generation:missed_cleavages` | `1` |
| `-variable_mods ...` | `-generation:variable_modifications` | none |
| `-observed none\|rt\|rt+im` | refinement write-in flags with `-no_filter` | `none` |
| `-effort fast\|default\|thorough` | `train:epochs`, `stop:*`, `cohort:*` presets | `default` |
| `-keep_models DIR` | `-tune_out_models` | temporary |

Why it is here rather than done:

- The report has to come from a first-pass search. DIALibGen does not run a
  search engine, and the report's precursors must match the library's digest
  and modification settings. The mode needs a compatibility check that fails
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

## `generate` is several times slower than DIA-NN's predictor

The same whole-proteome FASTA and digest (20,431 proteins and as many shuffled
entrapment twins; 2.65 M peptides, 5.75 M precursors at charge 2-4, 69 M
transitions), on one 192-core node with 96 threads each:

| Tool | Wall | CPU | Peak memory | Output |
|---|---|---|---|---|
| DIALibGen 0.10.1 `generate`, timsTOF, NCE 30 | 18 min 26 s | 3 h 47 min | 14.5 GB | TSV |
| DIA-NN 2.0 `--predictor --gen-spec-lib` | 2 min 43 s | | | `.speclib` |

About 7 x. Not yet a fair measurement: the two builds started together and
shared the node for DIA-NN's 2.7 minutes, DIALibGen wrote a 69 M-row TSV while
DIA-NN wrote its binary library, and the models differ (AlphaPeptDeep RT, MS2
and CCS against DIA-NN's own). The same predictor is also `-run`'s largest step
after extraction (1,464 s for 5.14 M targets and their decoys at 16 threads),
so it counts against the 60-minute resource gate as well.

Why it matters: `generate` is the first thing every user runs, and its cost
grows with every variable modification, charge state and missed cleavage.

Why it is here rather than done: where the 18 minutes go (digest, the three
models, writing) has not been measured.

## Benchmark the ONNX predictors on a GPU (H100)

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
   **Done on perf/run-speed** (with B2 decoy tokens and IonLadder, and B3:
   the searched pairs keep the prefilter's assays): the full timsTOF run
   1:26:29 -> 1:16:16, the prefilter 1,427 -> 1,028 s, mostly 0.11.2's
   predictor; report table identical (design document, "The prediction port
   and the decoys"). The OpenSWATH cap stays at 8 (B1 re-measured: 16 is
   slower on both runs, a passive wait policy changes nothing). Next for
   `-run`: the extraction straggler window (one window ends each timsTOF call
   alone), the calibration's OpenSWATH call (270 s on 1.8 cores), and the serial
   hashing and loading (about 100 s).
3. **CPU cost per prediction**, the remaining ~1.3x gap to DIA-NN at 8-16 threads
   (9,100 vs 6,386 CPU-s at 8 threads). First profile ONNX Runtime per operation
   (MS2 transformer layers, RT/CCS LSTMs). Anything that changes outputs
   (reduced precision; int8 was rejected) needs an owner decision.
4. **GPU re-benchmark** once 0.11.2's host-side changes are in (target: well
   under 100 s per whole proteome on one H100), and decide on bundling
   cuDNN 9 for CUDA 12 with a GPU package (about 1.5 GB of NVIDIA libraries,
   licence notices, no GPU runner in CI).
5. **Reuse predictions across runs** (`-run` predicts the same library for every run).
6. **Hardware cache counters** (`perf`) once IBMI allows them, to plot cache
   misses per phase. The 2026-10-11 analysis inferred cache behaviour from
   placement and batch-size experiments instead.

## Final performance comparison: four platform arms

The final performance report compares the same workloads on four arms:

| Arm | Machine | Compute used |
|---|---|---|
| Laptop CPU | MacBook Pro, Apple M4 Max (12 performance + 4 efficiency cores, 128 GB) | CPU only |
| Laptop CPU + GPU | same MacBook Pro, 40-core Apple GPU | CPU plus GPU (ONNX Runtime CoreML execution provider for prediction; LibTorch MPS for fine-tuning) |
| Server CPU | IBMI H100 node (2 x AMD EPYC 9654, 192 cores) | CPU only |
| Server CPU + GPU | same node, NVIDIA H100 (one and two GPUs) | CPU plus CUDA, cuDNN 9 for CUDA 12 |

Workloads: whole-proteome `generate`, `tune` of the RT (and CCS) model on one
run's identifications, and `refine`; the prediction step of `-run` once it is
on main. Compare each arm against DIA-NN 2.0 on the same machine (its predictor
is CPU-only) and report wall time, CPU time, peak memory, energy where the
platform reports it (macOS powermetrics, nvidia-smi power draw), and output
equivalence against the server-CPU library (byte- or content-identical on CPU;
per-column differences for GPU paths, TF32 off by default).

Rules as for the 2026-10-11 analysis: every DIALibGen run traced, at least
4 threads, one run per configuration plus a stated noise measurement for the
headline numbers, quiet machines (load recorded), the same FASTA, settings and
binaries per release.

Open work before this can run:
- **Apple GPU path.** DIALibGen has no CoreML/MPS execution path yet. ONNX
  Runtime's CoreML provider must support the PeptDeep graphs (check operator
  coverage; fall back loudly per model); LibTorch MPS for training. Equivalence
  against the CPU output must be measured.
- **A signed macOS build with the speed changes,** which 0.11.2 provides for the
  CPU arm.
- **Server arm:** spock or data. Both need cuDNN 9 for CUDA 12; see
  HANDOFF-cuDNN-GPU-nodes.md. Timing must avoid other users' load.
