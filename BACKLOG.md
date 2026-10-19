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

## Fragment intensities: all-charge b/y from raw spectra, from now on

Owner rule (2026-10-19): every fragment-intensity benchmark scores b **and**
y ions at all fragment charges, observed in the raw MS2 spectra. The method is
`scripts/eval/fragment_intensity_eval.py`, described in
[docs/testing.md](docs/testing.md#fragment-intensity-accuracy-standard-method).
A search engine's fragment export is no longer the observation; it is only a
consistency check.

**To do at the next report update:** replace the intensity section of
`docs/report/` (0.12.1: "Fragment intensities", and the summary-table row
"median cosine to observed y-ion intensities"). The current numbers come from
DIA-NN 1.7.12's six singly charged y ions per precursor. A draft paragraph and
table are in the evaluation notes outside the repository (project folder,
`reports/astral-ms2-params/ASTRAL-MS2.md`, Appendix A).

Measured on Astral run 01 (42,785 precursors, all b/y ions with z ≤ 2, unobserved
ions counted as 0):

| | DIALibGen | DIA-NN 2.0 | run 02 (replicate) |
|:--|--:|--:|--:|
| old score: singly charged y ions, 40,419 precursors (report 0.12.1) | 0.986 | 0.977 | |
| raw-spectrum score | **0.941** | 0.839 | |
| raw-spectrum score, 34,980 precursors identified in both runs | 0.953 | 0.855 | 0.974 |

Acceptance: the report's intensity numbers, figure and limits come from this
script, with the ion set, convention, tolerance and precursor count stated.
The replicate ceiling is given when a second run exists.

## A linear RT scale for untuned libraries: measured, not adopted

The timsTOF benchmark suggested that the stock RT model's raw output, being
less linear in run time than DIA-NN's iRT (straight-line map 2.83 against
1.81 min), was why DIA-NN 2.0 found 4.1 % fewer precursors with DIALibGen's
untuned library on timsTOF. The raw output saturates for late peptides: on
every run checked, raw 0.85-0.88 spans the last fifth of the gradient.

Measured (2026-10-19): a fixed monotone map (26 knots) fitted on six public
DIA runs of other studies (PXD016662, PXD022589 ×2, PXD029738, PXD047184,
PXD047901; DIA-NN 2.0 with untuned DIALibGen libraries), applied to the
untuned library before the search. It makes the RT more linear everywhere
(straight-line error on PXD029836 2.94 → 2.60 min, on Astral 0.372 → 0.321 min)
and changes DIA-NN 2.0's precursors at 1 % FDR by:

| data (runs) | 0.12.1 raw RT | linear map | change |
|:--|--:|--:|--:|
| Astral, ProteoBench (01 / 02 / 03) | 89,232 / 90,362 / 89,318 | 90,889 / 91,718 / 90,669 | +1.9 / +1.5 / +1.5 % |
| timsTOF Pro, PXD041421 (3197 / 9898) | 83,296 / 109,045 | 84,317 / 109,692 | +1.2 / +0.6 % |
| timsTOF Pro 2, PXD029836, charges 2-4 (01 / 02) | 130,064 / 129,908 | 127,723 / 127,527 | −1.8 / −1.8 % |

A map fitted on PXD029836 run 01 itself (the best possible straight
line for that gradient) also lost identifications there (−0.5 / −0.7 %), with
DIA-NN's RT window narrowed from 8.9 to 7.6 min. The non-linear RT scale is
therefore not why that benchmark trails DIA-NN untuned, and the map was not
made the default: it changes every library for a mean gain of about 0.2 % per dataset
over these datasets and loses on the one it was meant for. Candidate, if the
owner wants it: an opt-in `-generation:rt_scale linear` (the implementation,
the table and the fitting script are kept in the project's evaluation notes,
`release-validation/timstof-improvements/`). The remaining untuned timsTOF gap
is open; candidates are fragment intensities (the collision energy, below and
`timstof.yaml`) and the z=3 fragments.

## Fragment charge 3 for 3+ and 4+ precursors (needs a retrained MS2 head)

Measured on the timsTOF benchmark (PXD029836, raw diaPASEF spectra, all b/y
ions at z 1-3): z=3 fragments carry 24 % of the observed intensity of 4+
precursors (y3+ 17 %, b3+ 7 %) and 6 % of 3+, and doubly charged b ions 14 %
of 3+. DIALibGen predicts no z=3 fragment; DIA-NN 2.0 neither.

**No public model provides it** (checked 2026-10-19): AlphaPeptDeep's
`pretrained_models.zip`, `_v2` and `_v3` (release tag `pre-trained-models`)
have one MS2 model, `generic/ms2.pth`, with 8 output channels (b/y z1-2 and
their modloss variants), the channels the shipped ONNX model has; peptdeep
1.5.1 and alphaDIA 2.2.0 default to a maximum fragment charge of 2. The output
layer is fixed by the checkpoint: peptdeep re-initialises it (random weights)
when z=3 types are asked for, so a z=3 model must be trained. Outside the
family, Prosit_2023_intensity_timsTOF (Zenodo 8211811, CC-BY-4.0 weights)
predicts b/y z1-3, but it is a different network (GRU, length ≤ 30, only
carbamidomethylation and oxidation), a second predictor rather than a drop-in.

The smallest path: fine-tune `generic/ms2.pth` with
`charged_frag_types` b/y z1-3 (only the last linear layer changes, 64 → 6 or
→ 12 with modloss), instrument timsTOF, on public timsTOF spectra (PXD043844,
ProteomeTools on timsTOF Pro, and diaPASEF transfer libraries extracted with
z=3); export with OpenMS's `export_peptdeep_models_to_onnx.py` unchanged. In
DIALibGen: take the channel count from the model's output shape, allow
`max_fragment_charge` 3, pin the new model's hash, add the derived-weights
notice. Acceptance: the C++ output matches peptdeep within 1e-5; the
all-charge cosine of 3+/4+ precursors rises on a diaPASEF run the model did
not see, with no loss of DIA-NN identifications.

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
3. **CPU cost per prediction**, the remaining ~1.3x gap to DIA-NN at 8-16 threads
   (9,100 vs 6,386 CPU-s at 8 threads). First profile ONNX Runtime per operation
   (MS2 transformer layers, RT/CCS LSTMs). Anything that changes outputs
   (reduced precision; int8 was rejected) needs an owner decision.
4. **GPU re-benchmark** (done in 0.12.1: 77.7 s per whole proteome on one
   H100). cuDNN 9 is not bundled; the Linux packages ship ONNX Runtime's CUDA
   build, and CUDA 12 and cuDNN 9 come from the system.
5. **Reuse predictions across runs** (`-run` predicts the same library for every run).
6. **Hardware cache counters** (`perf`) once IBMI allows them, to plot cache
   misses per phase. The 2026-10-11 analysis inferred cache behaviour from
   placement and batch-size experiments instead.
