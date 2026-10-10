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

## Outputs that DIA-NN reads completely

Checked with DIA-NN 2.0 and a post-hardening 0.11.0 build:

- DIA-NN loads the TSV's precursors but reports "0 protein groups": it does
  not use the `Protein.Group` column of a library TSV. Users must add
  `--fasta ... --reannotate --met-excision`, and a first-pass report made
  without it is useless to `refine`/`tune`. In our benchmark pipeline, writing
  the column as `Protein.Ids` holding bare accessions made DIA-NN load protein
  groups; other header names have not been tested.
- DIA-NN refuses the Parquet output: it requires its own layout, including the
  `Flags` and `Protein.Ids` columns. A DIA-NN-2-layout Parquet export would let
  one file serve as both record and search library.
- The `-out` help text, and with it `docs/parameters.md`, the GUI help and the
  class documentation, still calls both formats a "DIA-NN spectral library".
  The source string should say that `.tsv` is DIA-NN's dialect and `.parquet`
  is DIALibGen's own layout; the generated reference then follows.

Not done because the first two change an output contract that existing
pipelines read; they need a versioned option and round-trip tests against
DIA-NN.

## Variable terminal modifications are silently ignored

`-generation:variable_modifications 'Acetyl (N-term)'` exits 0 and writes the
same library as a run without the option; nothing is logged. With
`-generation:max_variable_modifications 2` the same option yields `.(Acetyl)`
forms. `Amidated (C-term)` behaves the same way. `Acetyl (Protein N-term)` is
never applied, at any setting. Fixed terminal modifications work.

Cause: OpenMS `ModifiedPeptideGenerator::applyVariableModifications` matches
by residue in its one-modification path, and generation keeps no protein
position, so protein-terminal specificity cannot be honoured. To do: apply
residue-less terminal modifications at the default limit, track
protein-terminal peptides (including the Met-excised N-terminus), and until
then warn or refuse when a requested variable modification matched no peptide.
Add a regression test that counts modified forms, and check how DIA-NN and
OpenSWATH read the `.(Acetyl)PEPTIDE` notation.

## The SciexTOF warning is wrong

Generation warns that `SciexTOF` predictions "will be close to Lumos". At equal
NCE they are not: the median spectral cosine against `Lumos` is 0.85, the same
distance as the trained `timsTOF` label, while `ThermoTOF` is at 0.99. Split the
warning in `src/Generate.cpp`, and measure `SciexTOF` against `Lumos` on SCIEX
data before recommending either.

## OpenSWATH export

`TargetedFileConverter` does not read DIALibGen's TSV, which uses DIA-NN's
column names. [Usage](docs/usage.md#openswath) documents a column mapping. A
native `-out library.pqp` or an OpenSWATH-dialect TSV would remove that step;
it is not done because decoy semantics differ (see the note in usage) and
the export needs its own round-trip test against OpenMS.

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
host work. Decision (owner, 2026-10-10): keep TF32 ON on GPUs (ONNX
Runtime's default). GPU libraries then differ slightly from CPU ones
(~1 % of precursors pick a different fragment set, CCS up to ~1.9 A^2);
TF32 off would match the CPU but cost ~26 s per proteome.
