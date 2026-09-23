# Usage and migration to 0.11

Run `DIALibGen --helphelp` for the options compiled into your executable.
The [parameter reference](parameters.md) is generated from its TOPP INI schema.

## Choose a mode

| Task | Command |
|---|---|
| Predict from FASTA | `DIALibGen -in proteins.fasta -out predicted.tsv` |
| Apply observed values | `DIALibGen -mode refine -in predicted.tsv -ids report.parquet -out refined.tsv` |
| Learn RT/CCS and predict the full library | `DIALibGen -mode tune -in predicted.tsv -ids report.parquet -out tuned.tsv` |
| Tune, then apply observed values, in one call | `DIALibGen -mode refine -tune -no_filter -in predicted.tsv -ids report.parquet -out adapted.tsv` |

The output format follows the extension. Write `.tsv` for a file a search
engine reads; DIA-NN 2 does not read DIALibGen's `.parquet`. Write `.parquet`
for a record that carries its own recipe, or for a library that only a later
DIALibGen step reads. The examples below use `.parquet` in that second sense.

Generation is the default. The [desktop app](../gui/README.md) exposes all three
modes, including RT/CCS fine-tuning controls and optional model export.

Non-default settings belonging to a different mode are rejected. A complete
default INI remains usable, but generation settings do not configure refinement
or tuning, and refinement flags do not configure generation.

## Generation settings

Every editable generation setting has a `-generation:<key>` TOPP option. For example:

```bash
DIALibGen -in proteins.fasta -out library.parquet \
  -generation:enzyme 'Trypsin/P' \
  -generation:precursor_charges 2 3 \
  -generation:peptide_length 7 30 \
  -generation:fixed_modifications 'Carbamidomethyl (C)' \
  -generation:variable_modifications 'Oxidation (M)' \
  -generation:instrument timsTOF -generation:nce 30 \
  -generation:n_terminal_methionine_excision true
```

Quote each modification containing spaces. Boolean generation options take
`true`/`false`; numeric ranges take two values. In JSON, use `[]` to clear a
modification list, for example `{"fixed_modifications": []}` for a
non-alkylated preparation.

Variable terminal modifications need `-generation:max_variable_modifications 2`
or higher. With the default of `1`, `Acetyl (N-term)` or `Amidated (C-term)` is
accepted but silently not applied, and `Acetyl (Protein N-term)` is never
applied. Check the result, for example with `grep -c '(Acetyl)' library.tsv`.
Raising the limit also allows two residue-specific variable modifications per
peptide, so the library grows.

```bash
DIALibGen -write_config defaults.json
DIALibGen -config defaults.json -generation:nce 35 -write_config effective.json
DIALibGen -write_ini workflow.ini
```

For generation, precedence is explicit command-line settings, then supplied
INI settings, then JSON, then built-in defaults. An INI written by
`-write_ini` materializes settings: those entries count as supplied when the
INI is loaded. Remove entries that should instead come from JSON. Unknown
JSON keys and invalid values are errors.

The default instrument is `QE`. Choose the instrument used for acquisition;
accepted aliases include `Astral`/`Orbitrap Astral` (Lumos), `Exploris` (QE),
`timsTOF Pro` (timsTOF), and `ZenoTOF` (SciexTOF). Omitting NCE uses the
instrument default: 25 for Lumos and its aliases, 30 otherwise. The native
option `-generation:nce -1` means automatic, and generated INI files preserve
that automatic setting. The effective JSON config records the resolved NCE and
its source; reusing that numeric value fixes NCE until it is changed or the
native option is set back to -1.

Predicted RT is normalized model output by default. Set
`-generation:irt_rescale true` for iRT calibration. It uses the bundled
standards unless you also pass `-irt_standards standards.tsv`; that option has
no effect on its own. CCS and derived 1/K0 are distinct
quantities; `-generation:derive_ion_mobility false` omits the conversion.

## Refinement

```bash
DIALibGen -mode refine -in predicted.parquet -ids report.parquet \
  -out refined.parquet -write_im -out_report residuals.tsv
```

The default workflow filters the library to identified precursors and
replaces their RT with observed values. Precursor, peptide/global and protein
q-value gates default to 0.01. An enabled gate requires its corresponding
column; missing columns are errors. Reports containing multiple runs must be
split before use. When only some matched precursors have observed RT, RT writing
requires whole-library RT tuning first; otherwise predictions and observed minutes
could be mixed. Use `-no_write_rt` to retain the original values.

Useful options:

- `-write_im`: also use observed 1/K0, subject to charge and ramp limits.
- `-no_filter`: retain unidentified precursors. Pair it with `-no_write_rt` to
  keep library RT values, or with `-tune -tune_heads rt` (or `both`) to re-predict
  the whole library in reference-run minutes before writing observed RT. CCS-only
  tuning does not establish compatible RT units and is refused in this combination.
- `-empirical_library`: explicitly declare a pre-filtered empirical reference
  whose missing report gates are recorded as bypassed.
- `-out_report residuals.tsv`: save residuals measured before replacement.
- `-write_intensity`: use observed fragment intensities only when the report
  contains supported fragment annotations and values. DIA-NN 2 fragment-column
  layouts that cannot be interpreted are refused; this is not a universal
  report-to-library conversion.

A refined library is specific to the reference's RT domain and acquisition
method. Residuals after an observed-value overwrite would be zero by
construction, so they are not evidence of predictive accuracy. Filtering and
observed write-in also change the hypotheses available to the search engine;
validate downstream error rates with controls appropriate to that comparison.

## Model tuning

```bash
DIALibGen -mode tune -in predicted.parquet -ids report.parquet \
  -out tuned.parquet -tune_out_models tuned-models \
  -train:lr 0.0001 -machine:threads 4
```

Tuning adapts RT and CCS with LibTorch and re-predicts the complete input
library with ONNX Runtime. It preserves precursor keys, m/z and fragment
intensities; it does not apply observed-value refinement to that output.
Train on a representative reference run and evaluate transfer on separate
runs. Tuning needs at least 100 units in each of the protein-held-out
validation and test cohorts, in practice several hundred identified precursors;
smaller reports are refused.

The bundled models are the default starting point. `-tune_models DIR` selects
another model directory. `-tune_heads rt` or `ccs` selects one head;
`both` is the default. `-tune_out_models DIR` retains ONNX models and training
sidecars instead of using temporary model files. Tuned ONNX files also embed
provenance metadata.

Training options are grouped under `filter:`, `cohort:`, `train:`, `stop:` and
`machine:`. Protein-level cohorts and validation govern checkpoint selection;
the tool reports rejection and fallback decisions. Options such as
`-cohort:full_fit` or `-cohort:no_inner_val` change the hold-out contract and
must be reported when presenting accuracy results. CPU is the portable
release configuration. Source builds can request `-machine:device cuda:0`
when linked to a compatible CUDA LibTorch installation.

Tune mode uses training filters such as `-filter:q_value`. Refinement-only
settings such as `-q_precursor`, `-q_global`, `-q_protein` and observed-value
replacement flags are refused when changed from their defaults in tune mode.

To generate a new library using saved models:

```bash
DIALibGen -in proteins.fasta -out adapted.parquet \
  -generation:rt_model tuned-models/peptdeep_rt_dynamic.onnx \
  -generation:ccs_model tuned-models/peptdeep_ccs_dynamic.onnx \
  -generation:free_cysteine_rt_correction false
```

Here generation emits the tuned model's normalized RT. Tune mode instead
multiplies it by the recorded `rt_max_minutes` to return reference-run minutes.
The stock free-cysteine offset is disabled because the tuned model can learn
that effect itself. Keep the model's scale and acquisition context with the
result; these are not generic iRT models.

## Files and provenance

Input libraries can be DIALibGen Parquet or DIA-NN-dialect TSV. Output format
follows `-out`'s extension. The Parquet layout (one row per precursor, fragments
as lists, recipe in the metadata) is DIALibGen's own: DIA-NN 2 does not read it,
so write `.tsv` for a search and `.parquet` for the record. Generation embeds the effective recipe and FASTA/model
hashes in Parquet. Refinement and tuning additionally write JSON provenance
sidecars beside their output, including input hashes and the selected mode.
Each tuned head includes its full recipe, cohort counts, filtering, seed,
stopping policy and evaluation even when temporary model files are discarded.
Model digests are SHA-256; library, report and FASTA cache digests are explicitly
labelled FNV-1a64. These cache fingerprints are not cryptographic integrity checks.
Parquet embeds the same provenance under `odia.config_json`; plain TSV has no
schema metadata.

Existing library outputs are refused. Writers stage complete files in a private
directory beside the destination and rename them into place only after a
successful close. Refinement prepares its library, provenance and optional
report before publishing the library; a system interruption can leave an orphan
sidecar. A forced termination (including GUI Cancel) may also leave a hidden
`.dialibgen-tmp-*` directory beside the destination. After the process exits,
that abandoned directory can be removed. Atomic replacement prevents readers
from seeing partial writes; it does not promise durability after a power loss.

The training seed is reproducible within the same build and device setup;
numeric results can differ across platforms, runtimes and CUDA kernels. The
wall-clock training budget is checked at epoch boundaries and includes report
loading and validation; final evaluation and export can exceed that budget.

Keep the provenance and any training sidecars with the library. A successful
write does not establish that a library improves a downstream search; the
[historical benchmark](benchmark.md) gives the evidence available so far.

## Using a library with a search engine

### DIA-NN

Write the library as `.tsv`. The points below were checked with DIA-NN 2.0.

```bash
diann --f run01.d --lib predicted.tsv --fasta proteins.fasta --reannotate --met-excision \
      --out first_pass.parquet --threads 16
```

- DIA-NN 2.0 refuses DIALibGen's Parquet; it expects its own Parquet layout.
- DIA-NN loads the TSV's precursors but not its `Protein.Group` column ("0
  protein groups"). `--fasta` with `--reannotate` restores the annotation.
  Without it DIA-NN leaves `Protein.Group` empty and sets every protein q-value
  to 1. `tune` holds its validation and test cohorts out by protein, drops rows
  without a protein group and then stops with "no usable observations after
  filtering"; `refine` stops at its protein q-value gate.
- Reannotation matches library peptides against DIA-NN's own digest, so add
  `--met-excision`: DIALibGen removes the initiator methionine by default
  (`-generation:n_terminal_methionine_excision true`) and DIA-NN does not.
  Without the flag those N-terminal precursors (1,132 of 47,801 in a
  600-protein test library) stay without a protein group; `tune` counts them as
  "no protein group" and `refine` drops them at the protein gate. A non-tryptic
  `-generation:enzyme` additionally needs the matching `--cut`.
- DIA-NN warns that reannotation should not be combined with a raw-data search.
  The combined call is the form checked end to end here. The separate step is
  `diann --lib predicted.tsv --fasta proteins.fasta --reannotate --met-excision
  --gen-spec-lib --out-lib predicted.diann.parquet` (seconds, no raw data). Its
  output is DIA-NN's own Parquet; DIA-NN filters fragments when it writes it (the
  test library lost 9 precursors and 15 % of its transitions), and a search with
  it was not checked. Keep `predicted.tsv` as the `-in` of `tune`/`refine`.
- Pass one `--f` per run; DIA-NN does not take a shell glob after a single
  `--f`. For timsTOF data DIA-NN advises `--mass-acc 15 --mass-acc-ms1 15`.
- `refine` applies precursor, global and protein q-value gates of 0.01. If no
  observation passes, for example because a small library gives DIA-NN too few
  proteins for protein-level FDR, it stops with "no reference observation passed
  the gates" instead of writing an empty library.
- Libraries are generated without decoys by default; DIA-NN makes its own.

### OpenSWATH

OpenSWATH's `TargetedFileConverter` expects different column names. The mapping
is one-to-one:

| DIALibGen TSV | OpenSWATH TSV |
|---|---|
| `Precursor.Mz` | `PrecursorMz` |
| `Product.Mz` | `ProductMz` |
| `Relative.Intensity` | `LibraryIntensity` |
| `RT` | `NormalizedRetentionTime` |
| `IM` | `PrecursorIonMobility` |
| `Modified.Sequence` | `ModifiedPeptideSequence`, and without the bracketed names `PeptideSequence` |
| `Precursor.Charge` | `PrecursorCharge` |
| `Protein.Group` | `ProteinId` |
| `Precursor.Id` | `TransitionGroupId` |
| `Fragment.Type`, `Fragment.Charge`, `Fragment.Series.Number` | `FragmentType`, `FragmentCharge`, `FragmentSeriesNumber` |
| `Decoy` | `Decoy` |

[`scripts/to_openswath.py`](../scripts/to_openswath.py) applies it. The script
is a helper in the source repository, not part of the release archives; it
needs Python 3 with pandas:

```bash
DIALibGen -in proteins.fasta -out predicted.tsv -generation:irt_rescale true
python3 scripts/to_openswath.py predicted.tsv openswath.tsv
TargetedFileConverter -in openswath.tsv -out library.pqp
OpenSwathDecoyGenerator -in library.pqp -out library_decoys.pqp
```

Checked with OpenMS 3.5 on a 47,810-precursor library: the PQP and TraML keep
every precursor and all 497,895 transitions with RT and ion mobility, modification
names become UniMod accessions, and decoy generation succeeds. This was checked with Carbamidomethyl (C) and Oxidation (M). OpenMS writes the
mobility under the TraML term "ion mobility drift time" with unit millisecond;
the number is still 1/K0 in Vs/cm². A search with `OpenSwathWorkflow` has not
been tested; it also needs iRT assays (`-tr_irt`) for your sample, either the
spiked iRT peptides or confidently identified endogenous ones. Use `-generation:irt_rescale true` so
that `NormalizedRetentionTime` is on the iRT scale, and generate without decoys:
a DIALibGen decoy carries its target's sequence with shifted fragment m/z, which
is DIA-NN's convention and is misread by tools that recompute fragments from the
sequence. `tune` reads DIA-NN reports only; `refine` also accepts a pre-filtered library
with `-empirical_library`. Neither reads OpenSWATH result files.

## Planned: one-step workflow

Adapting a library to a run currently means `generate` followed by `tune`,
`refine`, or `refine -tune`, which already tunes before writing observed values.
A planned `-mode auto` takes a FASTA and a first-pass report and performs all
of it with a small set of top-level options
(`-instrument`, `-nce`, `-charges`, `-observed`, `-effort`). It is **not
implemented yet**; the design, the simplified parameter set and the open
questions are in the [backlog](../BACKLOG.md#one-step-mode-generate-tune-and-refine-in-one-invocation).
Until then, [the README's full-cycle example](../README.md#examples) shows the
explicit commands.

## Migration

- Existing `DIALibGen -in FASTA -config JSON -out FILE` commands remain valid.
  Generation now also exposes those settings directly to TOPP workflows.
- Replace `DIALibRefine ...` with `DIALibGen -mode refine ...` for observed-value
  reconstruction. Refinement option names are retained.
- Replace `DIALibRefine -tune ...` with `DIALibGen -mode tune ...` for whole-library
  adaptation. Tune mode does not first filter or overwrite observations; use a
  separate refine invocation when observed-value reconstruction is intended.
- `-threads` now defaults to 1, consistently with TOPP. Use `-threads 0` for
  automatic inference parallelism or an explicit count. `-machine:threads`
  controls training separately.
- Met excision now includes every allowed missed-cleavage peptide. Generated
  peptide space can increase; cache fingerprint v5 prevents reuse of libraries
  made with the old omission. Frozen 0.10.1 benchmark libraries are historical
  artifacts and should not be overwritten by regenerated outputs.
- Generation accepts distinct precursor charges 1–8 and fragment charges 1–2,
  matching the supported prediction/refinement contract. It rejects missing or
  non-finite model predictions and removes assays with no surviving fragments.
  Ambiguous-residue peptides are counted and skipped, while valid peptides from
  the same protein remain eligible; a terminal FASTA `*` is accepted.
- A source build without `DIALIBGEN_BUILD_FINETUNE` supports generation and
  observed-value refinement; requesting tune reports the missing capability.
