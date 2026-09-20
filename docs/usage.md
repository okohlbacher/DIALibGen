# Usage and migration to 0.11

Run `DIALibGen --helphelp` for the options compiled into your executable.
The [parameter reference](parameters.md) is generated from its TOPP INI schema.

## Choose a mode

| Task | Command |
|---|---|
| Predict from FASTA | `DIALibGen -in proteins.fasta -out predicted.parquet` |
| Apply observed values | `DIALibGen -mode refine -in predicted.parquet -ids report.parquet -out refined.parquet` |
| Learn RT/CCS and predict the full library | `DIALibGen -mode tune -in predicted.parquet -ids report.parquet -out tuned.parquet` |

Generation is the default. The desktop app currently provides generation;
use the command line for refinement and tuning.

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
`-generation:irt_rescale true` for iRT calibration with the bundled standards,
or supply `-irt_standards standards.tsv`. CCS and derived 1/K0 are distinct
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
split before use.

Useful options:

- `-write_im`: also use observed 1/K0, subject to charge and ramp limits.
- `-no_filter`: retain unidentified precursors. This creates a mixture of
  observed and predicted entries; interpret it accordingly.
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
runs.

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

Input libraries can be DIA-NN-dialect Parquet or TSV. Output format follows
`-out`'s extension. Generation embeds the effective recipe and FASTA/model
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
