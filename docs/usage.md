# Usage and migration to 0.11

Run `DIALibGen --helphelp` for the options compiled into your executable.
The [parameter reference](parameters.md) is generated from its TOPP INI schema.

## Choose a mode

| Task | Command |
|---|---|
| Predict from FASTA | `DIALibGen -in proteins.fasta -out predicted.parquet` |
| Apply observed values | `DIALibGen -mode refine -in predicted.parquet -ids report.parquet -out refined.parquet` |
| Learn RT/CCS and predict the full library | `DIALibGen -mode tune -in predicted.parquet -ids report.parquet -out tuned.parquet` |

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
`timsTOF Pro` (timsTOF), and `ZenoTOF` (SciexTOF). `SciexTOF` and `ThermoTOF`
have no trained weights in the shipped MS2 model and print a warning:
`ThermoTOF` predicts what `Lumos` predicts (median spectral cosine 0.993 at the
same NCE), `SciexTOF` does not (0.914, a tenth of precursors below 0.64), and
neither has been validated on SCIEX data. Omitting NCE uses the
instrument default: 25 for Lumos and its aliases, 30 otherwise. The native
option `-generation:nce -1` means automatic, and generated INI files preserve
that automatic setting. The effective JSON config records the resolved NCE and
its source; reusing that numeric value fixes NCE until it is changed or the
native option is set back to -1.

### Modifications

`fixed_modifications` and `variable_modifications` take OpenMS names such as
`Carbamidomethyl (C)`, `Oxidation (M)`, `Acetyl (N-term)`,
`Amidated (C-term)` or `Acetyl (Protein N-term)`;
`max_variable_modifications` (default 1) caps the variable ones per peptide.
Peptide-terminal modifications without a residue apply at any maximum (until
0.11.1 they were dropped at 1).

Protein-terminal modifications (`Acetyl (Protein N-term)`,
`Amidated (Protein C-term)`, ...) apply only to peptides that begin or end a
protein. With `n_terminal_methionine_excision` (the default) the residue after
an initiator Met counts as the protein N-terminus as well, so both
`(UniMod:1)MSAAAGK` and `(UniMod:1)SAAAGK` are made from a protein
`MSAAAGK...`. A variable one is added to every form of such a peptide that has
room for one more variable modification. A fixed one modifies the peptide
where it is protein-terminal; a peptide that also occurs inside a protein is
kept unmodified as well. This is DIA-NN's `--var-mod UniMod:1,42.010565,*n
--met-excision`: on a 40-protein FASTA both make the same 92 acetylated forms.
Until 0.11.2 protein-terminal modifications were refused.

Generation refuses, with exit code 6 (ILLEGAL_PARAMETERS) and no output, a
requested modification that leaves no modified precursor in the whole
library: a variable one with no site within the length and m/z windows (or
`max_variable_modifications` 0), and a protein-terminal one, fixed or
variable, when no protein-terminal peptide passes the windows.

The log lists the precursors carrying each variable modification. The TSV and
the DIA-NN Parquet write modifications as UniMod accessions, terminal ones as
DIA-NN does: `(UniMod:1)PEPTIDE` for an N-terminal and `PEPTIDE(UniMod:2)` for
a C-terminal one (0.11.2 wrote `PEPTIDE.(UniMod:2)`; DIA-NN 2.0 reads both and
reports the former, and OpenMS parses both; refine joins a report to the
library in either spelling). DIA-NN knows only some accessions
by default: `(UniMod:1)` loads as is, while for example C-terminal amidation
needs `--mod UniMod:2,-0.984016` in the DIA-NN command, or DIA-NN 2.0 stops on
an unknown modification. DIA-NN itself has no C-terminal modifications
(`--var-mod` with `c` or `*c` makes none). The OpenSWATH outputs use OpenMS's
notation, `.(UniMod:1)PEPTIDE` and `PEPTIDE.(UniMod:2)`.

Predicted RT is normalized model output by default. Set
`-generation:irt_rescale true` for iRT calibration with the bundled standards,
or supply `-irt_standards standards.tsv`. CCS and derived 1/K0 are distinct
quantities; `-generation:derive_ion_mobility false` omits the conversion.

## Speed, threads and memory

`-threads` follows the TOPP default of 1. `-threads 0` uses one thread per
**physical** core of the CPUs the process may run on (its affinity mask on
Linux): hardware threads of one core share its arithmetic units, and on a
2 x 96-core node 96 sessions on all 384 hardware threads were 12 % slower than
48 on the physical cores.

On the CPU, generation runs one ONNX Runtime session per thread, never more
than the physical cores (0.11.1 capped them at 16), with one shared copy of
each model's weights. It feeds them batches of 64 peptides of one length,
longest work first, so that every session stays busy to the end of a call;
the result does not depend on the batch size or the thread count. Host work
(digestion, enumeration, fragment selection, writing) runs on the same
threads. RT and CCS are predicted first and MS2 last, and the library is
written while MS2 runs: each complete Parquet row group, or every TSV row, as
soon as its MS2 results are in. Nothing appears under the output name before
the file is complete. The worker threads are kept for the whole run and, on a
machine with several NUMA nodes, placed on the nodes in turn, so each keeps
its memory local.

Measured on 2 x 96-core EPYC 9654 for a 5.8 M-precursor library (Parquet):
1,156 s and 4.2 GiB peak memory at 8 threads, 595 s and 4.7 GiB at 16, 158 s
and 9.5 GiB at 96, and 136 s and 15.1 GiB at `-threads 0` (192 physical
cores). Past one socket's cores the speed-up costs CPU: at 192 threads the run
used 56 % more CPU-seconds than at 96, because fully loaded cores run at a lower clock. On a
shared machine a thread count well below the core count is the efficient
choice.

OpenMS's LAPACK loads OpenBLAS, which starts a spinning thread per hardware
thread before DIALibGen's code runs: about 14 CPU-seconds on a 384-thread
machine, for BLAS calls generation never makes. Set `OPENBLAS_NUM_THREADS=1`
in the environment to avoid it; the desktop app does so for the runs it
starts.

[Tracing](tracing.md) records where a run spends its time (`ODIA_TRACE`).

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
  -train:lr 0.0001 -machine:threads 16
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
when linked to a compatible CUDA LibTorch installation (see
[building](building.md#cpu-and-cuda)); add `-tune_predict_gpu` to run the
re-prediction on the GPU as well. On one H100, tuning the RT head on an Astral
run (77,029 sequences, 59 epochs) and re-predicting a 3.5 M-precursor library
took 62 s, against 12.5 min on 64 CPU threads, with the same held-out accuracy
(TEST calibrated SD 0.350 min either way). CUDA training is not bit-identical to
CPU training.

`-machine:threads` (default 16) is the number of CPU threads LibTorch trains
with. Each training step is a chain of small operations, so the count stops
paying early: for the RT head on an Astral run (77,029 sequences) on a
2 x 96-core EPYC 9654, an epoch took 11.3 s at 4 threads, 6.7 s at 8, 4.4 s
at 16, 3.7 s at 32 and 3.5 s at 64, for 380, 450, 600, 990 and 1,850
CPU-seconds of training over 8 epochs. The tuned weights depend on the thread
count (it fixes the order of LibTorch's reductions): they are reproducible
with the same build at the same count, which `.tune.json` and the tuned ONNX
metadata record as `threads`. The default is a fixed number, not derived from
`-threads` or the machine, so the same command trains the same model
anywhere the build is the same. Up to 0.11.2 the default was 4: tuned models
made with default settings differ from 0.11.2's, and `-machine:threads 4`
reproduces them.

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

## GPU

Generation runs the RT, MS2 and CCS models on an NVIDIA GPU when the ONNX
Runtime it links has the CUDA execution provider and a GPU is visible;
otherwise on the CPU. Each model prints the device it ended up on:

```text
[RT] provider CUDA device 0 sessions 4
```

Requirements on Linux: an NVIDIA driver for CUDA 12 or later, the CUDA 12
runtime libraries (cudart, cuBLAS, cuBLASLt) and **cuDNN 9 for CUDA 12**,
including the unversioned `libcudnn.so`, which ONNX Runtime opens by that
name. Without cuDNN only the MS2 model (a transformer) runs on the GPU; RT and
CCS (convolution and LSTM layers) cannot create a CUDA session and run on the
CPU, which on one H100 makes a whole-proteome library about 2.6 times slower.
Install cuDNN system-wide (`apt install libcudnn9-cuda-12 libcudnn9-dev-cuda-12`
from NVIDIA's CUDA repository) or in user space (`pip install nvidia-cudnn-cu12`
or conda `cudnn=9 cuda-version=12`), create the `libcudnn.so -> libcudnn.so.9`
link if the package lacks it, and put that directory on `LD_LIBRARY_PATH`.

When a GPU exists but a model cannot use it, DIALibGen prints a `WARNING`
naming the model, the reason (for example `dlopen failed for libcudnn.so`) and
the fix, and records it in Parquet output as the `odia.warnings` metadata key
(tune: `device_warnings` in the provenance). A machine without a GPU, or a
CPU-only ONNX Runtime, runs on the CPU without a warning.
`-machine:require_gpu` turns every CPU fallback, that case included, into an
error.

Each model runs `-machine:gpu_sessions` CUDA sessions per GPU (default 4).
Each session has its own stream, so the host-side encoding of one batch
overlaps the GPU work of another; on an H100, 2-4 sessions make RT and CCS
30-45 % faster than one and 8 is slower again. `-machine:gpu_devices 0 1`
spreads the sessions over several GPUs, round robin; by default the first
visible GPU that works is used. Neither changes the output. CPU runs use one
session per thread (see [Speed, threads and memory](#speed-threads-and-memory)).

With the models on a GPU, most of the remaining time is host work: digestion,
precursor enumeration, fragment selection and writing. These run on `-threads`
threads (the output does not depend on the count); fragment selection overlaps
the prediction of the next block, and the library is written while MS2 runs.
On one H100 a whole human-sized proteome (5.8 M precursors) took 91 s at
`-threads 16`, on two 65 s.

`-machine:tf32` controls TensorFloat-32 math on the GPU, for prediction and
for `-machine:device cuda` training. It is **off** by default, so a GPU library
agrees with the CPU library to float rounding (measured on a whole proteome:
RT within 2e-6 normalised units, CCS within 6e-4 A^2, fragment intensities
within 1.2e-5; the same fragment set for all but 32 of 5.8 M precursors, ties
at the 12-fragment cap; every spectrum with a base peak of at least 0.01 has
cosine >= 0.985 against the CPU's). `-machine:tf32 true` rounds matmul and
convolution inputs to 10 mantissa bits: before the host-side speed-ups of
0.11.2 it made one H100 about 25 s faster (82 s instead of 106 s; two H100s
69 s instead of 76 s), but about 1 % of precursors then
select a different fragment set, the worst informative spectrum has cosine
0.955 against the CPU, and CCS moves by up to ~1.9 A^2 (median 0.02 A^2).
Either way, GPU output is deterministic for a given setting: it does not depend
on `-threads`, the number of sessions or the number of GPUs.

Since 0.11.2 CCS is predicted before fragment selection. On the CPU the
library is still exactly 0.11.1's; on a GPU, where the batch composition
shapes every value at float rounding, the CCS and 1/K0 of precursors kept
after selection dropped others can differ from 0.11.1's in the last bits.
`ODIA_STREAM=0` restores the old order (and the one-pass writer).

## Files and provenance

Input libraries can be DIALibGen Parquet, DIA-NN TSV or a DIA-NN 2.0 Parquet
library. The output format follows `-out`'s extension, or `-out_format` (all
modes):

| `-out_format` | `-out` | layout |
|---|---|---|
| `auto` (default) | `.parquet`, `.tsv`, `.pqp` | `dialibgen_parquet`, `diann_tsv`, `openswath_pqp` |
| `dialibgen_parquet` | `.parquet` | DIALibGen's compact layout with provenance; DIA-NN does not read it |
| `diann_tsv` | `.tsv` | DIA-NN's TSV library, read by DIA-NN 1.7 and 2.x (below) |
| `diann2_parquet` | `.parquet` | DIA-NN 2.0's own Parquet library layout (below) |
| `openswath_pqp` | `.pqp` | OpenSWATH's PQP, written by OpenMS ([OpenSWATH](#openswath)) |
| `openswath_tsv` | `.tsv` | OpenSWATH's assay TSV, written by OpenMS |

The default outputs are unchanged: a 0.11.3 `.parquet` or `.tsv` is byte for
byte what 0.11.2 writes. Generation embeds the effective recipe and FASTA/model
hashes in Parquet. Refinement and tuning additionally write JSON provenance
sidecars beside their output, including input hashes and the selected mode.
Each tuned head includes its full recipe, cohort counts, filtering, seed,
stopping policy and evaluation even when temporary model files are discarded.
Model digests are SHA-256; library, report and FASTA cache digests are explicitly
labelled FNV-1a64. These cache fingerprints are not cryptographic integrity checks.
Parquet embeds the same provenance under `odia.config_json`; plain TSV has no
schema metadata.

### The TSV library

`.parquet` is DIALibGen's own layout; DIA-NN does not read it (DIA-NN 2.0
stops: the file "is missing the 'Flags' column"). `-out_format diann2_parquet`
writes DIA-NN's layout instead (below). `.tsv` is DIA-NN's TSV library format, and DIA-NN 1.7 and 2.x both read the whole
file: precursors, fragments, RT, ion mobility, decoy flags, protein
accessions and entry names, and proteotypicity. No `--fasta --reannotate`
is needed to get protein groups. The columns:

| column | content | read by DIA-NN as |
|---|---|---|
| `Precursor.Id` | `<sequence><charge>`, `_decoy` appended for a decoy | not read (DIALibGen's own key) |
| `ModifiedPeptide` | sequence with modifications as UniMod accessions, e.g. `C(UniMod:4)`, `(UniMod:1)` for an N-terminal one | precursor sequence |
| `PrecursorCharge`, `PrecursorMz` | | charge, m/z |
| `Decoy` | 0 or 1; a decoy carries its target's sequence | decoy flag |
| `iRT` | library RT, as stored (see below) | reference RT |
| `IonMobility` | predicted or observed 1/K0 | ion mobility (2.x; 1.7 has none) |
| `ProductMz`, `LibraryIntensity` | | fragment m/z, relative intensity |
| `FragmentType`, `FragmentCharge`, `FragmentSeriesNumber`, `FragmentLossType` | `b`/`y`, charge, ordinal, `noloss`/`H2O`/`NH3` | fragment annotation |
| `UniprotID` | accessions of the protein group, `P1;P2` | protein ids (2.x: `Protein.Ids`) |
| `ProteinName` | entry names, `NAME1_HUMAN;NAME2_HUMAN` | protein names |
| `Proteotypic` | 1 if the group has exactly one protein | proteotypicity |
| `Protein.Group` | the group verbatim, `sp\|P1\|NAME1_HUMAN;...` | not read |
| `CCS` | collision cross-section | not read |

The undotted names are deliberate: DIA-NN 1.7 reads only those, and DIA-NN
2.x accepts them as well. DIA-NN reads neither `Protein.Group` nor the CCS;
they are there so that DIALibGen reads back exactly what it wrote.
DIALibGen still reads the dotted TSVs (`Modified.Sequence`, `Precursor.Mz`,
`Protein.Group` only, ...) that 0.11.1 and earlier wrote, and DIA-NN-style
TSVs without the two DIALibGen columns. The library records no gene names,
so there is no `Genes` column; for gene-level protein groups, search with
`--fasta` and `--reannotate`.

RT is written as stored. DIA-NN fits its own map from library RT to the
run's retention time, so the unit does not matter to it: a generated
library's RT is the model's normalised 0..1 output (or iRT with
`irt_rescale`), a tuned or refined library's is the reference run's
minutes. A DIA-NN 1.7.12 search of an Astral run with a generated library
found 84,076 precursors with RT scaled to 0..100 and 84,161 with it left
at 0..1.

`-mode convert -in lib.parquet -out lib.tsv` writes an existing library in
the other format, content unchanged, and works in either direction; with
`-out_format` it writes any of the layouts above.

### The DIA-NN 2.0 Parquet library

`-out_format diann2_parquet` (with `-out lib.parquet`) writes the layout DIA-NN
2.0 writes itself with `--out-lib lib.parquet`: one row per transition, these
27 columns in this order, none nullable:

| column | type | content |
|---|---|---|
| `Precursor.Id` | string | `<sequence><charge>`, `_decoy` appended for a decoy |
| `Modified.Sequence`, `Stripped.Sequence` | string | UniMod notation as in the TSV; residues only |
| `Precursor.Charge`, `Proteotypic`, `Decoy` | int64 | as in the TSV |
| `N.Term`, `C.Term` | int64 | 0: the library does not record protein-terminal peptides |
| `RT`, `IM` | float | as stored; 1/K0 0 when absent |
| `Q.Value`, `Peptidoform.Q.Value`, `PG.Q.Value` | float | 0, as in DIA-NN's predicted libraries |
| `PTM.Site.Confidence` | float | 1 |
| `Precursor.Mz`, `Product.Mz`, `Relative.Intensity` | float | |
| `Fragment.Type`, `Fragment.Charge`, `Fragment.Series.Number`, `Fragment.Loss.Type` | string, int64, int64, string | |
| `Exclude.From.Quant`, `Flags` | int64 | 0 |
| `Protein.Ids`, `Protein.Group` | string | accessions, `P1;P2` |
| `Protein.Names` | string | entry names |
| `Genes` | string | empty: the library records no gene names |

DIA-NN ignores the file's metadata, so the generation recipe travels in it as
in DIALibGen's own layout (`odia.config_json`); the file is not
cache-eligible. Checked with DIA-NN 2.0: a 3,000-protein library loads with
the same 2,997 protein isoforms, 3,379 protein groups and 503,945 precursors
from the TSV and from this file, and a search of an Astral run's first seven
minutes finds the same 350 precursors and 120 protein groups at 1 % FDR, with
identical reports. DIALibGen reads the file back (`-in`), with the protein
group as accessions. DIA-NN makes its own decoys; generate for DIA-NN with
`decoys` `none` (the default).

### OpenSWATH

`-out library.pqp` writes OpenSWATH's PQP, `-out library.tsv -out_format
openswath_tsv` its assay TSV. Both are written by OpenMS
(`TransitionPQPFile`, `TransitionTSVFile`) from an OpenMS
`TargetedExperiment` built from the library: one precursor per library
precursor (`TransitionGroupId` `<n>_<sequence>_<charge>`), the protein group's
members as proteins with their accessions, RT as `NormalizedRetentionTime`
(`LIBRARY_RT`), 1/K0 as `PrecursorIonMobility` (`LIBRARY_DRIFT_TIME`), every
fragment as a detecting and quantifying transition. No column mapping is
needed any more:

```bash
DIALibGen -in proteins.fasta -out library.pqp -generation:irt_rescale true
OpenSwathDecoyGenerator -in library.pqp -out library_decoys.pqp
OpenSwathWorkflow -in run.mzML -tr library_decoys.pqp -out_features run.osw
```

Decoys: the export writes what the library holds. Generate without decoys
(the default) and let `OpenSwathDecoyGenerator` make them. A DIALibGen decoy
carries its target's sequence with shifted fragment m/z (DIA-NN's
convention); OpenMS keys PQP peptides by sequence, so the target's peptide and
protein rows would be marked as decoys too, and DIALibGen warns when a library
with decoys is exported.

Checked with OpenMS 3.5: a 1,000-protein library (105,791 precursors,
1,235,554 transitions) round-trips through `TargetedFileConverter` (PQP to
TSV and TSV to PQP) with every precursor, transition, peptide and protein,
RT and 1/K0; `OpenSwathDecoyGenerator` and `OpenSwathWorkflow` run on the
exported PQP. On an Astral run's first seven minutes, the extracted fragment
traces of 62 of 96 precursors that DIA-NN identified there peak within 15 s
of DIA-NN's retention time (median 1.7 s). A calibrated, scored OpenSWATH
search (iRT, PyProphet) has not been benchmarked: on that uncalibrated slice
PyProphet reported no identifications at 1 % FDR. Use `-generation:irt_rescale true` so that the RT is on the iRT
scale, and give `OpenSwathWorkflow` iRT assays for your sample (`-tr_irt`) or
its `-auto_irt` sampling. OpenMS's writers take the whole library in memory
as a `TargetedExperiment`: measured, about 0.8 kB per transition for the PQP
and 2 kB for the OpenSWATH TSV, on top of the library. A whole human proteome
with charges 1-4 (about 70 M transitions) needs roughly 60 GB for the PQP;
export it on a machine that has the memory, or export a FASTA subset. The
PQP is deterministic: `-mode convert` of a generated Parquet writes the same
bytes as generating the PQP directly.

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
  one thread per physical core or an explicit count. `-machine:threads`
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
