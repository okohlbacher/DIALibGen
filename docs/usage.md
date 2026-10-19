# Usage and migration to 0.12

`DIALibGen --help` lists the main options (23), `DIALibGen --help -mode <mode>`
one mode's, `DIALibGen --helphelp` all of them. The
[parameter reference](parameters.md) is generated from the TOPP INI schema and
marks which options `--help` shows. An unknown option prints one line with the
closest name (`Unknown option '-instrument' (did you mean
-generation:instrument?)`) and exit code 6, not the whole help.

## Choose a mode

| Task | Command |
|---|---|
| Predict from FASTA | `DIALibGen -in proteins.fasta -out predicted.tsv` |
| Apply observed values | `DIALibGen -mode refine -in predicted.tsv -ids report.parquet -out refined.tsv` |
| Learn RT/CCS and predict the full library | `DIALibGen -mode tune -in predicted.tsv -ids report.parquet -out tuned.tsv` |
| Tune, then apply observed values, in one call | `DIALibGen -mode refine -tune -no_filter -in predicted.tsv -ids report.parquet -out adapted.tsv` |
| Change the format | `DIALibGen -in tuned.parquet -out tuned.tsv` |
| Tune with `-run` instead of `-ids` (**experimental**) | `DIALibGen -mode tune -in predicted.tsv -run run.mzML -out tuned.tsv` |

The output format follows the extension. Write `.tsv` for a file a search
engine reads; DIA-NN 2 does not read DIALibGen's `.parquet`. Write `.parquet`
for a record that carries its own recipe, or for a library that only a later
DIALibGen step reads. The examples below use `.parquet` in that second sense.

Without `-mode` (a stock INI's `mode = generate` does not count), a FASTA as
`-in` generates and a library (`.parquet`, `.tsv`, unless it starts with `>`)
converts; a library with `-ids` or `-run` is refused with "give -mode refine
or -mode tune". With `-write_config` the mode stays `generate`. The
[desktop app](../gui/README.md) exposes generate, refine and tune, including
RT/CCS fine-tuning controls and optional model export.

Non-default settings belonging to a different mode are rejected. A complete
default INI remains usable, but generation settings do not configure refinement
or tuning, and refinement flags do not configure generation. A setting that has
no effect in the run (`-irt_standards` without `-generation:irt_rescale true`,
`-generation:recompute_decoy_mz true` with `-generation:decoys none`,
`-generation:max_variable_modifications` without variable modifications,
`-generation:derive_ion_mobility false` without CCS, the `-intensity_*` options
or `-allow_mixed_intensity` without `-write_intensity`, the `-im_*` options
without `-write_im`, `-machine:coreml_units` off CoreML) prints a `WARNING`
that names it; from 0.13 these are refused (`-irt_standards` and `machine:*`
only ever warn: an INI written on a Mac must still run on Linux).

## Generation settings

Every editable generation setting has a `-generation:<key>` TOPP option. For example:

```bash
# The report's Astral library: the instrument implies the rest (instrument defaults, below).
DIALibGen -in proteins.fasta -out library.parquet \
  -generation:instrument Astral -generation:variable_modifications 'Oxidation (M)'

# Everything spelled out (0.11's form; the same library)
DIALibGen -in proteins.fasta -out library.parquet \
  -generation:instrument Astral -generation:nce 25 -generation:precursor_mz 380 980 \
  -generation:fragment_mz 150 2000 -generation:precursor_charges 1 2 3 4 \
  -generation:peptide_length 7 30 -generation:fixed_modifications 'Carbamidomethyl (C)' \
  -generation:variable_modifications 'Oxidation (M)' -generation:max_variable_modifications 1 \
  -generation:decoys none -generation:enzyme Trypsin/P -generation:missed_cleavages 1 \
  -generation:n_terminal_methionine_excision true
```

Quote each modification containing spaces. Boolean generation options take
`true`/`false`; numeric ranges take two values. In JSON, use `[]` to clear a
modification list, for example `{"fixed_modifications": []}` for a
non-alkylated preparation.

Variable terminal modifications apply at the default
`-generation:max_variable_modifications 1` (since 0.11.2; `Acetyl (N-term)`,
`Amidated (C-term)` and `Acetyl (Protein N-term)` included). The log counts the
precursors of each variable modification. Raising the limit allows two
variable modifications per peptide, so the library grows.

```bash
DIALibGen -write_config defaults.json
DIALibGen -config defaults.json -generation:nce 35 -write_config effective.json
DIALibGen -write_ini workflow.ini
```

For generation, precedence is explicit command-line settings, then supplied
INI settings, then JSON, then the instrument's defaults, then default.yaml
([instrument defaults](#instrument-defaults)). An INI written by
`-write_ini` materializes settings: those entries count as supplied when the
INI is loaded. Remove entries that should instead come from JSON. Unknown
JSON keys and invalid values are errors.

The default instrument is `QE` (default.yaml's `instrument`). Choose the
instrument used for acquisition; accepted aliases include `Astral`/`Orbitrap
Astral` (Lumos), `Exploris` (QE), `timsTOF Pro` (timsTOF), and `ZenoTOF`
(SciexTOF). Without `-generation:instrument` the log says so: `no
-generation:instrument: QE assumed, so CCS and 1/K0 are not predicted`.
`SciexTOF` and `ThermoTOF` have no trained weights in the shipped MS2 model and
print a warning: `ThermoTOF` predicts what `Lumos` predicts (median spectral
cosine 0.993 at the same NCE), `SciexTOF` does not (0.914, a tenth of
precursors below 0.64), and neither has been validated on SCIEX data. Omitting
NCE uses the instrument's default: 25 for Astral, Lumos and its aliases, 30
otherwise. The native option `-generation:nce -1` means automatic, and
generated INI files preserve that automatic setting. The effective JSON config
records the resolved NCE and its source; reusing that numeric value fixes NCE
until it is changed or the native option is set back to -1.

### Instrument defaults

What depends on the instrument lives in YAML files, one per instrument, and
`default.yaml` for everything an instrument file leaves out. The release ships
them in `share/DIALibGen/instruments/` (in the macOS app and the Windows and
Linux packages beside the models; a source tree has them in `data/instruments/`):

| file | instrument | what it sets |
|---|---|---|
| `default.yaml` | the default instrument (`QE`) and every other one | precursor m/z 350-1200, fragment m/z 200-1800, charges 1-4; per MS2 model slot the NCE (Lumos 25, others 30) and ion mobility (timsTOF only) |
| `astral.yaml` | Astral, Orbitrap Astral, Astral Zoom (MS2 model: Lumos) | precursor m/z 380-980, fragment m/z 150-2000, charges 1-4, NCE 25, no ion mobility |
| `timstof.yaml` | timsTOF and its Pro, SCP, HT, Ultra and Flex models | charges 2-4 (since 0.13; 0.12.1 had 1-4: `-generation:precursor_charges 1 2 3 4`), precursor m/z 350-1200, fragment m/z 200-1800, NCE 30, ion mobility with a 1/K0 correction (since 0.13; see [Ion mobility](#ion-mobility)) |

Each value carries its source in a comment. Precedence, per value: an
explicit option (command line, INI, `-config` JSON) > the instrument's file >
`default.yaml`. So `-generation:instrument Astral` alone gives the report's
Astral library, and `-generation:precursor_mz 350 1200` widens its window.
The log names the files, in one line, and the values an instrument file
implied:

```text
instrument defaults: astral.yaml (Astral) over default.yaml
Astral method defaults: precursor m/z 380-980, fragment m/z 150-2000, charges 1 2 3 4 (... to change them, e.g. -generation:precursor_mz 350 1200 for a wider window)
```

`-debug 1` adds the full paths and the directories searched. The recipe
(`-write_config`, the Parquet's `odia.config_json`) records the resolved values
with their source (`nce_source`, `precursor_mz_source`, `fragment_mz_source`,
`precursor_charges_source`: `instrument-default:Astral`, `default`, `TOPP` or
`config`) and the files under `instrument_defaults` (name and FNV-1a64 of
`default.yaml` and of the instrument's file). The resolved values are in the
library's fingerprint as before; the files themselves are not, so a library
built from the same values is the same library whichever file supplied them.

An instrument file has these keys (unknown keys, wrong types and values out of
range are errors that name the file and the key):

```yaml
name: Exploris 480                # required; matched ignoring case, ' ', '-', '_',
aliases: [E480, Orbitrap Exploris 480]  # a leading "Orbitrap" and a trailing model number
model_instrument: QE              # required: the MS2 model's slot (QE, Lumos, timsTOF, SciexTOF, ThermoTOF)
ion_mobility: false               # what -generation:predict_ion_mobility auto does
nce: 27                           # > 0, <= 100
nce_by_charge: {4: 32}            # NCE of some precursor charges, over nce (since 0.13)
precursor_mz: [400, 1000]         # [min, max]
fragment_mz: [200, 1800]
precursor_charges: [2, 3, 4]      # 1-8, no repeats
mobility_correction:              # ion-mobility instruments: see "Ion mobility" below
  ccs_model_fnv1a64: c26eb20995b991f5   # the CCS model the lines were fitted for
  intercept: 0.0237               # 1/K0 = intercept + slope x derived 1/K0, for charges
  slope: 0.9894                   #   without a line of their own
  by_charge:
    2: {intercept: 0.0370, slope: 0.9757}
```

`nce_by_charge` goes with the file's `nce`: the precursors of a listed charge
are predicted at its NCE, the others at `nce`. An explicit
`-generation:nce` (one value for every charge) replaces both. A `-config`
JSON may set `"nce_by_charge": {"2": 25, "4": 35}` itself, or `null` for one
NCE; the recipe records the resolved map and `nce_by_charge_source`, so a
library replays with the same energies.

**Collision energy on diaPASEF.** The best NCE for the timsTOF slot depends on
the method's collision-energy ramp, so `timstof.yaml` keeps 30. Median cosine
of all b/y ions at fragment charges 1-3 from the raw spectra, and DIA-NN 2.0
precursors at 1 % FDR (charges 2-4), against NCE 30:

| method | NCE 20 | NCE 25 | NCE 30 | best per charge |
|:--|--:|--:|--:|:--|
| PXD041421 (CE 26-46 eV), run 3197 | 0.875, +0.4 % IDs | 0.881, +1.0 % | **0.888** | 30 for 2+, 3+, 4+ |
| PXD041421, run 9898 | 0.870, −1.7 % | 0.877, −0.7 % | **0.885** | 30 for 2+, 3+, 4+ |
| PXD029836 (CE 20-57 eV), run 01 | **0.884**, +2.2 % | 0.881, +1.0 % | 0.868 | 20 for 2+ and 3+, 40 for 4+ |

For a method with a ramp like PXD029836's, copy `timstof.yaml`, set
`nce: 20` (or add `nce_by_charge: {4: 40}`) and pass it with
`-instrument_file`.

`default.yaml` has `instrument`, `precursor_mz`, `fragment_mz`,
`precursor_charges` (all required) and `models` (per MS2 model slot: `nce`,
`ion_mobility`). An instrument without a file of its own (Lumos, Exploris, QE
HF-X, ...) uses `default.yaml` with its model slot's NCE and ion mobility.

To add an instrument, or to change a shipped one, put a file in a directory of
your own; the first directory that has a file of a given name wins, so a user's
`astral.yaml` replaces the shipped one and a user's `default.yaml` the
defaults:

1. `$DIALIBGEN_INSTRUMENT_DIR` (a list, separated by `:`, on Windows `;`);
2. the user's directory: `~/.config/DIALibGen/instruments` (Linux, or
   `$XDG_CONFIG_HOME/DIALibGen/instruments`), `~/Library/Application
   Support/DIALibGen/instruments` (macOS), `%APPDATA%\DIALibGen\instruments`
   (Windows);
3. the shipped `share/DIALibGen/instruments`.

`-instrument_file my-instrument.yaml` reads one instrument's defaults from that
file instead; `-generation:instrument` may then be omitted, or must be one of
the file's names. Two files of one directory that claim the same name are an
error.

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
`-generation:irt_rescale true` for iRT calibration. It uses the bundled
standards unless you also pass `-irt_standards standards.tsv` (recorded in the
recipe as `irt_standards`); that option has no effect on its own and says so. CCS and derived 1/K0 are distinct
quantities; `-generation:derive_ion_mobility false` omits the conversion.

### Ion mobility

CCS and 1/K0 are predicted only for an instrument with ion mobility
(`ion_mobility: true` in its file, or for its MS2 model slot in default.yaml):
`timsTOF` and its aliases. For every other instrument, the default `QE`
included, the library carries neither (empty `IonMobility` and `CCS` in the
TSV, nulls in the Parquet, 0 in the DIA-NN 2.0 Parquet, no drift time in
OpenSWATH), and generation is faster by the CCS phase. The log says which, in
one line:

```
Astral has no ion mobility: CCS and 1/K0 not predicted (-generation:predict_ion_mobility true to predict them)
```

**1/K0 correction.** The 1/K0 that Mason-Schamp gives for the stock CCS
model's prediction (N2, 305 K; the constants of AlphaPeptDeep's own
conversion, alphabase `ccs_to_mobility_bruker`, within 5e-5) is low on
diaPASEF: by a median 0.014 1/K0 on the timsTOF benchmark (PXD029836), more
at low mobility and higher charge. `timstof.yaml` therefore carries a
`mobility_correction`: straight lines, per precursor charge, from the derived
to the observed 1/K0, fitted on two public timsTOF Pro runs of another study
(PXD041421, two instruments). On PXD029836, which the fit did not see, the
median |Δ1/K0| of held-out proteins falls from 0.0154 to 0.0104 and the bias
from −0.0137 to −0.0002 (DIA-NN 2.0's own prediction: 0.0102). The CCS column
is unchanged; only 1/K0 is corrected. The lines belong to the CCS model they
were fitted for (`ccs_model_fnv1a64`): a tuned CCS model, which learned the
run's own 1/K0, is not corrected, and the log says so.
`-generation:mobility_correction off` writes the uncorrected value (0.12.1's).
The recipe records `mobility_correction_applied` (the file and the lines).

`-generation:predict_ion_mobility` overrides the instrument: `auto` (default)
decides by it, `true` always predicts CCS, `false` never does. With CCS
predicted, `-generation:derive_ion_mobility` (default `true`) adds 1/K0;
without CCS it has no effect. Until 0.11.5 every library carried both. Apart
from these two columns and the recorded recipe (`predict_ion_mobility` in the
config, `im=none` in the fingerprint), an Orbitrap library is the one 0.11.5
writes, and a timsTOF library differs in the recorded
recipe and, since 0.13, in its default charges (2-4) and its corrected 1/K0
(both from timstof.yaml). DIA-NN
1.7.12 and 2.0 load such a TSV or DIA-NN 2.0 Parquet library (2.0 reports its
1/K0 as 0). `-tune_heads auto` (the default), `ccs` or `both` on such a
library with a diaPASEF report trains on the run's observed 1/K0 and then adds
CCS and 1/K0 to every precursor. A 0.11 recipe without `predict_ion_mobility`
follows this rule too: an Astral recipe from 0.11 now gives a library without
CCS.

## Speed, threads and memory

`-threads` defaults to one thread per **physical** core of the CPUs the process
may run on (its affinity mask and cgroup CPU quota on Linux) when neither the
command line nor the INI sets it; the log says `-threads 0: N threads (...)`.
Hardware threads of one core share its arithmetic units: on a 2 x 96-core node
96 sessions on all 384 hardware threads were 12 % slower than 48 on the
physical cores. An INI's value is respected, 1 included (`--help` shows TOPP's
registered default, 1, which DIALibGen does not apply); `-write_ini` writes
`threads 0`. Under SLURM without CPU binding or a cgroup the whole node counts:
pass `-threads $SLURM_CPUS_PER_TASK`.

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

Measured on 2 x 96-core EPYC 9654 for a 5.8 M-precursor library (Parquet),
one traced run each on the same node:

| version | 16 threads | 96 threads |
|---|---|---|
| 0.11.2 | 593.6 s | 157.9 s |
| 0.11.3 | 682.8 s | 146.2 s |

0.11.2 ran a pruned MS2 model, which 0.11.3 removed; without it MS2 needs
about 15-19 % more CPU, so 0.11.3 is 15 % slower at 16 threads and 7.5 %
faster at 96. In 0.11.2's measurements peak memory was 4.2 GiB at 8 threads,
4.7 GiB at 16, 9.5 GiB at 96 and 15.1 GiB at `-threads 0` (192 physical
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
q-value gates: only the precursor gate is on by default (`-q_precursor 0.01`);
the global (peptide) and protein gates (`-q_global`, `-q_protein`) are off
(1) since 0.12, as every recorded refine run and tuning's training set had it.
An enabled gate requires its corresponding column; missing columns are errors. Reports containing multiple runs must be
split before use. When only some matched precursors have observed RT, RT writing
requires whole-library RT tuning first; otherwise predictions and observed minutes
could be mixed. Use `-no_write_rt` to retain the original values.

Useful options:

- `-write_im`: also use observed 1/K0, subject to charge and ramp limits.
- `-no_filter`: retain unidentified precursors. Pair it with `-no_write_rt` to
  keep library RT values, or with `-tune` (any heads but `ccs`) to re-predict
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
DIALibGen -mode tune -in predicted.parquet -ids report.parquet -out tuned.parquet
DIALibGen -mode tune -in predicted.parquet -ids report.parquet \
  -out tuned.parquet -tune_out_models tuned-models \
  -train:lr 0.0001 -machine:train_threads 16
```

Tuning adapts RT and CCS with LibTorch and re-predicts the complete input
library with ONNX Runtime. It preserves precursor keys, m/z and fragment
intensities; it does not apply observed-value refinement to that output.
Train on a representative reference run and evaluate transfer on separate
runs. Tuning needs at least 100 units in each of the protein-held-out
validation and test cohorts, in practice several hundred identified precursors;
smaller reports are refused.

The bundled models are the default starting point. `-models DIR` selects
another model directory (for generate, tune and `-run`; `-tune_models` is its
deprecated name). `-tune_heads auto`, the default, tunes RT and adds CCS when
the report has the ion mobility for it: the same data selection and cohort
check as training decide, before any epoch. A report without 1/K0 (Astral,
Exploris) tunes RT only, with one line saying so; one with too few usable 1/K0
values drops CCS with a `WARNING`. `rt`, `ccs` or `both` choose explicitly
(`both` fails on a report without ion mobility). The provenance records the
requested and the tuned heads (`tune.heads`). `-tune_out_models DIR` retains ONNX models and training
sidecars instead of using temporary model files. Tuned ONNX files also embed
provenance metadata.

Training options are grouped under `filter:`, `cohort:`, `train:`, `stop:` and
`machine:` (`--helphelp`). `-write_config` in tune (or refine with `-tune`)
writes them as the recipe's `tuning` object, which `-config` reads back; options
on the command line or in the INI override it. Protein-level cohorts and validation govern checkpoint selection;
the tool reports rejection and fallback decisions. Options such as
`-cohort:full_fit` or `-cohort:no_inner_val` change the hold-out contract and
must be reported when presenting accuracy results. CPU is the portable
release configuration on Linux and Windows; on Apple silicon the release
build trains on the GPU by default (`-machine:device auto` resolves to `mps`
when LibTorch has Metal Performance Shaders, else to `cpu`; `-machine:device
cpu` trains on the CPU; see [Apple GPU](#apple-gpu)). Source builds can
request `-machine:device cuda:0` when linked to a compatible CUDA LibTorch
installation (see [building](building.md#cpu-and-cuda)). The re-prediction
follows `-machine:device` as generation does: CUDA when `auto` finds a usable
GPU with cuDNN 9 (or `cuda` forces it), else the CPU (`-tune_predict_gpu` is
deprecated). On one H100, tuning the RT head on an Astral
run (77,029 sequences, 59 epochs) and re-predicting a 3.5 M-precursor library
took 62 s, against 12.5 min on 64 CPU threads, with the same held-out accuracy
(TEST calibrated SD 0.350 min either way). CUDA training is not bit-identical to
CPU training, and two CUDA runs need not agree bit for bit either.

`-machine:train_threads` (deprecated name: `-machine:threads`) is the number of
CPU threads LibTorch trains with. Each
training step is a chain of small operations, so the count stops paying
early: for the RT head on an Astral run (77,029 sequences) on a 2 x 96-core
EPYC 9654, an epoch took 11.3 s at 4 threads, 6.7 s at 8, 4.4 s at 16, 3.7 s
at 32 and 3.5 s at 64, for 380, 450, 600, 990 and 1,850 CPU-seconds of
training over 8 epochs. On Linux the tuned weights depend on the thread
count: 4, 16 and 32 threads gave slightly different models (tuned RT a median
of 0.01-0.02 min apart, at most 0.26 min), and the same count gave the same model
bit for bit. On a Mac the count does not change the model (4, 8 and 16
threads gave identical weights). `.tune.json` and the tuned ONNX metadata
record the count as `threads`. The default (`0`) is a fixed number per platform, not
derived from `-threads` or the machine: **16 on Linux and Windows, 8 on
macOS**. On an M4 Max (12 performance and 4 efficiency cores) every operator
forks LibTorch's 16-thread OpenMP team and the join waits for the efficiency
cores' share: 12.4 s per epoch at 16 threads against 8.1 s at 8 (and 4), with
2,299 s of system time in an 847-s tune. A Mac and a Linux machine do not
train identical models at any count; Windows has not been measured. Up to
0.11.2 the default was 4: on Linux, tuned models made with default settings
differ from 0.11.2's, and `-machine:train_threads 4` reproduces them. On Apple silicon
the default device is the GPU anyway (`-machine:device auto` = `mps`, below),
where the count does not apply.

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

## Built-in identification (experimental)

`-run run.mzML` (or an mzPeak archive, `-run run.mzpeak`) replaces `-ids` in `refine` and `tune`: DIALibGen identifies
precursors in the run itself (fragment-evidence candidates, OpenSWATH
extraction, a semi-supervised discriminant, target-decoy FDR at 1 %), writes
the identifications as a report (`-out_ids`, default `<out>.ids.parquet`) and
continues exactly as with `-ids` on that report. Its settings are the
`-search:` options ([parameters](parameters.md)).

It is **experimental**. The entrapment (honesty) gate has not passed (held-out
test runs and screened entrapment are outstanding), and the resources gate,
16 GB and at most 60 minutes on 8 cores, is not met. `-ids` with the report of
an established search engine remains the validated route. `-write_intensity`
and `-min_fragments` are not available with `-run`; `-write_im` and the CCS
head need a diaPASEF run with a per-peak 1/K0 array. Gates, design and
measurements: [built-in identification](design/built-in-identification.md).

mzPeak runs (experimental as well): `-run` takes a `.mzpeak` archive or its
unpacked directory besides mzML (`.mzML`, `.mzML.gz`), chosen by the name. The
archive is read with mzpeak-openms (vendored, `third_party/mzpeak-openms`)
into the same OpenSWATH maps an mzML of the same data gives: every spectrum
through the same SwathFile consumers, in memory or in per-window cache files
(`search:readoptions`; `auto` decides by the decoded size, 3 GB). A diaPASEF
frame that the archive keeps whole (mzpeak-convert's default for a `.d`) is
split into one spectrum per isolation window by the windows' 1/K0 bands. RT
comes from the archive's spectrum time (double precision); an mzML export of
the same archive carries the float32 scan start time, so RTs differ from such
an mzML by up to about 1e-7 relative. Builds configured with
`-DDIALIBGEN_WITH_MZPEAK=OFF` refuse `.mzpeak` runs.

## GPU

With the default `-machine:device auto`, generation picks the fastest route
that works on the machine and says which, in one line, before predicting:

| machine | MS2 | RT and CCS |
|---|---|---|
| Linux/Windows, ONNX Runtime with CUDA, NVIDIA GPU visible, cuDNN 9 found | CUDA | CUDA |
| the same without cuDNN 9 (a `WARNING` with the fix) | CUDA | CPU, while MS2 runs |
| Mac (Apple silicon), library of at least 50,000 precursors per CPU thread | Apple GPU (CoreML) | CPU, while MS2 runs |
| anything else (no GPU, a CPU-only ONNX Runtime, a small library on a Mac) | CPU | CPU |

```text
device: prediction on CUDA (cuDNN 9.1.0)
device: MS2 on the Apple GPU (CoreML), RT and CCS on the CPU: 3.4 M precursors, CoreML pays from 800 k (-machine:device cpu to stay on the CPU)
device: prediction on the CPU: 13 k precursors, the Apple GPU pays only from 400 k (-machine:device coreml to use it anyway)
```

`-machine:device cpu`, `cuda` (`cuda:N` for device N), `coreml` or `mps` force
a route (until 0.11.5, `cpu` changed only the training device);
`-machine:device cpu` with `-machine:require_gpu` is refused. Each model
then prints the device it ended up on:

```text
[RT] provider CUDA device 0 sessions 4
```

The release packages for Linux (x64 and arm64, the tarball and the desktop
packages) ship ONNX Runtime's CUDA build (conda-forge `onnxruntime-cpp` 1.31.0
`cuda129`; on arm64 the SBSA build, for server ARM such as Grace: Jetson is not
supported). They do **not** ship CUDA or cuDNN: with an NVIDIA GPU, the CUDA 12
runtime libraries and cuDNN 9 installed, `auto` predicts on the GPU; without a
GPU it runs on the CPU without a warning; with a GPU but without CUDA 12 or
cuDNN 9 it falls back and the `WARNING` names what is missing. The arm64 GPU
route is built the same way but has not been run on aarch64 NVIDIA hardware.

Requirements on Linux: an NVIDIA driver for CUDA 12 or later, the CUDA 12
runtime libraries (cudart, cuBLAS, cuBLASLt) and **cuDNN 9 for CUDA 12**
(`libcudnn.so.9` and its `libcudnn_*.so.9` libraries; ONNX Runtime 1.30 and
1.31 need no unversioned `libcudnn.so` link: on an H100 RT and CCS ran on CUDA
with only the `.so.9` files). Before predicting, DIALibGen opens
`libcudnn.so.9` (else `libcudnn.so`) the way ONNX Runtime does (the loader's
search: `LD_LIBRARY_PATH`, the run path, the system cache, then ONNX
Runtime's own directory) and checks that it is cuDNN 9. A GPU counts as
visible with the NVIDIA driver loaded (`/proc/driver/nvidia/version`, or
WSL 2's `/dev/dxg`) unless `CUDA_VISIBLE_DEVICES` hides it. Without it only the MS2 model (a transformer)
runs on the GPU; RT and CCS (convolution and LSTM layers) run on the CPU, on
a thread of their own while the GPU predicts MS2, and a `WARNING` says why
and how to fix it, naming a directory where a cuDNN was found off the library
path (a conda environment, a pip `nvidia-cudnn-cu12` wheel, the CUDA toolkit
or a system directory: looked at, never loaded from). With RT and CCS on the
CPU, a whole-proteome library on one H100 took about 2.6 times as long (0.11.3,
phases in order).
Install cuDNN system-wide (`apt install libcudnn9-cuda-12 libcudnn9-dev-cuda-12`
from NVIDIA's CUDA repository) or in user space (`pip install nvidia-cudnn-cu12`
or conda `cudnn=9 cuda-version=12`), and put that directory on
`LD_LIBRARY_PATH`.

When a GPU exists but a model cannot use it, DIALibGen prints a `WARNING`
naming the model, the reason and the fix, and records it in Parquet output as
the `odia.warnings` metadata key (tune: `device_warnings` in the provenance).
A machine without a GPU, or a CPU-only ONNX Runtime, runs on the CPU without a
warning. `-machine:require_gpu` asks for the GPU for every model and turns
every CPU fallback, those cases included, into an error.

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
A whole human-sized proteome (5.8 M precursors) at `-threads 16`, same node
for both versions:

| version | 1 H100 | 2 H100 |
|---|---|---|
| 0.11.2 | 88.9 s | 65.2 s |
| 0.11.3 | 82.4 s | 47.3 s |

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

### Apple GPU

On Apple silicon, `-machine:device coreml` runs the MS2 model on the GPU
through ONNX Runtime's CoreML execution provider (the conda-forge and release
builds of ONNX Runtime for macOS arm64 include it); `-machine:device mps`
does the same for prediction and additionally trains on the GPU
(`-mode tune`, Metal Performance Shaders through LibTorch). Training on the
GPU is the default on Apple silicon: `-machine:device auto` (the default)
resolves to `mps` when LibTorch has the back end, and `-machine:device cpu`
is the way back. Prediction (`generate`) uses the GPU by default from 0.12
on, for a library large enough: loading the compiled shapes from the cache
costs 15-30 s (up to 60 s when macOS has to read them from disk again), and
for a small library the CPU is done sooner. Whole `generate` runs on an M4
Max (Astral settings, cache warm, one run each, 0.12 development build,
`benchmark/laptop/breakeven` in the project folder):

| precursors | threads | CPU | CoreML |
|---|---|---|---|
| 14 k | 1 | 6.8 s | 15.2 s |
| 111 k | 1 | 48.9 s | 24.6 s |
| 111 k | 4 | 15.6 s | 56.2 s (probably shapes read from disk again) |
| 482 k | 4 | 65.6 s | 31.7 s |
| 14 k | 16 | 1.4 s | 15.3 s |
| 111 k | 16 | 6.8 s | 39.6 s |
| 482 k | 16 | 33.7 s | 33.7 s |
| 1.64 M | 16 | 116.4 s | 67.9 s |
| 3.40 M | 16 | 251.1 s | 133.1 s |
| 6.96 M | 16 | 678.6 s | 274.6 s |

The break-even grows with the CPU threads, about 40 k precursors at 1
thread, 230 k at 4 and 480 k at 16. `-machine:device auto` therefore runs
MS2 on CoreML from 50,000 precursors per thread (`-threads 16` or `0`:
800 k; `-threads 4`: 200 k), a margin for a cold page cache
(`ODIA_COREML_MIN_PRECURSORS` sets the number), after the probe below, and
RT and CCS on the CPU without a warning. **A default library on Apple silicon then agrees with the CPU's to
float rounding, not byte for byte** (below; the `coreml_matches_cpu` test
holds RT, 1/K0 and CCS to 0 and intensities to 1e-4).
`-machine:device cpu` keeps the CPU library. The first CoreML run after a
restart or a long pause can take 40 s longer while macOS loads Core ML.

```text
device: MS2 on the Apple GPU (CoreML), RT and CCS on the CPU: 3.4 M precursors, CoreML pays from 800 k (-machine:device cpu to stay on the CPU)
[MS2] provider CoreML units CPUAndGPU concurrent 4
[RT] provider CPU device -1 sessions 8
```

With `-machine:device coreml` (or `mps`) every model asks for CoreML, whatever
the library's size, and RT and CCS fall back with a warning:

```text
WARNING: RT: runs on the CPU, not CoreML: CoreML has no operator for LSTM x2. Fix: none needed: ...
```

What runs where, and why:

- **MS2 on the GPU.** CoreML compiles a model per batch shape (rows x
  residues), so the predictor fixes both: every batch is padded to a
  compiled batch size (2048 rows, `PeptDeepPredictor::MAX_BATCH_ROWS`; from
  22 residues on also 512 and 1024, where the rarer lengths padded 15-38 %
  of their rows: on a 5.8 M-precursor library 8.8 % of the rows run are
  padding instead of 11.5 %, and those rows are the slow ones) and one
  session is compiled per residue count and batch size, with the whole
  graph on CoreML except the `OneHot` of the instrument index. The shapes
  the library needs (its residue counts, each with the batch sizes a chunk
  of that many precursors can take; a session holds its compiled model in
  memory, about 100 MB) are compiled, or loaded from the cache, in the
  background from the moment the predictor exists, in the order the chunks
  will ask for them, so the builds overlap the RT and CCS phases; a worker
  that needs a shape still being built waits for that one alone. Core ML compiles one model at a time whatever
  the number of builders (42 shapes for 7-30 residues: 105 s cold with 4, 8
  or 16 builders, about 2.5 s each), and cached loads gain a little from
  parallelism (23 s with 8 builders against 34 s one by one). (Before 0.11.4
  each shape was built by the first worker to need it while the others
  waited, one after the other, inside the MS2 phase: 66 s of a cold-cache
  whole-proteome run, 16 s warm.) (With dynamic shapes Core ML's compiler
  rejects the model's reshapes and the provider is slower than one CPU
  thread; the predictor also rewrites the models' scalar `PRelu` into
  `LeakyRelu`, which the provider's MLProgram builder gets wrong.) Core ML
  computes in float32: measured on `example/proteins.fasta`, intensities
  agree with the CPU's within 1.3e-6, no precursor selects a different
  fragment set, and two CoreML runs write the same bytes whatever
  `-machine:gpu_sessions` (here: concurrent batches; 2 gain about 10 % over
  1 on an M4 Max, 4 nothing more). On a whole human library (Astral
  settings, 3,488,457 precursors) the CoreML library has the CPU's
  precursors in the CPU's order, identical RT, mobility and CCS,
  intensities within 4.9e-6, 8 precursors with a different fragment set
  (ties at the 12-fragment cap) and a minimum spectrum cosine of 0.987
  against the CPU's (one precursor below 0.99, seven below 0.999) -- the
  same profile as CUDA with TF32 off.
- **The probe.** Before a model is compiled in the run itself, a child
  process (`DIALibGen --probe-coreml MODEL UNITS`) compiles its probe shape
  and runs one batch through it, once per process and model: a Core ML
  that aborts the process (a macOS virtual machine without a usable GPU,
  GitHub's macOS runners among them) then costs the child, not the run,
  and the model runs on the CPU with a `WARNING` that names the signal or
  exit code (`ODIA_COREML_PROBE=0` skips the child, `ODIA_COREML_PROBE_SECONDS`
  bounds it, 300). A shape Core ML refuses to compile inside the run falls
  back to a CPU session for its batches, also loudly. `-machine:require_gpu`
  turns either into an error. Measured on macOS 26 with an M4 Max; on macOS
  before 14.4 ONNX Runtime skips Core ML's optimisation hints and says so
  in the log (`ConfigureOptimizationHints`), which does not stop the model
  from running, and the probe decides the rest.
- **RT and CCS on the CPU, while the GPU predicts MS2.** CoreML has no LSTM
  operator, and these models are LSTMs: ONNX Runtime would keep the LSTM
  nodes on the CPU and hand the rest to CoreML, which is slower than the
  plain CPU path. Since MS2 does not depend on them, the RT phase (with the
  iRT calibration) and the CCS phase run on a thread of their own while the
  GPU predicts MS2, and the streaming writer waits for both. The GPU phase
  has a host side of its own (Core ML's threads, the encoding workers) and
  slows down when every core is busy, so the overlapped CPU phases run on
  half the cores (`ODIA_OVERLAP_RESERVE` sets how many they leave free).
  Measured on the Astral library on an M4 Max, untimed and under a
  changing background load: phases in order 238 s; overlapped on all 16
  cores 232 s (the GPU phase stretched from 153 to 229 s); on 12 cores
  208 s; on 8 cores 189 s. The CPU calls are the same calls whatever the
  session count, so the library is the same; `ODIA_OVERLAP=0` runs the
  phases in order. On CUDA all three models run on the GPU, and the phases
  stay in order: overlapped there (0.11.4), they contended for the GPU and
  a whole-proteome run took 89.2 s instead of 79.5 s on one H100. Running
  the LSTMs on CoreML was tried and rejected (`scripts/coreml_unroll_lstm.py`
  unrolls them per sequence length into MatMul/Sigmoid/Tanh steps): the
  unrolled graphs agree with the CPU's LSTM kernel to 6e-7 (RT) and
  4.3e-4 A^2 (CCS) on CoreML, within the tolerances the MS2 path is held
  to, but the GPU runs them at the speed of one CPU thread (16 tokens:
  16.2 k rows/s against 18.9 k per thread; 32 tokens: 9.5 k against
  10.2 k) and the CPU path has one session per core. The predictor
  probes each model's node placement and, when an arithmetic operator stays
  on the CPU, uses the CPU path and says so: a `WARNING` naming the model
  and the operator, recorded in the output's provenance as `odia.warnings`
  (tune: `device_warnings`). `-machine:require_gpu` turns that into an
  error, so on a Mac it refuses every library (RT always falls back); since 0.12
it does so before anything is predicted.
- **Compute units.** `-machine:coreml_units gpu` (CPU and GPU, the default)
  was the fastest on an M4 Max: MS2 at 2048 x 14 tokens ran 15-18 times
  faster than one CPU thread; `all` lets Core ML use the Neural Engine
  (about the same), `ane` (CPU and Neural Engine) was 4-5 times slower
  than the GPU for these batches, `cpu` is Core ML's own CPU path (for
  checks).
- **Compiled-model cache.** Compiling one shape takes 3-4 s and loading it
  back 0.8 s, so compiled models are kept under
  `~/Library/Caches/DIALibGen/coreml/<model>-<content hash>/<rows>x<tokens>`
  (about 65 MB per MS2 shape; 42 shapes, 2.7 GB, for 7-30 residues with the
  finer buckets). `DIALibGen -mode prepare` compiles the bundled MS2 model's
  shapes for `-generation:peptide_length` (default 7-30) into the cache
  once, for an installation or before a first run, and says where they are;
  without it the first `coreml` run on a proteome compiles them in the
  background while RT and CCS run. `ODIA_COREML_CACHE=<dir>` moves the
  cache, `ODIA_COREML_CACHE=0` keeps nothing, and `ODIA_COREML_CACHE_MAX`
  caps it (4G by default; bytes, or with K, M or G; `0`: no cap): above the
  cap the least recently used shapes of *other* models are removed when a
  run's predictor goes, never the running model's own. The cache is keyed
  by the model's content, so a tuned model never meets a stock model's
  compiled graph. Each run ends its MS2 phase with a line such as
  `[MS2] CoreML: 42 shapes built in 35.1 s (summed), 5808438 rows run,
  563914 padding rows (8.8 %)`.
- **Training on MPS** (`-mode tune`, the default `-machine:device auto` on
  Apple silicon, or `mps` by name): `auto` first allocates one small tensor
  on MPS and sums it; a LibTorch that reports MPS but cannot use it (a
  virtual machine) then trains on the CPU, with a log line and the reason
  under `device_note` in the provenance (`device_requested: auto`,
  `device: cpu`), while `-machine:device mps` by name fails loudly. The RT
  and CCS heads train on the GPU with no operator falling back to the CPU
  in LibTorch 2.10 (the log would say `not currently supported on the MPS
  backend`; set `PYTORCH_ENABLE_MPS_FALLBACK=1` in the environment **before
  starting DIALibGen** if a future LibTorch needs it, LibTorch reads it at
  load time). Measured on a 3,200-precursor fixture, 30 epochs: 18.6 s of
  training on MPS against 55 s on 16 CPU threads. MPS training is not
  bit-identical to CPU training (the provenance says so under
  `reproducibility`); the tuned weights differed by an L2 norm of 0.001 of a
  total change of 8.62 in that run. Two MPS runs on the same Mac differ from
  each other as well (tuned RT a median of 0.004 min apart, at most 0.37 min). The re-prediction after tuning runs on
  the CPU on a Mac (CoreML has no LSTM operator for RT and CCS).

## Messages and diagnostics

The log says what was done and, for anything that went wrong, the cause and
the fix, in a line or two. Warnings start with `WARNING:` (also those of the
TOPP log: `-log`). The `[RT]`, `[MS2]` and `[CCS]` provider lines keep their
form (`[MS2] provider CUDA device 0 sessions 4`) for scripts.

`-debug 1` (or the environment variable `ODIA_VERBOSE=1`, for the library
used without the tool) adds the details for diagnosis: the device checks
(ONNX Runtime's providers, whether an NVIDIA GPU is visible, the cuDNN
loader's result and the directories looked in), ONNX Runtime's own warnings, the CoreML node placement and per-run shape summary,
the directories searched for the models and the iRT standards, and the iRT
line's coefficients.

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
| `IonMobility` | predicted or observed 1/K0; empty without ion mobility | ion mobility (2.x; 1.7 has none) |
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

The training seed is reproducible within the same build and device setup
(on Linux, also the same `-machine:train_threads`; not on a GPU);
numeric results can differ across platforms, runtimes and CUDA kernels. The
wall-clock training budget is checked at epoch boundaries and includes report
loading and validation; final evaluation and export can exceed that budget.

Keep the provenance and any training sidecars with the library. A successful
write does not establish that a library improves a downstream search; the
[historical benchmark](benchmark.md) gives the evidence available so far.

## Using a library with a search engine

### DIA-NN

Write the library as `.tsv` ([the TSV library](#the-tsv-library)) or with
`-out_format diann2_parquet` ([DIA-NN 2.0 Parquet](#the-dia-nn-20-parquet-library)).

```bash
diann --f run01.d --lib predicted.tsv --fasta proteins.fasta --met-excision \
      --out first_pass.parquet --threads 16
```

- DIA-NN reads the protein groups from the TSV; `--fasta` with `--reannotate`
  adds gene-level grouping, which the library does not record.
- `--met-excision` matches DIALibGen's default removal of the initiator
  methionine (`-generation:n_terminal_methionine_excision true`), which DIA-NN
  does not do by default. A non-tryptic `-generation:enzyme` additionally needs
  the matching `--cut`.
- `tune` and `refine` read the report of **one** run and need its protein
  groups: `tune` holds its validation and test cohorts out by protein, and
  `refine` applies precursor, global and protein q-value gates of 0.01. If no
  observation passes, for example because a small library gives DIA-NN too few
  proteins for protein-level FDR, `refine` stops with "no reference observation
  passed the gates" instead of writing an empty library.
- Pass one `--f` per run; DIA-NN does not take a shell glob after a single
  `--f`. For timsTOF data DIA-NN advises `--mass-acc 15 --mass-acc-ms1 15`.
- Libraries are generated without decoys by default; DIA-NN makes its own.

For OpenSWATH, write the PQP directly ([OpenSWATH](#openswath)).

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

From 0.11 to 0.12 (0.12.1 is the first release of this line):

| 0.11 | 0.12 |
|---|---|
| `-threads` 1 unless given | all physical cores unless the command line or the INI sets it (an INI's 1 is respected) |
| `-tune_heads both` | `auto`: RT, plus CCS when the report has ion mobility |
| refine gates `-q_precursor`, `-q_global`, `-q_protein` at 0.01 | `-q_precursor 0.01` only; `-q_global 1 -q_protein 1` are the defaults |
| m/z windows and charges 350-1200 / 200-1800 / 1-4 for every instrument | from the [instrument files](#instrument-defaults): Astral 380-980 / 150-2000 / 1-4; timsTOF as before (350-1200 / 200-1800 / 1-4) |
| CCS and 1/K0 for every instrument | only with ion mobility (timsTOF) |
| `-mode generate` unless given | a library as `-in` converts |
| every option in `--help` | 23 in `--help`, the rest in `--helphelp` |
| unknown option: the error and the whole help | one line with the closest name |
| `convert` to Parquet drops the recipe | keeps the fingerprint and the recipe |

Deprecated in 0.12, removed in 0.13 (accepted with a `WARNING`; giving the old
and the new name with different values is an error; `-write_ini` writes only
the new names, and moves an old INI's values onto them): `-tune_models` (use
`-models`), `-machine:threads` (`-machine:train_threads`), `-tune_predict_gpu`
(re-prediction follows `-machine:device`), `-search:passes` (1 is its only
value) and `-search:allow_bootstrap` (`ODIA_SEARCH_ALLOW_BOOTSTRAP=1`).

An INI written by 0.11 loads in every mode: its values apply as given, and
because `-write_ini` wrote every default, it keeps the 0.11 defaults it lists
(`-tune_heads both`, the 0.01 gates, the 350-1200 window, `-threads 1`); the log
names them in one line, and deleting those lines gives the 0.12 defaults. The
mode check accepts them where they belong to another mode.

- Existing `DIALibGen -in FASTA -config JSON -out FILE` commands remain valid.
  Generation now also exposes those settings directly to TOPP workflows.
- Replace `DIALibRefine ...` with `DIALibGen -mode refine ...` for observed-value
  reconstruction. Refinement option names are retained.
- Replace `DIALibRefine -tune ...` with `DIALibGen -mode tune ...` for whole-library
  adaptation. Tune mode does not first filter or overwrite observations; use a
  separate refine invocation when observed-value reconstruction is intended.
- 0.11: `-threads` defaulted to 1, consistently with TOPP (0.12: all physical
  cores). On macOS, `-machine:cores performance` makes `-threads 0` count the
  performance cores only: on an M4 Max 12 threads instead of 16 needed 18 %
  less CPU time for the same library at -3 % to +16 % wall, so it saves
  energy, not time. `-machine:train_threads` controls training separately.
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
