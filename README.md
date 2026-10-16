# DIALibGen

Generate, refine and tune DIA spectral libraries with one TOPP-compatible
executable. DIALibGen uses OpenMS and AlphaPeptDeep in C++; it needs no Python
runtime.

## Overview

A DIA search needs a spectral library: for every peptide precursor, the
retention time, fragment intensities and ion mobility to look for. DIALibGen
predicts such a library directly from a protein FASTA with the open
AlphaPeptDeep models, as a scriptable TOPP tool that records how each library
was made. It fills the role of a search engine's built-in predictor, outside
the search engine. DIALibGen can then adapt the library to one of your runs.
It finds confidently identified peptides in the run itself, by targeted trace
extraction and target-decoy error control (`-run`, experimental), or takes them
from an external search's report (`-ids`). `refine` writes the observed values
into the library. `tune` retrains the RT and CCS predictors on the run and
re-predicts the entire library, including the peptides that have not been
identified.

```mermaid
flowchart LR
    fasta(["Protein FASTA"]) --> gen["<b>generate</b><br/>digest, then predict RT,<br/>fragments, CCS and 1/K0"]
    gen --> lib[("predicted<br/>library")]
    lib --> ident["<b>identify</b><br/>targeted trace extraction,<br/>scoring, target-decoy FDR"]
    run(["DIA run<br/>(mzML or mzPeak)"]) -->|"-run"| ident
    ident --> ids[("identifications<br/>of one run")]
    rep(["external search<br/>report"]) -->|"-ids"| ids
    lib --> tune["<b>tune</b><br/>retrain RT/CCS,<br/>re-predict every precursor"]
    ids --> tune
    lib --> refine["<b>refine</b><br/>keep identified precursors,<br/>write observed RT and 1/K0"]
    ids --> refine
    tune --> tl[("tuned library<br/>+ tuned models")]
    refine --> rl[("refined<br/>library")]
    tl --> search{{"DIA search<br/>DIA-NN, OpenSWATH"}}
    rl --> search
    classDef step fill:#dbeafe,stroke:#1d4ed8,color:#0b1b3f
    class gen,ident,tune,refine step
```

Blue boxes are the steps of the one DIALibGen executable; the rounded boxes are
your files and the hexagon is the search engine that uses the result.
Cylinders are what the steps produce. `identify` is not a mode of its own: it
runs inside `refine` and `tune` when they are given the raw run
(`-run run.mzML`), and its report is kept (`-out_ids`). Given an external
search's report instead (`-ids report.parquet`), they skip it. `-run` is
**experimental** ([below](#built-in-identification-experimental)); `-ids` with
an external search's report is the validated route. `tune` and `refine` can be used on
their own or combined: `-mode refine -tune` retrains and re-predicts the
library first and then writes the observed values, in one invocation; add
`-no_filter` to keep the precursors that were not identified. Models saved
with `-tune_out_models` can predict new libraries with `generate`. A predicted
library can also be searched as it is.

| Mode | Input | Result |
|---|---|---|
| `generate` (default) | Protein FASTA | Predicted RT, fragment intensities, CCS and 1/K0 derived from it |
| `refine` | Library, and one DIA run (`-run`, experimental) or that run's identification report (`-ids`) | Library filtered to identified precursors, with observed RT and optional mobility/intensities |
| `tune` | Library, and one DIA run (`-run`, experimental) or that run's identification report (`-ids`) | RT/CCS models adapted to that run, then predictions for the complete input library |

Version **0.12.1** incorporates library refinement and model training from
DIALibRefine. See [migration and usage](docs/usage.md),
[parameter reference](docs/parameters.md), [changes](CHANGELOG.md) and the
[backlog](BACKLOG.md).

## Install

Download the CLI archive or desktop installer for your platform from
[Releases](https://github.com/okohlbacher/DIALibGen/releases). Keep the archive's
`bin`, `lib` and `share` directories together. The CLI is one executable with
bundled runtime libraries, OpenMS data and the three prediction models.
Each platform's `DIALibGen-sources-<platform>.json` manifest lists its corresponding
dependency sources, recipes and patches. Large source archives are split into
numbered parts; [third-party notices](THIRD-PARTY-NOTICES.md#corresponding-sources)
explain how to reconstruct them. Runtime inventories and license texts are
included with the binaries.

```bash
curl -fsSLO https://github.com/okohlbacher/DIALibGen/releases/latest/download/DIALibGen-macos-arm64.tar.gz
tar xzf DIALibGen-macos-arm64.tar.gz
./bin/DIALibGen --help
```

On macOS, the [Homebrew tap](https://github.com/okohlbacher/homebrew-dialibrarygenerator)
also provides the CLI and desktop app:

```bash
brew install --cask okohlbacher/dialibrarygenerator/dialibgen-cli
brew install --cask okohlbacher/dialibrarygenerator/dialibgen
```

Signed 0.11.0 candidate `8cd1726` took between one and two seconds for its
first CLI launch on fresh hosted Macs (archive 2.0 s, Homebrew 0.8 s). Earlier
releases showed much longer delays. These observations do not guarantee startup
time on every Mac; the hashes and conditions are recorded in the
[macOS startup measurements](docs/testing.md#macos-startup).

## Usage

The commands below need DIALibGen **0.11.0 or later**; earlier releases have
neither `-mode` nor the `-generation:` options. Run `DIALibGen --helphelp` for
every option of your executable. The [usage guide](docs/usage.md) explains each
mode in depth; the [parameter reference](docs/parameters.md) lists all options.

| Task | Mode | Command |
|---|---|---|
| Predict a library from a FASTA | `generate` | `DIALibGen -in proteins.fasta -out predicted.tsv` |
| Write observed values into a library | `refine` | `DIALibGen -mode refine -in predicted.tsv -ids report.parquet -out refined.tsv` |
| Retrain RT/CCS on a run, re-predict the whole library | `tune` | `DIALibGen -mode tune -in predicted.tsv -ids report.parquet -out tuned.tsv` |
| Tune, then write observed values, in one call | `refine` with `-tune` | `DIALibGen -mode refine -tune -no_filter -in predicted.tsv -ids report.parquet -out adapted.tsv` |
| Tune on a raw run, no external search (experimental) | `tune` with `-run` | `DIALibGen -mode tune -in predicted.tsv -run run.mzML -out tuned.tsv` |

Existing output files are never overwritten. Options that belong to a different
mode are rejected when changed from their defaults.

### Examples

Generate a library for a timsTOF method (the instrument's defaults: charges
1-4, CCS and 1/K0; all physical cores), and one for the report's Astral
method (precursor m/z 380-980, fragment m/z 150-2000; see
[instrument defaults](docs/usage.md#instrument-defaults)):

```bash
DIALibGen -in proteins.fasta -out predicted.tsv -generation:instrument timsTOF
DIALibGen -in proteins.fasta -out astral.tsv -generation:instrument Astral \
  -generation:variable_modifications 'Oxidation (M)'
```

Match the digestion, modifications, charge and m/z ranges to the sample
preparation and acquisition. Defaults are Trypsin/P, one missed cleavage,
peptide lengths 7–30, precursor charges 1–4, precursor m/z 350–1200, fixed
Carbamidomethyl (C) and no decoys; the search engine creates its own decoys.

```bash
DIALibGen -in proteins.fasta -out predicted.parquet \
  -generation:precursor_charges 2 3 \
  -generation:variable_modifications 'Oxidation (M)' \
  -generation:missed_cleavages 2 \
  -generation:irt_rescale true
```

Keep a recipe and reuse it. Explicit command-line settings override the file:

```bash
DIALibGen -write_config generation.json
DIALibGen -in proteins.fasta -config generation.json -out predicted.tsv
```

A saved recipe fixes the NCE it resolved. When you change the instrument on
top of it, also pass `-generation:nce -1` (instrument default) or a value.
`.parquet` preserves the effective recipe and input/model hashes in metadata.
`.tsv` is DIA-NN's TSV library format, which DIA-NN 1.7 and 2.x read
completely, protein groups included (see [the column list](docs/usage.md#the-tsv-library));
retain the config separately. `-mode convert` turns an existing `.parquet` into a `.tsv`.
`-out_format diann2_parquet` writes DIA-NN 2.0's own Parquet library layout, and
`-out library.pqp` OpenSWATH's PQP ([output layouts](docs/usage.md#files-and-provenance)). With
`irt_rescale=false`, predicted RT is the model's normalized output, not iRT.

A complete cycle with DIA-NN 2: predict, search once, adapt to the run, search
again. `tune` and `refine` need a report of **one** run with protein groups;
DIA-NN takes them from the TSV, and `--fasta ... --reannotate` adds gene-level
grouping (see [DIA-NN](#dia-nn)).

```bash
# 1. predict; DIA-NN reads the TSV form
DIALibGen -in proteins.fasta -out predicted.tsv -generation:instrument timsTOF

# 2. first-pass search of one representative run
diann --f run01.d --lib predicted.tsv --fasta proteins.fasta --reannotate --met-excision \
      --out first_pass.parquet --threads 16

# 3a. adapt the predictors to the run and re-predict every precursor ...
DIALibGen -mode tune -in predicted.tsv -ids first_pass.parquet \
  -out tuned.tsv -tune_out_models tuned-models -machine:train_threads 8

# 3b. ... or keep only what was identified, with its observed RT and 1/K0 ...
DIALibGen -mode refine -in predicted.tsv -ids first_pass.parquet \
  -out refined.tsv -write_im -out_report residuals.tsv

# 3c. ... or both in one call: tune, re-predict every precursor, then write the observed values
DIALibGen -mode refine -tune -no_filter -in predicted.tsv -ids first_pass.parquet \
  -out adapted.tsv -write_im -machine:train_threads 8

# 4. final search: one --f per run (DIA-NN does not take a shell glob after a single --f)
diann --f run01.d --f run02.d --f run03.d --lib tuned.tsv \
      --fasta proteins.fasta --reannotate --met-excision --out final.parquet --threads 16
```

Which one? `tune` keeps the whole library and changes only its predicted RT and
ion mobility, so the final search tests the same hypotheses as the first.
`refine` without `-no_filter` shrinks the library to what the first search
identified. That makes the second search faster, but it also changes what its
error-rate estimate means, and precursors missed in the first pass stay missed.
Neither is guaranteed to improve a search; check on your own data. `refine`
applies precursor, global and protein q-value gates of 0.01 and stops when
nothing passes, which happens with very small libraries. `tune` needs several
hundred identified precursors. See the [usage guide](docs/usage.md) for
interpretation.

Tuned and observed retention times are in the reference run's minutes and are
specific to its gradient and method. Evaluate the transfer before using such a
library on a different method. The provenance sidecar (`*.refine.json`) records
quality gates, input fingerprints, residuals and, for tuning, each head's
recipe, cohorts, seed and model hashes.

Models saved with `-tune_out_models` can predict a new library later:

```bash
DIALibGen -in other.fasta -out other_adapted.tsv -generation:instrument timsTOF \
  -generation:rt_model tuned-models/peptdeep_rt_dynamic.onnx \
  -generation:ccs_model tuned-models/peptdeep_ccs_dynamic.onnx \
  -generation:free_cysteine_rt_correction false
```

This library's RT is the tuned model's normalized output, not minutes; only
`-mode tune` rescales to the reference run.

### Built-in identification (experimental)

**`-run` is experimental.** Two of its validation gates are still open: the
entrapment (honesty) gate has not passed on the test runs and a screened
entrapment design, and the resources gate (16 GB and at most 60 minutes on 8
cores) is not met: a whole timsTOF run takes about an hour on 96 cores. Use
`-ids` with an established search engine's report for libraries you rely on;
the [design document](docs/design/built-in-identification.md) has the gates
and the measurements.

With `-run`, `refine` and `tune` need no external search. DIALibGen finds
candidate precursors whose predicted fragments co-occur in the run's spectra,
extracts their traces with OpenSWATH (stock OpenMS), scores them with a
semi-supervised discriminant and keeps those that pass target-decoy error
control at 1 %. Decoys are built in memory; target and decoy are each
predicted from their own sequence, so the two compete on equal terms.

```bash
DIALibGen -mode tune -in predicted.parquet -run run.mzML -out tuned.tsv -threads 8
DIALibGen -mode refine -tune -no_filter -in predicted.parquet -run run.mzML -out adapted.tsv
```

- The run is a centroided DIA run as mzML or as an mzPeak archive (`.mzpeak`,
  or its unpacked directory; also experimental), chosen by the extension. For
  timsTOF diaPASEF, convert the `.d` with `mzpeak-convert` to `.mzpeak` (its
  default archive keeps each frame whole; DIALibGen splits it by the windows'
  1/K0 bands) or to mzML with a per-peak 1/K0 array (`--to mzml`). An mzPeak
  run is read with the vendored mzpeak-openms reader into the same OpenSWATH
  maps an mzML gives.
- The identifications are written to `<out>.ids.parquet` (`-out_ids`), in the
  columns an `-ids` report has, and can be reused with `-ids`.
- A diaPASEF run is searched with its ion mobility: the library's 1/K0 is
  calibrated to the run, each precursor is extracted from its one diaPASEF
  window within a 1/K0 window (`-search:im_window`, automatic by default), and
  the report carries each identification's observed 1/K0, measured at its
  elution apex over the window's whole 1/K0 range (a value read inside the
  1/K0 window would be pulled towards the prediction). `-tune_heads
  ccs|both` and `-write_im` therefore work with such a run; on a run without
  ion mobility (or with `-search:im_window -1`) they are refused before
  anything is searched (the default `-tune_heads auto` tunes RT only there). The CCS head is then tuned on the precursors the
  1/K0 window let through, which favours those the library already
  predicted well.
- Measured on an Orbitrap Astral and a timsTOF run: most identifications are
  also found by DIA-NN (95-97 %), and entrapment estimates the error at
  0.7-0.9 % at nominal 1 %. That is not yet the honesty gate, which needs the
  held-out test runs and screened entrapment.

### Instruments

`-generation:instrument` and `-generation:nce` select the fragment-intensity
prediction. They do not change predicted RT or CCS. Names are case-insensitive;
spaces, hyphens, underscores, a leading `Orbitrap` and a trailing model number
are ignored, so `-generation:instrument 'Orbitrap Exploris 480'` and `'timsTOF
Pro 2'` work (quote names with spaces). Vendor names and other affixes are not
stripped; if your instrument's full name is refused, use the value in the
"Pass" column. A name that is not recognised is an error, not a fallback.

| Your instrument | Pass | Model label | Default NCE | Note |
|---|---|---|---|---|
| Thermo Q Exactive, Plus, HF, HF-X; Orbitrap Exploris | `QE` | QE | 30 | The default. |
| Bruker timsTOF Pro, Pro 2, SCP, HT, Ultra, flex | `timsTOF` | timsTOF | 30 | |
| Thermo Orbitrap Astral; Tribrids (Fusion, Lumos, Eclipse); Velos, Elite | `Astral` or `Lumos` | Lumos | 25 | The model's reference label; AlphaPeptDeep groups the Astral here. |
| SCIEX ZenoTOF, TripleTOF | `ZenoTOF` or `SciexTOF` | SciexTOF | 30 | Accepted with a warning; see below. |
| ThermoTOF | `ThermoTOF` | ThermoTOF | 30 | Accepted with a warning; predictions are close to `Lumos`. |

NCE here is an input of the fragment model, not an instrument setting.
Orbitrap methods state a normalized collision energy, which you can pass on.
timsTOF and SCIEX methods set collision energies in eV; do **not** enter an eV
value. Keep the default, or compare a few values (25 to 40) on your own data:
on the one timsTOF diaPASEF method we measured, NCE 40 agreed better with the
observed fragments than the default 30, with little effect on identifications.
The effective configuration records the NCE that was used and whether it came
from you or from the instrument default.

The warning printed for `SciexTOF` says that predictions will be close to
`Lumos`. They are not: at equal NCE the predicted spectra differ from `Lumos`
about as much as the `timsTOF` label's do (median spectral cosine 0.85, against
0.99 for `ThermoTOF`). The label's weights in the shipped model look untrained,
and which of `SciexTOF` and `Lumos` fits SCIEX data better has not been
measured. If fragment agreement matters, generate both and compare.

1/K0 is derived from predicted CCS for every instrument unless
`-generation:derive_ion_mobility false` is set.

### Input and output formats

| File | Format | Notes |
|---|---|---|
| `-in` for `generate` | Protein FASTA | Peptides with ambiguous residues are skipped and counted. |
| `-out *.tsv` | DIA-NN's TSV library, one row per transition | Read completely by DIA-NN 1.7 and 2.x, protein groups included. Carries no recipe; keep the config next to it. |
| `-out *.parquet` | DIALibGen's compact Parquet, one row per precursor with fragment lists | From `generate` it embeds recipe, FASTA and model hashes; from `refine`/`tune` the same provenance as the sidecar. Input for `refine`/`tune`. Not read by DIA-NN 2; `-out_format diann2_parquet` writes DIA-NN 2.0's own Parquet layout. |
| `-out *.pqp` | OpenSWATH's PQP, written by OpenMS | `-out_format openswath_tsv` writes OpenSWATH's assay TSV instead ([output layouts](docs/usage.md#files-and-provenance)). |
| `-in` for `refine`/`tune` | DIALibGen Parquet, DIA-NN TSV, DIA-NN 2.0 Parquet library, or long-format Parquet (one row per transition, the TSV's columns) | |
| `-run` (experimental) | One centroided DIA run, mzML or mzPeak (`.mzpeak`) | Instead of `-ids`; see [built-in identification](#built-in-identification-experimental). |
| `-ids` | DIA-NN `report.parquet` of **one** run (`report.tsv` is not read) | Needs precursor q-values and protein annotation. A multi-run report is refused; filter it to one `Run` first. `refine` also accepts a pre-filtered long-format Parquet library with `-empirical_library`. |
| `-out_report` | TSV | Gate counts and residual summary (n, mean, sd, p95), measured before observed values are written. |
| `*.refine.json` | JSON sidecar of `refine`/`tune` | Provenance; keep it with the library. |

The TSV's columns are listed in the [usage guide](docs/usage.md#the-tsv-library).
`RT` is the model's normalized output unless `-generation:irt_rescale true` is
set (iRT) or the library was tuned or refined (minutes of the reference run;
`-rt_unit minmax` rescales to 0–100, `-no_write_rt` keeps the input RT). Ion
mobility is 1/K0 in Vs/cm².

#### DIA-NN

Use the **TSV** output, or `-out_format diann2_parquet`. DIA-NN reads the
TSV's protein annotation; `--met-excision` matches DIALibGen's default removal
of the initiator methionine, and a non-tryptic `-generation:enzyme` needs the
matching `--cut`. `tune` and `refine` need the report of **one** run with
protein groups. For timsTOF data DIA-NN advises fixed mass accuracies
(`--mass-acc 15 --mass-acc-ms1 15`). The default library has no decoys; let
DIA-NN generate them. Details: [DIA-NN](docs/usage.md#dia-nn).

#### OpenSWATH

`-out library.pqp` writes OpenSWATH's PQP through OpenMS's own writer; no
column mapping is needed:

```bash
DIALibGen -in proteins.fasta -out library.pqp -generation:irt_rescale true
OpenSwathDecoyGenerator -in library.pqp -out library_decoys.pqp
```

Generate without decoys (the default) and let OpenSWATH create its own.
`tune` reads DIA-NN reports only; neither mode reads OpenSWATH result files.
See [OpenSWATH](docs/usage.md#openswath) for what was checked.

## TOPP workflows and models

```bash
DIALibGen --helphelp
DIALibGen -write_ini DIALibGen.ini
DIALibGen -ini DIALibGen.ini -in proteins.fasta -out predicted.tsv
mkdir -p ctd && DIALibGen -write_ctd ctd/
```

`-threads` defaults to one thread per **physical** core unless the command line
or the INI sets it (TOPP's own default, 1, which `--help` shows, is not
applied); generation runs one inference session per thread. Training has its
own `-machine:train_threads` setting (default 16 on Linux and Windows, 8 on
macOS). Explicit generation CLI/INI settings override the optional JSON config,
which overrides the [instrument defaults](docs/usage.md#instrument-defaults);
CLI settings override INI settings. `--help` lists the main options,
`--help -mode <mode>` one mode's, `--helphelp` all.

Release builds include `peptdeep_{rt,ms2,ccs}_dynamic.onnx`. Override them with
`-models DIR`, `DIALIBGEN_MODEL_DIR` or the generation model-path options.
For a source installation, `scripts/fetch-models.sh --dir models` downloads and
verifies the pinned models. See [building](docs/building.md) and
[third-party notices](THIRD-PARTY-NOTICES.md).

## Desktop app

The [desktop app](gui/README.md) provides generation, observed-value refinement
and RT/CCS fine-tuning. Select the mode, input files and output, then adjust
the settings supplied by the same executable. Fine-tuning includes training
filters, cohorts, optimizer and stopping settings, and optional saved ONNX
models. The app requires unused output paths and keeps each mode's settings
separate.

## Limitations

- Generation emits b/y fragments without neutral losses. It does not tune the
  fragment-intensity model.
- Refinement requires the columns needed by enabled quality gates. Empirical
  libraries and fragment-intensity write-in have additional input contracts;
  empirical references require long-format Parquet, and compact references are
  refused. DIA-NN 2.x fragment parsing has synthetic tests and a regression
  using one real DIA-NN 2.0 `--export-quant` report row; other versions and
  export settings remain unverified.
- `-run` (built-in identification) is experimental: its entrapment and
  resource validation gates are open. It reads centroided mzML and mzPeak.
- Prediction supports one modification per residue. Refinement normalizes known
  modification names to UniMod accessions and reports unresolved names.
- A library can contain at most 4,294,967,295 transitions. Larger inputs must
  be split; oversized operations are refused before stored offsets can overflow.
- Portable release training uses CPU. CUDA training is supported by source
  builds with compatible CUDA LibTorch and runtime libraries.
- GPU prediction needs an ONNX Runtime with the CUDA provider (the Linux x64 and
  arm64 release packages include one; CUDA and cuDNN are not shipped), the
  CUDA 12 runtime and cuDNN 9 for CUDA 12; without cuDNN, RT and CCS run on the
  CPU with a warning; without a GPU everything runs on the CPU without one
  ([GPU](docs/usage.md#gpu)). The arm64 GPU route is untested on hardware.
- CWL/JSON descriptor export needs OpenMS built with `ENABLE_TDL=ON`; CTD and
  INI export do not.
- Windows releases target x64. The macOS desktop app and Homebrew casks
  require macOS 14 or later.

See [validation and measured coverage](docs/testing.md) for what the tests
establish and what remains untested. Whether a tuned or refined library improves
a search is an empirical question; the [usage guide](docs/usage.md#files-and-provenance)
points to the evidence available so far.

## Roadmap

A one-step mode that generates, tunes and refines in a single invocation, with
a small set of simplified options, is planned but not implemented. Tuning and
refining in one call already exists (`-mode refine -tune`); the planned mode
adds generation and the simplified options. The design and other open items
are in the [backlog](BACKLOG.md).

## Citation and licence

Cite the components and methods used:

- AlphaPeptDeep — Zeng et al., *Nature Communications* 13, 7238 (2022).
- OpenMS — Röst et al., *Nature Methods* 13, 741–748 (2016).
- iRT standards, when used — Escher et al., *Proteomics* 12, 1111–1121 (2012).
- Observed-value library reconstruction — Charkow et al., *Reference-Based
  Library Construction Improves Performance in low-input diaPASEF Workflows*,
  bioRxiv, DOI 10.64898/2026.04.29.721088.

BSD-3-Clause. See [LICENSE](LICENSE) and
[THIRD-PARTY-NOTICES.md](THIRD-PARTY-NOTICES.md). Bundled dependencies retain
their own terms, including Intel OpenMP's separate Windows redistribution terms.
