# DIALibGen

Generate, refine and tune DIA spectral libraries with one TOPP-compatible
executable. DIALibGen uses OpenMS and AlphaPeptDeep in C++; it needs no Python
runtime.

| Mode | Input | Result |
|---|---|---|
| `generate` (default) | Protein FASTA | Predicted RT, fragment intensities, CCS and optional 1/K0 |
| `refine` | Library and one run's identification report | Library filtered to identified precursors, with observed RT and optional mobility/intensities |
| `tune` | Library and one run's identification report | RT/CCS models adapted to that run, then predictions for the complete input library |

Version **0.11.0** incorporates library refinement and model training from
DIALibRefine. See [migration and usage](docs/usage.md),
[parameter reference](docs/parameters.md), and [changes](CHANGELOG.md).
Final release checks are tracked in [BACKLOG.md](BACKLOG.md).

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

The signed 0.11.0 candidate measured about one second on its first CLI launch
on fresh hosted Macs, through both archive and Homebrew delivery. Earlier
releases showed much longer delays. These observations do not guarantee startup
time on every Mac; the candidate, hashes and conditions are recorded in
[BACKLOG.md](BACKLOG.md).

## Generate a library

```bash
DIALibGen -in proteins.fasta -out predicted.parquet \
  -generation:instrument timsTOF -generation:nce 30 -threads 4
```

Match the digestion, modifications, charge and m/z ranges to the preparation
and acquisition. Defaults include Trypsin/P, one missed cleavage, peptide
lengths 7–30, charges 1–4 and fixed Carbamidomethyl (C). No decoys are added by
default; the consuming search engine can create them.

Every generation setting is available through a native TOPP option, INI or
JSON. Boolean options take `true` or `false`; lists take separate arguments.

```bash
DIALibGen -in proteins.fasta -out predicted.parquet \
  -generation:precursor_charges 2 3 \
  -generation:variable_modifications 'Oxidation (M)' \
  -generation:irt_rescale true
DIALibGen -write_config generation.json
DIALibGen -in proteins.fasta -config generation.json -out predicted.tsv
```

`.parquet` preserves the effective recipe and input/model hashes in metadata.
`.tsv` exports DIA-NN's library dialect; retain the config separately. With
`irt_rescale=false`, predicted RT is the model's normalized output, not iRT.

## Refine or tune a library

Refinement writes observed values from the reference run. It filters to
confident identifications by default; `-write_im` also replaces mobility.
To retain unidentified precursors with `-no_filter`, also use `-no_write_rt`,
or tune the RT head so the whole library first shares the reference run's RT units.

```bash
DIALibGen -mode refine -in predicted.parquet -ids report.parquet \
  -out refined.parquet -write_im -out_report residuals.tsv
```

Tuning learns RT/CCS from the report and predicts the complete library,
including unseen precursors. It preserves precursor keys and does not write
observations into the output. Saved models can be reused in later generation.

```bash
DIALibGen -mode tune -in predicted.parquet -ids report.parquet \
  -out tuned.parquet -tune_out_models tuned-models
```

Use a report containing one run. Provenance records quality gates, input
fingerprints, residuals and, when tuning, each head's complete recipe, cohorts,
seed and model SHA-256 hashes, even when models are temporary. Observed RT and tuned RT are
specific to the reference run's gradient; evaluate transfer to a different
method before using them there. See [usage and interpretation](docs/usage.md).

## TOPP workflows and models

```bash
DIALibGen --helphelp
DIALibGen -write_ini DIALibGen.ini
DIALibGen -ini DIALibGen.ini
DIALibGen -write_ctd ctd/
```

`-threads` follows the TOPP default of **1**. Set `-threads 0` to request all
available CPU inference sessions, capped at 16 to limit memory. Training has
its own `-machine:threads` setting. Explicit generation CLI/INI settings
override the optional JSON config; CLI settings override INI settings.

Release builds include `peptdeep_{rt,ms2,ccs}_dynamic.onnx`. Override them with
`DIALIBGEN_MODEL_DIR`, the generation model-path options, or `-tune_models`.
For a source installation, `scripts/fetch-models.sh --dir models` downloads and
verifies the pinned models. See [building](docs/building.md) and
[third-party notices](THIRD-PARTY-NOTICES.md).

## Desktop app

The [desktop app](gui/README.md) provides the generation workflow: select a
FASTA, output and settings, then generate a library. Its form uses the
executable's generation defaults. Refinement and tuning are available through
the CLI. The app requires an unused output path and validates thread counts
before launch.

## Evidence and limitations

A completed K562 benchmark of **DIALibGen 0.10.1 / DIALibRefine 0.3.0-dev**
found tuned DIALibGen within the registered whole-run margins of tuned DIA-NN,
with 0.79% / 0.34% fewer precursor identifications on two exposed technical
replicates. This is descriptive evidence, not equivalence or independent
confirmation. The 0.11.0 digestion fix changes generated peptide space; those
numbers do not validate unrestricted 0.11.0 outputs. Details and qualifications
are in [the historical benchmark summary](docs/benchmark.md).

- Generation emits b/y fragments without neutral losses. It does not tune the
  fragment-intensity model.
- Refinement requires the columns needed by enabled quality gates. Empirical
  libraries and fragment-intensity write-in have additional input contracts;
  empirical references require long-format Parquet, and compact references are
  refused. DIA-NN 2.x fragment-column parsing has synthetic contract tests;
  interoperability with a real `--export-quant` fixture remains to be checked.
- Prediction supports one modification per residue. Refinement normalizes known
  modification names to UniMod accessions and reports unresolved names.
- A library can contain at most 4,294,967,295 transitions. Larger inputs must
  be split; oversized operations are refused before stored offsets can overflow.
- Portable release training uses CPU. CUDA training is supported by source
  builds with compatible CUDA LibTorch and runtime libraries.
- CWL/JSON descriptor export needs OpenMS built with `ENABLE_TDL=ON`; CTD and
  INI export do not.
- Windows releases target x64. Platform-specific validation and macOS startup
  measurements are tracked in [BACKLOG.md](BACKLOG.md).

See [validation and measured coverage](docs/testing.md) for what the tests
establish and what remains untested.

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
