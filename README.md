# DIALibGen

Build an in-silico DIA spectral library from a FASTA: enzymatic digest →
precursor enumeration → AlphaPeptDeep predictions of retention time, MS2
fragment intensities and collision cross-section → Parquet or DIA-NN-dialect
TSV.

C++, no Python runtime. A TOPP-compatible tool built on
[OpenMS](https://github.com/OpenMS/OpenMS), which is consumed as an installed,
read-only dependency and never modified.

> **Status: pre-release (0.2.0).** Extracted from
> [OpenDIAlyzer](https://github.com/okohlbacher/ODIA), where it was developed.
> The output formats and the config schema are **not yet frozen** — see
> [Known limitations](#known-limitations) before depending on them.

There is a [command-line tool](#usage) and a [desktop app](#desktop-app) for
macOS, Windows and Linux. Both run the same binary and take the same config.

## Why this exists

Library generation is a different job from searching, on a different cadence: a
library is a cross-run artefact built once and reused, a search is per-run.

Two things make this tool worth having:

**It fills a hole in the OpenMS ecosystem.** Of the 151 tools in an OpenMS
install, none turns a FASTA into a predicted DIA library. `Digestor` digests,
`OpenSwathAssayGenerator` takes an *existing* transition list, `AssayGeneratorMetabo`
does metabolites. OpenMS `develop` has `PeptDeepRTInference`/`MS2Inference`/`CCSInference`,
but `PeptDeepInputBuilder` states plainly that *"modified peptides are
intentionally rejected for now"* and zero-fills the modification tensor — its
only entry points are `buildUnmodified*Batch()`. Carbamidomethylation is not
optional in practice: getting the alkylation state wrong cost this project's own
reference 12% of its identifications. **Modification-aware encoding is the piece
that is genuinely not available upstream**, and it is why this tool binds ONNX
Runtime directly rather than using OpenMS's own binding.

**The library carries the recipe that produced it.** Every content-affecting
parameter lives in one JSON file, and that JSON is embedded verbatim in the
Parquet output alongside the FASTA hash and content hashes of all three models.
A library states how it was built, without reference to anyone's defaults.

## What it is not

It is **not** more accurate than DIA-NN's predictor. Measured against a DIA-NN
*empirical* library on a Bruker timsTOF diaPASEF run (37,193 identified
precursors):

| Axis | Where this tool stands |
|---|---|
| Precursor coverage | Tied — 100.00% of observed precursors. That denominator is DIA-NN's own discovery set. |
| Precursor & fragment m/z | Tied (+4.60 µDa median, p99 0.046 mDa). |
| Fragment intensity | Within 0.015 spectral angle, on overlapping-but-not-identical fragment sets. |
| Fragment selection | **Behind.** Carries 86.5% of observed transitions; 2.53% of the misses are the observed base peak. Localised to doubly-charged y ions. |
| Retention time | **Behind** — 2.04× on a linear calibration, 1.37× after a monotone fit this tool does not ship. Elution *order* is near-tied (ρ 0.9936 vs 0.9956). |
| Ion mobility | **Behind** — 2.82% relative error against measured 1/K0, about 1.8× DIA-NN's. |

Caveats a reader should apply: these are from **one instrument**, and DIA-NN is
both the competitor and the measuring instrument — the "ground truth" library is
DIA-NN's own output, censored by its predictions and by its detection. Treat the
table as a statement of where this tool sits, not as a benchmark.

## Installing

Release builds for macOS, Windows and Linux — a CLI archive and a desktop
installer per platform — are attached to each
[release](https://github.com/okohlbacher/DIALibGen/releases). The CLI
archives are **self-contained**: unpack and run `bin/DIALibGen`. The
libraries it needs and OpenMS's own data travel with it, so there is nothing to
install and no OpenMS to set up.

```bash
curl -fsSLO https://github.com/okohlbacher/DIALibGen/releases/latest/download/DIALibGen-macos-arm64.tar.gz
tar xzf DIALibGen-macos-arm64.tar.gz
./bin/DIALibGen --help
```

There is also a [Homebrew tap](https://github.com/okohlbacher/homebrew-dialibrarygenerator)
with a cask for the app and a cask for the CLI:

```bash
brew install --cask okohlbacher/dialibrarygenerator/dialibgen       # desktop app
brew install --cask okohlbacher/dialibrarygenerator/dialibgen-cli   # CLI on PATH
```

Two casks because the app already carries its own copy of the CLI, so one cask
installing both would put the same tree on disk twice. Both require macOS 14:
the binaries are built for 13.3 (libc++ shipped `std::to_chars` there) and
Homebrew can only name whole releases, so the cask rounds up rather than promise
a machine it cannot load on.

> **The first run of the CLI takes about five minutes.** It is not stuck. macOS
> checks each of the 145 bundled libraries with Apple one at a time, and a
> `.tar.gz` cannot carry a stapled ticket that would answer for all of them at
> once. Every later run starts in about a second. See BACKLOG.md.

**Signed and notarized** since 0.9.0, with a Developer ID, so a Homebrew-installed
copy runs without a prompt. The `.dmg` and the `.app` inside it each carry a
stapled ticket and need no network; the CLI tarball cannot carry one — no archive
format can — so a quarantined copy is checked with Apple once, online, on first
run.

## Building

Requires an installed OpenMS, Apache Arrow/Parquet ≥ 19, ONNX Runtime and
nlohmann/json. All four are OpenMS dependencies already, except that ONNX
Runtime is behind OpenMS's `WITH_ONNX` option.

The environment CI builds and tests in, which is the combination known to work
on Linux x64/arm64 and macOS x64/arm64:

```bash
micromamba create -n dialibgen -c conda-forge -c bioconda \
  bzip2 cmake coin-or-cbc coin-or-utils cxx-compiler eigen glpk hdf5 \
  libarrow-acero libarrow-dataset libboost-devel libcurl libparquet \
  libsvm libzip ninja nlohmann_json numpy onnxruntime onnxruntime-cpp \
  openms=3.5.0 pyarrow python qt6-main xerces-c zlib
micromamba activate dialibgen
```

On macOS add `llvm-openmp` — Apple's clang ships no OpenMP runtime, and without
it the pragmas become no-ops and the build is slower, not broken. `python`,
`numpy`, `pyarrow` and `onnxruntime` are for the test suite, not the tool.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/your/prefix
cmake --build build -j
ctest --test-dir build --output-on-failure
cmake --install build
```

The prediction and end-to-end tests need the three `.onnx` models. Without them
those tests are not registered at all, so `ctest` reports a smaller suite that
passes — point `-DODIA_MODEL_DIR=` at a directory holding all three to run them.

The install is **relocatable**: the binary finds its data tables relative to its
own path (`<prefix>/share/DIALibGen`), so it works from a package, a
copied tree or a macOS `.app` and not only from the tree it was built in.

The sources compile at C++20, but OpenMS's imported CMake target declares
`INTERFACE_COMPILE_FEATURES cxx_std_23`, so consumers are raised to C++23
regardless. In practice the binding constraint is your OpenMS install's ABI —
in particular its Boost SONAMEs — not this code's language level.

## Models

The tool predicts with three ONNX exports of AlphaPeptDeep models:
`peptdeep_rt_dynamic.onnx`, `peptdeep_ms2_dynamic.onnx`,
`peptdeep_ccs_dynamic.onnx`.

**Every release from 0.10.1 carries them**, in `share/DIALibGen/models` next to
the binary — the second entry in the tool's own search order — so an installed
copy predicts straight away with nothing set. See
[THIRD-PARTY-NOTICES.md](THIRD-PARTY-NOTICES.md) for their provenance and the
paper to cite.

A build from source has no models unless you point `-DODIA_MODEL_DIR=` at some,
and no tagged OpenMS release ships them either: OpenMS downloads them only when
built from `develop` with `WITH_ONNX=ON`, which defaults off. For those cases,
and to refresh or verify a set:

```bash
dialibgen-fetch-models          # download, verify and install all three
dialibgen-fetch-models --check  # verify what is there, download nothing
```

It checks every file against a pinned SHA256 and installs where the tool looks.
`--dir DIR` puts them elsewhere, `--prefix DIR` targets another installation.
Re-running it is free: files that are already correct are left alone. In a
source checkout the same script is `scripts/fetch-models.sh`.

To use models from somewhere else entirely, name the directory and it wins over
the bundled ones:

```bash
export DIALIBGEN_MODEL_DIR=/path/to/models
```

…or name them individually in the config:

```json
{ "rt_model":  "/path/to/peptdeep_rt_dynamic.onnx",
  "ms2_model": "/path/to/peptdeep_ms2_dynamic.onnx",
  "ccs_model": "/path/to/peptdeep_ccs_dynamic.onnx" }
```

Relative model paths are resolved relative to the **config file**, so a config
plus a model directory is portable. The full search order is: the config,
`DIALIBGEN_MODEL_DIR`, `share/DIALibGen/models` beside the executable (the
bundled ones), then `share/OpenMS/models`. When it finds nothing it says which
file is missing and lists every directory it searched.

## Usage

```bash
DIALibGen -write_config effective.json   # see every default, materialised
DIALibGen -in proteins.fasta -config my.json -out library.parquet
```

It is a TOPP tool, so it also speaks the workflow dialect:

```bash
DIALibGen -write_ini  DIALibGen.ini   # OpenMS INI
DIALibGen -write_ctd  ./ctd/                    # KNIME/Galaxy descriptor
```

`-threads` defaults to **0 = all available cores**. `--help` shows OpenMS's own
line for it, which says `1`, with a correction printed underneath; the tool's
default is the one that applies.

`example/proteins.fasta` and `example/default.json` are a runnable starting
point. `-out` chooses the format by extension: `.parquet` carries the embedded
recipe, `.tsv` is the DIA-NN dialect and **cannot** — a TSV is not reproducible
from itself.

An unknown key in the config is rejected, not defaulted: a library built from a
typo looks exactly like one built correctly.

### Configuration

`-write_config` is the authoritative reference — it materialises every default.
The keys that most often need changing:

| Key | Default | Note |
|---|---|---|
| `enzyme` | `Trypsin/P` | Not `Trypsin`: cutting before proline regardless is what DIA-NN does, and keeping `Trypsin` cost 169,044 peptides on the human proteome. |
| `fixed_modifications` | `["Carbamidomethyl (C)"]` | Set to `[]` for a non-alkylated preparation. This is the single setting most likely to be wrong, and being wrong is expensive. |
| `precursor_charges` | `[1,2,3,4]` | `[2,3]` covers only 92.76% of observed precursors. |
| `decoys` | `"none"` | Deliberate: a library is an interchange artefact and the consumer decides its own null. DIA-NN searches shipped decoys *in addition* to its own. |
| `irt_rescale` | `false` | Off means the RT column is the model's raw 0..1 output, **not** iRT, and is not interchangeable with another tool's iRT library. Set true to export. |
| `derive_ion_mobility` | `true` | Emits 1/K0 alongside CCS. Off costs a diaPASEF consumer the entire mobility dimension. |

## Desktop app

`gui/` is a Tauri 2 + React front-end that runs the same binary: pick a FASTA,
pick an output, point it at the models, press go. See
[gui/README.md](gui/README.md).

The form is generated from the tool's own `-write_config` output, so the GUI
offers exactly the parameters the CLI has, with exactly its defaults. Settings
reach the CLI as a config file, so a library built from the GUI carries the same
embedded recipe as one built from the command line.

## Known limitations

Read these before treating output as authoritative:

- **The cache fingerprint does not capture everything.** It covers the digest
  parameters, the model contents, the decoy method and the RT domain, but not
  every field of the config. Two libraries that agree on all of those are
  treated as interchangeable.
- **`Fragment.Loss.Type` is always `"noloss"`.** The generator emits no neutral
  losses, so the column is only meaningful to a consumer that writes a
  loss-bearing library back out through this writer.
- **No model weights are shipped**, and no tagged OpenMS release contains them.
  The tool searches `DIALIBGEN_MODEL_DIR`, its own `share/` and
  `share/OpenMS/models`, and names what it could not find — but you have to
  supply the files.
- **`-write_cwl` and `-write_json` need an OpenMS built with `ENABLE_TDL=ON`.**
  Without it the tool refuses them and says so; `-write_ctd` works everywhere.
- **Windows is x64 only.** There is no arm64 Windows build; bioconda has no
  win-64 OpenMS either, so that leg builds OpenMS from source and is slower
  than the others on a cold cache.
- **The Windows builds are not signed.** The `.msi` and the setup `.exe` need
  SmartScreen's "More info → Run anyway" on first download. That certificate is
  a separate purchase and is not in place yet; macOS has been signed and
  notarized since 0.9.0.

## Provenance

Extracted from OpenDIAlyzer at commit `70a21a2` as a fresh history. The design
rationale, the full axis-by-axis library comparisons and the measurements quoted
above live in that repository's `doc/`.

## Citing

This tool is a wrapper around models and standards published by others. If you
use it, cite them:

- **AlphaPeptDeep** — Zeng *et al.*, *Nat. Commun.* **13**, 7238 (2022).
- **iRT standards** — Escher *et al.*, *Proteomics* **12**, 1111–1121 (2012).
- **OpenMS** — Röst *et al.*, *Nat. Methods* **13**, 741–748 (2016).

## Licence

BSD-3-Clause. See [LICENSE](LICENSE) and
[THIRD-PARTY-NOTICES.md](THIRD-PARTY-NOTICES.md).
