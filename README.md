# DIALibraryGenerator

Build an in-silico DIA spectral library from a FASTA: enzymatic digest →
precursor enumeration → AlphaPeptDeep predictions of retention time, MS2
fragment intensities and collision cross-section → Parquet or DIA-NN-dialect
TSV.

C++, no Python runtime. A TOPP-compatible tool built on
[OpenMS](https://github.com/OpenMS/OpenMS), which is consumed as an installed,
read-only dependency and never modified.

> **Status: pre-release (0.1.0).** Extracted from
> [OpenDIAlyzer](https://github.com/okohlbacher/ODIA), where it was developed.
> The output formats and the config schema are **not yet frozen** — see
> [Known limitations](#known-limitations) before depending on them.

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

## Building

Requires an installed OpenMS, Apache Arrow/Parquet ≥ 23, ONNX Runtime and
nlohmann/json. All four are OpenMS dependencies already, except that ONNX
Runtime is behind OpenMS's `WITH_ONNX` option.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```

The sources compile at C++20, but OpenMS's imported CMake target declares
`INTERFACE_COMPILE_FEATURES cxx_std_23`, so consumers are raised to C++23
regardless. In practice the binding constraint is your OpenMS install's ABI —
in particular its Boost SONAMEs — not this code's language level.

## Models

The tool predicts with three ONNX exports of AlphaPeptDeep models:
`peptdeep_rt_dynamic.onnx`, `peptdeep_ms2_dynamic.onnx`,
`peptdeep_ccs_dynamic.onnx`.

**They are not shipped here, and they are not in any tagged OpenMS release.**
OpenMS downloads them from `archive.openms.de` against pinned SHA256s, but only
when built from `develop` with `WITH_ONNX=ON`, which defaults off. On a released
OpenMS the files are simply absent and you must supply them yourself.

Point the tool at them from the config:

```json
{ "rt_model":  "/path/to/peptdeep_rt_dynamic.onnx",
  "ms2_model": "/path/to/peptdeep_ms2_dynamic.onnx",
  "ccs_model": "/path/to/peptdeep_ccs_dynamic.onnx" }
```

Relative model paths are resolved relative to the **config file**, so a config
plus a model directory is portable.

> The licence under which the published AlphaPeptDeep *weights* may be
> redistributed has not been established by this project. Nothing here
> redistributes them. If you package them, check first.

## Usage

```bash
DIALibraryGenerator -write_config effective.json   # see every default, materialised
DIALibraryGenerator -in proteins.fasta -config my.json -out library.parquet
```

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

## Known limitations

Carried over from extraction and **not** fixed here. Read these before treating
output as authoritative:

- **The fingerprint under-determines the library.** `min_relative_intensity` is
  formatted at 3 decimal places, so `1e-4` and `0.0` both render `0.000`; and no
  RT-domain token is recorded, so a raw-RT and an iRT library of the same content
  fingerprint identically.
- **An unknown `decoys` value is accepted** and silently treated as `mutate`,
  then written into the provenance as though it were real.
- **`-out` accepts any extension**, writing TSV for anything that is not
  `.parquet`.
- **`Fragment.Loss.Type` is written as a hardcoded `"noloss"`** by the Parquet
  writer. The generator itself never emits neutral losses, so this is reachable
  only by a consumer that writes a loss-bearing library back out.
- **No model is configured by default** and there is no environment fallback, so
  an out-of-the-box run fails inside the ONNX session constructor.
- **The data directory is a compile-time absolute path.** `-irt_standards`
  overrides it, and the default path (`irt_rescale=false`) reads no data file.
- **`schema_version` is accepted at any value and never read.**
- **CTD/CWL emission does not work out-of-tree.** `-write_ctd` and friends abort:
  OpenMS's `ToolHandler` consults a hard-coded tool list.

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
