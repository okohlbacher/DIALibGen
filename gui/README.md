# DIALibGen desktop app

The desktop app generates a predicted library from a protein FASTA using the
same `DIALibGen` executable as the command line. Release builds include the
three AlphaPeptDeep ONNX models. Refinement and RT/CCS tuning are available
through the [CLI](../docs/usage.md).

1. Select a FASTA and an output ending in `.parquet` or `.tsv`.
2. Match digestion, modifications and instrument settings to your experiment.
3. Leave the model directory empty to use bundled models, or select a directory
   containing all three ONNX files.
4. Select the thread count and generate. One thread is the default; zero
   requests automatic inference parallelism and uses more memory.

Parquet carries the effective generation recipe in its metadata. Save the
config separately when exporting TSV. Presets and last-used generation
settings are stored locally.

## How the form stays compatible

The backend calls `DIALibGen -mode generate -write_config` and builds the form
from that mode's effective settings. Refinement and tuning parameters are not
mixed into it. The layout adds labels, groups and choices; unknown generation
settings appear under Advanced. Model paths and provenance-only fields are
managed separately from editable parameters.

## Development

The app uses Tauri 2, React and the platform webview. Install Node.js 22 and Rust,
then run:

```bash
npm ci
DIALIBGEN_BIN=/path/to/build/DIALibGen npm run tauri dev
```

Binary lookup checks `DIALIBGEN_BIN`, the bundled
`src-tauri/resources/dialibgen/bin/DIALibGen`, then `PATH`. A local binary must
have access to its runtime libraries and models.

```bash
npm run typecheck
npm test
cd src-tauri
cargo test
cargo clippy --all-targets -- -D warnings
```

`npm run tauri build` produces the platform installer. Release packaging
supplies the complete CLI tree as a Tauri resource. CMake's version guard
checks that the CLI and app version files agree.
