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

NCE starts blank, meaning the CLI chooses the current instrument's default
(for example, QE 30 or Lumos 25). Entering a number makes it explicit; clear the
field to return to automatic selection. Imported configs and saved presets
with a numeric NCE retain that explicit value when the instrument changes.
Precursor charges must be unique integers from 1 to 8; fragment charges are
limited to the model's supported values, 1 and 2.

Cancellation stops the child process immediately. A cancelled or interrupted
run can leave a hidden `.dialibgen-tmp-*` staging directory beside the selected
output. After confirming that no DIALibGen process is using it, you may remove
that directory manually. The GUI does not delete matching directories because
another run may own them. Output replacement is atomic during normal operation;
it does not promise durability after a power loss.

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
For a development binary that finds models only through OpenMS's compiled-in
data path, select that model directory explicitly or set `DIALIBGEN_MODEL_DIR`;
the GUI cannot discover a different binary's compiled-in search paths.

```bash
npm run typecheck
npm run test:coverage
cd src-tauri
cargo test
cargo clippy --all-targets -- -D warnings
```

Frontend coverage includes every production TypeScript/TSX file, including the
native bridge and entry point. Reports are written to `coverage/`. The Rust
tests separately exercise configuration files, settings recovery, and the
Tauri invoke boundary with a real child process (success, failure, cancellation); frontend coverage does not measure the native backend or an
installed application's interaction with the platform webview.

`npm run tauri build` produces the platform installer. Release packaging
supplies the complete CLI tree as a Tauri resource. CMake's version guard
checks that the CLI and app version files agree.

The packaging hook generates `resources/third-party-licenses/` with dependency
license texts, exact-version source URLs and integrity checksums. It covers
production npm dependencies and a conservative Cargo dependency closure for the
build target. `npm run licenses -- <target-triple>` regenerates the inventory;
new dependencies with missing license text fail the build until their notices
are supplied. Bundled CLI and system runtime libraries have separate notices.
