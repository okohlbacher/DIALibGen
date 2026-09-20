# DIALibGen desktop app

The desktop app generates, refines and fine-tunes spectral libraries using the
same `DIALibGen` executable as the command line. Release builds include the
three AlphaPeptDeep ONNX models and CPU training.

## Choose a workflow

| Mode | Inputs | Result |
|---|---|---|
| Generate | Protein FASTA | A predicted spectral library |
| Refine | Library and one run's DIA-NN report | A library filtered and updated with accepted observations |
| Fine-tune | Library and one run's DIA-NN report | The complete library re-predicted with adapted RT/CCS models |

Choose an unused output ending in `.parquet` or `.tsv`. A library input may
use either format; identification reports use Parquet. The same report-column,
quality-gate and RT-unit rules apply as in the [CLI](../docs/usage.md).

Generation lets you set digestion, modifications, instrument and fragment
selection. Leave NCE blank for the instrument default, such as QE 30 or Lumos
25; a numeric value remains explicit when the instrument changes. Precursor
charges must be unique integers from 1 to 8; fragment charges are limited to
1 and 2. One inference thread is the default; zero requests automatic
parallelism and uses more memory.

Refinement exposes quality gates, filtering, observed RT/mobility and fragment
intensity replacement. It can also fine-tune the models before applying those
observations. Keeping unidentified precursors while writing observed RT
requires successful RT tuning; keeping the original RT values is another
option. The CLI refuses combinations that would mix RT units.

## Fine-tuning

1. Select Fine-tune, an input library, a single-run DIA-NN report and a new output.
2. Choose RT, CCS or both heads. Leave the model directory blank for discovered
   models (bundled in releases), or choose a directory containing the selected heads' ONNX files.
3. Set the training recipe. Controls include observation filters, protein
   cohorts, learning rate, epochs, batch size, stopping rules, device and seed.
   Training threads are separate from inference threads.
4. Optionally choose a directory to retain tuned ONNX models, provenance and
   training trajectories. Without it, models are temporary and the output
   library is retained.
5. Start the run and inspect its log and completion result.

Fine-tune preserves the whole input library and writes predictions; choose
Refine to filter it or replace values with observations. The training quality
filter is separate from refinement's output-quality gates. CPU is supported
by the portable packages; CUDA options require a compatible source build.

The training controls and defaults come from the executable's TOPP parameter
schema. Full-fit and no-inner-validation settings change which observations
are held out, as their descriptions explain. Evaluate transfer on a separate
run; a close fit to the training run is not evidence of improved search results.

## Settings and outputs

Each mode retains its own settings. Presets and saved desktop JSON include the
selected mode and training recipe. Reopen these files in the desktop app; they
are not the CLI's plain `-config` files. Plain CLI JSON can also be loaded for
the selected mode, and older generation presets remain usable.

The desktop uses the model directory shown in the form for every selected
head. Imported generation JSON's individual model paths are replaced by that
directory's models. Use the CLI when models must come from different directories.
When reusing a tuned RT model in Generate, turn off **free cysteine rt correction**
under Advanced; generation emits normalized RT. The [CLI example](../docs/usage.md#model-tuning)
explains the model scale and the stock correction.

Parquet includes provenance. Refinement and tuning also write a `.refine.json`
sidecar beside the library. Retained tuned models have `.tune.json` and
`.trajectory.tsv` sidecars. Existing outputs are refused; choose new names or
an unused model-output directory.

Cancellation stops the child process immediately. A cancelled or interrupted
run can leave a hidden `.dialibgen-tmp-*` staging directory beside an output.
After confirming that no DIALibGen process is using it, you may remove that
directory manually. Completed model files from an earlier head may also remain.
The GUI does not remove unrelated files. Output replacement is atomic during
normal operation; it does not promise durability after a power loss.

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
