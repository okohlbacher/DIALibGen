# DIALibraryGenerator desktop GUI

A cross-platform desktop front-end for the DIALibraryGenerator CLI, built with
[Tauri 2](https://tauri.app) (Rust backend + the OS's own webview) and React.
The CLI stays the source of truth: the GUI shells out to it, streams its
progress, and never reimplements a single parameter.

## The one idea worth knowing

**The form is built from the tool's own `-write_config` output.** On startup the
backend runs `DIALibraryGenerator -write_config`, which materialises every
default, and hands the JSON to the frontend. `src/paramLayout.ts` is an
*overlay* on that — prose, grouping, choice lists — and infers each widget's
type from the tool's own default value.

Two consequences, both deliberate:

- The GUI cannot offer a key the CLI does not have, nor default one differently.
- A key the CLI **gains** shows up in the GUI immediately, under Advanced,
  labelled as not yet described. A new parameter that the GUI silently hides is
  worse than one with a missing description, and this is what stops the two
  drifting apart without anyone noticing. There is a test for it.

The config reaches the CLI as a **file**, not as command-line values — which is
how the CLI wants it, and is what lets the Parquet writer embed the recipe
verbatim.

## Layout

```
gui/
  src/                  React frontend
    api.ts              the window.dialibgen bridge over Tauri invoke/events
    paramLayout.ts      the overlay: grouping, prose, choices, type inference
    ParamField.tsx      one widget per inferred kind
    App.tsx             inputs, model picker, run, log
    testing/mockBridge.ts  an in-memory window.dialibgen for the tests
  src-tauri/            Rust backend
    src/dialibgen.rs    resolve/probe the binary, ask it for its config, run it,
                        stream stderr as events, cancel
    src/settings.rs     named presets + last-used, atomic JSON in the config dir
    tauri.conf.json     window, bundle, icons
    resources/dialibgen/  the bundled CLI + share/ (dev: symlinks; release: real files)
```

## Models

The three AlphaPeptDeep ONNX exports are **not shipped** — no tagged OpenMS
release contains them and their redistribution licence is unestablished. The GUI
therefore has a model-directory picker and refuses to run until all three are
present, naming the ones it could not find. See the top-level README.

## Develop

```bash
npm install
npm run tauri dev
```

`resolve_binary` looks for the CLI in this order: `DIALIBGEN_BIN`, the bundled
`resources/dialibgen/bin/DIALibraryGenerator`, then `DIALibraryGenerator` on
`PATH`. For dev, point it at a local build:

```bash
DIALIBGEN_BIN=/path/to/build/DIALibraryGenerator npm run tauri dev
```

…or drop a symlink at `src-tauri/resources/dialibgen/bin/DIALibraryGenerator`.
The directory skeleton is tracked (with `.gitkeep`) because `tauri-build`
refuses to configure when a declared resource path is missing.

## Check

```bash
npm run typecheck && npm test           # frontend
cd src-tauri && cargo test && cargo clippy --all-targets -- -D warnings
```

## Build

```bash
npm run tauri build              # release .app/.dmg/.msi/.deb/.AppImage
npm run tauri build -- --debug   # faster, unsigned, for local checking
```

## Versions

Every file carrying the version is checked against `project(VERSION)` at CMake
configure time — see `cmake/version-guard.cmake` for what a bump has to touch.
