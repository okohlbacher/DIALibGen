# Validation

The release suite covers generation, observed-value refinement and native RT/CCS
training through the public API and the same TOPP executable users install.
Run it with all three bundled models and Python test dependencies available:

```bash
ctest --test-dir build --output-on-failure
cd gui
npm ci
npm run test:coverage
npm run typecheck
npm run build
cd src-tauri
cargo test
cargo clippy --all-targets -- -D warnings
```

Python, NumPy, PyArrow and Python ONNX Runtime are test tools; portable CLI
packages do not depend on them. Release CI checks the required test names and
rejects skipped tests. A source build without models can run a smaller suite,
which does not establish prediction or training correctness.

## What the checks establish

- Native encoder/predictor comparisons use an independent Python reference,
  pinned output values, real ONNX models and tolerances stated in each test.
- Generation checks digestion, modification uniqueness, decoy identities and
  masses, charge domains, calibration, intensity pruning and deterministic
  prediction across worker counts. Capacity tests check the 32-bit boundary
  without allocating billions of transitions.
- IO tests use independent Arrow fixtures, malformed nested arrays, missing
  values, multi-batch files, strict text parsing and injected write failures.
- Refinement tests exercise report gates, joins, duplicate selection, RT/IM and
  fragment-intensity replacement, target/decoy symmetry and output provenance.
  Regressions cover mixed RT-unit refusal, contextual modification-name joins,
  unresolved-token counts, compact empirical-reference refusal and missing
  decoy fragments with restriction both enabled and disabled.
- Training tests verify held-out protein cohorts, rejection paths, stopping
  rules, optimizer updates, selection of the best checkpoint, ONNX writeback,
  complete provenance and prediction of precursors absent from the report.
- GUI tests cover forms, effective JSON, native Tauri IPC and a real child
  process lifecycle using a controlled executable, including cancellation.
  A configuration fixture comes from the real CLI and is checked for drift;
  browser tests exercise decimal NCE/m/z fields and whole-number threads.
  Native tests verify that an occupied output path, including a dangling
  symlink, is refused without starting a child or changing the existing file.

## Release verification

Each release requires the complete CPU generation/refinement/training suite,
an installed CMake consumer and execution of all three modes after runtime
relocation on Linux x64/ARM64, macOS x64/ARM64 and Windows x64. The desktop
packages must contain the expected executable, models and runtime libraries,
and pass installation and launch checks. Successful compilation alone does
not establish runtime or installer correctness.

macOS gates require valid signatures, accepted notarization, stapled desktop
tickets and fresh-machine startup measurements. Distribution gates attribute
bundled libraries and the final AppImage payload, retain upstream notices,
and verify the corresponding-source manifests, archive parts and checksums.
The tagged draft must contain every expected asset before publication.
After publication, both Homebrew casks must pin the final CLI/DMG SHA-256
hashes and pass normal installation checks on Apple Silicon and Intel Macs.
A required workflow step is not evidence that the check passed.

The 0.11.0 candidates passed CPU generation/refinement/training, installed SDK
consumers and relocated execution on all five target platforms. The
[Unix packaging run at `8cd1726`](https://github.com/okohlbacher/DIALibGen/actions/runs/35534187789)
passed both Linux Debian/AppImage installation and native-window checks, both
Mac signing/notarization/stapling checks, corresponding-source verification
and fresh-machine startup measurements. The
[Windows run at `222bf79`](https://github.com/okohlbacher/DIALibGen/actions/runs/35537203627)
verified MSI administrative extraction, an actual NSIS installation, identical
embedded payloads and all three CLI modes. Native UI Automation confirmed a
rendered, enabled control after the configuration response, followed by clean
window closure and complete NSIS removal. MSI installation/upgrade behavior is
not established by extraction.

The [follow-up code checks at `4219a21`](https://github.com/okohlbacher/DIALibGen/actions/runs/35535916865)
passed all four Unix builds and the desktop suites on Linux, macOS and Windows.
They include the application-exit cleanup regression without a mocked Windows
application.

These are candidate results. Final tagged-build results, public asset checksums
and Homebrew installation verification belong to the
[0.11.0 release record](https://github.com/okohlbacher/DIALibGen/releases/tag/v0.11.0);
tagged-build and asset checks must pass before publication. Homebrew installation
checks then validate the public downloads before the tap is updated.

## Measured coverage and limits

The local Linux GCC coverage run on 20 September 2026 measured **91.3% of C++
lines (4,213/4,612)** and **95.7% of functions (287/300)**. Raw GCC branch
coverage was **48.2% (8,324/17,270)**, including compiler-generated exception
edges. The run passed all 45 CTest cases after the RT-unit, modification-join,
capacity, numeric-parser and overwrite corrections. A second build without
training passed all 41 applicable tests and its installed SDK consumer.

A Linux installed-package check on the same date linked and called all three
libraries from a moved install prefix. The exported CMake targets contained no
build-machine ONNX Runtime or LibTorch path. This local result does not stand
in for the required consumer checks on the other release platforms.

The frontend suite, including desktop refinement and fine-tuning, passed
**92 tests** and measured **100% lines (365/365)**, **100% functions (164/164)**,
**98.03% statements (498/508)** and **93.32% branches (573/614)**. Tests cover
mode-specific inputs and settings, all native tuning controls, selected-head
model checks, stale asynchronous results, configuration round-trips and errors.
Type checking and the frontend build also passed. The native
GUI passed **24 macOS Rust tests** and Clippy with warnings treated as errors.
Instrumented native tests measured **83.54% of production-library lines
(543/650)**, including the new training schema, command validation and model
directory selection, cancellation during startup and shutdown. Bootstrap and
several native error paths remain uncovered. Linux also checks isolation of the
bundled CLI from the AppImage launcher’s library environment.

These are local coverage measurements, not claims of complete coverage or
measurements of the final packaged artifacts.

Actual local macOS desktop interactions exercised both-head tuning, refinement,
combined tuning/refinement, selected-model generation, overwrite refusal and
cancellation followed by a successful restart. Tuning retained all 4,701 input
precursor keys in the synthetic fixture. The real-child application-exit
regression verifies process reaping and temporary configuration removal; a final
manual Command-Q retest was unavailable while the test Mac was locked. Linux
installer checks establish native-window startup and embedded CLI execution,
not full interactive form coverage.

Unexecuted code includes platform-specific bootstrap and error paths, GPU
provider selection and CUDA training. CPU release checks do not validate a
CUDA source build. Alongside synthetic DIA-NN 2.x fixtures, the native Parquet
loader test checks all 12 fragment identities, m/z, quantities and scores from
one real DIA-NN 2.0 `--export-quant` report row, including zero quantities and a
negative score. Its source hash and engine build are recorded in
`test/refinement/intensity_match.cpp`; peptide/protein identifiers are not copied.
This covers that row's token grammar and column types, not every DIA-NN version
or export setting.
Tests on synthetic reports prove software contracts; they do not establish
improved proteomics performance. The qualified
[benchmark results](benchmark.md) remain the scientific evidence.

## macOS startup

Older signed CLI archives showed very different first-launch times: 335 s for
a 0.9.0 tarball and about 33 s for a 0.10.0 Homebrew cask installation. Subsequent
launches were about one second. Those observations do not establish the
behavior of the new bundle.

The release measurement workflow compares archive and cask delivery on
separate fresh VMs, recording artifact hashes, quarantine attributes, first
and warm launch times. It does not launch or assess the candidate before the
first timed invocation, apart from Homebrew's own installation checks.

The archived 0.10.1 comparison completed on separate macOS 14.8.9 ARM VMs:
raw tarball first/warm launch **0.554 / 0.063 s**, Homebrew cask
**0.525 / 0.061 s**. Homebrew retained quarantine attributes; the curl archive
had none. The older long delay did not reproduce on these hosted machines.
[Measurement run](https://github.com/okohlbacher/DIALibGen/actions/runs/35521436604).

The signed, notarized **0.11.0 candidate `8cd1726`** was then measured on two
fresh macOS 14.8.9 ARM VMs with Homebrew 6.0.20:

| Delivery | Installation | First `--help` | Warm `--help` | Quarantine attributes |
|---|---:|---:|---:|---:|
| curl + tar | 3.806 s | 1.959 s | 0.391 s | 0 |
| Homebrew cask | 7.561 s | 0.792 s | 0.188 s | 4,263 |

Both routes used archive SHA-256
`6da8aa0a54573bd827a4e407685aa354157dc6d26cc68ddc9ac4a28fd027cedd`.
The cask used the published command-wrapper layout and unmodified Homebrew
quarantine behavior. A localhost staging proxy served the CI artifact before
publication, so installation timing excludes public GitHub download latency.
Each invocation exited successfully and reported version 0.11.0. Raw records,
signatures and policy logs are attached to the
[candidate measurement run](https://github.com/okohlbacher/DIALibGen/actions/runs/35534187789).

The several-minute delay did not reproduce in either comparison. These are
single observations per delivery route on hosted VMs, not a guarantee for
every Mac or network. These prior-candidate measurements do not replace checks
against the final release artifacts.
