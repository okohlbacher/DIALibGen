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

Final 0.11.0 verification is in progress. The final commit, platform job results
and public artifact/installer verification have not yet been recorded here;
the local measurements and earlier candidate evidence below do not replace them.

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
**86 tests** and measured **100% lines (351/351)**, **100% functions (164/164)**,
**98.34% statements (476/484)** and **93.82% branches (562/599)**. Tests cover
mode-specific inputs and settings, all native tuning controls, selected-head
model checks, stale asynchronous results, configuration round-trips and errors. Type checking and the frontend build also passed. The native
GUI passed **23 macOS Rust tests** and Clippy with warnings treated as errors.
The earlier generation-only production-library measurement was **76.58% of
lines (327/427)**; coverage must be remeasured for the new tuning backend before
final verification. Linux also checks isolation of the bundled CLI from the
AppImage launcher’s library environment.

These are local measurements, not final platform CI results or claims of
complete coverage. The final release record must identify the validated commit.

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

The signed, notarized **0.11.0 candidate `79932aa`** was then measured on two
fresh macOS 14.8.9 ARM VMs with Homebrew 6.0.20:

| Delivery | Installation | First `--help` | Warm `--help` | Quarantine attributes |
|---|---:|---:|---:|---:|
| curl + tar | 2.502 s | 0.987 s | 0.214 s | 0 |
| Homebrew cask | 7.844 s | 1.329 s | 0.245 s | 387 |

Both routes used archive SHA-256
`12e89c76a2ad0d221bd23b0cf925fcf1017df15bc88be56c5a5f089650e96ad4`
and executable code-directory hash `14f127b03fc40a1cc2a27ed24d06b3a1713646ed`.
The cask used the published command-wrapper layout and unmodified Homebrew
quarantine behavior. A localhost staging proxy served the CI artifact before
publication, so installation timing excludes public GitHub download latency.
Each invocation exited successfully and reported version 0.11.0. Raw records,
signatures and policy logs are attached to the
[candidate measurement run](https://github.com/okohlbacher/DIALibGen/actions/runs/35522179576).

The several-minute delay did not reproduce in either comparison. These are
single observations per delivery route on hosted VMs, not a guarantee for
every Mac or network. These prior-candidate measurements do not replace checks
against the final release artifacts.
