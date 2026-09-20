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

Release workflows require an installed CMake consumer and execution of all
three modes after runtime relocation on each platform. macOS packaging requires
signing, notarization and startup checks. Distribution gates attribute actual
bundled libraries and AppImage payload files, retain upstream notices and supply
checksum-verified corresponding-source assets. These gates describe the required
release evidence; their presence in the workflow is not a passing result.
Final candidate jobs and artifact verification remain tracked in
[BACKLOG.md](../BACKLOG.md).

## Measured coverage and limits

The local Linux GCC coverage run on 20 September 2026 measured **91.3% of C++
lines (4,213/4,612)** and **95.7% of functions (287/300)**. Raw GCC branch
coverage was **48.2% (8,324/17,270)**, including compiler-generated exception
edges. The run passed all 45 CTest cases after the RT-unit, modification-join,
capacity, numeric-parser and overwrite corrections. A second build without
training passed all 41 applicable tests and its installed SDK consumer.

The current frontend suite passed **59 tests** and measured **100% lines
(277/277)**, **100% functions (113/113)**, **98.43% statements (315/320)** and
**89.84% branches (239/266)** after the default-fixture, output-path and numeric
input corrections. Type checking and the frontend build also passed. The native
GUI passed **20 Rust tests** and Clippy with warnings treated as errors. Its
macOS production-library coverage is **76.58% of lines (327/427)** after the
output-path and automatic model-discovery regressions.

These are local measurements, not final platform CI results or claims of
complete coverage. The final release record must identify the validated commit.

Unexecuted code includes platform-specific bootstrap and error paths, GPU
provider selection and CUDA training. CPU release checks do not validate a
CUDA source build. The DIA-NN 2.x fragment-column tests use synthetic fixtures;
a real `--export-quant` fixture remains necessary to establish interoperability.
Tests on synthetic reports prove software contracts; they do not establish
improved proteomics performance. The qualified
[benchmark results](benchmark.md) remain the scientific evidence.
