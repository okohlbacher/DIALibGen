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
[0.12.1 release record](https://github.com/okohlbacher/DIALibGen/releases/tag/v0.12.1);
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

## Fragment-intensity accuracy (standard method)

From October 2026 every fragment-intensity number in reports, benchmarks
and release notes is measured with this method. It covers b **and** y ions at
every fragment charge, read from the raw MS2 spectra. Any other observation
source must be labelled as an exception.

The earlier method scored predictions against a search engine's fragment
export. DIA-NN 1.7.12's `--report-lib-info` gives six singly charged y ions
per precursor, chosen by ordinal rather than intensity, so that method could
not see b ions or doubly charged fragments.

The method has two steps, both in `scripts/eval/fragment_intensity_eval.py`
(numpy and pyarrow only):

```bash
# observed intensities of one run: identified precursors (DIA-NN report, Q.Value <= 0.01) + the run's mzML
python3 scripts/eval/fragment_intensity_eval.py extract \
  --report report.tsv --mzml run.mzML --out run.ions.parquet          # --ppm 10 --apex-halfwidth 1
# score one or more libraries (the first is the paired reference); a second run's extraction as the noise ceiling
python3 scripts/eval/fragment_intensity_eval.py score --ions run.ions.parquet \
  --library DIALibGen=lib.parquet --library DIANN=other.parquet \
  --replicate run02=run02.ions.parquet --out-prefix result
```

There are two extraction modes. Both write the same per-ion Parquet, so
`score` treats them alike:

| mode | input | apex | m/z match | default tolerance |
|:--|:--|:--|:--|:--|
| `fragment_intensity_eval.py extract` | mzML (Orbitrap, Astral; centroided) | the scan nearest the apex in the isolation window, ±1 scan | most intense centroid per scan, summed | 10 ppm |
| `tims_extract.py` | Bruker `.d` (diaPASEF), DIA-NN 2.x report with `IM` | the frame nearest the apex in the diaPASEF window holding the precursor's m/z and apex 1/K0, ±1 frame; scans within ±0.02 1/K0 | all raw TOF peaks within the tolerance, summed | 15 ppm |

```bash
python3 scripts/eval/tims_extract.py --report report.parquet --d run.d \
  --timsdata /path/to/libtimsdata.so --out run.ions.parquet           # --ppm 15 --im-halfwidth 0.02 --workers 16
```

`tims_extract.py` reads the `.d` through Bruker's timsdata library
(`libtimsdata.so` or `timsdata.dll`). That library comes from Bruker's TDF SDK,
and DIA-NN ships a copy in its installation directory. It is not part of this
repository and is loaded only when `--timsdata` names it; a missing library
stops with a message saying so. In the timsTOF mode the shared flag also
requires the other precursor's apex 1/K0 within 2 × the 1/K0 half-width.

**Extraction (mzML mode).** For each identified precursor:

- **Scans.** It takes the MS2 scan nearest the report's apex RT, among the
  scans whose isolation window contains the precursor m/z, plus one scan on
  either side in the same window.
- **Ions.** It considers every b and y ion with ordinal 1 to n−1, no neutral
  loss, and fragment charge 1 to min(3, precursor charge), inside the scan
  range.
- **Intensity.** Each ion's intensity is the most intense centroid within
  **10 ppm**, summed over those scans. 10 ppm suits Astral and Orbitrap data;
  use 20 ppm for TOF analysers and say so.
- **Unobserved ions** are kept with intensity 0.
- **Ambiguous ions,** where two theoretical ions of the precursor lie within
  2 × the tolerance of each other, are flagged and never scored.
- **Shared ions** are flagged when the matched peak is also a b or y ion
  (z ≤ 3) of another identified precursor whose elution range covers a used
  scan in the same window. Narrow-window DIA is chimeric, so the shared flag
  is the interference check.

**Scoring.** Every predictor is scored on the same precursors and the same
ions:

- **Precursors:** those that every predictor contains and that have at least
  three observed ions in the ion set.
- **Ion sets:**
  - `union`: the ions any of the given libraries stores.
  - `theo2`: every theoretical ion, z ≤ min(2, precursor z).
  - `theo3`: every theoretical ion, z ≤ min(3, precursor z).
- **Conventions:** `zero` counts unobserved ions as 0; `excluded` keeps only
  observed ions. Either way, an ion a library does not store is predicted 0. A
  predictor with no intensity on any scored ion gets cosine 0, never a dropped
  precursor, so a library missing a whole ion group (for example charge-3
  fragments) is scored as a miss on that group, not hidden.
- **Variants:** `all` and `unshared` (shared ions removed).

**Outputs:**

- the median cosine and spectral angle, with paired differences against the
  first library;
- a breakdown by ion type × fragment charge × precursor charge: theoretical
  and observed ions per precursor, share of observed and of predicted
  intensity, and cosine within the group;
- per-precursor scores.

**Report** `theo2`/`zero` as the headline figure, with `union`, `excluded`
and `unshared` as sensitivity rows. Add the replicate ceiling (run 2 scored
as a predictor of run 1 on the same ions) whenever a second run of the same
sample exists.

`test/check_fragment_intensity_eval.py` (ctest
`fragment_intensity_eval`) checks the method on a synthetic run: apex choice,
scan summing, the ppm cut, both conventions, the shared flag, the zero score
for a predictor with nothing on a group, and all three library layouts
(DIALibGen Parquet, DIA-NN Parquet and TSV). For `tims_extract.py` it checks
only the import and the refusals, because CI has no Bruker library. On one
diaPASEF run (300 precursors) the repository version gives the same
intensities, for all 17,771 ions, as the copy the timsTOF benchmark used.

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
