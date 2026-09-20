# Remaining validation and scope

The 0.11.0 integration combines the generation and refinement backlogs. This
file distinguishes implemented changes from checks that still need evidence.

## Final release gates — pending

The fixes are implemented, but a passing earlier candidate does not validate
the final release commit. Before publishing 0.11.0:

- Run the complete CPU generation/refinement/training suite, installed CMake
  consumer and relocated-bundle smoke checks on Linux x64/ARM64, macOS
  x64/ARM64 and Windows x64. Attach the exact commit and job results.
- Build the desktop installers on the same platforms. Check the packaged
  executable, bundled models, native dependencies and installed GUI behavior.
- Verify macOS signatures, accepted notarization, stapled desktop tickets and
  candidate startup measurements against the final artifacts.
- Verify that every binary has its required dependency inventory, notices and
  corresponding-source asset. Inspect the completed desktop source inventories,
  especially libraries added during AppImage packaging, and confirm that the
  source archives upload within the release-host size limit.
- Record final coverage and test counts after the review corrections. The
  existing measurements in [docs/testing.md](docs/testing.md) are qualified
  baselines, not a claim that every branch is covered.
- Publish only after the tagged draft release contains all expected assets.
  Update both Homebrew casks to 0.11.0 with the final ARM64/x64 CLI and DMG
  SHA-256 digests, then verify installation from the public downloads.

Linux ARM64 and Windows require actual training evidence: the earlier conda
LibTorch 2.10.0 ARM package crashed during LSTM execution, and Windows uses a
separate dependency build. A successful compile does not close these gates.

## macOS first-launch evidence — earlier candidates

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
every Mac or network. The measured result does not justify changing the
archive format. These measurements close the investigation of the historical
delay on the tested candidates; they do not replace the final artifact check.

## Integrated work

- **Cross-run transfer:** measured in the completed 0.10.1 / DIALibRefine
  0.3.0-dev K562 benchmark. Results and limitations are summarized in
  [docs/benchmark.md](docs/benchmark.md); untouched-data confirmation remains
  future scientific validation.
- **Embedded tuned-model provenance:** ONNX metadata accompanies the JSON
  training sidecars.
- **Obsolete Python tuning helpers:** not imported into the unified product;
  runtime training is C++/LibTorch.
- **Documentation drift:** the native parameter reference is generated from the
  executable and supports a `--check` gate.
- **Met-excised missed-cleavage peptides:** restored, covered by a regression,
  with cache fingerprint v5.
- **RT units:** retaining unidentified entries while writing observed RT requires
  whole-library RT re-prediction in reference-run minutes; otherwise use
  `-no_write_rt`. CCS-only tuning cannot establish that RT contract.
- **Modification joins:** additional names resolve through the OpenMS database
  with residue/terminal context. Unresolved library tokens are warned about and
  counted in provenance.
- **Decoy intensities:** replacement preserves matching target/decoy fragment
  sets, including the no-restriction path.
- **Capacity:** generation, IO and library transformations guard 32-bit transition
  offsets; libraries exceeding 4,294,967,295 transitions must be split.
- **Desktop defaults:** tests use a real CLI configuration fixture; decimal NCE
  and m/z values remain editable, thread counts are whole numbers, and existing
  output paths are refused before a child process starts.
- **Distribution:** runtime archives omit development files. Dependency notices,
  exact-source inventories and accompanying source archives are part of the
  packaging gates; final platform evidence remains required above.

## Input and scientific limits

- DIA-NN 2.x `Fr.N.Id` parsing is covered by synthetic fixtures. A real
  `--export-quant` fixture is still needed to establish interoperability with
  that exported grammar. Unsupported fragment layouts fail explicitly.
- Compact Parquet is supported for input libraries, but empirical reference
  libraries must use the long format. Compact empirical references are refused.
- The OpenMS sequence representation supports one modification per residue;
  prediction of stacked modifications on a single residue is unsupported.
- The prior cross-run benchmark used exposed technical replicates. Confirmation
  on untouched data remains scientific work, not a result implied by the software
  contract tests or the 0.11.0 release.

## Explicit product scope

Portable release archives provide CPU training. CUDA training is supported by
source builds with compatible CUDA LibTorch and runtime libraries; distributing
a separate CUDA/cuDNN bundle is outside this release's scope. A Python-on-GPU
comparison is a research control, not a prerequisite for any accuracy or
performance claim made in the current product documentation.

Installed CMake package relocation passed on Linux on 20 September 2026:
a consumer linked and called all three libraries from a moved install prefix.
The exported targets contain no build-machine ONNX Runtime or LibTorch path.
The same consumer check is required by the platform release workflows.
