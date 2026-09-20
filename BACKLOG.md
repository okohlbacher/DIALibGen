# Remaining validation and scope

The 0.11.0 integration combines the generation and refinement backlogs. This
file distinguishes implemented changes from checks that still need evidence.

## macOS first-launch measurement — pending

Older signed CLI archives showed very different first-launch times: 335 s for
a 0.9.0 tarball and about 33 s for a 0.10.0 Homebrew cask installation. Subsequent
launches were about one second. Those observations do not establish the
behavior of the new bundle.

The release measurement workflow must compare archive and cask delivery on
fresh machines that have not validated the tested code-directory hashes.
Record OS/architecture, artifact hashes, quarantine state, first and warm
launch times, and the delivery route. A warm retry on the same Mac is not a
second cold measurement. Choose a packaging change only after that comparison;
wrapping loose CLI files in a stapled container does not itself demonstrate
that the installed files avoid individual validation.

Evidence to attach before closing: workflow URL and both route measurements
for the released 0.11.0 artifact. **Pending.**

## Portable CPU training — release CI pending

- Linux ARM64: exercise the full training path with the selected runtime. The
  earlier conda LibTorch 2.10.0 aarch64 package crashed in LSTM execution; a
  successful compile is insufficient evidence that the replacement works.
- Windows x64: verify generation, refinement, tuning and relocated-bundle
  execution with the packaged LibTorch runtime.

Evidence to attach before closing: successful platform jobs and installed
bundle smoke-test logs. **Pending.**

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
  with cache fingerprint v4.

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
