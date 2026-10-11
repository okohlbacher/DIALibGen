# Build from source

DIALibGen requires CMake 3.21+, a C++23 compiler compatible with the installed
OpenMS, OpenMS 3.5, Arrow/Parquet 19+, ONNX Runtime and nlohmann/json. OpenMS's
public headers also require Boost, Eigen and its other exported dependencies.
LibTorch and Qt6 Core enable RT/CCS training and are enabled by default.

Point `CMAKE_PREFIX_PATH` at your installed dependency prefixes. On macOS,
install an OpenMP runtime such as `llvm-openmp`. Match LibTorch's C++ ABI to
the compiler and the other dependencies.

```bash
scripts/fetch-models.sh --dir models
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH='/path/to/dependencies;/path/to/libtorch' \
  -DODIA_MODEL_DIR="$PWD/models" -DCMAKE_INSTALL_PREFIX="$PWD/install"
cmake --build build -j 4
ctest --test-dir build --output-on-failure
cmake --install build
./install/bin/DIALibGen --helphelp
```

`fetch-models.sh` downloads the three stock AlphaPeptDeep models against pinned
SHA256 digests and copies DIALibGen's pruned MS2 graph,
`data/models/peptdeep_ms2_pruned.onnx`, beside them. The install carries the
pruned graph in `share/DIALibGen/models` in any case. It is derived from the
pinned MS2 model by `scripts/prune-ms2-modloss.py` (onnx and onnxruntime for
Python), which checks the graph against the original before writing it.

Use `-DTorch_DIR=/path/to/libtorch/share/cmake/Torch` if LibTorch is not found.
To omit training and its runtime dependency, configure with
`-DDIALIBGEN_BUILD_FINETUNE=OFF`. The resulting executable still generates
libraries and performs observed-value refinement.

The runtime is native C++; Python, NumPy, PyArrow and the Python ONNX Runtime
package are test dependencies only. Model-dependent tests require all three
ONNX models in `ODIA_MODEL_DIR`. A smaller passing suite without models is
not a prediction or end-to-end validation. See [validation and coverage](testing.md).

## CPU and CUDA

Portable releases bundle CPU training. A source build linked to CUDA-enabled
LibTorch can use `-machine:device cuda:0`; provide the matching CUDA and cuDNN
runtime libraries. GPU prediction additionally needs an ONNX Runtime build
with the CUDA execution provider (for example conda-forge `onnxruntime-cpp`
built for CUDA 12) and, at run time, cuDNN 9 for CUDA 12 with the unversioned
`libcudnn.so`: without it the RT and CCS models fall back to the CPU with a
warning (see [GPU](usage.md#gpu)). The standard CPU archive is not a CUDA
bundle, and no package bundles cuDNN.

A tested CUDA build on Linux x86_64 (H100, driver 580, CUDA 12.8 toolkit):
pytorch.org's `libtorch-shared-with-deps-2.10.0+cu128` (the same LibTorch API
as the CPU release) with

```bash
cmake ... -DTorch_DIR=$LIBTORCH/share/cmake/Torch \
  -DCMAKE_CUDA_COMPILER=/usr/local/cuda-12.8/bin/nvcc \
  -DTORCH_CUDA_ARCH_LIST=9.0 -DCMAKE_CUDA_STANDARD=17
```

`TORCH_CUDA_ARCH_LIST` skips LibTorch's GPU-detection test project, which
otherwise inherits this project's C++23 setting as a CUDA dialect nvcc does not
know. That LibTorch ships the complete cuDNN 9 for CUDA 12 (`lib/libcudnn*.so.9`),
so cuDNN stays on for training. ONNX Runtime can use the same copy: put a
`libcudnn.so -> $LIBTORCH/lib/libcudnn.so.9` link in a directory on
`LD_LIBRARY_PATH`, and both libraries share one loaded cuDNN.

Linux ARM64 and Windows use different native runtime packages. The
[release checks](testing.md#release-verification) require actual CPU training
on each platform; a successful compilation does not establish that it works.

## Check documentation

The native parameter table comes from the compiled TOPP schema:

```bash
python3 scripts/update-parameters.py ./build/DIALibGen
python3 scripts/update-parameters.py ./build/DIALibGen --check
```

Commit a refreshed [parameter reference](parameters.md) whenever the schema
changes. `--check` exits unsuccessfully if the checked-in table differs.
See [gui/README.md](../gui/README.md) for desktop development and checks.

## Release verification

Tagged builds attach packages to a draft release. The
[validation guide](testing.md#release-verification) defines the required
platform, installer, signing and distribution checks before publication.
