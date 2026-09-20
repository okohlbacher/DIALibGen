# Build from source

DIALibGen requires CMake 3.21+, a compiler compatible with the installed
OpenMS, OpenMS 3.5, Arrow/Parquet 19+, ONNX Runtime and nlohmann/json. OpenMS's
public headers also require Boost, Eigen and its other exported dependencies.
LibTorch enables RT/CCS training and is enabled by default.

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

Use `-DTorch_DIR=/path/to/libtorch/share/cmake/Torch` if LibTorch is not found.
To omit training and its runtime dependency, configure with
`-DDIALIBGEN_BUILD_FINETUNE=OFF`. The resulting executable still generates
libraries and performs observed-value refinement.

The runtime is native C++; Python, NumPy, PyArrow and the Python ONNX Runtime
package are test dependencies only. Model-dependent tests require all three
ONNX models in `ODIA_MODEL_DIR`. A smaller passing suite without models is
not a prediction or end-to-end validation.

## CPU and CUDA

Portable releases bundle CPU training. A source build linked to CUDA-enabled
LibTorch can use `-machine:device cuda:0`; provide the matching CUDA and cuDNN
runtime libraries. GPU prediction additionally needs an ONNX Runtime build
with the CUDA execution provider. The standard CPU archive is not a CUDA
bundle.

Linux ARM64 and Windows CPU training have dedicated release validation because
their LibTorch packaging differs. Consult [BACKLOG.md](../BACKLOG.md) for
unresolved release checks rather than assuming that a successful compilation
proves training works.

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

Tagged builds attach packages to a draft release. Publish only after every
platform passes the prediction, refinement, training, installed-library and
relocated-binary checks, and the expected assets are present. The macOS
archives and desktop images also undergo signature and notarization checks.
