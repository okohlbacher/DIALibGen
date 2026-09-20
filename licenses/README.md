# License text origins

These files reproduce upstream notices; line endings and trailing whitespace are
normalized. Release packaging additionally copies the exact notices from every
runtime package or SDK and records their versions in `runtime-dependencies.json`.

| Files | Origin |
|---|---|
| `OpenMS-LICENSE.txt`, `OpenMS-AUTHORS.txt` | OpenMS source distribution, https://github.com/OpenMS/OpenMS |
| `Apache-Arrow-LICENSE.txt` | Apache Arrow 25 conda package aggregate LICENSE, https://github.com/apache/arrow |
| `ONNX-Runtime-LICENSE.txt` | ONNX Runtime distribution LICENSE, https://github.com/microsoft/onnxruntime |
| `PyTorch-LICENSE`, `PyTorch-NOTICE` | Official PyTorch 2.10.0 CPU Windows wheel, https://download.pytorch.org/whl/cpu/torch-2.10.0%2Bcpu-cp311-cp311-win_amd64.whl (SHA256 `17a09465bab2aab8f0f273410297133d8d8fb6dd84dccbd252ca4a4f3a111847`) |
| `LLVM-OpenMP-LICENSE.txt` | LLVM OpenMP 22.1.8 conda package LICENSE, https://github.com/llvm/llvm-project |
| `Boost-LICENSE.txt` | Boost 1.91 package `LICENSE_1_0.txt`, https://www.boost.org/LICENSE_1_0.txt |
| `nlohmann-json-LICENSE.txt` | nlohmann/json 3.12 package LICENSE.MIT, https://github.com/nlohmann/json |
| `LGPL-3.0-only.txt` | Qt LGPL-3.0 license text, https://code.qt.io/cgit/qt/qtbase.git/tree/LICENSES/LGPL-3.0-only.txt |
| `GPL-3.0-or-later.txt` | GNU GPL version 3 text, https://www.gnu.org/licenses/gpl-3.0.txt; each component's grant determines “only” versus “or later” |
| `GCC-Runtime-Exception-3.1.txt` | GCC source distribution `COPYING.RUNTIME`, https://gcc.gnu.org/git/?p=gcc.git;a=blob_plain;f=COPYING.RUNTIME;hb=HEAD |

The GCC runtime exception does not change DIALibGen's BSD license or the licenses
of other dependencies. The generated component inventory records the applicable
grant for each distributed library. AlphaPeptDeep's Apache-2.0 notice and the
attribution of derived constants and model weights are in
`THIRD-PARTY-NOTICES.md`.
