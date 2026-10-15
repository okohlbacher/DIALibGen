// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace ODIA
{
  /// What this machine and this ONNX Runtime offer for prediction, found
  /// without building a session (probeDeviceCaps).
  struct DeviceCaps
  {
    bool cuda_provider = false;     ///< ONNX Runtime has the CUDA execution provider
    bool nvidia_gpu = false;        ///< an NVIDIA driver is loaded and CUDA_VISIBLE_DEVICES hides nothing
    bool cudnn = false;             ///< cuDNN 9 loads (libcudnn.so.9, else libcudnn.so)
    std::string cudnn_detail;       ///< its version, or why not (one line)
    std::string cudnn_elsewhere;    ///< a directory off the library path that holds it (a hint)
    bool coreml_provider = false;   ///< ONNX Runtime has the CoreML execution provider (macOS)
    std::vector<std::string> diagnostics;   ///< what was checked, for the verbose log
  };

  /// Finds the DeviceCaps. Linux, with the CUDA provider and an NVIDIA GPU:
  /// opens libcudnn.so.9 (else libcudnn.so) by name, through the loader's own
  /// search (RPATH, LD_LIBRARY_PATH, the ld.so cache) as ONNX Runtime does,
  /// and keeps it loaded; when that fails, only looks (never loads) in the
  /// usual places a pip, conda or system cuDNN sits, so that the warning can
  /// say where.
  DeviceCaps probeDeviceCaps();

  /// Where generate predicts: the provider of the MS2 model, and of the RT
  /// and CCS models (LSTMs: they need cuDNN on CUDA and do not run on
  /// CoreML at all).
  struct DeviceRoute
  {
    bool ms2_gpu = false;      ///< try a GPU for MS2 (prefer_gpu)
    bool ms2_coreml = false;   ///< ... through CoreML rather than CUDA
    bool lstm_gpu = false;     ///< try a GPU for RT and CCS
    bool lstm_coreml = false;  ///< ... through CoreML (only when forced: it falls back, with a WARNING)
    std::string summary;       ///< one line: which route and why
    std::string warning;       ///< a soft fallback: reason and fix, one line; empty if none
  };

  /// Library size from which -machine:device auto runs MS2 on CoreML, for
  /// @p threads CPU threads: below it, loading the compiled shapes (15-30 s,
  /// up to 60 s when they are not in the page cache) costs more than the GPU
  /// saves over the CPU. 50,000 precursors per thread: measured on an M4 Max
  /// the warm break-even was about 40 k at 1 thread, 230 k at 4 and 480 k
  /// at 16 (docs/usage.md#apple-gpu);
  /// ODIA_COREML_MIN_PRECURSORS overrides it.
  std::size_t coremlMinPrecursors(unsigned threads);

  /// The route for -machine:device @p device (auto, cpu, cuda, cuda:N,
  /// coreml, mps) on @p caps, for a library of @p precursors. auto: CoreML
  /// for MS2 on a Mac from @p coreml_min precursors on; CUDA for MS2 when the
  /// CUDA provider and a GPU are there, and for RT and CCS only with cuDNN;
  /// otherwise the CPU. @p require (-machine:require_gpu) asks for the GPU
  /// everywhere, so that every way of not getting it fails.
  /// @p ccs: whether CCS is predicted at all (only for wording).
  DeviceRoute chooseDeviceRoute(const std::string& device, bool require, const DeviceCaps& caps,
                                std::size_t precursors, std::size_t coreml_min, bool ccs = true);

  /// Whether an NVIDIA GPU can be visible to this process at all: not when
  /// CUDA_VISIBLE_DEVICES hides every device, and on Linux only with the
  /// NVIDIA kernel driver loaded. Elsewhere assumed, so a failure there still
  /// warns.
  bool nvidiaGpuVisible();

  /// Verbosity of the diagnostics (provider probes, loader errors, search
  /// paths): TOPP -debug, or ODIA_VERBOSE=1 for the library's own users.
  int diagnosticsLevel();
  void setDiagnosticsLevel(int level);
}
