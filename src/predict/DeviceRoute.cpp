// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <odia/DeviceRoute.h>

#include <onnxruntime_cxx_api.h>

#include <atomic>
#include <cstdlib>
#include <filesystem>
#include <string_view>
#include <system_error>
#include <vector>

#if defined(__linux__)
#include <dlfcn.h>
#endif

namespace ODIA
{
  namespace
  {
    std::atomic<int> g_diagnostics{-1};

    /// 5.8 M, 240 k, 950.
    std::string count(std::size_t n)
    {
      if (n >= 1000000)
      { return std::to_string(n / 1000000) + "." + std::to_string(n / 100000 % 10) + " M"; }
      if (n >= 1000) { return std::to_string(n / 1000) + " k"; }
      return std::to_string(n);
    }

#if defined(__linux__)
    /// Directories, off the loader's search path, where a cuDNN 9 for CUDA 12
    /// is commonly installed: a conda environment, a Python environment's
    /// nvidia-cudnn-cu12 wheel, the CUDA toolkit and the system library
    /// directories. Only looked at, never loaded from.
    std::vector<std::filesystem::path> cudnnPlaces()
    {
      namespace fs = std::filesystem;
      std::vector<fs::path> out;
      std::error_code ec;
      const auto wheels = [&](const fs::path& lib) {
        // <prefix>/lib/python3.X/site-packages/nvidia/cudnn/lib
        for (fs::directory_iterator it(lib, ec), end; !ec && it != end; it.increment(ec))
        {
          if (it->path().filename().string().rfind("python3", 0) == 0)
          { out.push_back(it->path() / "site-packages" / "nvidia" / "cudnn" / "lib"); }
        }
        ec.clear();
      };
      if (const char* conda = std::getenv("CONDA_PREFIX"); conda != nullptr && *conda != '\0')
      { out.push_back(fs::path(conda) / "lib"); wheels(fs::path(conda) / "lib"); }
      if (const char* venv = std::getenv("VIRTUAL_ENV"); venv != nullptr && *venv != '\0')
      { wheels(fs::path(venv) / "lib"); }
      if (const char* home = std::getenv("HOME"); home != nullptr && *home != '\0')
      { wheels(fs::path(home) / ".local" / "lib"); }
      for (const char* dir : {"/usr/local/cuda/lib64", "/usr/local/cuda/targets/x86_64-linux/lib",
                              "/usr/lib/x86_64-linux-gnu", "/usr/lib/aarch64-linux-gnu", "/usr/lib64", "/usr/lib"})
      { out.emplace_back(dir); }
      return out;
    }
#endif
  }

  int diagnosticsLevel()
  {
    int level = g_diagnostics.load();
    if (level < 0)
    {
      const char* v = std::getenv("ODIA_VERBOSE");
      level = v != nullptr && *v != '\0' && std::string_view(v) != "0" ? 1 : 0;
      g_diagnostics.store(level);
    }
    return level;
  }

  void setDiagnosticsLevel(int level)
  {
    const char* v = std::getenv("ODIA_VERBOSE");
    const int env = v != nullptr && *v != '\0' && std::string_view(v) != "0" ? 1 : 0;
    g_diagnostics.store(level > env ? level : env);
  }

  bool nvidiaGpuVisible()
  {
    if (const char* visible = std::getenv("CUDA_VISIBLE_DEVICES"); visible != nullptr)
    {
      const std::string_view v(visible);
      if (v.empty() || v == "-1" || v == "none" || v == "NoDevFiles") { return false; }
    }
#if defined(__linux__)
    std::error_code ec;
    return std::filesystem::exists("/proc/driver/nvidia/version", ec);
#else
    return true;
#endif
  }

  DeviceCaps probeDeviceCaps()
  {
    DeviceCaps c;
    std::string providers;
    for (const auto& name : Ort::GetAvailableProviders())
    {
      providers += (providers.empty() ? "" : ", ") + name;
      if (name == "CUDAExecutionProvider") { c.cuda_provider = true; }
      if (name == "CoreMLExecutionProvider") { c.coreml_provider = true; }
    }
    c.diagnostics.push_back("ONNX Runtime providers: " + providers);
#if defined(__APPLE__)
    // CUDA_VISIBLE_DEVICES aside, nvidiaGpuVisible() assumes a GPU off Linux
    // (so that a failure there still warns); a Mac has none.
    c.nvidia_gpu = false;
#else
    c.nvidia_gpu = nvidiaGpuVisible();
#endif
    c.diagnostics.push_back(std::string("NVIDIA GPU visible: ") + (c.nvidia_gpu ? "yes" : "no") +
#if defined(__linux__)
                            " (/proc/driver/nvidia/version, CUDA_VISIBLE_DEVICES)"
#else
                            ""
#endif
    );
#if !defined(__linux__)
    // Elsewhere cuDNN is not looked for: RT and CCS try CUDA, and fall back
    // with a WARNING if it is missing, as before 0.12.
    if (c.cuda_provider && c.nvidia_gpu) { c.cudnn = true; c.cudnn_detail = "cuDNN not checked here"; }
#else
    if (c.cuda_provider && c.nvidia_gpu)
    {
      // cuDNN 9 by its soname (what ONNX Runtime 1.30-1.31 needs: measured
      // on an H100, RT and CCS run on CUDA with only the .so.9 files), else
      // the unversioned name, by the loader's standard search. Kept loaded
      // when it is cuDNN 9: ONNX Runtime's own dlopen then finds this one.
      void* handle = nullptr;
      std::string opened;
      for (const char* name : {"libcudnn.so.9", "libcudnn.so"})
      {
        dlerror();
        handle = dlopen(name, RTLD_NOW | RTLD_LOCAL);
        const char* error = handle == nullptr ? dlerror() : nullptr;
        c.diagnostics.push_back(std::string("dlopen ") + name + ": " +
                                (handle != nullptr ? "loaded" : error != nullptr ? error : "failed"));
        if (handle != nullptr) { opened = name; break; }
      }
      if (handle != nullptr)
      {
        using Version = std::size_t (*)();
        const auto version = reinterpret_cast<Version>(dlsym(handle, "cudnnGetVersion"));
        const std::size_t v = version != nullptr ? version() : 0;
        Dl_info info{};
        const std::string where = version != nullptr && dladdr(reinterpret_cast<void*>(version), &info) != 0 &&
                                  info.dli_fname != nullptr ? info.dli_fname : opened;
        // cuDNN 9 numbers its version major*10000 + minor*100 + patch, 8 major*1000 + ...
        if (v >= 90000 && v < 100000)
        {
          c.cudnn = true;
          c.cudnn_detail = "cuDNN " + std::to_string(v / 10000) + "." + std::to_string(v / 100 % 100) + "." +
                           std::to_string(v % 100);
        }
        else
        {
          c.cudnn_detail = opened + " is not cuDNN 9 (version " + std::to_string(v) + ")";
          dlclose(handle);
        }
        c.diagnostics.push_back(where + ": cudnnGetVersion " + std::to_string(v));
      }
      else
      {
        c.cudnn_detail = "no cuDNN 9 (libcudnn.so.9) on the library path";
        std::error_code ec;
        for (const auto& dir : cudnnPlaces())
        {
          const bool here = std::filesystem::exists(dir / "libcudnn.so.9", ec) || std::filesystem::exists(dir / "libcudnn.so", ec);
          c.diagnostics.push_back("looked for cuDNN in " + dir.string() + ": " + (here ? "found" : "no"));
          if (here && c.cudnn_elsewhere.empty()) { c.cudnn_elsewhere = dir.string(); }
        }
      }
    }
#endif
    return c;
  }

  std::size_t coremlMinPrecursors(unsigned threads)
  {
    if (const char* v = std::getenv("ODIA_COREML_MIN_PRECURSORS"); v != nullptr && *v != '\0')
    {
      char* end = nullptr;
      const long long x = std::strtoll(v, &end, 10);
      if (end != v && x >= 0) { return static_cast<std::size_t>(x); }
    }
    return std::size_t{50000} * (threads > 0 ? threads : 1);
  }

  DeviceRoute chooseDeviceRoute(const std::string& device, bool require, const DeviceCaps& caps,
                                std::size_t precursors, std::size_t coreml_min, bool ccs)
  {
    DeviceRoute r;
    const std::string lstm = ccs ? "RT and CCS" : "RT";
    const bool coreml = device == "coreml" || device == "mps";
    if (device == "cpu")
    {
      r.summary = "prediction on the CPU (-machine:device cpu)";
      return r;
    }
    if (coreml || (require && caps.coreml_provider && device == "auto"))
    {
      // Forced (or required): every model asks for CoreML, and RT and CCS
      // fall back with a WARNING (no LSTM operator), as before 0.12.
      r.ms2_gpu = r.ms2_coreml = r.lstm_gpu = r.lstm_coreml = true;
      r.summary = "MS2 on the Apple GPU (CoreML), " + lstm + " on the CPU (" +
                  (coreml ? "-machine:device " + device : std::string("-machine:require_gpu")) + ")";
      return r;
    }
    if (device != "auto" || require)
    {
      // cuda, cuda:N, or the GPU required: CUDA for every model; a model
      // that cannot have it warns (or, required, fails) as before 0.12.
      r.ms2_gpu = r.lstm_gpu = true;
      r.summary = std::string("prediction on CUDA (") +
                  (device != "auto" ? "-machine:device " + device : std::string("-machine:require_gpu")) + ")";
      return r;
    }
    if (caps.coreml_provider)
    {
      if (precursors >= coreml_min)
      {
        r.ms2_gpu = r.ms2_coreml = true;
        r.summary = "MS2 on the Apple GPU (CoreML), " + lstm + " on the CPU: " + count(precursors) +
                    " precursors, CoreML pays from " + count(coreml_min) + " (-machine:device cpu to stay on the CPU)";
      }
      else
      {
        r.summary = "prediction on the CPU: " + count(precursors) + " precursors, the Apple GPU pays only from " +
                    count(coreml_min) + " (-machine:device coreml to use it anyway)";
      }
      return r;
    }
    if (caps.cuda_provider && caps.nvidia_gpu)
    {
      r.ms2_gpu = true;
      if (caps.cudnn)
      {
        r.lstm_gpu = true;
        r.summary = "prediction on CUDA (" + caps.cudnn_detail + ")";
      }
      else
      {
        r.summary = "MS2 on CUDA, " + lstm + " on the CPU (no cuDNN 9)";
        r.warning = lstm + (ccs ? " run" : " runs") + " on the CPU, not the GPU: " + caps.cudnn_detail +
                    ". Fix: put cuDNN 9 for CUDA 12 (libcudnn.so.9) on LD_LIBRARY_PATH" +
                    (caps.cudnn_elsewhere.empty() ? std::string() : " (one is in " + caps.cudnn_elsewhere + ")") +
                    "; docs/usage.md#gpu";
      }
      return r;
    }
    if (caps.cuda_provider)
    { r.summary = "prediction on the CPU: no NVIDIA GPU is visible"; }
    else if (caps.nvidia_gpu)
    {
      r.summary = "prediction on the CPU: this ONNX Runtime has no CUDA provider, so the NVIDIA GPU is not used "
                  "(docs/usage.md#gpu)";
    }
    else
    { r.summary = "prediction on the CPU: this ONNX Runtime has no GPU provider"; }
    return r;
  }
}
