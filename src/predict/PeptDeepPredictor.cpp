// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <odia/PeptDeepPredictor.h>

#include <odia/Hardware.h>
#include <odia/PeptDeepEncoder.h>
#include <odia/Trace.h>

#include "../library/Parallel.h"
#include "../library/ThreadPool.h"

#include <OpenMS/CHEMISTRY/AASequence.h>

#include <onnxruntime_cxx_api.h>

#include <algorithm>
#include <atomic>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <mutex>
#include <string_view>
#include <system_error>
#include <thread>
#include <limits>
#include <stdexcept>

#if defined(_WIN32)
#include <process.h>
#else
#include <unistd.h>
#endif

namespace
{
  /// The model path in the character type ONNX Runtime's API actually takes.
  ///
  /// Ort::Session takes ORTCHAR_T*, which is wchar_t* on Windows and char*
  /// everywhere else, so passing std::string::c_str() compiles on POSIX and
  /// fails on MSVC with "no overloaded function could convert all the argument
  /// types". std::filesystem::path::c_str() is ORTCHAR_T* on both, and it also
  /// carries a non-ASCII path across correctly, which a naive widening would
  /// not.
  ///
  /// The returned path is a temporary at every call site; its c_str() stays
  /// valid for the whole full-expression, which is exactly as long as the
  /// Ort::Session constructor needs it. Do not "fix" it into a dangling
  /// pointer by hoisting the c_str() out.
  std::filesystem::path ortPath(const std::string& p) { return std::filesystem::path(p); }

  /// True for the CUDA errors that mean "there is no such GPU" rather than
  /// "this GPU cannot run the model": probing devices 0..7 on a one-GPU box,
  /// or any device with CUDA_VISIBLE_DEVICES empty, is the normal CPU case.
  bool noSuchDevice(const std::string& what)
  {
    return what.find("no CUDA-capable device") != std::string::npos ||
           what.find("invalid device ordinal") != std::string::npos ||
           what.find("CUDA driver version is insufficient") != std::string::npos;
  }

  /// Whether an NVIDIA GPU can be visible to this process at all: not when
  /// CUDA_VISIBLE_DEVICES hides every device, and on Linux only with the
  /// NVIDIA kernel driver loaded. Elsewhere assumed, so a failure there still
  /// warns.
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

  /// ONNX Runtime's message without the source path, function signature and
  /// status code it prefixes: "... [ONNXRuntimeError] : 9 : NOT_IMPLEMENTED :
  /// cuDNN is unavailable ...: dlopen failed for libcudnn.so: ..." keeps the
  /// part after the last " : ", which is what a user can act on.
  std::string conciseOrtError(const std::string& what)
  {
    if (what.find("[ONNXRuntimeError]") == std::string::npos) { return what; }
    const auto last = what.rfind(" : ");
    return last == std::string::npos ? what : what.substr(last + 3);
  }

  /// Options every session shares. Generation runs one intra-op thread per
  /// session (sessions, not intra-op threads, scale these models); where a
  /// caller asks for more (0 = ONNX Runtime's default, all cores), the pool's
  /// threads must not spin while they wait for work: spinning cost 29-43 %
  /// more CPU at the same speed (16 sessions x 6 threads, perf report). The
  /// result does not depend on it.
  void configureSession(Ort::SessionOptions& options, int intra_op_threads)
  {
    options.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_ALL);
    if (intra_op_threads > 0) { options.SetIntraOpNumThreads(intra_op_threads); }
    options.AddConfigEntry("session.intra_op.allow_spinning", "0");
    options.AddConfigEntry("session.inter_op.allow_spinning", "0");
  }

  /// Registers the CUDA provider for one device. Through the V2 options,
  /// because the legacy struct cannot switch TF32 off: ONNX Runtime's CUDA
  /// provider enables it by default.
  void appendCuda(Ort::SessionOptions& options, int device, bool tf32)
  {
    const OrtApi& api = Ort::GetApi();
    OrtCUDAProviderOptionsV2* raw = nullptr;
    Ort::ThrowOnError(api.CreateCUDAProviderOptions(&raw));
    const std::unique_ptr<OrtCUDAProviderOptionsV2, void (*)(OrtCUDAProviderOptionsV2*)> cuda(
      raw, [](OrtCUDAProviderOptionsV2* o) { Ort::GetApi().ReleaseCUDAProviderOptions(o); });
    const std::string id = std::to_string(device);
    const char* keys[] = {"device_id", "use_tf32"};
    const char* values[] = {id.c_str(), tf32 ? "1" : "0"};
    Ort::ThrowOnError(api.UpdateCUDAProviderOptions(cuda.get(), keys, values, 2));
    options.AppendExecutionProvider_CUDA_V2(*cuda);
    // Every CUDA session (probes included) otherwise logs, at warning level,
    // that ONNX Runtime keeps the shape operators on the CPU -- by design --
    // and with four sessions per model those lines bury a real WARNING.
    // Errors still print.
    options.SetLogSeverityLevel(ORT_LOGGING_LEVEL_ERROR);
  }

  std::string fallbackFix(const std::string& reason)
  {
    std::string lower = reason;
    std::transform(lower.begin(), lower.end(), lower.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    const char* docs = "https://github.com/okohlbacher/DIALibGen/blob/main/docs/usage.md#gpu";
    if (lower.find("cudnn") != std::string::npos)
    {
      return std::string("install cuDNN 9 for CUDA 12, including the unversioned libcudnn.so "
                         "(apt: libcudnn9-cuda-12 libcudnn9-dev-cuda-12; or pip nvidia-cudnn-cu12 "
                         "plus a libcudnn.so symlink, on LD_LIBRARY_PATH); see ") + docs;
    }
    return std::string("check the NVIDIA driver and the CUDA 12 runtime libraries (cuBLAS, cuDNN 9); see ") + docs;
  }

  /// Numbers predictor calls and predictor instances in the trace (tracing only).
  std::atomic<std::int64_t> g_trace_calls{0};
  std::atomic<std::int64_t> g_trace_instances{0};
  /// Names the temporary optimised-model files apart within one process.
  std::atomic<std::int64_t> g_optimized_files{0};

  long processId()
  {
#if defined(_WIN32)
    return static_cast<long>(_getpid());
#else
    return static_cast<long>(::getpid());
#endif
  }

  /// One float initializer of an ONNX model: name, shape and its raw bytes.
  struct FloatInitializer
  {
    std::string name;
    std::vector<std::int64_t> dims;
    std::string_view raw;
  };

  /// The protobuf wire format, just enough of it to walk an ONNX file.
  struct WireReader
  {
    const char* p;
    const char* end;

    bool varint(std::uint64_t& v)
    {
      v = 0;
      for (int shift = 0; shift < 64 && p < end; shift += 7)
      {
        const auto byte = static_cast<std::uint8_t>(*p++);
        v |= static_cast<std::uint64_t>(byte & 0x7f) << shift;
        if (!(byte & 0x80)) { return true; }
      }
      return false;
    }
    /// Next field: number, wire type and, for length-delimited fields, the payload.
    bool next(std::uint32_t& field, std::uint32_t& type, std::string_view& payload, std::uint64_t& value)
    {
      std::uint64_t key = 0;
      if (p >= end || !varint(key)) { return false; }
      field = static_cast<std::uint32_t>(key >> 3);
      type = static_cast<std::uint32_t>(key & 7);
      switch (type)
      {
        case 0: return varint(value);
        case 1: if (end - p < 8) { return false; } p += 8; return true;
        case 5: if (end - p < 4) { return false; } p += 4; return true;
        case 2:
        {
          std::uint64_t n = 0;
          if (!varint(n) || n > static_cast<std::uint64_t>(end - p)) { return false; }
          payload = std::string_view(p, static_cast<std::size_t>(n));
          p += n;
          return true;
        }
        default: return false;   // groups are not used by ONNX
      }
    }
  };

  /// The float32 initializers of an ONNX model held in raw_data, which is
  /// how every weight matrix of these models is stored. Anything else -- int64
  /// shape constants, external data -- is left to ONNX Runtime. False if the
  /// bytes do not parse as a model.
  bool floatInitializers(std::string_view model, std::vector<FloatInitializer>& out)
  {
    std::uint32_t field = 0, type = 0;
    std::string_view payload;
    std::uint64_t value = 0;
    WireReader m{model.data(), model.data() + model.size()};
    std::string_view graph;
    while (m.p < m.end)
    {
      if (!m.next(field, type, payload, value)) { return false; }
      if (field == 7 && type == 2) { graph = payload; }           // ModelProto.graph
    }
    WireReader g{graph.data(), graph.data() + graph.size()};
    while (g.p < g.end)
    {
      if (!g.next(field, type, payload, value)) { return false; }
      if (field != 5 || type != 2) { continue; }                  // GraphProto.initializer
      FloatInitializer t;
      std::uint64_t data_type = 0, location = 0;
      WireReader r{payload.data(), payload.data() + payload.size()};
      std::string_view inner;
      while (r.p < r.end)
      {
        if (!r.next(field, type, inner, value)) { return false; }
        if (field == 1 && type == 0) { t.dims.push_back(static_cast<std::int64_t>(value)); }
        else if (field == 1 && type == 2)                         // packed dims
        {
          WireReader d{inner.data(), inner.data() + inner.size()};
          while (d.p < d.end)
          {
            if (!d.varint(value)) { return false; }
            t.dims.push_back(static_cast<std::int64_t>(value));
          }
        }
        else if (field == 2 && type == 0) { data_type = value; }
        else if (field == 8 && type == 2) { t.name = std::string(inner); }
        else if (field == 9 && type == 2) { t.raw = inner; }
        else if (field == 14 && type == 0) { location = value; }
      }
      std::uint64_t elements = 1;
      for (const auto d : t.dims) { elements *= static_cast<std::uint64_t>(d); }
      if (data_type == 1 && location == 0 && !t.name.empty() &&
          t.raw.size() == elements * sizeof(float))
      { out.push_back(std::move(t)); }
    }
    return true;
  }
}



namespace ODIA
{
  namespace
  {
#if ORT_API_VERSION < 23
    /// ONNX Runtime's pre-packed-weights container through the C API: before
    /// 1.23 the C++ API has no wrapper class (Ort::PrepackedWeightsContainer),
    /// and conda-forge's macOS x64 build is still 1.22.
    class PrepackedWeights
    {
    public:
      PrepackedWeights() = default;
      PrepackedWeights(const PrepackedWeights&) = delete;
      PrepackedWeights& operator=(const PrepackedWeights&) = delete;
      ~PrepackedWeights() { reset(); }
      void create()
      {
        reset();
        Ort::ThrowOnError(Ort::GetApi().CreatePrepackedWeightsContainer(&container_));
      }
      void reset()
      {
        if (container_ != nullptr) { Ort::GetApi().ReleasePrepackedWeightsContainer(container_); }
        container_ = nullptr;
      }
      OrtPrepackedWeightsContainer* get() const { return container_; }
      explicit operator bool() const { return container_ != nullptr; }
    private:
      OrtPrepackedWeightsContainer* container_ = nullptr;
    };
#endif

    /// Stamps a chunk-level trace scope with its model, session and rows.
    void tag(trace::Scope& scope, const char* kind, std::size_t session, std::size_t rows)
    {
      if (!scope.active()) { return; }
      scope.fields().model = kind;
      scope.fields().session = static_cast<std::int64_t>(session);
      scope.fields().rows = static_cast<std::int64_t>(rows);
    }
  }

  struct PeptDeepPredictor::Impl
  {
    /// ONNX Runtime keeps ONE environment per process and reference-counts
    /// it, so every predictor's Env is the same object underneath.
    Ort::Env env{ORT_LOGGING_LEVEL_WARNING, "odia"};
    /// The weights every CPU session reads, held once (shareWeights_).
    /// Declared before the options and sessions that point into them, so they
    /// are destroyed after.
    std::vector<std::vector<float>> shared_buffers;
    std::vector<Ort::Value> shared_values;
#if ORT_API_VERSION >= 23
    Ort::PrepackedWeightsContainer prepacked{nullptr};
#else
    PrepackedWeights prepacked;
#endif
    std::size_t shared_bytes = 0;
    Ort::SessionOptions options;
    std::unique_ptr<Ort::Session> session;
    /// Sessions 1..n. Session 0 is `session` above, so a single-session
    /// predictor allocates nothing extra and behaves exactly as before.
    std::vector<std::unique_ptr<Ort::Session>> replicas;

    Ort::Session& at(std::size_t i)
    {
      return i == 0 ? *session : *replicas[i - 1];
    }
    Ort::MemoryInfo memory{Ort::MemoryInfo::CreateCpu(OrtArenaAllocator, OrtMemTypeDefault)};
    Provider provider = Provider::CPU;
    int device = -1;
    std::vector<std::string> input_names;
    std::vector<std::string> output_names;
    std::string fallback;   ///< see fallbackReason()
    std::vector<int> session_devices;   ///< CUDA device per session; empty on CPU
    /// The names as Run takes them, built once: every batch passes the same.
    std::vector<const char*> input_ptrs;
    const char* output_ptr = nullptr;
    /// RT, MS2 or CCS, from the inputs; names the model in trace events and
    /// picks the CPU chunking (forEachChunk_).
    const char* kind = "-";
    /// ODIA_ORT_PROFILE: every session writes ONNX Runtime's per-op profile.
    bool ort_profiling = false;
  };

  PeptDeepPredictor::PeptDeepPredictor(const std::string& model_path,
                                       bool prefer_gpu, int intra_op_threads,
                                       int sessions, int gpu_device,
                                       const GpuOptions& gpu)
    : impl_(std::make_unique<Impl>())
  {
    configureSession(impl_->options, intra_op_threads);

    // ODIA_ORT_PROFILE=<prefix>: ONNX Runtime's own per-operator profile, one
    // JSON per session. Each session needs its own prefix -- ONNX Runtime
    // appends only a time stamp to the second, so sessions created in the same
    // second would overwrite each other's file. Opt-in and slow; use it on
    // small inputs.
    const char* ort_profile = std::getenv("ODIA_ORT_PROFILE");
    impl_->ort_profiling = ort_profile != nullptr && *ort_profile != '\0';
    if (impl_->ort_profiling)
    {
      // ONNX Runtime does not create the directory and reports no error when
      // it cannot open the file.
      std::error_code ec;
      const auto parent = std::filesystem::path(ort_profile).parent_path();
      if (!parent.empty()) { std::filesystem::create_directories(parent, ec); }
    }
    const std::int64_t instance = impl_->ort_profiling ? g_trace_instances++ : 0;
    // @p optimized: where ONNX Runtime should write the optimised graph
    // (session 0 only, when the weights are to be shared; see shareWeights_).
    const auto makeSession = [&](std::size_t index, const Ort::SessionOptions& options,
                                 const std::filesystem::path& optimized = {}) {
      trace::Scope created("session_create", 3);
      created.fields().session = static_cast<std::int64_t>(index);
      const auto open = [&](const Ort::SessionOptions& o) {
        return impl_->prepacked
#if ORT_API_VERSION >= 23
          ? std::make_unique<Ort::Session>(impl_->env, ortPath(model_path).c_str(), o, impl_->prepacked)
#else
          ? std::make_unique<Ort::Session>(impl_->env, ortPath(model_path).c_str(), o, impl_->prepacked.get())
#endif
          : std::make_unique<Ort::Session>(impl_->env, ortPath(model_path).c_str(), o);
      };
      if (!impl_->ort_profiling && optimized.empty()) { return open(options); }
      Ort::SessionOptions custom = options.Clone();
      if (!optimized.empty())
      {
        custom.SetOptimizedModelFilePath(optimized.c_str());
        // ONNX Runtime warns that an ENABLE_ALL graph saved as ONNX may be
        // hardware specific. True, and irrelevant: the file is read back by
        // this process only, for its weights, and then deleted.
        custom.SetLogSeverityLevel(ORT_LOGGING_LEVEL_ERROR);
      }
      if (impl_->ort_profiling)
      {
        const std::string prefix = std::string(ort_profile) + "_" +
          std::filesystem::path(model_path).stem().string() + "_i" + std::to_string(instance) +
          "_s" + std::to_string(index);
        custom.EnableProfiling(ortPath(prefix).c_str());
      }
      return open(custom);
    };

    // Several CPU sessions share one copy of the weights (shareWeights_), so
    // session 0 records the graph ONNX Runtime actually runs.
    std::filesystem::path optimized;
    const auto firstCpuSession = [&]() {
      if (sessions > 1)
      {
        try
        {
          // The pid keeps two concurrent runs on one machine apart, the
          // counter two predictors of one run.
          optimized = std::filesystem::temp_directory_path() /
            ("odia-" + std::filesystem::path(model_path).stem().string() + "-" +
             std::to_string(processId()) + "-" + std::to_string(g_optimized_files++) + ".onnx");
          return makeSession(0, impl_->options, optimized);
        }
        catch (const std::exception&)
        {
          // No writable temporary directory, say: share nothing, and let the
          // plain attempt below report a model that cannot be opened at all.
          std::error_code ec;
          if (!optimized.empty()) { std::filesystem::remove(optimized, ec); }
          optimized.clear();
        }
      }
      return makeSession(0, impl_->options);
    };

    // Why CUDA was not used, and whether that deserves a warning: a GPU that
    // exists but cannot run the model (cuDNN missing, out of memory) does; a
    // machine without one does not -- the same binary serves both.
    std::string reason;
    bool gpu_present = false;
    // CUDA: the session of each device that worked, in device order. The
    // probe is a real session built with the final options, so it is kept:
    // the k-th one is session k, the first session on that device.
    std::vector<std::unique_ptr<Ort::Session>> probed;
    if (prefer_gpu)
    {
      // Ask before attempting. Registering CUDA blindly on a CPU-only build
      // still falls back correctly, but ONNX Runtime logs a failure to load
      // libonnxruntime_providers_shared.so at error level first -- alarming
      // output for the expected case on a machine with no GPU, which is the
      // machine this is developed on.
      bool cuda_available = false;
      for (const auto& name : Ort::GetAvailableProviders())
      {
        if (name == "CUDAExecutionProvider") { cuda_available = true; break; }
      }

      if (!cuda_available) { reason = "this ONNX Runtime build has no CUDA provider"; }
      else
      {
        // PROBE the devices rather than trusting device 0. On this project's
        // GPU nodes device 0 is dead, and AppendExecutionProvider_CUDA happily
        // accepts it -- the failure only surfaces when the session is built.
        // A list (GpuOptions::devices) is probed whole and every device that
        // works is used; otherwise the first device that works.
        std::vector<int> candidates = gpu.devices;
        const bool listed = !candidates.empty();
        if (!listed && gpu_device >= 0) { candidates.push_back(gpu_device); }
        if (!listed && gpu_device < 0) { for (int dev = 0; dev <= 7; ++dev) { candidates.push_back(dev); } }
        for (const int dev : candidates)
        {
          trace::Scope probe_step("session_probe", 3);
          probe_step.fields().idx = dev;
          try
          {
            Ort::SessionOptions options;
            configureSession(options, intra_op_threads);
            appendCuda(options, dev, gpu.tf32);
            probed.push_back(makeSession(probed.size(), options));   // the real check
            impl_->session_devices.push_back(dev);
            if (!listed) { break; }
          }
          catch (const Ort::Exception& e)
          {
            // Keep the first real failure: it is the one the user can fix. A
            // device the user listed by hand is never "absent by design".
            // A CUDA build of ONNX Runtime on a machine without an NVIDIA
            // driver fails to load its CUDA libraries: no GPU, the expected
            // CPU case, not a GPU that cannot be used.
            if ((listed || (!noSuchDevice(e.what()) && nvidiaGpuVisible())) && !gpu_present)
            {
              gpu_present = true;
              reason = "CUDA device " + std::to_string(dev) + ": " + conciseOrtError(e.what());
            }
            else if (reason.empty())
            { reason = "CUDA device " + std::to_string(dev) + ": " + conciseOrtError(e.what()); }
          }
        }
        if (!impl_->session_devices.empty())
        {
          appendCuda(impl_->options, impl_->session_devices.front(), gpu.tf32);
          impl_->provider = Provider::CUDA;
          impl_->device = impl_->session_devices.front();
          impl_->session = std::move(probed.front());
        }
        if (impl_->provider != Provider::CUDA && reason.empty())
        { reason = "no usable CUDA device is visible"; }
      }
    }

    // On CUDA session 0 is the first probe, above: registering the provider
    // can succeed while creating a session with it fails, and the probe
    // covers both.
    if (impl_->provider != Provider::CUDA)
    {
      try { impl_->session = firstCpuSession(); }
      catch (const Ort::Exception& e)
      { throw std::runtime_error("cannot open model " + model_path + ": " + e.what()); }
    }

    if (prefer_gpu && impl_->provider == Provider::CPU)
    {
      const std::string model = std::filesystem::path(model_path).filename().string();
      // Requiring the GPU turns every way of ending up on the CPU into an
      // error, the expected no-GPU case included: whoever asked knows.
      if (gpu.require)
      {
        throw std::runtime_error("the GPU is required (-machine:require_gpu) but " + model +
                                 " cannot run on CUDA: " + reason + ". Fix: " + fallbackFix(reason));
      }
      // A silent fall back to CPU turns a 2.5-minute library into a
      // 6.5-minute one (H100, whole proteome) with nothing in the log to
      // say why: the caller prints and records this.
      if (gpu_present)
      { impl_->fallback = model + " runs on the CPU although CUDA was requested: " + reason + ". Fix: " + fallbackFix(reason); }
    }

    if (impl_->provider == Provider::CUDA && gpu_present && !gpu.devices.empty())
    {
      // Some listed devices work, some do not: run on the good ones, but say
      // so. (Probing without a list skips a dead device silently, as before.)
      const std::string model = std::filesystem::path(model_path).filename().string();
      if (gpu.require)
      { throw std::runtime_error("the GPU is required (-machine:require_gpu) but " + model + " cannot use " + reason); }
      impl_->fallback = model + " skips an unusable GPU: " + reason;
    }

    Ort::AllocatorWithDefaultOptions allocator;
    for (std::size_t i = 0; i < impl_->session->GetInputCount(); ++i)
    {
      impl_->input_names.emplace_back(
        impl_->session->GetInputNameAllocated(i, allocator).get());
    }
    for (std::size_t i = 0; i < impl_->session->GetOutputCount(); ++i)
    {
      impl_->output_names.emplace_back(
        impl_->session->GetOutputNameAllocated(i, allocator).get());
    }

    // Validate by name and count, not just count. Inputs are bound positionally,
    // and the CCS and MS2 models take three and five inputs respectively -- both
    // previously constructed cleanly here and then failed inside Run with a
    // red "Missing Input" from ONNX Runtime, which is exactly the alarming
    // output the provider probing was added to avoid.
    // Two inputs is the RT model, five is MS2. Anything else is neither.
    const bool looks_like_rt =
      impl_->input_names.size() == 2 &&
      (impl_->input_names[0] == "input_sequences" || impl_->input_names[0] == "aa_indices") &&
      impl_->input_names[1] == "mod_x";
    const bool looks_like_ms2 =
      impl_->input_names.size() == 5 &&
      (impl_->input_names[0] == "input_sequences" || impl_->input_names[0] == "aa_indices") &&
      impl_->input_names[1] == "mod_x" && impl_->input_names[2] == "charges" &&
      impl_->input_names[3] == "nce" && impl_->input_names[4] == "instrument_indices";
    const bool looks_like_ccs =
      impl_->input_names.size() == 3 &&
      (impl_->input_names[0] == "input_sequences" || impl_->input_names[0] == "aa_indices") &&
      impl_->input_names[1] == "mod_x" && impl_->input_names[2] == "charges";
    if (!looks_like_rt && !looks_like_ms2 && !looks_like_ccs)
    {
      std::error_code ec;
      if (!optimized.empty()) { std::filesystem::remove(optimized, ec); }
      std::string names;
      for (const auto& n : impl_->input_names) { names += (names.empty() ? "" : ", ") + n; }
      throw std::runtime_error(
        "model " + model_path + " is not a PeptDeep retention-time, CCS or MS2 model: "
        "found inputs (" + names + ")");
    }
    impl_->kind = looks_like_rt ? "RT" : looks_like_ms2 ? "MS2" : "CCS";
    for (const auto& n : impl_->input_names) { impl_->input_ptrs.push_back(n.c_str()); }
    impl_->output_ptr = impl_->output_names.at(0).c_str();

    // Replicas only after the model has been accepted, so a bad path or a
    // wrong model reports once rather than N times.
    //
    // On CUDA, sessions_per_device sessions per device, assigned round robin
    // (session k on device k mod n). Each session has its own CUDA stream, so
    // one session's host-side encoding and copies overlap another's kernels:
    // on an H100, 2-4 sessions gave RT and CCS +30-45 % and 8 was slower. The
    // output does not depend on the count: chunks are cut before dispatch.
    //
    // The probes are sessions 0 to n-1. The others are built concurrently,
    // one builder per device creating that device's sessions in order: CUDA
    // session creation (cuDNN, the device arena) costs about a second each,
    // and one device's builds do not wait for another's.
    if (impl_->provider == Provider::CUDA)
    {
      const std::vector<int> devices = impl_->session_devices;
      const std::size_t n = devices.size();
      impl_->session_devices.clear();
      const std::size_t total = n * static_cast<std::size_t>(std::max(1, gpu.sessions_per_device));
      for (std::size_t k = 0; k < total; ++k) { impl_->session_devices.push_back(devices[k % n]); }
      impl_->replicas.resize(total - 1);
      for (std::size_t k = 1; k < n; ++k) { impl_->replicas[k - 1] = std::move(probed[k]); }
      probed.clear();
      if (total > n)
      {
        std::atomic<std::size_t> next{0};
        std::mutex error_mutex;
        std::exception_ptr error;
        std::size_t error_session = total;
        const auto build = [&]() {
          for (std::size_t d = next++; d < n; d = next++)
          {
            for (std::size_t k = n + d; k < total; k += n)
            {
              try
              {
                Ort::SessionOptions options;
                configureSession(options, intra_op_threads);
                appendCuda(options, devices[d], gpu.tf32);
                impl_->replicas[k - 1] = makeSession(k, options);
              }
              catch (...)
              {
                // The lowest failing session is reported, as one builder would.
                const std::lock_guard<std::mutex> lock(error_mutex);
                if (k < error_session) { error_session = k; error = std::current_exception(); }
                break;
              }
            }
          }
        };
        parallel::team(n, [&](std::size_t) { build(); });
        if (error) { std::rethrow_exception(error); }
      }
    }
    else
    {
      // Replicas exist on CPU to work around intra-op parallelism scaling
      // badly with these small recurrent networks.
      try
      {
        // Session 0 again, now reading the shared weights, so that no session
        // keeps a private copy.
        if (!optimized.empty() && shareWeights_(optimized)) { impl_->session = makeSession(0, impl_->options); }
      }
      catch (const Ort::Exception&)
      {
        // Not sharing is slower to build and larger, never wrong: start the
        // options over without the shared weights.
#if ORT_API_VERSION >= 23
        impl_->prepacked = Ort::PrepackedWeightsContainer{nullptr};
#else
        impl_->prepacked.reset();
#endif
        impl_->options = Ort::SessionOptions{};
        configureSession(impl_->options, intra_op_threads);
        impl_->shared_values.clear();
        impl_->shared_buffers.clear();
        impl_->shared_bytes = 0;
      }
      if (sessions > 1)
      {
        // Built concurrently from the same read-only options and shared
        // weights (ONNX Runtime's prepacked-weights container takes its own
        // lock): 96 MS2 sessions one after another took 3.6 s, with up to 32
        // builders 1.8 s. Each slot is written by one builder only, and the
        // first failure is rethrown.
        impl_->replicas.resize(static_cast<std::size_t>(sessions - 1));
        const std::size_t builders = std::min<std::size_t>(impl_->replicas.size(), 32);
        std::atomic<std::size_t> next{0};
        std::mutex error_mutex;
        std::exception_ptr error;
        const auto build = [&]() {
          for (std::size_t i = next++; i < impl_->replicas.size(); i = next++)
          {
            try { impl_->replicas[i] = makeSession(i + 1, impl_->options); }
            catch (...)
            {
              const std::lock_guard<std::mutex> lock(error_mutex);
              if (!error) { error = std::current_exception(); }
              return;
            }
          }
        };
        parallel::team(builders, [&](std::size_t) { build(); });
        if (error)
        {
          std::error_code ec;
          if (!optimized.empty()) { std::filesystem::remove(optimized, ec); }
          std::rethrow_exception(error);
        }
      }
    }
    if (!optimized.empty())
    {
      std::error_code ec;
      std::filesystem::remove(optimized, ec);
    }
    if (trace::enabled())
    {
      trace::note(std::string("ort_version"), Ort::GetVersionString());
      trace::note(std::string("predictor.") + impl_->kind + ".model", model_path);
      trace::note(std::string("predictor.") + impl_->kind + ".sessions", std::to_string(sessionCount()));
      trace::note(std::string("predictor.") + impl_->kind + ".shared_weights",
                  std::to_string(impl_->shared_values.size()) + " tensors, " +
                  std::to_string(impl_->shared_bytes) + " bytes");
    }
  }

  bool PeptDeepPredictor::shareWeights_(const std::filesystem::path& optimized)
  {
    // Every session would otherwise hold its own copy of the weights and of
    // their pre-packed (GEMM-blocked) form: sixteen identical copies of
    // read-only data, competing for the same caches. ONNX Runtime can share
    // both, but only for initializers handed to it with AddInitializer, and it
    // substitutes those BY NAME after optimising the graph. The weights are
    // therefore taken from the graph session 0 actually optimised, not from the
    // file on disk -- an optimiser that rewrote a weight under its own name
    // would otherwise be silently undone. Any failure leaves the sessions as
    // they were: unshared, and correct.
    std::string bytes;
    {
      std::ifstream in(optimized, std::ios::binary);
      if (!in) { return false; }
      bytes.assign(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
    }
    std::vector<FloatInitializer> weights;
    if (!floatInitializers(bytes, weights) || weights.empty()) { return false; }

    const auto cpu = Ort::MemoryInfo::CreateCpu(OrtDeviceAllocator, OrtMemTypeDefault);
    impl_->shared_buffers.reserve(weights.size());
    impl_->shared_values.reserve(weights.size());
    for (const auto& w : weights)
    {
      // Copied out of the file buffer: raw_data has no alignment guarantee.
      auto& buffer = impl_->shared_buffers.emplace_back(w.raw.size() / sizeof(float));
      std::memcpy(buffer.data(), w.raw.data(), w.raw.size());
      impl_->shared_values.push_back(Ort::Value::CreateTensor<float>(
        cpu, buffer.data(), buffer.size(), w.dims.data(), w.dims.size()));
      impl_->options.AddInitializer(w.name.c_str(), impl_->shared_values.back());
      impl_->shared_bytes += w.raw.size();
    }
#if ORT_API_VERSION >= 23
    impl_->prepacked = Ort::PrepackedWeightsContainer();
#else
    impl_->prepacked.create();
#endif
    return true;
  }



  int PeptDeepPredictor::device() const { return impl_->device; }

  const std::string& PeptDeepPredictor::fallbackReason() const { return impl_->fallback; }

  const std::vector<int>& PeptDeepPredictor::sessionDevices() const { return impl_->session_devices; }

  std::size_t PeptDeepPredictor::sessionCount() const
  {
    return impl_->replicas.size() + 1;
  }

  void PeptDeepPredictor::forEachChunk_(
    const std::vector<OpenMS::AASequence>& peptides,
    const std::function<void(std::size_t, const std::vector<std::size_t>&)>& body)
  {
    // The chunk list is built up front and in the serial order, so which
    // session runs a chunk cannot change how the input was split -- and the
    // result can depend on the batch size at the last bit (see MAX_BATCH_ROWS).
    trace::Scope chunking("chunking", 3);
    std::vector<std::vector<std::size_t>> chunks;
    // On the CPU each MAX_BATCH_ROWS chunk is cut further, into pieces of
    // about CPU_CHUNK_ROWS. One call holds a few dozen length groups, and a
    // call ends only when its largest chunk does: with 2048-row chunks, 96
    // sessions ran 13 at a time on average. The cut cannot change a result:
    // on the CPU the MS2 model is bit-identical for every batch size, a single
    // row included, so MS2 is simply cut; the LSTM models (RT, CCS) are
    // bit-identical for every batch of two or more rows but not for a single
    // row, so their pieces never leave a one-row remainder -- it joins the
    // piece before it -- and a chunk that was a single row stays one. The GPU
    // keeps whole chunks: it wants large batches, and its result was measured
    // with them.
    const bool cpu = impl_->provider == Provider::CPU;
    const bool is_ms2 = std::string_view(impl_->kind) == "MS2";
    for (const auto& whole_group : PeptDeepEncoder::groupByLength(peptides))
    {
      for (std::size_t offset = 0; offset < whole_group.size(); offset += MAX_BATCH_ROWS)
      {
        const std::size_t end = std::min(offset + MAX_BATCH_ROWS, whole_group.size());
        if (!cpu || end - offset <= CPU_CHUNK_ROWS)
        {
          chunks.emplace_back(whole_group.begin() + static_cast<std::ptrdiff_t>(offset),
                              whole_group.begin() + static_cast<std::ptrdiff_t>(end));
          continue;
        }
        for (std::size_t b = offset; b < end;)
        {
          std::size_t e = std::min(b + CPU_CHUNK_ROWS, end);
          if (!is_ms2 && end - e == 1) { e = end; }
          chunks.emplace_back(whole_group.begin() + static_cast<std::ptrdiff_t>(b),
                              whole_group.begin() + static_cast<std::ptrdiff_t>(e));
          b = e;
        }
      }
    }
    // Costliest chunks first (longest-processing-time order): a chunk that
    // starts last is what leaves sessions idle at the end of a call. The cost
    // model is rows x (6 + residues), from the measured per-row MS2 cost of
    // 0.82-2.90 ms for 7-30 residues. Order only -- every chunk keeps its rows,
    // so no result can change. Longest-peptide-first was measured 11 % slower.
    if (cpu)
    {
      const auto cost = [&](const std::vector<std::size_t>& c) {
        return c.size() * (6 + peptides[c.front()].size());
      };
      std::stable_sort(chunks.begin(), chunks.end(),
                       [&](const auto& a, const auto& b) { return cost(a) > cost(b); });
    }
    // Tracing: one "dispatch" event per call, one "worker" per thread and one
    // "chunk" per batch, each stamped with the call number, so that queue
    // wait, tail idle and per-session load can be read off afterwards.
    const bool tracing = chunking.active();
    const std::int64_t call = tracing ? g_trace_calls++ : -1;
    if (tracing)
    {
      chunking.fields().model = impl_->kind;
      chunking.fields().call = call;
      chunking.fields().rows = static_cast<std::int64_t>(peptides.size());
      chunking.fields().idx = static_cast<std::int64_t>(chunks.size());
      chunking.end();
    }
    if (chunks.empty()) { return; }

    const std::size_t workers = std::min(sessionCount(), chunks.size());
    const double dispatch_t0 = trace::now();
    const double dispatch_cpu0 = tracing ? trace::threadCpu() : 0.0;
    const auto runChunk = [&](std::size_t session, std::size_t i, std::int64_t& rows_done)
    {
      if (!tracing) { body(session, chunks[i]); return; }
      const double t0 = trace::now();
      const double cpu0 = trace::threadCpu();
      body(session, chunks[i]);
      trace::Fields f;
      f.model = impl_->kind;
      f.call = call;
      f.idx = static_cast<std::int64_t>(i);
      f.session = static_cast<std::int64_t>(session);
      f.rows = static_cast<std::int64_t>(chunks[i].size());
      f.len = static_cast<std::int64_t>(peptides[chunks[i].front()].size());
      f.wait_s = t0 - dispatch_t0;
      trace::record("chunk", 3, t0, trace::now(), trace::threadCpu() - cpu0, f);
      rows_done += f.rows;
    };
    const auto recordWorker = [&](std::size_t session, double t0, double cpu0,
                                  std::int64_t chunks_done, std::int64_t rows_done)
    {
      trace::Fields f;
      f.model = impl_->kind;
      f.call = call;
      f.idx = chunks_done;
      f.session = static_cast<std::int64_t>(session);
      f.rows = rows_done;
      f.wait_s = t0 - dispatch_t0;
      trace::record("worker", 3, t0, trace::now(), trace::threadCpu() - cpu0, f);
    };
    const auto recordDispatch = [&]()
    {
      trace::Fields f;
      f.model = impl_->kind;
      f.call = call;
      f.idx = static_cast<std::int64_t>(chunks.size());
      f.session = static_cast<std::int64_t>(workers);
      f.rows = static_cast<std::int64_t>(peptides.size());
      trace::record("dispatch", 3, dispatch_t0, trace::now(), trace::threadCpu() - dispatch_cpu0, f);
    };

    if (workers <= 1)
    {
      std::int64_t rows_done = 0;
      for (std::size_t i = 0; i < chunks.size(); ++i) { runChunk(0, i, rows_done); }
      if (tracing)
      {
        recordWorker(0, dispatch_t0, dispatch_cpu0, static_cast<std::int64_t>(chunks.size()), rows_done);
        recordDispatch();
      }
      return;
    }

    std::atomic<std::size_t> next{0};
    // Checked once per chunk by every worker: a flag, not a lock.
    std::atomic<bool> failed{false};
    std::mutex first_error_mutex;
    std::exception_ptr first_error;

    const auto worker = [&](std::size_t session)
    {
      const double w0 = trace::now();
      const double wcpu0 = tracing ? trace::threadCpu() : 0.0;
      if (tracing && session > 0) { trace::nameWorker(impl_->kind, session); }
      std::int64_t chunks_done = 0, rows_done = 0;
      for (std::size_t i = next++; i < chunks.size(); i = next++)
      {
        // One thread's failure must stop the others rather than let them run
        // on for the rest of a two-million-peptide proteome. The first
        // exception is kept and rethrown on the calling thread.
        if (failed.load(std::memory_order_relaxed)) { break; }
        try
        {
          runChunk(session, i, rows_done);
          ++chunks_done;
        }
        catch (...)
        {
          const std::lock_guard<std::mutex> lock(first_error_mutex);
          if (!first_error) { first_error = std::current_exception(); }
          failed = true;
          break;
        }
      }
      if (tracing) { recordWorker(session, w0, wcpu0, chunks_done, rows_done); }
    };

    // Each worker owns one session; a worker that cannot be started leaves
    // its session idle and the others take its chunks. Session k runs on
    // the same pool worker call after call (the lowest free ones), so its
    // arena stays on that worker's NUMA node.
    parallel::team(workers, worker, true);
    if (tracing) { recordDispatch(); }

    if (first_error) { std::rethrow_exception(first_error); }
  }

  PeptDeepPredictor::~PeptDeepPredictor()
  {
    if (!impl_ || !impl_->ort_profiling) { return; }
    // End explicitly rather than at session destruction, so the file names
    // can be reported.
    Ort::AllocatorWithDefaultOptions allocator;
    for (std::size_t i = 0; i < sessionCount(); ++i)
    {
      try
      {
        const auto file = impl_->at(i).EndProfilingAllocated(allocator);
        const std::string name = file.get();
        std::cerr << "[trace] ONNX Runtime profile: " << name << "\n";
        trace::note(std::string("ort_profile.") + impl_->kind + ".s" + std::to_string(i), name);
      }
      catch (const std::exception&) {}
    }
  }

  PeptDeepPredictor::Provider PeptDeepPredictor::provider() const
  {
    return impl_->provider;
  }

  void PeptDeepPredictor::runBatch_(std::size_t session,
                                    const PeptDeepEncoder::Batch& batch,
                                    const std::vector<std::size_t>& group,
                                    std::vector<float>& out)
  {
    const std::int64_t rows = static_cast<std::int64_t>(batch.rows);
    const std::int64_t length = static_cast<std::int64_t>(batch.sequence_length);
    const std::int64_t width = static_cast<std::int64_t>(PEPTDEEP_MOD_ELEMENTS.size());

    std::array<std::int64_t, 2> aa_shape{rows, length};
    std::array<std::int64_t, 3> mod_shape{rows, length, width};

    std::vector<Ort::Value> inputs;
    inputs.push_back(Ort::Value::CreateTensor<std::int64_t>(
      impl_->memory, const_cast<std::int64_t*>(batch.aa_indices.data()),
      batch.aa_indices.size(), aa_shape.data(), aa_shape.size()));
    inputs.push_back(Ort::Value::CreateTensor<float>(
      impl_->memory, const_cast<float*>(batch.mod_x.data()),
      batch.mod_x.size(), mod_shape.data(), mod_shape.size()));

    trace::Scope running("ort_run", 4);
    tag(running, impl_->kind, session, group.size());
    auto results = impl_->at(session).Run(Ort::RunOptions{nullptr}, impl_->input_ptrs.data(),
                                       inputs.data(), inputs.size(),
                                       &impl_->output_ptr, 1);
    running.end();
    trace::Scope unpacking("unpack", 4);
    tag(unpacking, impl_->kind, session, group.size());

    // Validate what came back rather than trusting the shape. The MS2 model
    // returns [batch, seq_len-3, 8]; read blindly, the same loop would hand back
    // the first few floats of its b_z1 channel as retention times.
    const auto info = results.front().GetTensorTypeAndShapeInfo();
    if (info.GetElementType() != ONNX_TENSOR_ELEMENT_DATA_TYPE_FLOAT)
    {
      throw std::runtime_error("model output is not float32");
    }
    if (info.GetElementCount() != group.size() || info.GetShape().size() != 1)
    {
      throw std::runtime_error(
        "model returned " + std::to_string(info.GetElementCount()) + " values for " +
        std::to_string(group.size()) + " peptides; this does not look like a "
        "retention-time model");
    }

    const float* values = results.front().GetTensorData<float>();
    for (std::size_t i = 0; i < group.size(); ++i)
    {
      // Back to the caller's order: grouping reorders, and a caller that got
      // predictions silently permuted would have no way to notice.
      if (!std::isfinite(values[i]))
      { throw std::runtime_error("model returned a non-finite prediction for peptide index " + std::to_string(group[i])); }
      out[group[i]] = values[i];
    }
  }

  std::vector<PeptDeepPredictor::Spectrum>
  PeptDeepPredictor::predictMS2(const std::vector<OpenMS::AASequence>& peptides,
                                const std::vector<int>& charges, float nce,
                                const std::string& instrument,
                                std::vector<Failure>* failures)
  {
    if (charges.size() != peptides.size())
    {
      throw std::invalid_argument("one charge per peptide is required");
    }
    if (impl_->input_names.size() != 5)
    {
      throw std::runtime_error("this is not a PeptDeep MS2 model: it takes " +
                               std::to_string(impl_->input_names.size()) +
                               " inputs, not 5");
    }

    // NCE is one value for the whole call, so a bad one is not a per-peptide
    // failure and must not be reported as 2048 of them. It is also silent: the
    // model returns a plausible spectrum for any float it is given.
    if (!std::isfinite(nce) || nce <= 0.0f || nce > 100.0f)
    {
      throw std::invalid_argument("collision energy " + std::to_string(nce) +
                                  " is not a usable NCE");
    }

    const auto canonical_instrument = PeptDeepEncoder::canonicalInstrument(instrument);
    if (canonical_instrument.empty()) { throw std::invalid_argument("unknown instrument: " + instrument); }
    std::vector<Spectrum> out(peptides.size());
    std::mutex failure_mutex;

    forEachChunk_(peptides, [&](std::size_t session, const std::vector<std::size_t>& group)
    {
      trace::Scope encoding("encode", 4);
      tag(encoding, impl_->kind, session, group.size());
      PeptDeepEncoder::Batch batch;
      try
      {
        // Encoded where the rows lie, without copying them out first.
        batch = PeptDeepEncoder::encode(peptides, group, charges, nce, canonical_instrument);
      }
      catch (const std::exception&)
      {
        encoding.end();
        // As in predictRT: one unencodable peptide, or one impossible charge,
        // must not discard its whole chunk and take 2047 good peptides with it.
        // Retry singly so only the offending rows are recorded -- reporting the
        // whole chunk also gave every peptide a reason naming another one's
        // modification.
        if (failures == nullptr) { throw; }
        for (std::size_t i = 0; i < group.size(); ++i)
        {
          try
          {
            const auto single = PeptDeepEncoder::encode(
              peptides, std::vector<std::size_t>{group[i]}, charges, nce, canonical_instrument);
            runMS2Batch_(session, single, {group[i]}, out);
          }
          catch (const std::exception& inner)
          {
            const std::lock_guard<std::mutex> lock(failure_mutex);
            failures->push_back({group[i], inner.what()});
          }
        }
        return;
      }
      encoding.end();

      runMS2Batch_(session, batch, group, out);
    });

    // Sorted so the report does not depend on which session hit which chunk
    // first. Callers print these, and an order that changes run to run makes
    // two identical runs look different.
    if (failures != nullptr)
    {
      std::sort(failures->begin(), failures->end(),
                [](const Failure& a, const Failure& b) { return a.index < b.index; });
    }
    return out;
  }

  void PeptDeepPredictor::runMS2Batch_(std::size_t session,
                                       const PeptDeepEncoder::Batch& batch,
                                       const std::vector<std::size_t>& group,
                                       std::vector<Spectrum>& out)
  {
    {
      {
      const std::int64_t rows = static_cast<std::int64_t>(batch.rows);
      const std::int64_t length = static_cast<std::int64_t>(batch.sequence_length);
      const std::int64_t width = static_cast<std::int64_t>(PEPTDEEP_MOD_ELEMENTS.size());

      std::array<std::int64_t, 2> aa_shape{rows, length};
      std::array<std::int64_t, 3> mod_shape{rows, length, width};
      std::array<std::int64_t, 2> meta_shape{rows, 1};
      // Rank 1, unlike charges and nce. Feeding [rows, 1] is rejected outright.
      std::array<std::int64_t, 1> instrument_shape{rows};

      std::vector<Ort::Value> inputs;
      inputs.push_back(Ort::Value::CreateTensor<std::int64_t>(
        impl_->memory, const_cast<std::int64_t*>(batch.aa_indices.data()),
        batch.aa_indices.size(), aa_shape.data(), aa_shape.size()));
      inputs.push_back(Ort::Value::CreateTensor<float>(
        impl_->memory, const_cast<float*>(batch.mod_x.data()),
        batch.mod_x.size(), mod_shape.data(), mod_shape.size()));
      inputs.push_back(Ort::Value::CreateTensor<float>(
        impl_->memory, const_cast<float*>(batch.charges.data()),
        batch.charges.size(), meta_shape.data(), meta_shape.size()));
      inputs.push_back(Ort::Value::CreateTensor<float>(
        impl_->memory, const_cast<float*>(batch.nces.data()),
        batch.nces.size(), meta_shape.data(), meta_shape.size()));
      inputs.push_back(Ort::Value::CreateTensor<std::int64_t>(
        impl_->memory, const_cast<std::int64_t*>(batch.instrument_indices.data()),
        batch.instrument_indices.size(), instrument_shape.data(),
        instrument_shape.size()));

      trace::Scope running("ort_run", 4);
      tag(running, impl_->kind, session, group.size());
      auto results = impl_->at(session).Run(Ort::RunOptions{nullptr}, impl_->input_ptrs.data(),
                                         inputs.data(), inputs.size(),
                                         &impl_->output_ptr, 1);
      running.end();
      trace::Scope unpacking("unpack", 4);
      tag(unpacking, impl_->kind, session, group.size());

      const auto info = results.front().GetTensorTypeAndShapeInfo();
      const auto shape = info.GetShape();
      // shape[1] is checked as well as shape[0] and shape[2]. The comparison
      // against the reference cannot detect a wrong positions count -- both
      // sides take it from the same model -- and upstream OpenMS does not
      // trust it either, slicing to size()-1 rather than using what came back.
      if (info.GetElementType() != ONNX_TENSOR_ELEMENT_DATA_TYPE_FLOAT ||
          shape.size() != 3 || shape[0] != rows ||
          shape[1] != length - 3 ||
          shape[2] != static_cast<std::int64_t>(Spectrum::CHANNELS))
      {
        throw std::runtime_error(
          "MS2 output is not [batch, nAA - 1, 8]: the model returned a tensor of "
          "rank " + std::to_string(shape.size()) + " for " + std::to_string(rows) +
          " peptides of " + std::to_string(length - 2) + " residues");
      }

      const auto positions = static_cast<std::size_t>(shape[1]);
      const float* values = results.front().GetTensorData<float>();
      for (std::size_t i = 0; i < group.size(); ++i)
      {
        auto& spectrum = out[group[i]];
        spectrum.positions = positions;
        // Indexed by group[i], not i: the peptides were reordered into length
        // groups, and writing them back in batch order would silently permute
        // whole spectra between peptides of the same length.
        const auto* first = values + i * positions * Spectrum::CHANNELS;
        const auto* last = values + (i + 1) * positions * Spectrum::CHANNELS;
        if (!std::all_of(first, last, [](float v) { return std::isfinite(v); }))
        { throw std::runtime_error("MS2 model returned a non-finite spectrum for peptide index " + std::to_string(group[i])); }
        spectrum.intensities.assign(first, last);
      }
     }
    }
  }

  std::vector<float>
  PeptDeepPredictor::predictRT(const std::vector<OpenMS::AASequence>& peptides,
                               std::vector<Failure>* failures)
  {
    // The constructor accepts both model kinds, so the per-call guard has to be
    // here: feeding an MS2 model two of its five inputs fails inside Run with a
    // red "Missing Input" instead of a sentence a user can act on.
    if (impl_->input_names.size() != 2)
    {
      throw std::runtime_error("predictRT was given a model with " +
                               std::to_string(impl_->input_names.size()) +
                               " inputs; the retention-time model takes 2");
    }

    // NaN, not 0.0: an unfilled entry must be recognisable as "no prediction"
    // rather than pass for a real retention time.
    std::vector<float> out(peptides.size(), std::numeric_limits<float>::quiet_NaN());

    std::mutex failure_mutex;

    // Submitted in capped chunks; see MAX_BATCH_ROWS.
    forEachChunk_(peptides, [&](std::size_t session, const std::vector<std::size_t>& group)
    {
      trace::Scope encoding("encode", 4);
      tag(encoding, impl_->kind, session, group.size());
      PeptDeepEncoder::Batch batch;
      try
      {
        batch = PeptDeepEncoder::encode(peptides, group);
      }
      catch (const std::exception&)
      {
        encoding.end();
        // One unencodable peptide must not discard the whole run -- at proteome
        // scale that is minutes of inference lost to a single mass-shift
        // modification -- nor its whole chunk, which would take 2047 good
        // peptides with it. Retry the chunk one peptide at a time so only the
        // offending ones are recorded.
        if (failures == nullptr) { throw; }
        for (std::size_t i = 0; i < group.size(); ++i)
        {
          try
          {
            const auto single = PeptDeepEncoder::encode(peptides, std::vector<std::size_t>{group[i]});
            runBatch_(session, single, {group[i]}, out);
          }
          catch (const std::exception& inner)
          {
            const std::lock_guard<std::mutex> lock(failure_mutex);
            failures->push_back({group[i], inner.what()});
          }
        }
        return;
      }
      encoding.end();

      runBatch_(session, batch, group, out);
    });

    if (failures != nullptr)
    {
      std::sort(failures->begin(), failures->end(),
                [](const Failure& a, const Failure& b) { return a.index < b.index; });
    }
    return out;
  }

  std::vector<float>
  PeptDeepPredictor::predictCCS(const std::vector<OpenMS::AASequence>& peptides,
                                const std::vector<int>& charges,
                                std::vector<Failure>* failures)
  {
    if (charges.size() != peptides.size())
    {
      throw std::invalid_argument("one charge per peptide is required");
    }
    if (impl_->input_names.size() != 3)
    {
      throw std::runtime_error("predictCCS was given a model with " +
                               std::to_string(impl_->input_names.size()) +
                               " inputs; the CCS model takes 3");
    }

    // NaN, as in predictRT: an unfilled entry must be recognisable as "no
    // prediction" rather than pass for a collision cross-section of zero.
    std::vector<float> out(peptides.size(), std::numeric_limits<float>::quiet_NaN());

    std::mutex failure_mutex;

    forEachChunk_(peptides, [&](std::size_t session, const std::vector<std::size_t>& group)
    {
      trace::Scope encoding("encode", 4);
      tag(encoding, impl_->kind, session, group.size());
      PeptDeepEncoder::Batch batch;
      try
      {
        // The CCS model takes no NCE and no instrument, but the encoder's meta
        // overload fills both; they are simply not bound below. Reusing it
        // keeps the charge validation and scaling in one place.
        batch = PeptDeepEncoder::encode(peptides, group, charges, 30.0f, "QE");
      }
      catch (const std::exception&)
      {
        encoding.end();
        // One unencodable peptide, or one impossible charge, must not discard
        // its chunk. As in predictRT and predictMS2.
        if (failures == nullptr) { throw; }
        for (std::size_t i = 0; i < group.size(); ++i)
        {
          try
          {
            const auto single = PeptDeepEncoder::encode(
              peptides, std::vector<std::size_t>{group[i]}, charges, 30.0f, "QE");
            runCCSBatch_(session, single, {group[i]}, out);
          }
          catch (const std::exception& inner)
          {
            const std::lock_guard<std::mutex> lock(failure_mutex);
            failures->push_back({group[i], inner.what()});
          }
        }
        return;
      }
      encoding.end();

      runCCSBatch_(session, batch, group, out);
    });

    if (failures != nullptr)
    {
      std::sort(failures->begin(), failures->end(),
                [](const Failure& a, const Failure& b) { return a.index < b.index; });
    }
    return out;
  }

  void PeptDeepPredictor::runCCSBatch_(std::size_t session,
                                       const PeptDeepEncoder::Batch& batch,
                                       const std::vector<std::size_t>& group,
                                       std::vector<float>& out)
  {
    const std::int64_t rows = static_cast<std::int64_t>(batch.rows);
    const std::int64_t length = static_cast<std::int64_t>(batch.sequence_length);
    const std::int64_t width = static_cast<std::int64_t>(PEPTDEEP_MOD_ELEMENTS.size());

    std::array<std::int64_t, 2> aa_shape{rows, length};
    std::array<std::int64_t, 3> mod_shape{rows, length, width};
    std::array<std::int64_t, 2> charge_shape{rows, 1};

    std::vector<Ort::Value> inputs;
    inputs.push_back(Ort::Value::CreateTensor<std::int64_t>(
      impl_->memory, const_cast<std::int64_t*>(batch.aa_indices.data()),
      batch.aa_indices.size(), aa_shape.data(), aa_shape.size()));
    inputs.push_back(Ort::Value::CreateTensor<float>(
      impl_->memory, const_cast<float*>(batch.mod_x.data()),
      batch.mod_x.size(), mod_shape.data(), mod_shape.size()));
    inputs.push_back(Ort::Value::CreateTensor<float>(
      impl_->memory, const_cast<float*>(batch.charges.data()),
      batch.charges.size(), charge_shape.data(), charge_shape.size()));

    trace::Scope running("ort_run", 4);
    tag(running, impl_->kind, session, group.size());
    auto results = impl_->at(session).Run(Ort::RunOptions{nullptr}, impl_->input_ptrs.data(),
                                       inputs.data(), inputs.size(),
                                       &impl_->output_ptr, 1);
    running.end();
    trace::Scope unpacking("unpack", 4);
    tag(unpacking, impl_->kind, session, group.size());

    const auto info = results.front().GetTensorTypeAndShapeInfo();
    const auto shape = info.GetShape();
    // Rank 1, like the RT model and unlike what the format notes claim. A
    // [batch, 1] output would still have the right element count, so checking
    // the count alone would not notice.
    if (info.GetElementType() != ONNX_TENSOR_ELEMENT_DATA_TYPE_FLOAT ||
        shape.size() != 1 || shape[0] != rows)
    {
      throw std::runtime_error("CCS output is not [batch]");
    }

    const float* values = results.front().GetTensorData<float>();
    for (std::size_t i = 0; i < group.size(); ++i)
    {
      // Indexed by group[i]: the peptides were reordered into length groups.
      if (!std::isfinite(values[i]))
      { throw std::runtime_error("model returned a non-finite prediction for peptide index " + std::to_string(group[i])); }
      out[group[i]] = values[i];
    }
  }

} // namespace ODIA
