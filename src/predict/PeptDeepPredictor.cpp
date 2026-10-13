// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <odia/PeptDeepPredictor.h>

#include <odia/Hardware.h>
#include <odia/PeptDeepEncoder.h>
#include <odia/Trace.h>

#include "../library/Parallel.h"
#include "../library/ThreadPool.h"
#include "CoreMLSupport.h"
#include "OnnxWire.h"
#include "OrtEnv.h"

#include <OpenMS/CHEMISTRY/AASequence.h>

#include <onnxruntime_cxx_api.h>

#include <algorithm>
#include <atomic>
#include <cctype>
#include <chrono>
#include <iomanip>
#include <locale>
#include <sstream>
#include <condition_variable>
#include <map>
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
#include <unordered_map>
#include <limits>
#include <stdexcept>

#if defined(_WIN32)
#include <process.h>
#else
#include <unistd.h>
#endif
#if defined(__APPLE__)
#include <fcntl.h>
#include <mach-o/dyld.h>
#include <signal.h>
#include <spawn.h>
#include <sys/wait.h>
extern char** environ;
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

  /// Registers the CoreML provider: an MLProgram (the current Core ML model
  /// format; the legacy NeuralNetwork format has fewer operators), compiled
  /// for static shapes only -- with dynamic ones Core ML's compiler rejects
  /// the models' reshapes ("unbounded dimension") and the provider runs them
  /// slower than one CPU thread -- and cached on disk, since compiling one
  /// shape takes 3-4 s and loading it back 0.8 s.
  void appendCoreML(Ort::SessionOptions& options, const std::string& units,
                    const std::filesystem::path& cache, const char* log_id)
  {
    std::unordered_map<std::string, std::string> provider;
    provider["ModelFormat"] = "MLProgram";
    provider["MLComputeUnits"] = units;
    provider["RequireStaticInputShapes"] = "1";
    provider["EnableOnSubgraphs"] = "0";
    if (!cache.empty()) { provider["ModelCacheDirectory"] = cache.string(); }
    // ODIA_COREML_PROFILE=1: Core ML's compute plan (which unit runs each
    // operator) in the log, for diagnosis.
    if (const char* profile = std::getenv("ODIA_COREML_PROFILE"); profile != nullptr && *profile == '1')
    { provider["ProfileComputePlan"] = "1"; }
    options.AppendExecutionProvider("CoreML", provider);
    options.SetLogId(log_id);
  }

  /// The batch sizes CoreML sessions are compiled for, per residue count
  /// (@p tokens = residues + 2). A chunk is padded up to the smallest bucket
  /// that holds it; every bucket is another compiled shape per residue count
  /// (0.8-4 s each). Up to 21 residues a block has several full chunks per
  /// length and one partial one, 5-11 % padding: MAX_BATCH_ROWS alone. From
  /// 22 residues on the lengths are rarer and padded 15-38 %, so 512 and
  /// 1024 are added (whole human library, 200 k blocks: 11.5 % padded rows
  /// over the phase against 8.8 %; the saved rows are the slow ones, 2-3 k
  /// rows/s per session at 25-30 residues against 7-8 k at 7-8).
  /// ODIA_COREML_BUCKETS=512,2048 overrides it for every length and
  /// ODIA_COREML_BUCKETS=0 is MAX_BATCH_ROWS alone, for measurements.
  constexpr std::int64_t FINE_BUCKETS_FROM_RESIDUES = 22;
  std::vector<std::int64_t> coremlBucketsFor(std::int64_t tokens)
  {
    static const std::vector<std::int64_t> override = [] {
      std::vector<std::int64_t> buckets;
      if (const char* v = std::getenv("ODIA_COREML_BUCKETS"); v != nullptr && *v != '\0')
      {
        std::string_view rest(v);
        while (!rest.empty())
        {
          const auto comma = rest.find(',');
          const std::string token(rest.substr(0, comma));
          rest = comma == std::string_view::npos ? std::string_view() : rest.substr(comma + 1);
          char* end = nullptr;
          const long long x = std::strtoll(token.c_str(), &end, 10);
          if (end != token.c_str() && x > 0) { buckets.push_back(static_cast<std::int64_t>(x)); }
        }
        std::sort(buckets.begin(), buckets.end());
        if (buckets.empty() || buckets.back() < static_cast<std::int64_t>(ODIA::PeptDeepPredictor::MAX_BATCH_ROWS))
        { buckets.push_back(static_cast<std::int64_t>(ODIA::PeptDeepPredictor::MAX_BATCH_ROWS)); }
      }
      return buckets;
    }();
    if (!override.empty()) { return override; }
    constexpr auto max_rows = static_cast<std::int64_t>(ODIA::PeptDeepPredictor::MAX_BATCH_ROWS);
    if (tokens - 2 >= FINE_BUCKETS_FROM_RESIDUES) { return {512, 1024, max_rows}; }
    return {max_rows};
  }

  /// Where compiled CoreML models live between runs: ODIA_COREML_CACHE, or
  /// ~/Library/Caches/DIALibGen/coreml, else the temporary directory.
  /// ODIA_COREML_CACHE=0 keeps nothing (every session compiles its shape).
  std::filesystem::path coremlCacheRoot()
  {
    if (const char* v = std::getenv("ODIA_COREML_CACHE"); v != nullptr && *v != '\0')
    {
      if (std::string_view(v) == "0") { return {}; }
      return std::filesystem::path(v);
    }
    if (const char* home = std::getenv("HOME"); home != nullptr && *home != '\0')
    { return std::filesystem::path(home) / "Library" / "Caches" / "DIALibGen" / "coreml"; }
    std::error_code ec;
    const auto tmp = std::filesystem::temp_directory_path(ec);
    return ec ? std::filesystem::path() : tmp / "dialibgen-coreml";
  }

  /// The compiled-model cache's size cap: ODIA_COREML_CACHE_MAX in bytes,
  /// with an optional K, M or G suffix (4G by default; 0: no cap). Above it
  /// the least recently used shapes of other models are removed when a
  /// predictor goes (pruneCoremlCache).
  std::uintmax_t coremlCacheMax()
  {
    constexpr std::uintmax_t DEFAULT = std::uintmax_t{4} << 30;
    const char* v = std::getenv("ODIA_COREML_CACHE_MAX");
    if (v == nullptr || *v == '\0') { return DEFAULT; }
    char* end = nullptr;
    const auto x = std::strtoull(v, &end, 10);
    if (end == v) { return DEFAULT; }
    std::uintmax_t bytes = x;
    switch (std::tolower(static_cast<unsigned char>(*end)))
    {
      case 'k': bytes <<= 10; break;
      case 'm': bytes <<= 20; break;
      case 'g': bytes <<= 30; break;
      default: break;
    }
    return bytes;
  }

  std::uintmax_t directorySize(const std::filesystem::path& dir)
  {
    std::uintmax_t total = 0;
    std::error_code ec;
    for (std::filesystem::recursive_directory_iterator it(dir, ec), end; !ec && it != end; it.increment(ec))
    {
      if (it->is_regular_file(ec)) { total += it->file_size(ec); }
    }
    return total;
  }

  /// Keeps the cache under coremlCacheMax(): the shape directories
  /// (<root>/<model>-<hash>/<rows>x<tokens>) are removed least recently
  /// used first (their modification time, which buildCoreML refreshes on
  /// every use) until the total fits, except the shapes of @p keep, the
  /// model that is running: they are what the next run needs, whatever
  /// their size. Returns how many were removed.
  std::size_t pruneCoremlCache(const std::filesystem::path& root, const std::filesystem::path& keep)
  {
    const std::uintmax_t cap = coremlCacheMax();
    if (cap == 0 || root.empty()) { return 0; }
    struct Shape { std::filesystem::path dir; std::filesystem::file_time_type used; std::uintmax_t bytes; bool keep; };
    std::vector<Shape> shapes;
    std::uintmax_t total = 0;
    std::error_code ec;
    for (std::filesystem::directory_iterator model(root, ec), end; !ec && model != end; model.increment(ec))
    {
      if (!model->is_directory(ec)) { continue; }
      for (std::filesystem::directory_iterator shape(model->path(), ec), send; !ec && shape != send; shape.increment(ec))
      {
        if (!shape->is_directory(ec)) { continue; }
        Shape s{shape->path(), shape->last_write_time(ec), directorySize(shape->path()), model->path() == keep};
        total += s.bytes;
        shapes.push_back(std::move(s));
      }
    }
    if (total <= cap) { return 0; }
    std::sort(shapes.begin(), shapes.end(), [](const Shape& a, const Shape& b) { return a.used < b.used; });
    std::size_t removed = 0;
    for (const auto& s : shapes)
    {
      if (total <= cap) { break; }
      if (s.keep) { continue; }
      std::filesystem::remove_all(s.dir, ec);
      if (ec) { continue; }
      total -= s.bytes;
      ++removed;
      std::filesystem::remove(s.dir.parent_path(), ec);   // an emptied model directory; fails, harmlessly, when not empty
    }
    return removed;
  }

  /// Core ML's verdict on a model, taken in a child process (`DIALibGen
  /// --probe-coreml`): the compile of the probe shape and one Run happen
  /// there first, so that a Core ML that aborts the process (macOS virtual
  /// machines without a usable GPU, GitHub's macOS runners among them)
  /// costs the child, not the run. Once per process and model.
  struct ChildProbe
  {
    bool ran = false;      ///< a child could be started at all
    int exit_code = -1;    ///< 0 usable, 2 model stays on the CPU, else Core ML failed
    int signal = 0;        ///< the child died of this signal
    std::string reason;    ///< the child's one line
  };

  ChildProbe probeCoremlInChild(const std::string& model_path, const std::string& units)
  {
    ChildProbe probe;
#if defined(__APPLE__)
    // The child must not probe in turn: the environment says it is one.
    if (const char* in_child = std::getenv("ODIA_COREML_PROBE_CHILD"); in_child != nullptr && *in_child == '1') { return probe; }
    if (const char* off = std::getenv("ODIA_COREML_PROBE"); off != nullptr && std::string_view(off) == "0") { return probe; }
    static std::mutex mutex;
    static std::map<std::string, ChildProbe> done;
    const std::string key = model_path + "\n" + units;
    const std::lock_guard<std::mutex> lock(mutex);
    if (const auto it = done.find(key); it != done.end()) { return it->second; }
    std::vector<char> exe(4096);
    std::uint32_t size = static_cast<std::uint32_t>(exe.size());
    if (_NSGetExecutablePath(exe.data(), &size) != 0) { exe.resize(size); if (_NSGetExecutablePath(exe.data(), &size) != 0) { return probe; } }
    int pipefd[2];
    if (::pipe(pipefd) != 0) { return probe; }
    posix_spawn_file_actions_t actions;
    posix_spawn_file_actions_init(&actions);
    posix_spawn_file_actions_adddup2(&actions, pipefd[1], 1);
    posix_spawn_file_actions_addclose(&actions, pipefd[0]);
    posix_spawn_file_actions_addclose(&actions, pipefd[1]);
    std::vector<std::string> env_strings{"ODIA_COREML_PROBE_CHILD=1"};
    for (char** e = environ; e != nullptr && *e != nullptr; ++e)
    { if (std::string_view(*e).rfind("ODIA_COREML_PROBE_CHILD=", 0) != 0) { env_strings.emplace_back(*e); } }
    std::vector<char*> envp;
    for (auto& s : env_strings) { envp.push_back(s.data()); }
    envp.push_back(nullptr);
    std::string flag = "--probe-coreml", model = model_path, u = units;
    char* argv[] = {exe.data(), flag.data(), model.data(), u.data(), nullptr};
    pid_t pid = 0;
    const int rc = posix_spawn(&pid, exe.data(), &actions, nullptr, argv, envp.data());
    posix_spawn_file_actions_destroy(&actions);
    ::close(pipefd[1]);
    if (rc != 0) { ::close(pipefd[0]); return probe; }
    probe.ran = true;
    // The child's output and its exit, polled together; a child that does
    // not exit within the limit (ODIA_COREML_PROBE_SECONDS, 300) is killed.
    ::fcntl(pipefd[0], F_SETFL, ::fcntl(pipefd[0], F_GETFL) | O_NONBLOCK);
    std::string out;
    const auto drain = [&]() {
      char buffer[256];
      for (;;)
      {
        const auto n = ::read(pipefd[0], buffer, sizeof buffer);
        if (n <= 0) { break; }
        out.append(buffer, static_cast<std::size_t>(n));
      }
    };
    int limit = 300;
    if (const char* v = std::getenv("ODIA_COREML_PROBE_SECONDS"); v != nullptr && *v != '\0') { limit = std::max(1, std::atoi(v)); }
    int status = 0;
    const auto t0 = std::chrono::steady_clock::now();
    for (;;)
    {
      drain();
      const pid_t w = ::waitpid(pid, &status, WNOHANG);
      if (w == pid) { break; }
      if (w < 0) { ::close(pipefd[0]); probe.ran = false; return probe; }
      if (std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() > limit)
      {
        ::kill(pid, SIGKILL);
        ::waitpid(pid, &status, 0);
        ::close(pipefd[0]);
        probe.reason = "the probe process did not finish within " + std::to_string(limit) + " s";
        probe.exit_code = -1;
        done[key] = probe;
        return probe;
      }
      std::this_thread::sleep_for(std::chrono::milliseconds(20));
    }
    drain();
    ::close(pipefd[0]);
    if (const auto nl = out.find('\n'); nl != std::string::npos) { out.resize(nl); }
    probe.reason = out;
    if (WIFEXITED(status)) { probe.exit_code = WEXITSTATUS(status); }
    else if (WIFSIGNALED(status)) { probe.exit_code = -1; probe.signal = WTERMSIG(status); probe.reason = "the probe process died of signal " + std::to_string(probe.signal) + (out.empty() ? "" : " (" + out + ")"); }
    done[key] = probe;
#else
    (void)model_path; (void)units;
#endif
    return probe;
  }

  std::string fallbackFix(const std::string& reason)
  {
    std::string lower = reason;
    std::transform(lower.begin(), lower.end(), lower.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    const char* docs = "https://github.com/okohlbacher/DIALibGen/blob/main/docs/usage.md#gpu";
    if (lower.rfind("coreml", 0) == 0)
    {
      return std::string("this model runs on the CPU path (no slower than without -machine:device coreml); "
                         "see https://github.com/okohlbacher/DIALibGen/blob/main/docs/usage.md#apple-gpu");
    }
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
  /// Names the CoreML probe sessions' log ids apart within one process.
  std::atomic<std::int64_t> g_probes{0};

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

  using WireReader = ODIA::onnxwire::Reader;

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
    /// it, so every predictor's Env is the same object underneath (OrtEnv.h).
    Ort::Env& env = ort::env();
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
    std::string note;       ///< see providerNote()
    std::vector<int> session_devices;   ///< CUDA device per session; empty on CPU

    // CoreML: the model is compiled per batch shape (rows x residues + 2),
    // so the sessions are a map built as shapes turn up, from the model
    // bytes with PRelu rewritten (CoreMLSupport.h). All workers share them:
    // Run is thread-safe, and the GPU is saturated by one or two streams.
    std::string model_bytes;
    std::string dim_rows, dim_length;     ///< the free dimensions' names
    std::string coreml_units;
    std::filesystem::path coreml_cache;   ///< per model content; empty = no cache
    std::map<std::pair<std::int64_t, std::int64_t>, std::unique_ptr<Ort::Session>> coreml_sessions;
    std::mutex coreml_mutex;
    std::condition_variable coreml_built;
    std::size_t coreml_workers = 1;
    double coreml_build_seconds = 0;      ///< for the trace: compile and load time
    std::size_t coreml_builds = 0;
    /// The shapes GpuOptions::prebuild_* asked for, built in the background
    /// by these threads (joined by waitForShapes and the destructor).
    std::vector<std::pair<std::int64_t, std::int64_t>> prebuild;
    bool prebuild_discard = false;
    /// Shapes Core ML refused in-process: their batches run on `session`,
    /// a CPU session built on first need (sessionFor_).
    std::map<std::pair<std::int64_t, std::int64_t>, bool> coreml_cpu_fallback;
    std::vector<std::thread> builders;
    /// Rows run and rows padded (the difference is wasted GPU time), for
    /// the summary line.
    std::atomic<std::uint64_t> coreml_rows{0}, coreml_padded_rows{0};

    /// A CoreML session for one shape. Builds take 0.8 s (cached) to 4 s
    /// (compiled); a worker that needs a shape another worker is building
    /// waits for it rather than building a second copy.
    std::unique_ptr<Ort::Session> buildCoreML(std::int64_t rows, std::int64_t length, const char* log_id,
                                              bool verbose)
    {
      Ort::SessionOptions options;
      configureSession(options, 1);
      // Through the C API: the C++ wrapper of conda-forge's ONNX Runtime
      // 1.22 for macOS x64 has no AddFreeDimensionOverrideByName.
      Ort::ThrowOnError(Ort::GetApi().AddFreeDimensionOverrideByName(options, dim_rows.c_str(), rows));
      Ort::ThrowOnError(Ort::GetApi().AddFreeDimensionOverrideByName(options, dim_length.c_str(), length));
      std::filesystem::path cache;
      if (!coreml_cache.empty())
      {
        cache = coreml_cache / (std::to_string(rows) + "x" + std::to_string(length));
        std::error_code ec;
        std::filesystem::create_directories(cache, ec);
        if (ec) { cache.clear(); }
      }
      appendCoreML(options, coreml_units, cache, log_id);
      options.SetLogSeverityLevel(verbose ? ORT_LOGGING_LEVEL_VERBOSE : ORT_LOGGING_LEVEL_ERROR);
      const auto t0 = std::chrono::steady_clock::now();
      auto session = std::make_unique<Ort::Session>(env, model_bytes.data(), model_bytes.size(), options);
      if (!cache.empty())
      {
        // Used now: what the cache's least-recently-used pruning goes by.
        std::error_code ec;
        std::filesystem::last_write_time(cache, std::filesystem::file_time_type::clock::now(), ec);
      }
      // Workers build concurrently (sessionFor_): the counters are shared.
      const std::lock_guard<std::mutex> lock(coreml_mutex);
      coreml_build_seconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
      ++coreml_builds;
      return session;
    }

    /// Builds the claimed shapes of `prebuild` in order on this thread,
    /// shared with the other builders through the index. A shape that fails
    /// is released (sessionFor_ then builds it itself and reports the error)
    /// and the builder goes on with the next.
    void buildAhead(std::atomic<std::size_t>& next)
    {
      for (std::size_t i = next++; i < prebuild.size(); i = next++)
      {
        const auto key = prebuild[i];
        std::unique_ptr<Ort::Session> built;
        try
        {
          trace::Scope building("session_create", 3);
          building.fields().rows = key.first;
          building.fields().len = key.second;
          built = buildCoreML(key.first, key.second, "odia-coreml", false);
        }
        catch (...) {}
        const std::lock_guard<std::mutex> lock(coreml_mutex);
        if (built && !prebuild_discard) { coreml_sessions[key] = std::move(built); }
        else { coreml_sessions.erase(key); }   // failed, or compiled for the cache only
        coreml_built.notify_all();
      }
    }
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
    if (prefer_gpu && gpu.coreml)
    {
      // Apple. The provider is asked for explicitly (-machine:device coreml
      // or mps), so every way of not getting it is worth a WARNING.
      gpu_present = true;
      bool available = false;
      for (const auto& name : Ort::GetAvailableProviders())
      {
        if (name == "CoreMLExecutionProvider") { available = true; break; }
      }
      if (!available) { reason = "CoreML: this ONNX Runtime build has no CoreML provider"; }
      else
      {
        trace::Scope probe_step("session_probe", 3);
        // Core ML first in a process of its own: a compile or Run that
        // takes the process down (a macOS virtual machine without a usable
        // GPU) must not take the run with it.
        // Exit 0: it runs; 2: the child's placement check kept the model on
        // the CPU, and the one below says why (the compile itself survived);
        // anything else: Core ML failed, and the run stays off it.
        const ChildProbe child = probeCoremlInChild(model_path, gpu.coreml_units.empty() ? "CPUAndGPU" : gpu.coreml_units);
        if (child.ran && child.exit_code != 0 && child.exit_code != 2)
        {
          reason = "CoreML: Core ML failed in the probe process (" +
                   (child.reason.empty() ? "exit code " + std::to_string(child.exit_code) : child.reason) + ")";
        }
        else
        try
        {
          {
            std::ifstream in(model_path, std::ios::binary);
            if (!in) { throw std::runtime_error("cannot read " + model_path); }
            impl_->model_bytes.assign(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
          }
          if (!coreml::freeDimensionNames(impl_->model_bytes, impl_->dim_rows, impl_->dim_length))
          { throw std::runtime_error("the model's first input is not [batch, length] with named dimensions"); }
          const std::size_t rewritten = coreml::rewritePRelu(impl_->model_bytes);
          impl_->coreml_units = gpu.coreml_units.empty() ? "CPUAndGPU" : gpu.coreml_units;
          if (const auto root = coremlCacheRoot(); !root.empty())
          {
            impl_->coreml_cache = root / (std::filesystem::path(model_path).stem().string() + "-" +
                                          std::to_string(coreml::fnv64(impl_->model_bytes)));
          }
          // The probe is a real session for a real shape -- the smallest
          // bucket, 16 tokens (14 residues, the commonest tryptic length) --
          // whose verbose log says where ONNX Runtime put the nodes. It is
          // kept: it is the session that shape will run on.
          const std::int64_t probe_length = 16, probe_rows = coremlBucketsFor(probe_length).back();
          const std::string log_id = "odia-probe-" + std::to_string(g_probes++);
          ort::LogCapture capture(log_id);
          auto probe = impl_->buildCoreML(probe_rows, probe_length, log_id.c_str(), true);
          const auto placement = coreml::parsePlacement(capture.text());
          if (!placement.found) { throw std::runtime_error("ONNX Runtime did not report the node placement"); }
          std::map<std::string, std::size_t> on_cpu;
          for (const auto& op : placement.cpu_ops) { ++on_cpu[op]; }
          std::string heavy, light;
          for (const auto& [op, n] : on_cpu)
          {
            auto& where = coreml::heavyOp(op) ? heavy : light;
            where += (where.empty() ? "" : ", ") + op + (n > 1 ? " x" + std::to_string(n) : "");
          }
          if (!heavy.empty())
          {
            throw std::runtime_error("CoreML has no operator for " + heavy + ": ONNX Runtime keeps " +
                                     std::to_string(placement.cpu_ops.size()) + " of the model's nodes on the CPU (" +
                                     std::to_string(placement.coreml_nodes) + " CoreML partition" +
                                     (placement.coreml_nodes == 1 ? "" : "s") + ")");
          }
          if (!light.empty()) { impl_->note = "on the CPU: " + light; }
          impl_->provider = Provider::CoreML;
          impl_->coreml_sessions[{probe_rows, probe_length}] = std::move(probe);
          if (trace::enabled())
          {
            trace::note(std::string("predictor.") + std::filesystem::path(model_path).stem().string() + ".coreml",
                        std::to_string(placement.coreml_nodes) + " partitions, cpu ops: " + (light.empty() ? "none" : light) +
                        ", prelu rewritten: " + std::to_string(rewritten) + ", units " + impl_->coreml_units);
          }
        }
        catch (const Ort::Exception& e) { reason = "CoreML: " + conciseOrtError(e.what()); }
        catch (const std::exception& e) { reason = std::string(e.what()).rfind("CoreML", 0) == 0 ? e.what() : "CoreML: " + std::string(e.what()); }
        if (impl_->provider != Provider::CoreML)
        {
          impl_->model_bytes.clear();
          impl_->model_bytes.shrink_to_fit();
          impl_->coreml_sessions.clear();
        }
      }
    }
    else if (prefer_gpu)
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
    // covers both. On CoreML the probe is the first shape's session.
    if (impl_->provider == Provider::CPU)
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
      const char* wanted = gpu.coreml ? "CoreML" : "CUDA";
      if (gpu.require)
      {
        throw std::runtime_error("the GPU is required (-machine:require_gpu) but " + model +
                                 " cannot run on " + wanted + ": " + reason + ". Fix: " + fallbackFix(reason));
      }
      // A silent fall back to CPU turns a 2.5-minute library into a
      // 6.5-minute one (H100, whole proteome) with nothing in the log to
      // say why: the caller prints and records this.
      if (gpu_present)
      { impl_->fallback = model + " runs on the CPU although " + wanted + " was requested: " + reason + ". Fix: " + fallbackFix(reason); }
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
    Ort::Session& first = impl_->provider == Provider::CoreML ? *impl_->coreml_sessions.begin()->second : *impl_->session;
    for (std::size_t i = 0; i < first.GetInputCount(); ++i)
    {
      impl_->input_names.emplace_back(first.GetInputNameAllocated(i, allocator).get());
    }
    for (std::size_t i = 0; i < first.GetOutputCount(); ++i)
    {
      impl_->output_names.emplace_back(first.GetOutputNameAllocated(i, allocator).get());
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
    if (impl_->provider == Provider::CoreML)
    {
      // Workers, not sessions: the shape sessions are shared. Two concurrent
      // batches gain about 10 % over one on an M4 Max (the host-side
      // encoding of one overlaps the GPU work of the other); four gain
      // nothing more, and cost nothing but the option's default.
      impl_->coreml_workers = static_cast<std::size_t>(std::max(1, gpu.sessions_per_device));
      // The shapes of the peptide-length window, built in the background
      // from now on, several at a time (ODIA_COREML_BUILDERS, 8 by default),
      // in the order the chunks will ask for them: by length, the full
      // chunks' shape before the remainders'. Core ML compiles one model at
      // a time whatever the builder count (42 shapes cold: 105 s with 4, 8
      // or 16 builders, about 2.5 s each); cached loads gain a little (23 s
      // with 8 builders against 34 s one by one). The gain is the start:
      // the builds overlap the RT and CCS phases. Without this every shape
      // was built by the first worker to need it while the others waited,
      // one after the other, because the chunks of a block come in length
      // order: 66 s of a cold-cache whole-proteome run, 16 s warm.
      // The lengths to build: the library's, with the buckets a chunk of
      // that many rows can need (a chunk has at most MAX_BATCH_ROWS rows,
      // and no more than the length's precursors), or every bucket of the
      // window (prepare).
      std::vector<std::pair<int, std::uint64_t>> lengths = gpu.prebuild_lengths;
      if (lengths.empty() && gpu.prebuild_min_residues > 0 && gpu.prebuild_max_residues >= gpu.prebuild_min_residues)
      {
        for (int residues = gpu.prebuild_min_residues; residues <= gpu.prebuild_max_residues; ++residues)
        { lengths.emplace_back(residues, std::numeric_limits<std::uint64_t>::max()); }
      }
      std::sort(lengths.begin(), lengths.end());
      impl_->prebuild_discard = gpu.prebuild_discard;
      if (!lengths.empty())
      {
        for (const auto& [residues, count] : lengths)
        {
          if (residues <= 0 || count == 0) { continue; }
          const auto tokens = static_cast<std::int64_t>(residues) + 2;
          const auto buckets = coremlBucketsFor(tokens);
          const auto largest_chunk = static_cast<std::int64_t>(std::min<std::uint64_t>(count, MAX_BATCH_ROWS));
          std::int64_t top = buckets.back();
          for (const auto b : buckets) { if (b >= largest_chunk) { top = b; break; } }
          for (auto it = buckets.rbegin(); it != buckets.rend(); ++it)
          {
            if (*it > top) { continue; }
            const std::pair<std::int64_t, std::int64_t> key{*it, tokens};
            if (impl_->coreml_sessions.count(key) != 0) { continue; }
            impl_->coreml_sessions[key] = nullptr;   // claimed: being built
            impl_->prebuild.push_back(key);
          }
        }
        std::size_t builders = 8;
        if (const char* v = std::getenv("ODIA_COREML_BUILDERS"); v != nullptr && *v != '\0')
        {
          char* end = nullptr;
          const long long x = std::strtoll(v, &end, 10);
          if (end != v && x > 0) { builders = static_cast<std::size_t>(x); }
        }
        builders = std::min(builders, impl_->prebuild.size());
        // Shared by the builders; outlives them (joined in waitForShapes).
        impl_->builders.reserve(builders);
        auto next = std::make_shared<std::atomic<std::size_t>>(0);
        for (std::size_t b = 0; b < builders; ++b)
        {
          try { impl_->builders.emplace_back([impl = impl_.get(), next] { impl->buildAhead(*next); }); }
          catch (const std::system_error&) { break; }   // no thread: the workers build what is left
        }
        if (impl_->builders.empty())
        {
          const std::lock_guard<std::mutex> lock(impl_->coreml_mutex);
          for (const auto& key : impl_->prebuild) { impl_->coreml_sessions.erase(key); }
          impl_->prebuild.clear();
        }
      }
    }
    else if (impl_->provider == Provider::CUDA)
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

  const std::string& PeptDeepPredictor::providerNote() const { return impl_->note; }

  const std::vector<int>& PeptDeepPredictor::sessionDevices() const { return impl_->session_devices; }

  std::size_t PeptDeepPredictor::sessionCount() const
  {
    if (impl_->provider == Provider::CoreML) { return impl_->coreml_workers; }
    return impl_->replicas.size() + 1;
  }

  std::pair<std::size_t, double> PeptDeepPredictor::waitForShapes()
  {
    if (impl_->provider != Provider::CoreML) { return {0, 0.0}; }
    for (auto& t : impl_->builders) { if (t.joinable()) { t.join(); } }
    impl_->builders.clear();
    const std::lock_guard<std::mutex> lock(impl_->coreml_mutex);
    return {impl_->coreml_builds, impl_->coreml_build_seconds};
  }

  int PeptDeepPredictor::probeCoremlMain(const std::string& model_path, const std::string& units)
  {
    try
    {
      GpuOptions gpu;
      gpu.coreml = true;
      gpu.coreml_units = units;
      PeptDeepPredictor predictor(model_path, true, 1, 1, -1, gpu);
      if (predictor.provider() != Provider::CoreML)
      {
        std::cout << (predictor.fallbackReason().empty() ? std::string("CoreML: the model stays on the CPU") : predictor.fallbackReason()) << std::endl;
        return 2;
      }
      // One real batch through the probe shape: the Run is where a broken
      // Core ML runtime fails, not the compile.
      const std::vector<OpenMS::AASequence> peptides{OpenMS::AASequence::fromString("ELVISLIVESKING")};   // 14 residues: the probe shape
      const std::string kind = predictor.impl_->kind;
      if (kind == "MS2") { predictor.predictMS2(peptides, {2}); }
      else if (kind == "CCS") { predictor.predictCCS(peptides, {2}); }
      else { predictor.predictRT(peptides); }
      std::cout << "CoreML runs " << std::filesystem::path(model_path).filename().string() << std::endl;
      return 0;
    }
    catch (const std::exception& e)
    {
      std::cout << "CoreML: " << e.what() << std::endl;
      return 3;
    }
  }

  bool PeptDeepPredictor::coremlProviderAvailable()
  {
    for (const auto& name : Ort::GetAvailableProviders()) { if (name == "CoreMLExecutionProvider") { return true; } }
    return false;
  }

  std::filesystem::path PeptDeepPredictor::coremlCacheDir() const
  {
    return impl_->provider == Provider::CoreML ? impl_->coreml_cache : std::filesystem::path();
  }

  Ort::Session& PeptDeepPredictor::sessionFor_(std::size_t worker, std::int64_t rows, std::int64_t length)
  {
    if (impl_->provider != Provider::CoreML) { return impl_->at(worker); }
    const std::pair<std::int64_t, std::int64_t> key{rows, length};
    std::unique_lock<std::mutex> lock(impl_->coreml_mutex);
    for (;;)
    {
      if (impl_->coreml_cpu_fallback.count(key) != 0) { return *impl_->session; }
      const auto it = impl_->coreml_sessions.find(key);
      if (it == impl_->coreml_sessions.end())
      {
        impl_->coreml_sessions[key] = nullptr;   // claimed: being built
        break;
      }
      if (it->second) { return *it->second; }
      impl_->coreml_built.wait(lock);
    }
    lock.unlock();
    std::unique_ptr<Ort::Session> built;
    try
    {
      trace::Scope building("session_create", 3);
      building.fields().rows = rows;
      building.fields().len = length;
      built = impl_->buildCoreML(rows, length, "odia-coreml", false);
    }
    catch (const std::exception& e)
    {
      // A shape Core ML will not compile runs on the CPU, this batch and
      // every later one of that shape, and says so once: the padded rows
      // cost CPU time but the result is the CPU's.
      lock.lock();
      impl_->coreml_sessions.erase(key);
      impl_->coreml_built.notify_all();
      if (!impl_->session)
      {
        Ort::SessionOptions cpu;
        configureSession(cpu, 1);
        impl_->session = std::make_unique<Ort::Session>(impl_->env, impl_->model_bytes.data(), impl_->model_bytes.size(), cpu);
      }
      std::cerr << "WARNING: " << impl_->kind << ": CoreML could not build the " << rows << "x" << length
                << " shape (" << conciseOrtError(e.what()) << "); batches of that shape run on the CPU\n";
      impl_->coreml_cpu_fallback[key] = true;
      return *impl_->session;
    }
    lock.lock();
    Ort::Session& session = *(impl_->coreml_sessions[key] = std::move(built));
    impl_->coreml_built.notify_all();
    return session;
  }

  const PeptDeepEncoder::Batch& PeptDeepPredictor::padded_(const PeptDeepEncoder::Batch& batch,
                                                           PeptDeepEncoder::Batch& storage) const
  {
    if (impl_->provider != Provider::CoreML) { return batch; }
    const auto buckets = coremlBucketsFor(static_cast<std::int64_t>(batch.sequence_length));
    std::int64_t bucket = buckets.back();
    for (const auto b : buckets) { if (b >= static_cast<std::int64_t>(batch.rows)) { bucket = b; break; } }
    impl_->coreml_rows += batch.rows;
    impl_->coreml_padded_rows += static_cast<std::uint64_t>(bucket) - batch.rows;
    if (bucket == static_cast<std::int64_t>(batch.rows)) { return batch; }
    // Padded with copies of the last row: valid input, whose predictions are
    // dropped. The real rows are not affected by it (the models apply no
    // cross-row operation), and the CoreML library is compared with the
    // CPU's by the tests on that basis.
    const auto rows = static_cast<std::size_t>(bucket);
    storage = batch;
    storage.rows = rows;
    const auto grow = [&](auto& values, std::size_t per_row) {
      if (values.empty()) { return; }
      const auto last = (batch.rows - 1) * per_row;
      values.resize(rows * per_row);
      for (std::size_t r = batch.rows; r < rows; ++r)
      { std::copy_n(values.begin() + static_cast<std::ptrdiff_t>(last), per_row, values.begin() + static_cast<std::ptrdiff_t>(r * per_row)); }
    };
    grow(storage.aa_indices, batch.sequence_length);
    grow(storage.mod_x, batch.sequence_length * PEPTDEEP_MOD_ELEMENTS.size());
    grow(storage.charges, 1);
    grow(storage.nces, 1);
    grow(storage.instrument_indices, 1);
    return storage;
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
    if (impl_ && impl_->provider == Provider::CoreML)
    {
      // The builders refer to impl_: joined before anything of it goes.
      waitForShapes();
      const std::uint64_t rows = impl_->coreml_rows.load(), padded = impl_->coreml_padded_rows.load();
      std::ostringstream summary;
      summary.imbue(std::locale::classic());
      summary << std::fixed << std::setprecision(1) << impl_->coreml_builds << " shapes built in "
              << impl_->coreml_build_seconds << " s (summed), " << rows << " rows run, " << padded << " padding rows ("
              << (rows + padded > 0 ? 100.0 * static_cast<double>(padded) / static_cast<double>(rows + padded) : 0.0) << " %)";
      if (rows > 0) { std::cerr << "[" << impl_->kind << "] CoreML: " << summary.str() << "\n"; }
      if (trace::enabled()) { trace::note(std::string("predictor.") + impl_->kind + ".coreml_sessions", summary.str()); }
      // The sessions and the builders are done with the cache: the least
      // recently used shapes of other models go when it is over the cap.
      if (!impl_->coreml_cache.empty())
      {
        const std::size_t removed = pruneCoremlCache(impl_->coreml_cache.parent_path(), impl_->coreml_cache);
        if (removed > 0)
        { std::cerr << "[" << impl_->kind << "] CoreML: removed " << removed << " least recently used compiled shape(s) of other models from "
                    << impl_->coreml_cache.parent_path().string() << " (cap ODIA_COREML_CACHE_MAX)\n"; }
      }
    }
    if (!impl_ || !impl_->ort_profiling) { return; }
    // End explicitly rather than at session destruction, so the file names
    // can be reported.
    Ort::AllocatorWithDefaultOptions allocator;
    for (std::size_t i = 0; i < sessionCount() && impl_->provider != Provider::CoreML; ++i)
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
                                    const PeptDeepEncoder::Batch& unpadded,
                                    const std::vector<std::size_t>& group,
                                    std::vector<float>& out)
  {
    PeptDeepEncoder::Batch storage;
    const PeptDeepEncoder::Batch& batch = padded_(unpadded, storage);
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
    auto results = sessionFor_(session, rows, length).Run(Ort::RunOptions{nullptr}, impl_->input_ptrs.data(),
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
    if (info.GetElementCount() != static_cast<std::size_t>(rows) || info.GetShape().size() != 1)
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
                                       const PeptDeepEncoder::Batch& unpadded,
                                       const std::vector<std::size_t>& group,
                                       std::vector<Spectrum>& out)
  {
    {
      {
      PeptDeepEncoder::Batch storage;
      const PeptDeepEncoder::Batch& batch = padded_(unpadded, storage);
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
      auto results = sessionFor_(session, rows, length).Run(Ort::RunOptions{nullptr}, impl_->input_ptrs.data(),
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
                                       const PeptDeepEncoder::Batch& unpadded,
                                       const std::vector<std::size_t>& group,
                                       std::vector<float>& out)
  {
    PeptDeepEncoder::Batch storage;
    const PeptDeepEncoder::Batch& batch = padded_(unpadded, storage);
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
    auto results = sessionFor_(session, rows, length).Run(Ort::RunOptions{nullptr}, impl_->input_ptrs.data(),
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
