// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <odia/PeptDeepPredictor.h>

#include <odia/PeptDeepEncoder.h>
#include <odia/Trace.h>

#include <OpenMS/CHEMISTRY/AASequence.h>

#include <onnxruntime_cxx_api.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <iostream>
#include <mutex>
#include <thread>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#if defined(__linux__)
#include <pthread.h>
#include <sched.h>
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

  /// Numbers predictor calls and predictor instances in the trace (tracing only).
  std::atomic<std::int64_t> g_trace_calls{0};
  std::atomic<std::int64_t> g_trace_instances{0};

  /// Experiment knobs (perf/memory-cache). Each is read once; unset means the
  /// unchanged default behaviour, so a run without them is the baseline.
  ///   ODIA_BATCH_ROWS=<n>     rows per ONNX call instead of MAX_BATCH_ROWS;
  ///                           a 1-row tail is merged into the previous chunk
  ///                           (RT/CCS differ at the last bit for 1-row batches)
  ///   ODIA_SHARE_PREPACK=1    one PrepackedWeightsContainer per predictor, so
  ///                           the sessions share one copy of prepacked weights
  ///   ODIA_ARENA_SHRINK=1     RunOptions memory.enable_memory_arena_shrinkage
  ///   ODIA_PIN_SESSIONS=1     session s runs on the s-th CPU of the process
  ///                           affinity mask (numactl/taskset choose the CPUs)
  long knob(const char* name, long fallback)
  {
    const char* v = std::getenv(name);
    if (v == nullptr || *v == '\0') { return fallback; }
    char* end = nullptr;
    const long x = std::strtol(v, &end, 10);
    return (end != v && x > 0) ? x : fallback;
  }
  std::size_t batchRowsKnob()
  {
    static const long rows = knob("ODIA_BATCH_ROWS", 0);
    return static_cast<std::size_t>(rows);
  }
  bool sharePrepackKnob() { static const bool on = knob("ODIA_SHARE_PREPACK", 0) > 0; return on; }
  bool arenaShrinkKnob() { static const bool on = knob("ODIA_ARENA_SHRINK", 0) > 0; return on; }
  bool pinSessionsKnob() { static const bool on = knob("ODIA_PIN_SESSIONS", 0) > 0; return on; }

  /// The CPUs of the process affinity mask at first use, in ascending order.
  const std::vector<int>& allowedCpus()
  {
    static const std::vector<int> cpus = [] {
      std::vector<int> out;
#if defined(__linux__)
      cpu_set_t set;
      CPU_ZERO(&set);
      if (sched_getaffinity(0, sizeof(set), &set) == 0)
      {
        for (int c = 0; c < CPU_SETSIZE; ++c) { if (CPU_ISSET(c, &set)) { out.push_back(c); } }
      }
#endif
      return out;
    }();
    return cpus;
  }

  /// Pins the calling thread to one CPU; returns false where unsupported.
  bool pinTo(int cpu)
  {
#if defined(__linux__)
    cpu_set_t set;
    CPU_ZERO(&set);
    CPU_SET(cpu, &set);
    return pthread_setaffinity_np(pthread_self(), sizeof(set), &set) == 0;
#else
    (void)cpu;
    return false;
#endif
  }

  Ort::RunOptions makeRunOptions()
  {
    if (!arenaShrinkKnob()) { return Ort::RunOptions{nullptr}; }
    Ort::RunOptions options;
    options.AddConfigEntry("memory.enable_memory_arena_shrinkage", "cpu:0");
    return options;
  }
}

namespace ODIA
{
  namespace
  {
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
    Ort::Env env{ORT_LOGGING_LEVEL_WARNING, "odia"};
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
    /// RT, MS2 or CCS, from the inputs; names the model in trace events.
    const char* kind = "-";
    /// ODIA_ORT_PROFILE: every session writes ONNX Runtime's per-op profile.
    bool ort_profiling = false;
    /// ODIA_SHARE_PREPACK: one container of prepacked weights for all sessions.
    std::unique_ptr<Ort::PrepackedWeightsContainer> prepacked;
  };

  PeptDeepPredictor::PeptDeepPredictor(const std::string& model_path,
                                       bool prefer_gpu, int intra_op_threads,
                                       int sessions, int gpu_device)
    : impl_(std::make_unique<Impl>())
  {
    if (intra_op_threads > 0) { impl_->options.SetIntraOpNumThreads(intra_op_threads); }
    impl_->options.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_ALL);

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
    if (sharePrepackKnob() && !prefer_gpu) { impl_->prepacked = std::make_unique<Ort::PrepackedWeightsContainer>(); }
    const auto makeSession = [&](std::size_t index) {
      trace::Scope created("session_create", 3);
      created.fields().session = static_cast<std::int64_t>(index);
      if (!impl_->ort_profiling)
      {
        if (impl_->prepacked)
        {
          return std::make_unique<Ort::Session>(impl_->env, ortPath(model_path).c_str(), impl_->options,
                                                *impl_->prepacked);
        }
        return std::make_unique<Ort::Session>(impl_->env, ortPath(model_path).c_str(), impl_->options);
      }
      Ort::SessionOptions profiled = impl_->options.Clone();
      const std::string prefix = std::string(ort_profile) + "_" +
        std::filesystem::path(model_path).stem().string() + "_i" + std::to_string(instance) +
        "_s" + std::to_string(index);
      profiled.EnableProfiling(ortPath(prefix).c_str());
      return std::make_unique<Ort::Session>(impl_->env, ortPath(model_path).c_str(), profiled);
    };

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

      if (cuda_available)
      {
        // PROBE the devices rather than trusting device 0. On this project's
        // GPU nodes device 0 is dead, and AppendExecutionProvider_CUDA happily
        // accepts it -- the failure only surfaces when the session is built.
        // A silent fall back to CPU turns a 4-minute library into a 25-minute
        // one with nothing in the log to say why.
        const int first = gpu_device >= 0 ? gpu_device : 0;
        const int last  = gpu_device >= 0 ? gpu_device : 7;
        for (int dev = first; dev <= last; ++dev)
        {
          try
          {
            Ort::SessionOptions probe;
            probe.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_ALL);
            OrtCUDAProviderOptions cuda{};
            cuda.device_id = dev;
            probe.AppendExecutionProvider_CUDA(cuda);
            Ort::Session test(impl_->env, ortPath(model_path).c_str(), probe);   // the real check
            OrtCUDAProviderOptions keep{};
            keep.device_id = dev;
            impl_->options.AppendExecutionProvider_CUDA(keep);
            impl_->provider = Provider::CUDA;
            impl_->device = dev;
            break;
          }
          catch (const Ort::Exception&) { /* dead or absent: try the next */ }
        }
      }
    }

    try
    {
      impl_->session = makeSession(0);
    }
    catch (const Ort::Exception& e)
    {
      if (impl_->provider == Provider::CUDA)
      {
        // Registering the provider can succeed while creating the session with
        // it fails, so the fallback has to cover both.
        impl_->provider = Provider::CPU;
        impl_->device = -1;
        impl_->options = Ort::SessionOptions{};
        impl_->options.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_ALL);
        if (intra_op_threads > 0) { impl_->options.SetIntraOpNumThreads(intra_op_threads); }
        try
        {
          impl_->session = makeSession(0);
        }
        catch (const Ort::Exception& retry)
        { throw std::runtime_error("cannot open model " + model_path + ": " + retry.what()); }
      }
      else
      {
        throw std::runtime_error("cannot open model " + model_path + ": " + e.what());
      }
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
      std::string names;
      for (const auto& n : impl_->input_names) { names += (names.empty() ? "" : ", ") + n; }
      throw std::runtime_error(
        "model " + model_path + " is not a PeptDeep retention-time, CCS or MS2 model: "
        "found inputs (" + names + ")");
    }
    impl_->kind = looks_like_rt ? "RT" : looks_like_ms2 ? "MS2" : "CCS";

    // Replicas only after the model has been accepted, so a bad path or a
    // wrong model reports once rather than N times.
    //
    // CUDA is deliberately excluded: replicas exist to work around CPU
    // intra-op parallelism scaling badly, and on a GPU they would contend for
    // one device and multiply the weight copies in its memory instead.
    if (impl_->provider == Provider::CPU)
    {
      for (int i = 1; i < sessions; ++i)
      {
        impl_->replicas.push_back(makeSession(static_cast<std::size_t>(i)));
      }
    }
    if (trace::enabled())
    {
      trace::note(std::string("ort_version"), Ort::GetVersionString());
      trace::note(std::string("predictor.") + impl_->kind + ".sessions", std::to_string(sessionCount()));
      trace::note(std::string("knob.") + impl_->kind,
                  "batch_rows=" + std::to_string(batchRowsKnob() ? batchRowsKnob() : MAX_BATCH_ROWS) +
                  " share_prepack=" + std::to_string(impl_->prepacked ? 1 : 0) +
                  " arena_shrink=" + std::to_string(arenaShrinkKnob() ? 1 : 0) +
                  " pin_sessions=" + std::to_string(pinSessionsKnob() ? 1 : 0));
    }
  }

  int PeptDeepPredictor::device() const { return impl_->device; }

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
    // result depends on the batch size at the last bit (see MAX_BATCH_ROWS).
    trace::Scope chunking("chunking", 3);
    std::vector<std::vector<std::size_t>> chunks;
    const std::size_t knob_rows = batchRowsKnob();
    const std::size_t cap = knob_rows > 0 ? knob_rows : MAX_BATCH_ROWS;
    for (const auto& whole_group : PeptDeepEncoder::groupByLength(peptides))
    {
      for (std::size_t offset = 0; offset < whole_group.size(); offset += cap)
      {
        std::size_t stop = std::min(offset + cap, whole_group.size());
        // Knob only: never leave a 1-row tail behind (see ODIA_BATCH_ROWS).
        if (knob_rows > 0 && whole_group.size() - stop == 1) { stop = whole_group.size(); }
        chunks.emplace_back(
          whole_group.begin() + static_cast<std::ptrdiff_t>(offset),
          whole_group.begin() + static_cast<std::ptrdiff_t>(stop));
        if (stop == whole_group.size()) { break; }
      }
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
    std::mutex first_error_mutex;
    std::exception_ptr first_error;

    const auto worker = [&](std::size_t session)
    {
#if defined(__linux__)
      // ODIA_PIN_SESSIONS: session s always on the same CPU, so its arena
      // pages stay where they were first touched. The caller's mask is
      // restored after session 0's share of the work.
      cpu_set_t saved;
      bool restore = false;
      if (pinSessionsKnob() && !allowedCpus().empty())
      {
        if (session == 0) { restore = pthread_getaffinity_np(pthread_self(), sizeof(saved), &saved) == 0; }
        pinTo(allowedCpus()[session % allowedCpus().size()]);
      }
      struct Restore
      {
        bool on; cpu_set_t* set;
        ~Restore() { if (on) { pthread_setaffinity_np(pthread_self(), sizeof(*set), set); } }
      } restore_guard{restore, &saved};
#endif
      const double w0 = trace::now();
      const double wcpu0 = tracing ? trace::threadCpu() : 0.0;
      if (tracing && session > 0) { trace::nameWorker(impl_->kind, session); }
      std::int64_t chunks_done = 0, rows_done = 0;
      for (std::size_t i = next++; i < chunks.size(); i = next++)
      {
        // One thread's failure must stop the others rather than let them run
        // on for the rest of a two-million-peptide proteome. The first
        // exception is kept and rethrown on the calling thread.
        {
          const std::lock_guard<std::mutex> lock(first_error_mutex);
          if (first_error) { break; }
        }
        try
        {
          runChunk(session, i, rows_done);
          ++chunks_done;
        }
        catch (...)
        {
          const std::lock_guard<std::mutex> lock(first_error_mutex);
          if (!first_error) { first_error = std::current_exception(); }
          break;
        }
      }
      if (tracing) { recordWorker(session, w0, wcpu0, chunks_done, rows_done); }
    };

    std::vector<std::thread> threads;
    threads.reserve(workers - 1);
    for (std::size_t s = 1; s < workers; ++s) { threads.emplace_back(worker, s); }
    worker(0);
    for (auto& t : threads) { t.join(); }
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

    std::vector<const char*> in_names{impl_->input_names[0].c_str(),
                                      impl_->input_names[1].c_str()};
    std::vector<const char*> out_names{impl_->output_names[0].c_str()};

    trace::Scope running("ort_run", 4);
    tag(running, impl_->kind, session, group.size());
    auto results = impl_->at(session).Run(makeRunOptions(), in_names.data(),
                                       inputs.data(), inputs.size(),
                                       out_names.data(), out_names.size());
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
      std::vector<OpenMS::AASequence> subset;
      std::vector<int> subset_charges;
      for (const auto index : group)
      {
        subset.push_back(peptides[index]);
        subset_charges.push_back(charges[index]);
      }

      PeptDeepEncoder::Batch batch;
      try
      {
        batch = PeptDeepEncoder::encode(subset, subset_charges, nce, canonical_instrument);
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
              std::vector<OpenMS::AASequence>{subset[i]},
              std::vector<int>{subset_charges[i]}, nce, canonical_instrument);
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

      std::vector<const char*> in_names;
      for (const auto& n : impl_->input_names) { in_names.push_back(n.c_str()); }
      std::vector<const char*> out_names{impl_->output_names[0].c_str()};

      trace::Scope running("ort_run", 4);
      tag(running, impl_->kind, session, group.size());
      auto results = impl_->at(session).Run(makeRunOptions(), in_names.data(),
                                         inputs.data(), inputs.size(),
                                         out_names.data(), out_names.size());
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
      std::vector<OpenMS::AASequence> subset;
      subset.reserve(group.size());
      for (const auto index : group) { subset.push_back(peptides[index]); }

      PeptDeepEncoder::Batch batch;
      try
      {
        batch = PeptDeepEncoder::encode(subset);
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
            const auto single = PeptDeepEncoder::encode(
              std::vector<OpenMS::AASequence>{subset[i]});
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
      std::vector<OpenMS::AASequence> subset;
      std::vector<int> subset_charges;
      subset.reserve(group.size());
      subset_charges.reserve(group.size());
      for (const auto index : group)
      {
        subset.push_back(peptides[index]);
        subset_charges.push_back(charges[index]);
      }

      PeptDeepEncoder::Batch batch;
      try
      {
        // The CCS model takes no NCE and no instrument, but the encoder's meta
        // overload fills both; they are simply not bound below. Reusing it
        // keeps the charge validation and scaling in one place.
        batch = PeptDeepEncoder::encode(subset, subset_charges, 30.0f, "QE");
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
              std::vector<OpenMS::AASequence>{subset[i]},
              std::vector<int>{subset_charges[i]}, 30.0f, "QE");
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

    std::vector<const char*> in_names;
    for (const auto& n : impl_->input_names) { in_names.push_back(n.c_str()); }
    std::vector<const char*> out_names{impl_->output_names[0].c_str()};

    trace::Scope running("ort_run", 4);
    tag(running, impl_->kind, session, group.size());
    auto results = impl_->at(session).Run(makeRunOptions(), in_names.data(),
                                       inputs.data(), inputs.size(),
                                       out_names.data(), out_names.size());
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
