// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <cstddef>
#include <cstdint>
#include <string>

/// Opt-in performance tracing, for finding where the wall and CPU time go.
///
/// Off unless the environment variable ODIA_TRACE names an output directory.
/// When it is unset every entry point below returns after one test of a
/// boolean that is read once per process, nothing is allocated, no thread is
/// started and no file is touched -- the library a run writes does not depend
/// on whether it was traced. docs/tracing.md documents the variables and the
/// columns of every file written.
///
/// When it is set:
///   - Scope objects record events (begin, end, thread, thread CPU, and a
///     few integer fields such as rows and session) into a per-thread buffer,
///     which is handed over once, when the thread exits;
///   - Session starts one sampler thread that reads cheap, lock-free /proc
///     counters (Linux only) every ODIA_TRACE_INTERVAL_MS (default 100 ms),
///     and writes everything as TSV when it ends;
///   - the expensive views -- glibc malloc_info, /proc/self/numa_maps and
///     /proc/self/smaps_rollup, which hold the allocator's arena locks or
///     the mmap lock while they walk -- are read only when a phase (level 0)
///     begins, at most once per ODIA_TRACE_SNAPSHOT_MIN_S (default 30 s), on
///     the thread that begins it, and each read is timed.
namespace ODIA::trace
{
  namespace detail
  {
    bool readEnabled() noexcept;
  }

  /// Whether ODIA_TRACE is set. Read once; constant for the life of the process.
  inline bool enabled() noexcept
  {
    static const bool on = detail::readEnabled();
    return on;
  }

  /// Seconds since the trace epoch (steady clock). 0 when tracing is off.
  double now() noexcept;

  /// Optional per-event fields. -1 (and "-") mean "not applicable".
  struct Fields
  {
    std::int64_t call = -1;      ///< predictor call number (one per predictRT/MS2/CCS call)
    std::int64_t idx = -1;       ///< block or chunk index within its parent
    std::int64_t session = -1;   ///< ONNX session (= worker) index
    std::int64_t rows = -1;      ///< peptides / precursors / entries handled
    std::int64_t len = -1;       ///< residues per peptide (chunks are length-homogeneous)
    double wait_s = -1.0;        ///< queue wait: dispatch -> dequeue (chunks) or -> start (workers)
    const char* model = "-";     ///< RT, MS2 or CCS for predictor events
  };

  /// A timed region. Level 0 is a phase of the run: while it is open, events
  /// recorded on ANY thread are filed under its name. Higher levels nest.
  class Scope
  {
  public:
    /// Not started; see start().
    Scope() noexcept = default;
    explicit Scope(const char* name, int level = 1) noexcept
    {
      if (enabled()) { begin_(name, level); }
    }
    /// Start a default-constructed scope later, e.g. one declared early so
    /// that it outlives the locals whose destruction it should time.
    void start(const char* name, int level = 1) noexcept
    {
      if (!active_ && enabled()) { begin_(name, level); }
    }
    ~Scope() { if (active_) { end_(); } }
    Scope(const Scope&) = delete;
    Scope& operator=(const Scope&) = delete;

    bool active() const noexcept { return active_; }
    /// Fields recorded with the event when it ends.
    Fields& fields() noexcept { return fields_; }
    /// End now rather than at destruction.
    void end() noexcept { if (active_) { end_(); } }

  private:
    void begin_(const char* name, int level) noexcept;
    void end_() noexcept;

    bool active_ = false;
    int level_ = 0;
    const char* name_ = nullptr;
    const char* outer_phase_ = nullptr;
    double t0_ = 0.0;
    double cpu0_ = 0.0;
    Fields fields_;
  };

  /// Record an event whose times were taken by the caller (seconds from now()).
  void record(const char* name, int level, double t0, double t1, double thread_cpu_s,
              const Fields& fields) noexcept;

  /// CPU seconds consumed by the calling thread so far (NaN if unavailable).
  double threadCpu() noexcept;

  /// Hand the calling thread's events over now rather than at its exit: for
  /// threads that outlive the trace session (the worker pool).
  void flushThread() noexcept;

  /// Name the calling thread "odia-<model>-s<session>" (Linux; tracing only),
  /// so the sampler's per-thread rows say which worker a thread was.
  void nameWorker(const char* model, std::size_t session) noexcept;

  /// Add a key/value row to meta.tsv (e.g. which ONNX profile file a session wrote).
  void note(const std::string& key, const std::string& value);

  /// The trace directory, or "" when tracing is off.
  const std::string& directory();

  /// Starts the samplers and, on destruction, stops them and writes the files.
  /// One per process; does nothing when tracing is off.
  class Session
  {
  public:
    Session(int argc, const char** argv);
    ~Session();
    Session(const Session&) = delete;
    Session& operator=(const Session&) = delete;

  private:
    bool owner_ = false;
  };

} // namespace ODIA::trace
