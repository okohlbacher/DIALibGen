// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

// Switches for the output-neutral speed-ups prototyped on perf/generate-solutions.
//
// Every switch is read once per process from the environment. Unset (or 0),
// each one leaves the behaviour exactly as it was, so one binary can measure
// every change alone and in combination. ODIA_SOLUTIONS=1 turns on the
// recommended combination; an explicit variable still overrides it.
//
//   ODIA_MAX_SESSIONS=N  session cap instead of 16 (generate and -run);
//                        0 means "as many as -threads".
//   ODIA_MS2_ROWS=R      MS2 chunk rows instead of 2048 (MS2 is bit-identical
//                        for any split, including a single row).
//   ODIA_LSTM_ROWS=R     RT/CCS chunks of about R rows, cut out of the 2048-row
//                        chunks so that no chunk of one row appears that was
//                        not already one row (RT/CCS are bit-identical for any
//                        batch of 2 or more rows; a 1-row batch differs).
//   ODIA_MS2_BLOCK=N     precursors per MS2 call instead of 20000.
//   ODIA_LPT=1           dispatch chunks longest peptide first (order only).
//   ODIA_INTERLEAVE=1    interleave all later allocations over the NUMA nodes
//                        (what numactl --interleave=all does); under
//                        ODIA_SOLUTIONS on only above 32 threads.
//   ODIA_PAR_POST=1      MS2 post-processing (annotation, floor, ranking) on
//                        all worker threads, concatenated in order.
//   ODIA_PAR_ENUM=1      digest and precursor enumeration on all worker
//                        threads, merged in input order.
//   ODIA_PAR_WRITE=1     Parquet arrays built per column in parallel and
//                        column chunks encoded in parallel; TSV rows formatted
//                        in parallel and written in order.
//   ODIA_TRIM=1          malloc_trim(0) after every prediction phase.
//   ODIA_PAR_CREATE=1    ONNX sessions 1..n-1 created concurrently.
//   ODIA_MS2_PRUNED=path run this MS2 graph (modloss branch removed, channels
//                        4-7 zero) in place of the configured one; the
//                        configured file is still the one hashed into the
//                        library's provenance (prototype only).

#include <cstddef>
#include <cstdlib>
#include <functional>
#include <string>
#include <thread>
#include <vector>

namespace ODIA::perf
{
  /// Integer value of @p name, or @p fallback when unset, empty or not a number.
  inline long knob(const char* name, long fallback)
  {
    const char* v = std::getenv(name);
    if (v == nullptr || *v == '\0') { return fallback; }
    char* end = nullptr;
    const long x = std::strtol(v, &end, 10);
    return end != v ? x : fallback;
  }

  inline bool solutions() { static const bool v = knob("ODIA_SOLUTIONS", 0) > 0; return v; }

  /// A switch: explicit setting wins, otherwise on when ODIA_SOLUTIONS is set.
  inline bool flag(const char* name) { return knob(name, solutions() ? 1 : 0) > 0; }

  /// Session cap: 16 by default; 0 (or ODIA_SOLUTIONS) means no cap beyond -threads.
  inline unsigned sessionCap()
  {
    static const unsigned v = [] {
      const long x = knob("ODIA_MAX_SESSIONS", solutions() ? 0 : 16);
      return x <= 0 ? 0u : static_cast<unsigned>(x);
    }();
    return v;
  }
  inline std::size_t ms2Rows()   { static const long v = knob("ODIA_MS2_ROWS", solutions() ? 64 : 0); return v > 0 ? std::size_t(v) : 0; }
  inline std::size_t lstmRows()  { static const long v = knob("ODIA_LSTM_ROWS", solutions() ? 64 : 0); return v > 0 ? std::size_t(v) : 0; }
  inline std::size_t ms2Block()  { static const long v = knob("ODIA_MS2_BLOCK", solutions() ? 100000 : 20000); return v > 0 ? std::size_t(v) : 20000; }
  inline bool lpt()              { static const bool v = flag("ODIA_LPT"); return v; }
  /// Interleave: explicit 0/1 wins; under ODIA_SOLUTIONS it is on only above
  /// 32 worker threads. Measured on sub10: +6 % at 16 sessions (2048 rows), but
  /// at 96 sessions it spreads the busy threads over both sockets.
  inline bool interleave(unsigned threads)
  {
    static const long v = knob("ODIA_INTERLEAVE", -1);
    if (v >= 0) { return v > 0; }
    return solutions() && (threads == 0 || threads > 32);
  }
  inline bool parPost()          { static const bool v = flag("ODIA_PAR_POST"); return v; }
  inline bool parEnum()          { static const bool v = flag("ODIA_PAR_ENUM"); return v; }
  inline bool parWrite()         { static const bool v = flag("ODIA_PAR_WRITE"); return v; }
  inline bool parCreate()        { static const bool v = flag("ODIA_PAR_CREATE"); return v; }
  inline bool trim()             { static const bool v = flag("ODIA_TRIM"); return v; }
  inline std::string ms2Pruned() { const char* v = std::getenv("ODIA_MS2_PRUNED"); return v ? v : ""; }

  /// Worker threads for the CPU-side parallel loops (-threads; set by the tool).
  unsigned workerThreads();
  void setWorkerThreads(unsigned n);

  /// Runs fn(slice, begin, end) over @p n items cut into @p slices contiguous
  /// slices, one thread per slice (the caller runs slice 0). The first
  /// exception is rethrown on the caller.
  void parallelSlices(std::size_t n, std::size_t slices,
                      const std::function<void(std::size_t, std::size_t, std::size_t)>& fn);

  /// Releases free heap memory back to the system (glibc malloc_trim; no-op elsewhere).
  void trimHeap();

  /// Interleaves later page allocations over all NUMA nodes (Linux; no-op elsewhere).
  bool interleaveMemory();
}
