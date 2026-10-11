// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <odia/Hardware.h>

#include "ThreadPool.h"

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <exception>
#include <limits>
#include <mutex>
#include <system_error>
#include <thread>
#include <vector>

namespace ODIA::parallel
{
  /// Worker threads for host-side work: `requested`, or every physical core
  /// for 0 (the meaning of -threads 0; see physicalCores()).
  inline unsigned workers(unsigned requested)
  {
    return requested != 0 ? requested : physicalCores();
  }

  /// Calls body(i) for every i in [0, n) on up to `threads` threads, the
  /// calling thread included.
  ///
  /// This is how the host-side stages stay byte-identical at any thread count:
  /// body(i) must write only to slot i of a result the caller sized up front,
  /// and the caller then consumes the slots in index order. Nothing depends on
  /// which thread ran which index.
  ///
  /// If bodies throw, the exception of the LOWEST failing index is rethrown
  /// once all threads have stopped -- the one a serial loop would have hit
  /// first -- so the error message does not depend on scheduling either.
  ///
  /// Indices are handed out @p grain at a time, in order within a grain: one
  /// shared counter increment per index costs more than a body of a few
  /// hundred nanoseconds once dozens of threads contend for its cache line.
  template <class Body>
  void forEach(std::size_t n, unsigned threads, Body&& body, std::size_t grain = 1)
  {
    grain = std::max<std::size_t>(1, grain);
    const std::size_t grains = (n + grain - 1) / grain;
    const std::size_t count = std::min<std::size_t>(std::max(1u, threads), grains);
    if (count <= 1)
    {
      for (std::size_t i = 0; i < n; ++i) { body(i); }
      return;
    }

    std::atomic<std::size_t> next{0};
    std::atomic<std::size_t> stop_at{std::numeric_limits<std::size_t>::max()};
    std::mutex error_mutex;
    std::exception_ptr error;
    std::size_t error_index = std::numeric_limits<std::size_t>::max();

    const auto run = [&]
    {
      for (std::size_t g = next.fetch_add(1, std::memory_order_relaxed); g < grains;
           g = next.fetch_add(1, std::memory_order_relaxed))
      {
        const std::size_t end = std::min(n, (g + 1) * grain);
        for (std::size_t i = g * grain; i < end && i < stop_at.load(std::memory_order_relaxed); ++i)
        {
          try { body(i); }
          catch (...)
          {
            const std::lock_guard<std::mutex> lock(error_mutex);
            if (i < error_index) { error_index = i; error = std::current_exception(); }
            // Indices below i are still run, so a lower failure is still found.
            std::size_t current = stop_at.load();
            while (i < current && !stop_at.compare_exchange_weak(current, i)) {}
          }
        }
        if (g * grain >= stop_at.load(std::memory_order_relaxed)) { break; }
      }
    };

    // On the calling thread and up to count - 1 pool workers; fewer when no
    // more can be started, which changes only the time.
    team(count, [&](std::size_t) { run(); });
    if (error) { std::rethrow_exception(error); }
  }
} // namespace ODIA::parallel
