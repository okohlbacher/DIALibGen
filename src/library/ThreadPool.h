// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <cstddef>
#include <functional>
#include <future>

namespace ODIA::parallel
{
  /// The process's worker threads, created on first need and then kept.
  ///
  /// Until 0.11.2 every parallel step started threads of its own -- about
  /// 5,600 per whole-proteome run -- and each new thread could land on any
  /// core and allocate from a new glibc arena. Pool workers stay: worker k
  /// is bound to NUMA node k mod n (n = the nodes this process may run on;
  /// Linux, unless ODIA_POOL_PIN=0), so consecutive workers alternate
  /// between the sockets, and the memory a worker first touches -- an ONNX
  /// session's arena, say -- stays on its node.
  ///
  /// team() hands out the lowest free workers, so the sessions of a
  /// predictor (which pass `low`) meet the same workers call after call;
  /// other teams take the highest free ones. A worker that cannot be
  /// started (a process or cgroup thread limit) leaves the team smaller;
  /// callers distribute work dynamically, so that changes only the time.

  /// Calls fn(i) for i in [0, m), where m <= n: i = 0 on the calling thread,
  /// the others on pool workers. Returns m (at least 1). If calls throw, the
  /// first exception is rethrown once all have returned.
  std::size_t team(std::size_t n, const std::function<void(std::size_t)>& fn, bool low = false);

  /// Runs @p task on a pool worker; on this thread, before returning, when
  /// no worker can be started.
  std::future<void> async(std::function<void()> task);

  /// Workers created so far (for the trace).
  std::size_t poolSize();
} // namespace ODIA::parallel
