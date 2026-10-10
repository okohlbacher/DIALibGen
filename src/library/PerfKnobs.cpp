// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <odia/PerfKnobs.h>

#include <algorithm>
#include <atomic>
#include <exception>
#include <mutex>

#if defined(__linux__)
#include <malloc.h>
#include <sys/syscall.h>
#include <unistd.h>
#include <fstream>
#endif

namespace ODIA::perf
{
  namespace
  {
    std::atomic<unsigned> g_workers{0};
  }

  unsigned workerThreads()
  {
    const unsigned n = g_workers.load();
    return n > 0 ? n : std::max(1u, std::thread::hardware_concurrency());
  }

  void setWorkerThreads(unsigned n) { g_workers = n; }

  void parallelSlices(std::size_t n, std::size_t slices,
                      const std::function<void(std::size_t, std::size_t, std::size_t)>& fn)
  {
    slices = std::max<std::size_t>(1, std::min(slices, n));
    if (slices <= 1) { if (n > 0) { fn(0, 0, n); } return; }
    std::exception_ptr first;
    std::mutex m;
    const auto run = [&](std::size_t s) {
      const std::size_t b = n * s / slices, e = n * (s + 1) / slices;
      try { fn(s, b, e); }
      catch (...) { const std::lock_guard<std::mutex> lock(m); if (!first) { first = std::current_exception(); } }
    };
    std::vector<std::thread> threads;
    threads.reserve(slices - 1);
    for (std::size_t s = 1; s < slices; ++s) { threads.emplace_back(run, s); }
    run(0);
    for (auto& t : threads) { t.join(); }
    if (first) { std::rethrow_exception(first); }
  }

  void trimHeap()
  {
#if defined(__linux__) && defined(__GLIBC__)
    malloc_trim(0);
#endif
  }

  bool interleaveMemory()
  {
#if defined(__linux__) && defined(SYS_set_mempolicy)
    // Nodes from /sys, so that a machine with holes in its node numbering works.
    unsigned long mask[16] = {};
    std::ifstream in("/sys/devices/system/node/online");
    std::string list;
    if (!(in >> list)) { return false; }
    std::size_t pos = 0;
    while (pos < list.size())
    {
      const std::size_t comma = list.find(',', pos);
      const std::string part = list.substr(pos, comma == std::string::npos ? std::string::npos : comma - pos);
      const std::size_t dash = part.find('-');
      const unsigned lo = static_cast<unsigned>(std::stoul(part.substr(0, dash)));
      const unsigned hi = dash == std::string::npos ? lo : static_cast<unsigned>(std::stoul(part.substr(dash + 1)));
      for (unsigned node = lo; node <= hi && node < 1024; ++node) { mask[node / 64] |= 1ul << (node % 64); }
      if (comma == std::string::npos) { break; }
      pos = comma + 1;
    }
    constexpr int MPOL_INTERLEAVE_ = 3;
    return syscall(SYS_set_mempolicy, MPOL_INTERLEAVE_, mask, 1024ul + 1) == 0;
#else
    return false;
#endif
  }
}
