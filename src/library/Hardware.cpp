// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <odia/Hardware.h>

#include <algorithm>
#include <thread>

#if defined(__linux__)
#include <fstream>
#include <set>
#include <stdexcept>
#include <string>
#include <sched.h>
#include <sys/syscall.h>
#include <unistd.h>
#elif defined(__APPLE__)
#include <sys/sysctl.h>
#include <sys/types.h>
#elif defined(_WIN32)
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#include <vector>
#endif

namespace ODIA
{
  namespace
  {
    unsigned countPhysicalCores()
    {
#if defined(__linux__)
      cpu_set_t mask;
      CPU_ZERO(&mask);
      if (sched_getaffinity(0, sizeof mask, &mask) != 0) { return 0; }
      // One entry per core: its siblings list ("0,192" or "0-1") is the same
      // string for every hardware thread of the core.
      std::set<std::string> cores;
      for (int cpu = 0; cpu < CPU_SETSIZE; ++cpu)
      {
        if (!CPU_ISSET(cpu, &mask)) { continue; }
        std::ifstream in("/sys/devices/system/cpu/cpu" + std::to_string(cpu) + "/topology/thread_siblings_list");
        std::string siblings;
        if (!(in >> siblings)) { return 0; }   // no topology: let the caller fall back
        cores.insert(siblings);
      }
      return static_cast<unsigned>(cores.size());
#elif defined(__APPLE__)
      int cores = 0;
      std::size_t size = sizeof cores;
      if (sysctlbyname("hw.physicalcpu", &cores, &size, nullptr, 0) != 0 || cores < 1) { return 0; }
      return static_cast<unsigned>(cores);
#elif defined(_WIN32)
      DWORD bytes = 0;
      GetLogicalProcessorInformation(nullptr, &bytes);
      if (bytes == 0) { return 0; }
      std::vector<SYSTEM_LOGICAL_PROCESSOR_INFORMATION> info(bytes / sizeof(SYSTEM_LOGICAL_PROCESSOR_INFORMATION));
      if (!GetLogicalProcessorInformation(info.data(), &bytes)) { return 0; }
      unsigned cores = 0;
      for (const auto& entry : info) { if (entry.Relationship == RelationProcessorCore) { ++cores; } }
      return cores;
#else
      return 0;
#endif
    }
  } // namespace

  unsigned physicalCores()
  {
    static const unsigned cores = [] {
      const unsigned counted = countPhysicalCores();
      return counted > 0 ? counted : std::max(1u, std::thread::hardware_concurrency());
    }();
    return cores;
  }

  namespace
  {
#if defined(__linux__)
    /// The node numbers in a sysfs list such as "0-1,3".
    std::set<unsigned> onlineNodes()
    {
      std::set<unsigned> nodes;
      std::ifstream in("/sys/devices/system/node/has_memory");
      std::string list;
      if (!(in >> list))
      {
        std::ifstream online("/sys/devices/system/node/online");
        if (!(online >> list)) { return nodes; }
      }
      std::size_t pos = 0;
      try
      {
        while (pos < list.size())
        {
          const std::size_t comma = list.find(',', pos);
          const std::string part = list.substr(pos, comma == std::string::npos ? std::string::npos : comma - pos);
          const std::size_t dash = part.find('-');
          const unsigned lo = static_cast<unsigned>(std::stoul(part.substr(0, dash)));
          const unsigned hi = dash == std::string::npos ? lo : static_cast<unsigned>(std::stoul(part.substr(dash + 1)));
          for (unsigned node = lo; node <= hi && node < 1024; ++node) { nodes.insert(node); }
          if (comma == std::string::npos) { break; }
          pos = comma + 1;
        }
      }
      catch (const std::exception&) { nodes.clear(); }
      return nodes;
    }
#endif
  } // namespace

  unsigned numaNodes()
  {
#if defined(__linux__)
    const auto nodes = onlineNodes();
    return nodes.empty() ? 1u : static_cast<unsigned>(nodes.size());
#else
    return 1;
#endif
  }

  bool interleaveMemory()
  {
#if defined(__linux__) && defined(SYS_set_mempolicy)
    const auto nodes = onlineNodes();
    if (nodes.size() < 2) { return false; }
    unsigned long mask[1024 / (8 * sizeof(unsigned long))] = {};
    for (const unsigned node : nodes) { mask[node / (8 * sizeof(unsigned long))] |= 1ul << (node % (8 * sizeof(unsigned long))); }
    constexpr int MPOL_INTERLEAVE_ = 3;   // <linux/mempolicy.h>, without libnuma
    return ::syscall(SYS_set_mempolicy, MPOL_INTERLEAVE_, mask, 1024ul + 1) == 0;
#else
    return false;
#endif
  }
} // namespace ODIA
