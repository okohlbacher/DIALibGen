// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <odia/Hardware.h>

#include <algorithm>
#include <thread>

#include <string>
#include <vector>

#if defined(__linux__)
#include <cerrno>
#include <fstream>
#include <set>
#include <sstream>
#include <stdexcept>
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
#endif

namespace ODIA
{
  namespace
  {
    struct CoreCount
    {
      unsigned cores = 0;       ///< 0: could not be determined
      std::string reason;
    };

#if defined(__linux__)
    /// CPUs the process may run on, as a dynamically sized set: a fixed
    /// cpu_set_t holds 1,024 CPUs, and sched_getaffinity fails with EINVAL
    /// on a kernel configured for more.
    std::vector<int> affinityCpus()
    {
      for (int n = 1024; n <= (1 << 16); n *= 2)
      {
        cpu_set_t* set = CPU_ALLOC(n);
        if (set == nullptr) { return {}; }
        const std::size_t size = CPU_ALLOC_SIZE(n);
        CPU_ZERO_S(size, set);
        if (sched_getaffinity(0, size, set) == 0)
        {
          std::vector<int> cpus;
          for (int cpu = 0; cpu < n; ++cpu) { if (CPU_ISSET_S(cpu, size, set)) { cpus.push_back(cpu); } }
          CPU_FREE(set);
          return cpus;
        }
        CPU_FREE(set);
        if (errno != EINVAL) { break; }
      }
      return {};
    }

    std::string readLine(const std::string& path)
    {
      std::ifstream in(path);
      std::string line;
      std::getline(in, line);
      return line;
    }

    /// CPUs' worth of time the process's cgroups allow per period (cgroup v2
    /// cpu.max on the process's cgroup and its ancestors, or v1
    /// cpu.cfs_quota_us / cpu.cfs_period_us), rounded down and at least 1;
    /// 0 when there is no limit. A container limited with --cpus=4 sees every
    /// CPU of the host in its affinity mask but gets four CPUs' worth of time.
    unsigned cgroupCpuQuota()
    {
      unsigned limit = 0;
      const auto take = [&](long long quota, long long period) {
        if (quota <= 0 || period <= 0) { return; }
        const auto cpus = static_cast<unsigned>(std::max<long long>(1, quota / period));
        limit = limit == 0 ? cpus : std::min(limit, cpus);
      };
      std::ifstream groups("/proc/self/cgroup");
      std::string line;
      while (std::getline(groups, line))
      {
        // "0::/user.slice/..." (v2) or "4:cpu,cpuacct:/..." (v1)
        const auto first = line.find(':');
        const auto second = first == std::string::npos ? first : line.find(':', first + 1);
        if (second == std::string::npos) { continue; }
        const std::string controllers = line.substr(first + 1, second - first - 1);
        std::string path = line.substr(second + 1);
        if (controllers.empty())
        {
          // v2: the limit of any ancestor applies too.
          for (;;)
          {
            std::istringstream fields(readLine("/sys/fs/cgroup" + path + "/cpu.max"));
            std::string quota;
            long long period = 0;
            if ((fields >> quota >> period) && quota != "max")
            {
              try { take(std::stoll(quota), period); } catch (const std::exception&) {}
            }
            if (path.empty() || path == "/") { break; }
            const auto slash = path.find_last_of('/');
            path = slash == 0 || slash == std::string::npos ? std::string("/") : path.substr(0, slash);
          }
        }
        else if ((',' + controllers + ',').find(",cpu,") != std::string::npos)
        {
          for (const char* mount : {"/sys/fs/cgroup/cpu,cpuacct", "/sys/fs/cgroup/cpu"})
          {
            const std::string quota = readLine(mount + path + "/cpu.cfs_quota_us");
            const std::string period = readLine(mount + path + "/cpu.cfs_period_us");
            if (quota.empty() || period.empty()) { continue; }
            try { take(std::stoll(quota), std::stoll(period)); } catch (const std::exception&) {}
            break;
          }
        }
      }
      return limit;
    }
#endif

    CoreCount countPhysicalCores()
    {
      CoreCount c;
#if defined(__linux__)
      const auto cpus = affinityCpus();
      if (cpus.empty()) { return c; }
      // One entry per core: its siblings list ("0,192" or "0-1") is the same
      // string for every hardware thread of the core.
      std::set<std::string> cores;
      bool topology = true;
      for (const int cpu : cpus)
      {
        std::ifstream in("/sys/devices/system/cpu/cpu" + std::to_string(cpu) + "/topology/thread_siblings_list");
        std::string siblings;
        if (!(in >> siblings)) { topology = false; break; }
        cores.insert(siblings);
      }
      if (topology)
      {
        c.cores = static_cast<unsigned>(cores.size());
        c.reason = "physical cores of the " + std::to_string(cpus.size()) + " CPUs in the affinity mask";
      }
      else
      {
        c.cores = static_cast<unsigned>(cpus.size());
        c.reason = "CPUs in the affinity mask (no core topology in sysfs)";
      }
      if (const unsigned quota = cgroupCpuQuota(); quota > 0 && quota < c.cores)
      {
        c.reason += ", limited to " + std::to_string(quota) + " by the cgroup CPU quota";
        c.cores = quota;
      }
#elif defined(__APPLE__)
      int cores = 0;
      std::size_t size = sizeof cores;
      if (sysctlbyname("hw.physicalcpu", &cores, &size, nullptr, 0) != 0 || cores < 1) { return c; }
      c.cores = static_cast<unsigned>(cores);
      c.reason = "hw.physicalcpu";
      // Performance and efficiency cores alike: chunks are handed out one at
      // a time, so a slower core simply takes fewer of them.
      int performance = 0;
      size = sizeof performance;
      if (sysctlbyname("hw.perflevel0.physicalcpu", &performance, &size, nullptr, 0) == 0 && performance > 0 &&
          performance < cores)
      { c.reason += " (" + std::to_string(performance) + " performance and " + std::to_string(cores - performance) +
                    " efficiency cores)"; }
#elif defined(_WIN32)
      // ...Ex and RelationProcessorCore: every processor group, not only the
      // calling thread's (at most 64 logical processors).
      DWORD bytes = 0;
      GetLogicalProcessorInformationEx(RelationProcessorCore, nullptr, &bytes);
      if (bytes == 0) { return c; }
      std::vector<char> buffer(bytes);
      auto* first = reinterpret_cast<SYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX*>(buffer.data());
      if (!GetLogicalProcessorInformationEx(RelationProcessorCore, first, &bytes)) { return c; }
      unsigned cores = 0;
      for (DWORD offset = 0; offset < bytes;)
      {
        const auto* entry = reinterpret_cast<const SYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX*>(buffer.data() + offset);
        if (entry->Relationship == RelationProcessorCore) { ++cores; }
        if (entry->Size == 0) { break; }
        offset += entry->Size;
      }
      c.cores = cores;
      c.reason = "processor cores of all processor groups";
#endif
      return c;
    }

    const CoreCount& counted()
    {
      static const CoreCount c = [] {
        CoreCount r = countPhysicalCores();
        if (r.cores == 0)
        {
          r.cores = std::max(1u, std::thread::hardware_concurrency());
          r.reason = "hardware threads (the core count could not be read)";
        }
        return r;
      }();
      return c;
    }
  } // namespace

  unsigned physicalCores() { return counted().cores; }

  const std::string& physicalCoresReason() { return counted().reason; }

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
