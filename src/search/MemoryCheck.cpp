// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <odia/search/MemoryCheck.h>

#include <algorithm>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>

#if defined(_WIN32)
#  ifndef NOMINMAX
#    define NOMINMAX
#  endif
#  ifndef WIN32_LEAN_AND_MEAN
#    define WIN32_LEAN_AND_MEAN
#  endif
#  include <windows.h>
#elif defined(__APPLE__)
#  include <mach/mach.h>
#  include <mach/mach_host.h>
#  include <unistd.h>
#elif defined(__linux__)
#  include <malloc.h>   // malloc_trim (glibc)
#  include <unistd.h>
#endif

namespace ODIA::search
{
  std::string memorySize(std::uint64_t bytes)
  {
    char buffer[32];
    if (bytes < 1000000000ull) { std::snprintf(buffer, sizeof buffer, "%.1f MB", static_cast<double>(bytes) / 1e6); }
    else { std::snprintf(buffer, sizeof buffer, "%.1f GB", static_cast<double>(bytes) / 1e9); }
    return buffer;
  }

  namespace
  {
    std::string gigabytes(std::uint64_t bytes) { return memorySize(bytes); }

#if defined(__linux__)
    /// Resident bytes of this process (/proc/self/statm), 0 when unreadable.
    std::uint64_t residentBytes()
    {
      std::ifstream statm("/proc/self/statm");
      std::uint64_t size = 0, resident = 0;
      if (!(statm >> size >> resident)) { return 0; }
      const long page = sysconf(_SC_PAGESIZE);
      return page > 0 ? resident * static_cast<std::uint64_t>(page) : 0;
    }
#endif

#if defined(__linux__)
    /// The cgroup v2 limit's room (memory.max - memory.current), or 0 when
    /// there is no such limit (none set, cgroup v1, or not readable).
    std::uint64_t cgroupRoom()
    {
      std::ifstream self("/proc/self/cgroup");
      std::string line, path;
      while (std::getline(self, line))
      {
        if (line.rfind("0::", 0) == 0) { path = line.substr(3); break; }
      }
      if (path.empty()) { return 0; }
      const std::string dir = "/sys/fs/cgroup" + (path == "/" ? std::string() : path);
      std::ifstream max_file(dir + "/memory.max"), current_file(dir + "/memory.current");
      std::string max_text;
      std::uint64_t max = 0, current = 0;
      if (!(max_file >> max_text) || max_text == "max") { return 0; }
      try { max = std::stoull(max_text); }
      catch (...) { return 0; }
      if (!(current_file >> current)) { return 0; }
      return max > current ? max - current : 1;   // 1: a limit that is used up is still a limit
    }
#endif
  }

  std::uint64_t availableMemoryBytes()
  {
#if defined(_WIN32)
    MEMORYSTATUSEX status;
    status.dwLength = sizeof status;
    if (!GlobalMemoryStatusEx(&status)) { return 0; }
    return static_cast<std::uint64_t>(status.ullAvailPhys);
#elif defined(__APPLE__)
    vm_statistics64_data_t vm;
    mach_msg_type_number_t count = HOST_VM_INFO64_COUNT;
    if (host_statistics64(mach_host_self(), HOST_VM_INFO64, reinterpret_cast<host_info64_t>(&vm), &count) != KERN_SUCCESS)
    { return 0; }
    const long page = sysconf(_SC_PAGESIZE);
    if (page <= 0) { return 0; }
    // free_count includes the speculative pages (vm_stat lists them apart).
    return (static_cast<std::uint64_t>(vm.free_count) + static_cast<std::uint64_t>(vm.inactive_count)) *
           static_cast<std::uint64_t>(page);
#elif defined(__linux__)
    std::ifstream meminfo("/proc/meminfo");
    std::string key, unit;
    std::uint64_t value = 0, available = 0;
    while (meminfo >> key >> value)
    {
      std::getline(meminfo, unit);
      if (key == "MemAvailable:") { available = value * 1024; break; }
    }
    if (available == 0) { return 0; }
    const std::uint64_t room = cgroupRoom();
    return room > 0 ? std::min(available, room) : available;
#else
    return 0;
#endif
  }

  std::uint64_t releaseFreedMemory()
  {
#if defined(__linux__)
    const std::uint64_t before = residentBytes();
#  if defined(__GLIBC__)
    malloc_trim(0);
#  endif
    const std::uint64_t after = residentBytes();
    return before > after ? before - after : 0;
#else
    return 0;
#endif
  }

  Ms1MemoryDecision decideMs1InMemory(Ms1InMemory mode, std::uint64_t estimate, std::uint64_t available)
  {
    const std::string figures = "estimated " + gigabytes(estimate) + " in memory, " +
                                (available > 0 ? gigabytes(available) + " free" : std::string("free memory unknown"));
    if (estimate == 0) { return {false, "no MS1 spectra to hold"}; }
    if (mode == Ms1InMemory::On) { return {true, "search:ms1_in_memory on (" + figures + ")"}; }
    if (mode == Ms1InMemory::Off) { return {false, "search:ms1_in_memory off (" + figures + ")"}; }
    char factor[16];
    std::snprintf(factor, sizeof factor, "%g", ms1_in_memory_headroom);
    const double needed = ms1_in_memory_headroom * static_cast<double>(estimate);
    if (available == 0)
    { return {false, "search:ms1_in_memory auto: " + figures + ", so streamed; search:ms1_in_memory on holds it anyway"}; }
    if (static_cast<double>(available) >= needed)
    { return {true, "search:ms1_in_memory auto: " + figures + ", at least " + factor + "x the map"}; }
    return {false, "search:ms1_in_memory auto: " + figures + ", less than " + factor + "x the map (" +
                     gigabytes(static_cast<std::uint64_t>(needed)) + "), so streamed; search:ms1_in_memory on holds it anyway"};
  }
}
