// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <cstdint>
#include <string>

namespace ODIA
{
  /// Physical CPU cores this process may run on: hardware threads, not
  /// counted twice for SMT siblings.
  ///
  /// What `-threads 0` means. Two inference sessions on the two hardware
  /// threads of one core share its arithmetic units and run at about half
  /// speed each (0.486 of one instance alone, measured on Zen 4), and 96
  /// sessions on all 384 hardware threads of a 192-core machine were 12 %
  /// slower than 48 on the physical cores.
  ///
  /// The effective budget, not the machine. Linux: the CPUs in the affinity
  /// mask (taskset; cgroup cpusets, i.e. cpuset.cpus.effective, show up
  /// there), counted once per core by their sysfs `thread_siblings_list`
  /// (without sysfs topology: the CPUs of the mask), and no more than the
  /// cgroup CPU quota allows (v2 cpu.max on the cgroup and its ancestors, v1
  /// cfs quota; rounded down, at least 1). macOS: `hw.physicalcpu`,
  /// performance and efficiency cores alike (batches are handed out one at a
  /// time, so a slower core takes fewer). Windows: the processor cores of all
  /// processor groups (GetLogicalProcessorInformationEx). If nothing can be
  /// read, std::thread::hardware_concurrency(). Computed once; at least 1.
  unsigned physicalCores();

  /// How physicalCores() was arrived at, for the log.
  const std::string& physicalCoresReason();

  /// macOS: count only the performance cores (`hw.perflevel0.physicalcpu`)
  /// in physicalCores(), what `-machine:cores performance` asks for. Saves
  /// CPU energy, not time (measured on an M4 Max, 12 + 4 cores: 12 threads
  /// against 16 were -3 % to +16 % wall and 18 % less CPU time). Must be
  /// called before the first physicalCores(); no effect on other platforms
  /// or when the performance-core count cannot be read. Does not change the
  /// output.
  void preferPerformanceCores(bool only);

  /// Bytes of last-level (L3) cache serving the CPUs this process may run
  /// on: each L3 counted once, however many of those CPUs share it. Linux:
  /// sysfs cache `index*` of level 3 for the CPUs in the affinity mask (a
  /// process pinned to two of a socket's twelve core complexes sees their
  /// two L3s); macOS: `hw.l3cachesize` (Intel; Apple silicon has none);
  /// Windows: every level-3 cache (GetLogicalProcessorInformationEx). 0 when
  /// unknown. Computed once.
  std::uint64_t lastLevelCacheBytes();
  /// NUMA nodes with memory that this machine has online (Linux sysfs); 1
  /// elsewhere or when it cannot be read.
  unsigned numaNodes();

  /// Interleave this thread's later page allocations over all NUMA nodes, as
  /// `numactl --interleave=all` does; threads created afterwards inherit it.
  /// Linux only; false where unsupported or refused.
  bool interleaveMemory();
}
