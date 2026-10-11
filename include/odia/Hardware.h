// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

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
  /// Linux: the CPUs in the affinity mask (taskset, cgroup cpusets), counted
  /// once per core by their sysfs `thread_siblings_list`. macOS:
  /// `hw.physicalcpu`. Windows: the processor cores reported by
  /// GetLogicalProcessorInformation. Anything that cannot be read falls back
  /// to std::thread::hardware_concurrency(). Computed once; at least 1.
  unsigned physicalCores();

  /// NUMA nodes with memory that this machine has online (Linux sysfs); 1
  /// elsewhere or when it cannot be read.
  unsigned numaNodes();

  /// Interleave this thread's later page allocations over all NUMA nodes, as
  /// `numactl --interleave=all` does; threads created afterwards inherit it.
  /// Linux only; false where unsupported or refused.
  bool interleaveMemory();
}
