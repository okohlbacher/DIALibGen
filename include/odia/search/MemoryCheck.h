// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

/// How much memory the machine can give, and whether extraction holds the MS1
/// map in memory (`search:ms1_in_memory`). See
/// docs/design/built-in-identification.md, Resources.
#pragma once

#include <odia/search/SearchParams.h>

#include <cstdint>
#include <string>

namespace ODIA::search
{
  /// Bytes of physical memory this process can still get without swapping, or
  /// 0 when the system does not say:
  /// - Linux: MemAvailable of /proc/meminfo, and no more than the cgroup v2
  ///   limit leaves (memory.max - memory.current) when one is set;
  /// - macOS: free plus inactive pages (vm_stat's "free", "speculative" and
  ///   "inactive"), times the page size;
  /// - Windows: ullAvailPhys of GlobalMemoryStatusEx.
  std::uint64_t availableMemoryBytes();

  /// Gives memory that was freed but that the allocator keeps for reuse back
  /// to the system, and returns how many resident bytes that released (0
  /// where the system does not say). glibc keeps freed pages in the arena
  /// that allocated them (malloc_trim releases them); macOS and Windows
  /// return large blocks on free, and this does nothing there.
  std::uint64_t releaseFreedMemory();

  /// @p bytes as "12.3 GB", or "6.3 MB" below 1 GB (decimal units).
  std::string memorySize(std::uint64_t bytes);

  /// What auto asks of the free memory: this many times the estimated size of
  /// the MS1 map in memory (owner decision D7: room for the map and as much
  /// again for everything else that grows while it is held).
  inline constexpr double ms1_in_memory_headroom = 2.0;

  struct Ms1MemoryDecision
  {
    bool hold = false;
    /// Why, in words, for the log.
    std::string reason;
  };

  /// The rule: on and off are obeyed; auto holds the map when @p available is
  /// known (> 0) and at least ms1_in_memory_headroom times @p estimate. With
  /// no MS1 spectra (@p estimate 0) there is nothing to hold.
  Ms1MemoryDecision decideMs1InMemory(Ms1InMemory mode, std::uint64_t estimate, std::uint64_t available);
}
