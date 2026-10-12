/*

This file is part of the mzpeak.h project.  It is subject to the
license specified in the LICENSE file which can be found in the
top-level directory of this repository.

*/

#pragma once

#include <cstdint>
#include <memory>
#include <optional>
#include <utility>
#include <vector>

#include "mzpeak/util/planner.h"
#include "mzpeak/util/projection.h"
#include "mzpeak/util/slice.h"

namespace parquet::arrow {
class FileReader;
}

namespace MzPeak::Util {

class Parquet;

/**
 * Query executor.
 */
class Executor final {
public:
  /// Clean up signatures.
  using Slice = MzPeak::Util::Slice;

  /// Destructor.
  ~Executor();

  /**
   * Execute a query and extract the requested projection.
   */
  std::unique_ptr<Slice> execute(const Planner::Plan&);

  /**
   * The rows of every key in [@p first, @p last], one slice per key in
   * ascending order, taken from the planned row groups' KeyRuns: no filter,
   * no key column.  nullopt when a planned group holds no runs (its key is not
   * declared sorted, or the projection includes the key) or the keys do not
   * ascend across groups; the caller then filters instead.
   */
  std::optional<std::vector<std::pair<uint64_t, std::unique_ptr<Slice>>>>
  execute_runs(const Planner::Plan&, uint64_t first, uint64_t last);

private:
  friend class Parquet;

  /// Constructor.  `source` supplies decoded row groups and outlives this.
  Executor(Parquet& source, const Projection&);

  struct Impl;
  std::unique_ptr<Impl> impl_;
};

} // namespace MzPeak::Util
