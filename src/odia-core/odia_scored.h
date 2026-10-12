// Copyright (c) 2026, Oliver Kohlbacher and the OpenDIAlyzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
//
// odia_scored.h -- the score rows a classifier consumes, one per candidate peak group.
//
// Adapted from OpenDIAlyzer src/odia_scored.h (class ScoreRows): the sub-scores of every row
// live in ONE row-major block whose width is fixed before the first row, so a row is a slice,
// not an allocation. What changed is recorded in src/odia-core/MANIFEST.json. Dependency-free
// C++17: no OpenMS, no Eigen.
//
// Why a flat block: an OpenSWATH FeatureMap holds each feature's scores as string-keyed meta
// values behind per-feature heap allocations, many times the size of the numbers themselves.
// The classifier needs the numbers and four ids per row, nothing else, so a pipeline flattens
// each chunk's features into this table and frees them.

#ifndef ODIA_CORE_SCORED_H
#define ODIA_CORE_SCORED_H

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace odia::core
{

/// Candidate peak groups to classify. Public members, so a caller may also fill the vectors
/// directly; validate() checks that they describe one table.
struct ScoreTable
{
  /// Sub-score names in column order (OpenSWATH meta value names such as var_library_corr).
  /// Fixed before the first row.
  std::vector<std::string> feature_names;
  /// rows() x width() sub-scores, row-major. NaN (or any non-finite value) means missing.
  std::vector<float> values;
  /// Per row: the precursor this candidate peak group belongs to.
  std::vector<std::int64_t> group;
  /// Per row: the pair id a target precursor shares with its decoy. Every row of a precursor
  /// carries the same id; a precursor without a partner uses an id of its own (or a negative one).
  std::vector<std::int64_t> pair;
  /// Per row: 1 = decoy, 0 = target. Every row of a precursor carries the same value.
  std::vector<std::uint8_t> is_decoy;
  /// Per row: an id unique within the precursor (e.g. OpenSWATH's feature id). It fixes the
  /// canonical row order, so it must not depend on thread scheduling.
  std::vector<std::int64_t> feature_id;

  void setColumns(std::vector<std::string> names)
  {
    if (!group.empty()) { throw std::logic_error("odia::core::ScoreTable: columns fixed after the first row"); }
    feature_names = std::move(names);
  }
  std::size_t width() const { return feature_names.size(); }
  std::size_t rows() const { return group.size(); }

  void reserve(std::size_t n)
  {
    values.reserve(n * width());
    group.reserve(n);
    pair.reserve(n);
    is_decoy.reserve(n);
    feature_id.reserve(n);
  }

  /// One candidate peak group; @p row holds width() values in column order. Returns its index.
  template <typename T>
  std::size_t append(std::int64_t group_id, std::int64_t pair_id, bool decoy, std::int64_t feature,
                     const T* row)
  {
    for (std::size_t j = 0; j < width(); ++j) { values.push_back(static_cast<float>(row[j])); }
    group.push_back(group_id);
    pair.push_back(pair_id);
    is_decoy.push_back(decoy ? 1 : 0);
    feature_id.push_back(feature);
    return group.size() - 1;
  }

  const float* row(std::size_t i) const { return values.data() + i * width(); }

  /// Throws std::invalid_argument, with the counts, unless the vectors describe one table.
  void validate() const
  {
    const std::size_t n = rows();
    if (pair.size() != n || is_decoy.size() != n || feature_id.size() != n ||
        values.size() != n * width())
    {
      throw std::invalid_argument(
        "odia::core::ScoreTable: " + std::to_string(n) + " group ids, " +
        std::to_string(pair.size()) + " pair ids, " + std::to_string(is_decoy.size()) +
        " decoy flags, " + std::to_string(feature_id.size()) + " feature ids and " +
        std::to_string(values.size()) + " values for " + std::to_string(width()) + " columns");
    }
  }

  std::size_t bytes() const
  {
    std::size_t s = values.capacity() * sizeof(float) +
                    (group.capacity() + pair.capacity() + feature_id.capacity()) * sizeof(std::int64_t) +
                    is_decoy.capacity();
    for (const auto& name : feature_names) { s += sizeof(std::string) + name.capacity(); }
    return s;
  }

  /// Release everything, capacity included.
  void clear()
  {
    std::vector<std::string>().swap(feature_names);
    std::vector<float>().swap(values);
    std::vector<std::int64_t>().swap(group);
    std::vector<std::int64_t>().swap(pair);
    std::vector<std::uint8_t>().swap(is_decoy);
    std::vector<std::int64_t>().swap(feature_id);
  }
};

} // namespace odia::core

#endif // ODIA_CORE_SCORED_H
