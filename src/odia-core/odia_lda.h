// Copyright (c) 2026, Oliver Kohlbacher and the OpenDIAlyzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
//
// odia_lda.h -- semi-supervised LDA scorer (mProphet / pyProphet lineage).
//
// Vendored from OpenDIAlyzer src/odia_lda.h, LDA only: the gradient-boosted-tree and network
// learners, bagging, anchor training and their OpenMP nesting are not part of this copy. Origin
// commit, file hashes and every change against the original are recorded in
// src/odia-core/MANIFEST.json. Dependency-free C++17: no OpenMS, no Eigen.
//
// ALGORITHM
//   Inputs: features (n rows x m sub-scores, row-major), labels (1 = target, 0 = decoy) and
//   group (precursor id; a precursor has several candidate peak groups = rows). FDR is a
//   per-precursor quantity, so it is computed on the best-scoring row per group.
//   1. z-standardise each feature column over all rows; a missing (non-finite) cell becomes 0,
//      i.e. it is imputed at the column mean.
//   2. k-fold cross-validation BY PAIR: fold = hash(pair id, seed) % k for every group of the
//      pair, so a row is never scored by a model trained on its own precursor or on its own
//      target-decoy partner, and the assignment never looks at a label or at row order.
//   3. Per fold, train on the other folds:
//        - seed: the single most target/decoy-separating feature (max |Welch t|) ranks each
//          group; an LDA of the best target rows against the best decoy rows then gives a real
//          multivariate direction to start from (a single feature is often too weak to ignite
//          the loop -- a flat q for every group and zero trained iterations);
//        - n_iter times: positives = best row of each target group at q <= train_fdr (0.15 on
//          the first iteration, then 0.05); negatives = best row of each decoy group (top
//          decoys only, as pyProphet, so the classes do not also differ in peak rank);
//          Fisher LDA w = Sw^-1 (mu_pos - mu_neg), Sw the pooled within-class covariance,
//          solved by Cholesky with a ridge on the diagonal.
//      Apply the fold's final w to its held-out rows, then rescale the fold's scores so that
//      folds are commensurable before they are pooled: median and MAD of the best score of
//      every held-out group, targets and decoys alike, so the transform is label-blind.
//   4. q-values on the best d-score per group (odia_fdr.h), broadcast to the group's rows.
//   Deterministic given `seed`; folds are trained in parallel and never share state, so the
//   result does not depend on the thread count.

#ifndef ODIA_CORE_LDA_H
#define ODIA_CORE_LDA_H

#include "odia_fdr.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace odia::core
{

struct LdaParams
{
  int n_folds = 3;                  ///< cross-validation folds (by group); clamped to [2, groups]
  int n_iter = 3;                   ///< semi-supervised iterations per fold
  double train_fdr_initial = 0.15;  ///< FDR for the FIRST training-set selection: the loop starts
                                    ///< from the seed, so a strict cut can select too few positives
                                    ///< to fit anything (pyProphet uses 0.15 for the same reason)
  double train_fdr = 0.05;          ///< FDR for later iterations, once a real discriminant exists
  double ridge = 1e-6;              ///< diagonal regularisation of the within-class covariance
  unsigned seed = 42;               ///< fold assignment only (foldOfPair)
  bool use_pi0 = false;             ///< Storey pi0. false = honest/conservative; true inflates IDs
  bool top_decoys_only = true;      ///< negatives = each decoy group's BEST row, not all rows
  bool normalize_folds = true;      ///< rescale each fold's held-out scores (median/MAD, label-
                                    ///< blind) before pooling
  int threads = 0;                  ///< fold-level threads; 0 = OpenMP default. Result-neutral.
};

struct LdaResult
{
  std::vector<double> dscore;   ///< per input row (peak group)
  std::vector<double> qvalue;   ///< per input row: its group's q-value (pooled estimator)
  std::vector<double> pvalue;   ///< per input row: tail probability under the decoy null
  std::vector<double> pep;      ///< per input row: local FDR (posterior error probability)
  std::vector<int> fold;        ///< per input row: the cross-validation fold that scored it
  int n_folds = 0;              ///< folds actually used (n_folds clamped to [2, groups])
  /// Per fold: the location and scale its held-out scores were rescaled by (0 and 1 when the
  /// fold was left unscaled), and how many folds were left unscaled for want of a spread.
  std::vector<double> fold_location;
  std::vector<double> fold_scale;
  int folds_unscaled = 0;
  /// Semi-supervised iterations that FITTED a discriminant, and those that did not (too few
  /// confident positives, a failed solve, or a fold with no training groups; a failed seed fit
  /// also counts as skipped). If trained == 0 the scores come from the single-feature seed, not
  /// from a learned model.
  int n_iterations_trained = 0;
  int n_iterations_skipped = 0;
};

/// The fold of a target-decoy pair: a hash of the pair id and the seed, reduced modulo @p folds.
///
/// A function of the pair alone, so both members of a pair always land in one fold, the
/// assignment cannot see a label, and a pair keeps its fold when other pairs join or leave the
/// candidate set (chunking, subsetting). Fold sizes are balanced in expectation, not exactly.
/// SplitMix64 finaliser: the same value on every platform, unlike std::shuffle.
inline int foldOfPair(std::int64_t pair, unsigned seed, int folds)
{
  const auto mix = [](std::uint64_t v) {
    v += 0x9e3779b97f4a7c15ULL;
    v = (v ^ (v >> 30)) * 0xbf58476d1ce4e5b9ULL;
    v = (v ^ (v >> 27)) * 0x94d049bb133111ebULL;
    return v ^ (v >> 31);
  };
  const std::uint64_t h = mix(static_cast<std::uint64_t>(pair) ^ mix(static_cast<std::uint64_t>(seed)));
  return static_cast<int>(h % static_cast<std::uint64_t>(folds < 1 ? 1 : folds));
}

/// Label-blind location and scale of a set of scores: the median and 1.4826 x the median
/// absolute deviation (the normal-consistent MAD). `valid` is false for fewer than two values or
/// no spread, in which case nothing should be rescaled.
struct RobustScale
{
  double location = 0.0;
  double scale = 1.0;
  bool valid = false;
};

inline RobustScale robustScale(std::vector<double> values)
{
  RobustScale out;
  const std::size_t n = values.size();
  if (n < 2) { return out; }
  const auto median = [n](std::vector<double>& v) {
    const std::size_t mid = n / 2;
    std::nth_element(v.begin(), v.begin() + static_cast<std::ptrdiff_t>(mid), v.end());
    const double upper = v[mid];
    if (n % 2 == 1) { return upper; }
    const double lower = *std::max_element(v.begin(), v.begin() + static_cast<std::ptrdiff_t>(mid));
    return 0.5 * (lower + upper);
  };
  const double location = median(values);
  for (double& v : values) { v = std::abs(v - location); }
  const double scale = 1.4826 * median(values);
  if (!std::isfinite(location) || !std::isfinite(scale) ||
      !(scale > std::numeric_limits<double>::epsilon()))
  {
    return out;
  }
  out.location = location;
  out.scale = scale;
  out.valid = true;
  return out;
}

namespace detail
{

inline double dot(const double* a, const double* b, std::size_t m)
{
  double result = 0.0;
  for (std::size_t j = 0; j < m; ++j) { result += a[j] * b[j]; }
  return result;
}

// Solve A x = b for symmetric positive-definite A (row-major m x m). The copy of A is replaced
// by its lower-triangular Cholesky factor.
inline bool choleskySolve(std::vector<double> a, const std::vector<double>& b,
                          std::vector<double>& x)
{
  const std::size_t m = b.size();
  for (std::size_t i = 0; i < m; ++i)
  {
    for (std::size_t j = 0; j <= i; ++j)
    {
      double value = a[i * m + j];
      for (std::size_t k = 0; k < j; ++k) { value -= a[i * m + k] * a[j * m + k]; }
      if (i == j)
      {
        if (!(value > 0.0) || !std::isfinite(value)) { return false; }
        a[i * m + i] = std::sqrt(value);
      }
      else
      {
        a[i * m + j] = value / a[j * m + j];
      }
    }
  }

  std::vector<double> y(m, 0.0);
  for (std::size_t i = 0; i < m; ++i)
  {
    double value = b[i];
    for (std::size_t j = 0; j < i; ++j) { value -= a[i * m + j] * y[j]; }
    y[i] = value / a[i * m + i];
  }

  x.assign(m, 0.0);
  for (std::size_t ii = m; ii > 0; --ii)
  {
    const std::size_t i = ii - 1;
    double value = y[i];
    for (std::size_t j = i + 1; j < m; ++j) { value -= a[j * m + i] * x[j]; }
    x[i] = value / a[i * m + i];
    if (!std::isfinite(x[i])) { return false; }
  }
  return true;
}

} // namespace detail

/// Semi-supervised LDA with cross-validation and target-decoy q-values.
///
/// @p x is row-major, labels.size() rows by @p m columns; non-finite cells are missing. It is
/// taken by value and standardised in place: std::move a large matrix in.
/// labels[i] is 1 (target) or 0 (decoy); every row of a group must carry the same label.
/// group[i] is the precursor id shared by that precursor's candidate peak groups.
/// pair[i] is the id a target precursor shares with its decoy, supplied by the caller; every row
/// of a group must carry the same pair id. A precursor without a partner has a pair id of its own.
/// Throws std::invalid_argument on inconsistent sizes, or a group with mixed labels or pair ids.
inline LdaResult scoreSemiSupervisedLDA(std::vector<double> x, std::size_t m,
                                        const std::vector<int>& labels,
                                        const std::vector<std::int64_t>& group,
                                        const std::vector<std::int64_t>& pair,
                                        const LdaParams& params = LdaParams())
{
  const std::size_t n = labels.size();
  if (group.size() != n || pair.size() != n || x.size() != n * m)
  {
    throw std::invalid_argument("odia::core::scoreSemiSupervisedLDA: " + std::to_string(n) +
                                " labels, " + std::to_string(group.size()) + " group ids, " +
                                std::to_string(pair.size()) + " pair ids and " +
                                std::to_string(x.size()) + " cells for " + std::to_string(m) +
                                " columns do not describe one table");
  }
  LdaResult result;
  result.dscore.assign(n, 0.0);
  result.qvalue.assign(n, 1.0);
  result.pvalue.assign(n, 1.0);
  result.pep.assign(n, 1.0);
  result.fold.assign(n, 0);
  if (n == 0 || m == 0) { return result; }

  // Global z-standardisation is unsupervised. Missing cells become 0 = the column mean.
  std::vector<double> mean(m, 0.0);
  std::vector<double> count(m, 0.0);
  for (std::size_t i = 0; i < n; ++i)
  {
    for (std::size_t j = 0; j < m; ++j)
    {
      const double v = x[i * m + j];
      if (std::isfinite(v)) { mean[j] += v; count[j] += 1.0; }
    }
  }
  for (std::size_t j = 0; j < m; ++j)
  {
    if (count[j] > 0.0) { mean[j] /= count[j]; }
  }
  std::vector<double> sum_squared_deviation(m, 0.0);
  for (std::size_t i = 0; i < n; ++i)
  {
    for (std::size_t j = 0; j < m; ++j)
    {
      const double v = x[i * m + j];
      if (std::isfinite(v)) { const double d = v - mean[j]; sum_squared_deviation[j] += d * d; }
    }
  }
  std::vector<double> sd(m, 1.0);
  for (std::size_t j = 0; j < m; ++j)
  {
    sd[j] = count[j] > 1.0 ? std::sqrt(sum_squared_deviation[j] / (count[j] - 1.0)) : 0.0;
    if (!(sd[j] > std::numeric_limits<double>::epsilon()) || !std::isfinite(sd[j])) { sd[j] = 1.0; }
  }
  // Standardised in place: x is taken by value, so a caller that moves its matrix in pays for
  // one copy of it, not two.
  std::vector<double>& z = x;
  for (std::size_t i = 0; i < n; ++i)
  {
    for (std::size_t j = 0; j < m; ++j)
    {
      double& v = z[i * m + j];
      v = std::isfinite(v) ? (v - mean[j]) / sd[j] : 0.0;
    }
  }
  const auto zrow = [&](std::size_t row) { return z.data() + row * m; };

  // Groups in first-occurrence order, rows of each group in row order (CSR).
  std::unordered_map<std::int64_t, std::size_t> group_lookup;
  group_lookup.reserve(n);
  std::vector<std::size_t> group_of_row(n);
  std::vector<int> group_label;
  std::vector<std::int64_t> group_id;
  std::vector<std::int64_t> group_pair;
  for (std::size_t i = 0; i < n; ++i)
  {
    const auto inserted = group_lookup.emplace(group[i], group_id.size());
    const int label = labels[i] == 1 ? 1 : 0;
    if (inserted.second)
    {
      group_label.push_back(label);
      group_id.push_back(group[i]);
      group_pair.push_back(pair[i]);
    }
    else if (group_label[inserted.first->second] != label)
    {
      throw std::invalid_argument("odia::core::scoreSemiSupervisedLDA: group " +
                                  std::to_string(group[i]) + " mixes target and decoy rows");
    }
    else if (group_pair[inserted.first->second] != pair[i])
    {
      throw std::invalid_argument("odia::core::scoreSemiSupervisedLDA: group " +
                                  std::to_string(group[i]) + " has rows in pairs " +
                                  std::to_string(group_pair[inserted.first->second]) + " and " +
                                  std::to_string(pair[i]));
    }
    group_of_row[i] = inserted.first->second;
  }
  const std::size_t group_count = group_id.size();
  std::vector<std::size_t> first(group_count + 1, 0);
  for (std::size_t i = 0; i < n; ++i) { ++first[group_of_row[i] + 1]; }
  for (std::size_t g = 0; g < group_count; ++g) { first[g + 1] += first[g]; }
  std::vector<std::size_t> members(n);
  {
    std::vector<std::size_t> fill(first.begin(), first.end() - 1);
    for (std::size_t i = 0; i < n; ++i) { members[fill[group_of_row[i]]++] = i; }
  }

  if (group_count < 2) { return result; } // leakage-free training is impossible

  int folds = params.n_folds;
  if (folds < 2) { folds = 2; }
  if (folds > static_cast<int>(group_count)) { folds = static_cast<int>(group_count); }

  // Fold assignment by PAIR. ODIA shuffled target groups and decoy groups independently, so a
  // target and its own decoy usually fell into different folds: the partner of a held-out target
  // then trained the model that scored it. It also used std::shuffle, whose output differs
  // between standard libraries. A hash of the pair id has neither problem.
  std::vector<int> group_fold(group_count, 0);
  for (std::size_t g = 0; g < group_count; ++g)
  {
    group_fold[g] = foldOfPair(group_pair[g], params.seed, folds);
  }
  for (std::size_t i = 0; i < n; ++i) { result.fold[i] = group_fold[group_of_row[i]]; }
  result.n_folds = folds;

  // Folds are independent BY CONSTRUCTION: fold f trains on the groups not assigned to f and
  // writes result.dscore only for the rows of groups assigned to f. Nothing else is shared except
  // the two counters, which are reduced, so the result does not depend on completion order or on
  // the number of threads. Inside a fold everything runs serially in a fixed order.
  result.fold_location.assign(static_cast<std::size_t>(folds), 0.0);
  result.fold_scale.assign(static_cast<std::size_t>(folds), 1.0);
  int n_trained = 0, n_skipped = 0, n_unscaled = 0;
#ifdef _OPENMP
  const int team = params.threads > 0 ? params.threads : omp_get_max_threads();
#pragma omp parallel for schedule(dynamic, 1) num_threads(team) reduction(+ : n_trained, n_skipped, n_unscaled)
#endif
  for (int fold = 0; fold < folds; ++fold)
  {
    std::vector<std::size_t> train_groups;
    train_groups.reserve(group_count);
    for (std::size_t g = 0; g < group_count; ++g)
    {
      if (group_fold[g] != fold) { train_groups.push_back(g); }
    }
    if (train_groups.empty())
    {
      n_skipped += 1 + std::max(0, params.n_iter);   // held-out rows keep d-score 0
      if (params.normalize_folds) { ++n_unscaled; }
      continue;
    }

    std::vector<double> w(m, 0.0);
    const auto score_row = [&](std::size_t row) { return detail::dot(w.data(), zrow(row), m); };

    // Best-scoring row of each group; the FIRST maximum in the group's row order wins.
    const auto best_rows_of = [&](const std::vector<std::size_t>& groups,
                                  std::vector<std::size_t>& out, std::vector<double>& out_score) {
      out.assign(groups.size(), 0);
      out_score.assign(groups.size(), 0.0);
      for (std::size_t i = 0; i < groups.size(); ++i)
      {
        const std::size_t g = groups[i];
        std::size_t best_row = members[first[g]];
        double best = score_row(best_row);
        for (std::size_t k = first[g] + 1; k < first[g + 1]; ++k)
        {
          const double sc = score_row(members[k]);
          if (sc > best) { best = sc; best_row = members[k]; }
        }
        out[i] = best_row;
        out_score[i] = best;
      }
    };

    // w = Sw^-1 (mu_pos - mu_neg): Fisher's discriminant on two explicit row sets.
    const auto fit_lda = [&](const std::vector<std::size_t>& positive_rows,
                             const std::vector<std::size_t>& negative_rows,
                             std::vector<double>& out) -> bool {
      if (positive_rows.size() < 2 || negative_rows.size() < 2) { return false; }
      std::vector<double> positive_mean(m, 0.0);
      std::vector<double> negative_mean(m, 0.0);
      for (const std::size_t row : positive_rows)
      {
        const double* zr = zrow(row);
        for (std::size_t j = 0; j < m; ++j) { positive_mean[j] += zr[j]; }
      }
      for (const std::size_t row : negative_rows)
      {
        const double* zr = zrow(row);
        for (std::size_t j = 0; j < m; ++j) { negative_mean[j] += zr[j]; }
      }
      for (std::size_t j = 0; j < m; ++j)
      {
        positive_mean[j] /= static_cast<double>(positive_rows.size());
        negative_mean[j] /= static_cast<double>(negative_rows.size());
      }

      std::vector<double> covariance(m * m, 0.0);
      const auto add_within_class = [&](const std::vector<std::size_t>& rows,
                                        const std::vector<double>& class_mean) {
        for (const std::size_t row : rows)
        {
          const double* zr = zrow(row);
          for (std::size_t j = 0; j < m; ++j)
          {
            const double dj = zr[j] - class_mean[j];
            for (std::size_t k = 0; k <= j; ++k)
            {
              covariance[j * m + k] += dj * (zr[k] - class_mean[k]);
            }
          }
        }
      };
      add_within_class(positive_rows, positive_mean);
      add_within_class(negative_rows, negative_mean);
      const double dof = static_cast<double>(positive_rows.size() + negative_rows.size() - 2);
      for (std::size_t j = 0; j < m; ++j)
      {
        for (std::size_t k = 0; k <= j; ++k)
        {
          covariance[j * m + k] /= dof;
          covariance[k * m + j] = covariance[j * m + k];
        }
      }

      std::vector<double> difference(m);
      for (std::size_t j = 0; j < m; ++j) { difference[j] = positive_mean[j] - negative_mean[j]; }

      // A positive ridge is tried exactly; a zero/negative one starts at a tiny numerical ridge.
      // The ridge grows geometrically on failure.
      double ridge = params.ridge > 0.0 && std::isfinite(params.ridge) ? params.ridge : 1e-12;
      bool solved = false;
      for (int attempt = 0; attempt < 10 && !solved; ++attempt)
      {
        std::vector<double> regularised = covariance;
        for (std::size_t j = 0; j < m; ++j) { regularised[j * m + j] += ridge; }
        solved = detail::choleskySolve(regularised, difference, out);
        ridge *= 10.0;
      }
      double norm_squared = 0.0;
      for (const double value : out) { norm_squared += value * value; }
      return solved && norm_squared > std::numeric_limits<double>::epsilon();
    };

    // Single-feature bootstrap: the signed feature with the largest Welch |t| between all target
    // and all decoy training rows.
    std::size_t best_feature = 0;
    double best_abs_t = -1.0;
    double best_difference = 1.0;
    for (std::size_t j = 0; j < m; ++j)
    {
      double sum[2] = {0.0, 0.0};
      double sum_sq[2] = {0.0, 0.0};
      std::size_t class_n[2] = {0, 0};
      for (const std::size_t g : train_groups)
      {
        const int cls = group_label[g];
        for (std::size_t k = first[g]; k < first[g + 1]; ++k)
        {
          const double v = zrow(members[k])[j];
          sum[cls] += v;
          sum_sq[cls] += v * v;
          ++class_n[cls];
        }
      }
      if (class_n[0] == 0 || class_n[1] == 0) { continue; }
      const double class_mean0 = sum[0] / static_cast<double>(class_n[0]);
      const double class_mean1 = sum[1] / static_cast<double>(class_n[1]);
      const double var0 = class_n[0] > 1
                            ? std::max(0.0, (sum_sq[0] - sum[0] * class_mean0) /
                                              static_cast<double>(class_n[0] - 1))
                            : 0.0;
      const double var1 = class_n[1] > 1
                            ? std::max(0.0, (sum_sq[1] - sum[1] * class_mean1) /
                                              static_cast<double>(class_n[1] - 1))
                            : 0.0;
      const double difference = class_mean1 - class_mean0;
      const double standard_error = std::sqrt(var0 / static_cast<double>(class_n[0]) +
                                              var1 / static_cast<double>(class_n[1]));
      const double abs_t = standard_error > 0.0
                             ? std::abs(difference) / standard_error
                             : (difference == 0.0 ? 0.0 : std::numeric_limits<double>::infinity());
      if (abs_t > best_abs_t)
      {
        best_abs_t = abs_t;
        best_feature = j;
        best_difference = difference;
      }
    }
    w[best_feature] = best_difference < 0.0 ? -1.0 : 1.0;

    std::vector<std::size_t> bg_row;
    std::vector<double> bg_score;

    // Seed: an LDA of all best target rows against all best decoy rows under the single-feature
    // ranking. The target class is contaminated (most target groups are false), which costs the
    // direction magnitude rather than orientation; it only has to ignite the loop.
    {
      std::vector<std::size_t> seed_pos, seed_neg;
      best_rows_of(train_groups, bg_row, bg_score);
      for (std::size_t i = 0; i < train_groups.size(); ++i)
      {
        (group_label[train_groups[i]] == 1 ? seed_pos : seed_neg).push_back(bg_row[i]);
      }
      std::vector<double> next_w;
      if (fit_lda(seed_pos, seed_neg, next_w)) { w.swap(next_w); }
      else { ++n_skipped; }   // a failed seed is recorded, never silent
    }

    for (int iteration = 0; iteration < std::max(0, params.n_iter); ++iteration)
    {
      std::vector<RankedGroup> ranked;
      ranked.reserve(train_groups.size());
      best_rows_of(train_groups, bg_row, bg_score);
      for (std::size_t i = 0; i < train_groups.size(); ++i)
      {
        RankedGroup r;
        r.group_index = train_groups[i];
        r.best_row = bg_row[i];
        r.label = group_label[train_groups[i]];
        r.score = bg_score[i];
        ranked.push_back(r);
      }
      assignQValues(ranked, params.use_pi0);

      const double raw_fdr = (iteration == 0) ? params.train_fdr_initial : params.train_fdr;
      const double train_fdr = std::max(0.0, std::min(1.0, raw_fdr));
      std::vector<std::size_t> positive_rows;
      for (const auto& candidate : ranked)
      {
        if (candidate.label == 1 && candidate.qvalue <= train_fdr) { positive_rows.push_back(candidate.best_row); }
      }
      // Too few confident positives to fit an m-dimensional discriminant: skip, and count it.
      if (positive_rows.size() < m + 2)
      {
        ++n_skipped;
        continue;
      }

      std::vector<std::size_t> negative_rows;
      if (!params.top_decoys_only)
      {
        for (const std::size_t g : train_groups)
        {
          if (group_label[g] == 1) { continue; }
          for (std::size_t k = first[g]; k < first[g + 1]; ++k) { negative_rows.push_back(members[k]); }
        }
      }
      else
      {
        // The ranking scan above already found every training group's best row with this model.
        for (std::size_t i = 0; i < train_groups.size(); ++i)
        {
          if (group_label[train_groups[i]] == 1) { continue; }
          negative_rows.push_back(bg_row[i]);
        }
      }
      if (negative_rows.size() < 2) { continue; }

      // A failed fit keeps the previous model and is counted as skipped.
      std::vector<double> next_w;
      if (fit_lda(positive_rows, negative_rows, next_w)) { w.swap(next_w); ++n_trained; }
      else { ++n_skipped; }
    }

    // Held-out scoring: this model has seen no row of the groups it scores.
    for (std::size_t g = 0; g < group_count; ++g)
    {
      if (group_fold[g] != fold) { continue; }
      for (std::size_t k = first[g]; k < first[g + 1]; ++k)
      {
        result.dscore[members[k]] = score_row(members[k]);
      }
    }

    // Every fold has its own weight vector, defined only up to scale and offset, so raw scores are
    // not comparable across folds. Rescale this fold's held-out scores by the median and MAD of
    // the best-per-group score over ALL of its groups, targets and decoys alike.
    //
    // ODIA standardised each fold to the mean and sd of its DECOYS' best scores. Those statistics
    // are computed from the labels of the very rows the null is then built from: with few decoys
    // in a fold, two decoys are forced to +/-1/sqrt(2) while an exchangeable null target is not,
    // so target and decoy stop being exchangeable -- the assumption target-decoy FDR rests on.
    // This transform never looks at a label. Both members of a pair share the fold (above), so
    // the rescale is also the same monotone map for both and cannot change which member wins.
    if (params.normalize_folds)
    {
      std::vector<double> best_scores;
      for (std::size_t g = 0; g < group_count; ++g)
      {
        if (group_fold[g] != fold) { continue; }
        double best = result.dscore[members[first[g]]];
        for (std::size_t k = first[g] + 1; k < first[g + 1]; ++k)
        {
          best = std::max(best, result.dscore[members[k]]);
        }
        best_scores.push_back(best);
      }
      const RobustScale rs = robustScale(std::move(best_scores));
      if (rs.valid)
      {
        for (std::size_t g = 0; g < group_count; ++g)
        {
          if (group_fold[g] != fold) { continue; }
          for (std::size_t k = first[g]; k < first[g + 1]; ++k)
          {
            result.dscore[members[k]] = (result.dscore[members[k]] - rs.location) / rs.scale;
          }
        }
        result.fold_location[static_cast<std::size_t>(fold)] = rs.location;
        result.fold_scale[static_cast<std::size_t>(fold)] = rs.scale;
      }
      else
      {
        ++n_unscaled;   // no spread carries no scale information; the scores stay as they are
      }
    }
  }

  result.n_iterations_trained = n_trained;
  result.n_iterations_skipped = n_skipped;
  result.folds_unscaled = n_unscaled;

  std::vector<RankedGroup> final_ranked;
  final_ranked.reserve(group_count);
  for (std::size_t g = 0; g < group_count; ++g)
  {
    std::size_t best_row = members[first[g]];
    for (std::size_t k = first[g] + 1; k < first[g + 1]; ++k)
    {
      if (result.dscore[members[k]] > result.dscore[best_row]) { best_row = members[k]; }
    }
    RankedGroup r;
    r.group_index = g;
    r.best_row = best_row;
    r.label = group_label[g];
    r.score = result.dscore[best_row];
    final_ranked.push_back(r);
  }
  assignQValues(final_ranked, params.use_pi0);
  for (const auto& ranked_group : final_ranked)
  {
    const std::size_t g = ranked_group.group_index;
    for (std::size_t k = first[g]; k < first[g + 1]; ++k)
    {
      result.qvalue[members[k]] = ranked_group.qvalue;
      result.pvalue[members[k]] = ranked_group.pvalue;
      result.pep[members[k]] = ranked_group.pep;
    }
  }
  return result;
}

} // namespace odia::core

#endif // ODIA_CORE_LDA_H
