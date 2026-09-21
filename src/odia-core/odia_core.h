// Copyright (c) 2026, Oliver Kohlbacher, the OpenDIAlyzer contributors and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause
//
// odia_core.h -- score candidate peak groups and control the FDR, in one call.
//
// The entry point a search pipeline uses. It wraps the vendored LDA (odia_lda.h) and FDR
// (odia_fdr.h) code with the preparation OpenDIAlyzer did in src/opendialyzer.cpp before
// scoring (OswRows::canonicalize, dropUninformativeColumns_); what changed is recorded in
// src/odia-core/MANIFEST.json. Dependency-free C++17: no OpenMS, no Eigen.
//
//   ScoreTable t;                                   // odia_scored.h
//   t.setColumns({"var_library_corr", ...});
//   t.append(precursor, pair, is_decoy, feature_id, scores);   // per candidate peak group
//   odia::core::ScoredResult r = odia::core::scoreAndControl(t, odia::core::Options());
//   // r.groups[g].qvalue: precursor q after concatenated target-decoy competition
//
// Every step is label-blind except training and the q-values themselves: the column filter,
// the imputation, the fold assignment and the fold rescaling never read a label. The result is
// a function of the table's CONTENT: rows are put into a canonical order (precursor, then
// feature id) first, so neither the order the rows arrive in nor the thread count changes it.

#ifndef ODIA_CORE_CORE_H
#define ODIA_CORE_CORE_H

#include "odia_fdr.h"
#include "odia_lda.h"
#include "odia_scored.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace odia::core
{

struct Options
{
  /// The classifier: 3 folds by pair, 3 iterations, train FDR 0.15 then 0.05, ridge 1e-6,
  /// pi0 off, top decoys only, label-blind fold rescaling. See LdaParams.
  LdaParams lda;
  /// Precursor-level q over the pair winners; Ratio is (D_win+1)/N_dec_win / (T_win/N_tar_win).
  QEstimator estimator = QEstimator::Ratio;
  /// The cut the diagnostic counts are taken at.
  double report_q = 0.01;
  /// Sub-scores the fit must not see at all, by exact column name. For a discriminant that
  /// selects retention-time or ion-mobility calibration anchors, the RT and 1/K0 deviation scores
  /// (e.g. var_norm_rt_score, var_im_delta_score): anchors picked by the very agreement they are
  /// meant to correct would bias the correction toward the uncorrected state. Names that match
  /// no column are reported in Diagnostics::exclusions_unmatched, never silently ignored.
  std::vector<std::string> exclude_features;
};

/// One precursor.
struct GroupResult
{
  std::int64_t group = 0;
  std::int64_t pair = 0;
  bool is_decoy = false;
  std::size_t rows = 0;          ///< its candidate peak groups
  std::size_t best_row = 0;      ///< input row index of the best-scoring one
  double score = 0.0;            ///< that row's d-score
  bool winner = false;           ///< won its pair, or had no competitor present
  double qvalue = 1.0;           ///< after concatenated competition; 1 for a loser
  double pep = 1.0;              ///< after concatenated competition; 1 for a loser
  double pooled_qvalue = 1.0;    ///< DIAGNOSTIC: ODIA's pooled estimator, no competition
};

struct Diagnostics
{
  std::size_t rows = 0;
  std::size_t groups = 0;
  std::size_t target_groups = 0;
  std::size_t decoy_groups = 0;
  std::size_t pairs_complete = 0;      ///< pairs with both precursors present
  std::size_t targets_unpaired = 0;    ///< target precursors whose decoy has no row
  std::size_t decoys_unpaired = 0;     ///< decoy precursors whose target has no row
  std::vector<std::string> features_used;
  std::vector<std::string> features_dropped;   ///< all missing, or constant (label-blind)
  std::vector<std::string> features_excluded;  ///< removed by Options::exclude_features
  std::vector<std::string> exclusions_unmatched;   ///< requested exclusions that named no column
  std::size_t cells_imputed = 0;       ///< missing cells of used columns, imputed at the column mean
  int n_folds = 0;
  int iterations_trained = 0;          ///< 0 means no discriminant was ever learned
  int iterations_skipped = 0;
  int folds_unscaled = 0;
  std::size_t target_winners = 0;
  std::size_t decoy_winners = 0;
  std::size_t targets_at_q = 0;        ///< winning targets at q <= report_q: the identifications
  std::size_t decoys_at_q = 0;         ///< winning decoys at q <= report_q
  std::size_t pooled_targets_at_q = 0; ///< targets at pooled q <= report_q
  /// (pooled_targets_at_q + 1) / (targets_at_q + 1): about 1 when pairs behave as exchangeable
  /// target-decoy pairs. Smoothed so it is always finite.
  double pooled_vs_paired = 1.0;
};

struct ScoredResult
{
  std::vector<double> dscore;          ///< per input row, in input order
  std::vector<GroupResult> groups;     ///< one per precursor, ascending group id
  Diagnostics diagnostics;
};

/// Rows in canonical order: by (group, feature_id). Throws std::invalid_argument if two rows
/// share both, since their order -- and therefore the result -- would then depend on input order.
inline std::vector<std::size_t> canonicalOrder(const ScoreTable& table)
{
  std::vector<std::size_t> order(table.rows());
  std::iota(order.begin(), order.end(), std::size_t{0});
  std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
    if (table.group[a] != table.group[b]) { return table.group[a] < table.group[b]; }
    if (table.feature_id[a] != table.feature_id[b]) { return table.feature_id[a] < table.feature_id[b]; }
    return a < b;
  });
  for (std::size_t k = 1; k < order.size(); ++k)
  {
    if (table.group[order[k]] == table.group[order[k - 1]] &&
        table.feature_id[order[k]] == table.feature_id[order[k - 1]])
    {
      throw std::invalid_argument("odia::core::canonicalOrder: rows " + std::to_string(order[k - 1]) +
                                  " and " + std::to_string(order[k]) + " share group " +
                                  std::to_string(table.group[order[k]]) + " and feature id " +
                                  std::to_string(table.feature_id[order[k]]));
    }
  }
  return order;
}

/// Columns that carry signal: at least one finite value, and not constant among the finite
/// values. Label-blind (it reads no label, so it cannot leak into the cross-validation); keeps
/// an all-missing column from turning every d-score into NaN. Columns flagged in @p skip are
/// not considered at all.
inline std::vector<std::size_t> informativeColumns(const ScoreTable& table,
                                                   const std::vector<char>& skip = {})
{
  std::vector<std::size_t> keep;
  const std::size_t m = table.width();
  for (std::size_t j = 0; j < m; ++j)
  {
    if (j < skip.size() && skip[j]) { continue; }
    bool have_first = false, varies = false;
    float first = 0.0f;
    for (std::size_t i = 0; i < table.rows() && !varies; ++i)
    {
      const float v = table.values[i * m + j];
      if (!std::isfinite(v)) { continue; }
      if (!have_first) { first = v; have_first = true; }
      else if (v != first) { varies = true; }
    }
    if (varies) { keep.push_back(j); }
  }
  return keep;
}

/// Score every candidate peak group with the semi-supervised LDA and control the FDR at the
/// precursor level by concatenated target-decoy competition within each pair.
///
/// Throws std::invalid_argument, with counts, when the table cannot be scored honestly: sizes
/// that do not describe one table, duplicate (group, feature id), a precursor whose rows mix
/// labels or pair ids, a pair with two targets or two decoys, no target or no decoy precursor,
/// or no informative column.
inline ScoredResult scoreAndControl(const ScoreTable& table, const Options& options = Options())
{
  table.validate();
  const std::size_t n = table.rows();
  const std::size_t width = table.width();
  if (n == 0) { throw std::invalid_argument("odia::core::scoreAndControl: the score table is empty"); }

  ScoredResult out;
  Diagnostics& diag = out.diagnostics;
  diag.rows = n;

  const std::vector<std::size_t> order = canonicalOrder(table);

  // Excluded columns are removed before anything else looks at them, so the fit -- seed,
  // iterations and held-out scoring alike -- never sees them.
  std::vector<char> excluded(width, 0);
  for (const std::string& name : options.exclude_features)
  {
    bool matched = false;
    for (std::size_t j = 0; j < width; ++j)
    {
      if (table.feature_names[j] == name) { excluded[j] = 1; matched = true; }
    }
    if (!matched) { diag.exclusions_unmatched.push_back(name); }
  }
  const std::vector<std::size_t> kept = informativeColumns(table, excluded);
  {
    std::vector<char> is_kept(width, 0);
    for (const std::size_t j : kept) { is_kept[j] = 1; diag.features_used.push_back(table.feature_names[j]); }
    for (std::size_t j = 0; j < width; ++j)
    {
      if (excluded[j]) { diag.features_excluded.push_back(table.feature_names[j]); }
      else if (!is_kept[j]) { diag.features_dropped.push_back(table.feature_names[j]); }
    }
  }
  if (kept.empty())
  {
    throw std::invalid_argument("odia::core::scoreAndControl: none of the " + std::to_string(width) +
                                " sub-score columns varies (" + std::to_string(diag.features_excluded.size()) +
                                " excluded); nothing to classify on");
  }

  // The classifier's inputs, in canonical order. Missing stays NaN here; the LDA imputes it at
  // the column mean when it standardises.
  const std::size_t m = kept.size();
  std::vector<double> x(n * m);
  std::vector<int> label(n);
  std::vector<std::int64_t> group(n), pair(n);
  for (std::size_t k = 0; k < n; ++k)
  {
    const std::size_t i = order[k];
    const float* row = table.row(i);
    for (std::size_t c = 0; c < m; ++c)
    {
      const float v = row[kept[c]];
      if (std::isfinite(v)) { x[k * m + c] = static_cast<double>(v); }
      else { x[k * m + c] = std::numeric_limits<double>::quiet_NaN(); ++diag.cells_imputed; }
    }
    label[k] = table.is_decoy[i] ? 0 : 1;
    group[k] = table.group[i];
    pair[k] = table.pair[i];
  }

  // Precursors: contiguous and ascending in canonical order.
  std::vector<std::size_t> group_first;
  for (std::size_t k = 0; k < n; ++k)
  {
    if (k == 0 || group[k] != group[k - 1]) { group_first.push_back(k); }
  }
  group_first.push_back(n);
  const std::size_t groups = group_first.size() - 1;
  for (std::size_t g = 0; g < groups; ++g)
  {
    (label[group_first[g]] == 1 ? diag.target_groups : diag.decoy_groups)++;
  }
  diag.groups = groups;
  if (diag.target_groups == 0 || diag.decoy_groups == 0)
  {
    throw std::invalid_argument("odia::core::scoreAndControl: " + std::to_string(diag.target_groups) +
                                " target and " + std::to_string(diag.decoy_groups) +
                                " decoy precursors; target-decoy FDR needs both");
  }
  {
    std::unordered_map<std::int64_t, std::int64_t> member[2];   // pair -> precursor, by label
    for (std::size_t g = 0; g < groups; ++g)
    {
      const std::size_t k = group_first[g];
      if (pair[k] < 0) { continue; }
      const auto ins = member[label[k]].emplace(pair[k], group[k]);
      if (!ins.second)
      {
        throw std::invalid_argument("odia::core::scoreAndControl: pair " + std::to_string(pair[k]) +
                                    " has two " + (label[k] == 1 ? "target" : "decoy") +
                                    " precursors, " + std::to_string(ins.first->second) + " and " +
                                    std::to_string(group[k]));
      }
    }
  }

  // Throws on a precursor with mixed labels or pair ids.
  const LdaResult lda = scoreSemiSupervisedLDA(std::move(x), m, label, group, pair, options.lda);
  diag.n_folds = lda.n_folds;
  diag.iterations_trained = lda.n_iterations_trained;
  diag.iterations_skipped = lda.n_iterations_skipped;
  diag.folds_unscaled = lda.folds_unscaled;

  out.dscore.assign(n, 0.0);
  for (std::size_t k = 0; k < n; ++k) { out.dscore[order[k]] = lda.dscore[k]; }

  // Best row per precursor: the first maximum in canonical order (lowest feature id on a tie).
  out.groups.resize(groups);
  std::vector<std::int64_t> group_pair(groups);
  std::vector<int> group_label(groups);
  std::vector<double> group_score(groups);
  for (std::size_t g = 0; g < groups; ++g)
  {
    std::size_t best = group_first[g];
    for (std::size_t k = group_first[g] + 1; k < group_first[g + 1]; ++k)
    {
      if (lda.dscore[k] > lda.dscore[best]) { best = k; }
    }
    GroupResult& gr = out.groups[g];
    gr.group = group[best];
    gr.pair = pair[best];
    gr.is_decoy = label[best] == 0;
    gr.rows = group_first[g + 1] - group_first[g];
    gr.best_row = order[best];
    gr.score = lda.dscore[best];
    group_pair[g] = gr.pair;
    group_label[g] = label[best];
    group_score[g] = gr.score;
  }

  const Competition comp = concatenatedCompetition(group_pair, group_label, group_score, options.estimator);
  const std::vector<double> pooled = pooledQValues(group_label, group_score);
  for (std::size_t g = 0; g < groups; ++g)
  {
    GroupResult& gr = out.groups[g];
    gr.winner = comp.winner[g] != 0;
    gr.qvalue = comp.qvalue[g];
    gr.pep = comp.pep[g];
    gr.pooled_qvalue = pooled[g];
    if (gr.winner && gr.qvalue <= options.report_q) { (gr.is_decoy ? diag.decoys_at_q : diag.targets_at_q)++; }
    if (!gr.is_decoy && gr.pooled_qvalue <= options.report_q) { ++diag.pooled_targets_at_q; }
  }
  diag.pairs_complete = comp.pairs_complete;
  diag.targets_unpaired = comp.targets_unpaired;
  diag.decoys_unpaired = comp.decoys_unpaired;
  diag.target_winners = comp.target_winners;
  diag.decoy_winners = comp.decoy_winners;
  diag.pooled_vs_paired = (static_cast<double>(diag.pooled_targets_at_q) + 1.0) /
                          (static_cast<double>(diag.targets_at_q) + 1.0);
  return out;
}

/// Entity-level (peptide or protein) q-values from the precursor results.
struct EntityResult
{
  std::vector<Entity> entities;      ///< one per (key, label), with q, p and PEP
  std::vector<double> group_qvalue;  ///< per ScoredResult::groups entry: its entity's q
  std::vector<double> group_pep;     ///< per ScoredResult::groups entry: its entity's PEP
  std::size_t pairs = 0;             ///< complete target/decoy entity pairs
};

/// Roll each precursor's best d-score up to `key_of_group[g]` (aligned with result.groups; an
/// empty key leaves that precursor out, q 1) and control the FDR at that level: the entity
/// takes its best member's score (context FDR).
///
/// A DECOY precursor must be given its TARGET's key (the target's modified sequence, the
/// target's protein group string): the target entity and the decoy entity of one key are then a
/// structural pair. `picked` = true collapses each such pair to its winner first, ties to the
/// decoy (picked protein-group FDR); entities that lose get q = PEP = 1. `picked` = false ranks
/// every entity (the pooled context estimator, as for peptides).
inline EntityResult entityQValues(const ScoredResult& result, const std::vector<std::string>& key_of_group,
                                  bool picked, QEstimator estimator = QEstimator::Ratio)
{
  if (key_of_group.size() != result.groups.size())
  {
    throw std::invalid_argument("odia::core::entityQValues: " + std::to_string(key_of_group.size()) +
                                " keys for " + std::to_string(result.groups.size()) + " precursors");
  }
  std::vector<double> score(result.groups.size());
  std::vector<int> label(result.groups.size());
  for (std::size_t g = 0; g < result.groups.size(); ++g)
  {
    score[g] = result.groups[g].score;
    label[g] = result.groups[g].is_decoy ? 0 : 1;
  }
  EntityResult out;
  out.entities = rollUp(key_of_group, score, label);
  std::vector<Entity> ranked = out.entities;
  if (picked) { ranked = pickedCompetition(ranked, &out.pairs); }
  else
  {
    std::unordered_map<std::int64_t, int> members;
    for (const auto& e : ranked) { ++members[e.pair]; }
    for (const auto& kv : members) { if (kv.second == 2) { ++out.pairs; } }
  }
  assignQValues(ranked, false, estimator);

  // Back to every entity (losers keep q = PEP = 1), then to every precursor.
  std::unordered_map<std::string, std::size_t> index[2];
  for (std::size_t e = 0; e < out.entities.size(); ++e)
  {
    index[out.entities[e].label].emplace(out.entities[e].id, e);
    out.entities[e].qvalue = 1.0;
    out.entities[e].pvalue = 1.0;
    out.entities[e].pep = 1.0;
  }
  for (const auto& r : ranked) { out.entities[index[r.label].at(r.id)] = r; }
  out.group_qvalue.assign(result.groups.size(), 1.0);
  out.group_pep.assign(result.groups.size(), 1.0);
  for (std::size_t g = 0; g < result.groups.size(); ++g)
  {
    if (key_of_group[g].empty()) { continue; }
    const Entity& e = out.entities[index[label[g]].at(key_of_group[g])];
    out.group_qvalue[g] = e.qvalue;
    out.group_pep[g] = e.pep;
  }
  return out;
}

} // namespace odia::core

#endif // ODIA_CORE_CORE_H
