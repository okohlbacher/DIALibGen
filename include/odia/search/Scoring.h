// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

/// Candidate peak groups in, q-values out: the odia-core classifier and FDR
/// (src/odia-core/) applied to one search, plus the run-level guards.
///
/// Precursor level: semi-supervised LDA over the sub-scores, the best peak
/// group per precursor, concatenated target-decoy competition within each
/// structural pair (ties to the decoy), q = (D + 1) / T. Peptide level
/// (Global.Q.Value): best precursor per canonical modified sequence; protein
/// level (PG.Q.Value): best precursor per Protein.Group string; both by picked
/// competition with the same estimator. A decoy is keyed by its target's
/// sequence and protein group, so target and decoy entities pair structurally
/// at every level.
///
/// Diagnostics, never gates: the pooled precursor estimate (no competition;
/// conservative when decoys of present targets light up with them), and the
/// null-pair balance -- targets vs decoys among the lowest-scoring pair
/// winners, which must be about even when decoys are exchangeable with null
/// targets.
#pragma once

#include <odia/search/CandidateSelector.h>
#include <odia/search/SearchParams.h>

#include <odia-core/odia_core.h>

#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace ODIA::search
{
  /// What extraction hands to scoring: one row per candidate peak group.
  ///
  /// `scores` carries the classifier's view (sub-scores, precursor = SearchSet
  /// index, pair id, decoy flag, feature id); the parallel vectors carry what
  /// the report needs from the chosen row. Fill it through add(), which takes
  /// the ids and labels from the SearchSet so an extractor cannot get them wrong.
  struct PeakGroups
  {
    odia::core::ScoreTable scores;
    std::vector<float> apex_rt;    ///< seconds, run time
    std::vector<float> rt_start;   ///< seconds (OpenSWATH leftWidth)
    std::vector<float> rt_stop;    ///< seconds (OpenSWATH rightWidth)
    std::vector<float> im;         ///< 1/K0 of the peak group; NaN when none

    /// Sub-score names, fixed before the first row (OpenSWATH meta value names).
    void setColumns(std::vector<std::string> names) { scores.setColumns(std::move(names)); }
    std::size_t rows() const { return scores.rows(); }

    /// One peak group of SearchSet precursor @p precursor. @p feature_id must be
    /// unique within the precursor and independent of thread scheduling and
    /// chunking (e.g. the peak group's rank by apex RT); it fixes the canonical
    /// row order and with it the result. @p values holds scores.width() values.
    std::size_t add(const SearchSet& set, std::size_t precursor, std::int64_t feature_id, const float* values,
                    float apex_rt_s, float rt_start_s, float rt_stop_s, float one_over_k0);

    /// Throws std::invalid_argument unless the table and the parallel vectors
    /// agree and every row's ids and label match the SearchSet.
    void validate(const SearchSet& set) const;

    /// Release everything.
    void clear();
  };

  /// The sub-scores search:rt_im_scores false removes from the fit.
  const std::vector<std::string>& rtImScoreNames();

  struct ScoringOutcome
  {
    odia::core::ScoredResult scored;     ///< per precursor with peak groups; groups[g].group is the SearchSet index
    std::vector<double> peptide_q;       ///< per scored.groups entry: Global.Q.Value
    std::vector<double> protein_q;       ///< per scored.groups entry: PG.Q.Value
    std::size_t peptide_pairs = 0;       ///< complete target/decoy peptide pairs
    std::size_t protein_pairs = 0;       ///< complete target/decoy protein-group pairs
    std::size_t peptides_at_q = 0;       ///< target peptides at q <= 0.01
    std::size_t proteins_at_q = 0;       ///< target protein groups at q <= 0.01

    bool entrapment = false;             ///< search:entrapment_tag was set
    odia::core::EntrapmentEstimate entrapment_estimate;
    std::string entrapment_db_basis;     ///< where the estimator's database ratio was counted
    std::size_t entrapment_shared = 0;   ///< identified targets whose group mixes tagged and untagged proteins (left out)

    /// Null-pair balance: winners of complete pairs in the lowest
    /// SearchParams::null_balance_fraction of winning scores, by class.
    std::size_t null_targets = 0;
    std::size_t null_decoys = 0;
    double null_ratio = 0.0;             ///< null_targets / null_decoys (NaN without decoys)
    double null_z = 0.0;                 ///< (T - D) / sqrt(T + D); > 0: decoys weaker than null targets

    bool selftest = false;               ///< search:selftest ran
    std::size_t selftest_label_swap_ids = 0;
    std::size_t selftest_random_label_ids = 0;
    std::size_t selftest_limit = 0;      ///< the most either check may report

    std::vector<std::string> warnings;   ///< non-fatal findings, also for the provenance
  };

  /// Guards on the extraction itself, before scoring: something was extracted,
  /// and the decoy:target ratio of precursors with peak groups is inside
  /// [decoy_ratio_low, decoy_ratio_high]. Throws SearchAbort with the counts.
  void checkExtraction(const SearchSet& set, const PeakGroups& groups, const SearchParams& params);

  /// Score, compete, roll up; the self-checks when search:selftest is set.
  /// Throws std::invalid_argument (from odia-core) on a table that cannot be
  /// scored honestly.
  ScoringOutcome scorePeakGroups(const SearchSet& set, const PeakGroups& groups, const SearchParams& params);

  /// The run-level guards after scoring. Each throws SearchAbort with its counts:
  /// no discriminant learned (iterations_trained == 0); more than
  /// max_target_fraction of the scored targets identified; fewer than min_ids
  /// identified (only with @p min_ids: Identifier writes the report of such a
  /// search for inspection first, then applies it); a failed self-check.
  void checkGuards(const ScoringOutcome& outcome, const SearchParams& params, bool min_ids = true);

  /// Identified target precursors (winners at q <= 0.01).
  std::size_t identifications(const ScoringOutcome& outcome);

  /// A protein group's entrapment class: Trap when every member starts with
  /// @p tag, Real when none does, Shared otherwise.
  enum class Entrapment { Real, Trap, Shared };
  Entrapment entrapmentClass(std::string_view group, const std::string& tag);
}
