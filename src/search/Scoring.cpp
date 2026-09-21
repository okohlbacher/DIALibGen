// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <odia/search/Scoring.h>

#include <odia/LibraryRefiner.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

namespace ODIA::search
{
  namespace
  {
    std::string percent(double part, double whole)
    {
      std::ostringstream s;
      s.setf(std::ios::fixed);
      s.precision(1);
      s << (whole > 0 ? 100.0 * part / whole : 0.0) << "%";
      return s.str();
    }

    std::string fixed(double v, int digits)
    {
      std::ostringstream s;
      s.setf(std::ios::fixed);
      s.precision(digits);
      s << v;
      return s.str();
    }

    odia::core::Options classifierOptions(const SearchParams& params)
    {
      odia::core::Options o;
      o.lda.threads = params.threads;
      o.estimator = odia::core::QEstimator::Count;
      o.report_q = SearchParams::identification_q;
      if (!params.rt_im_scores) { o.exclude_features = rtImScoreNames(); }
      return o;
    }

    /// A fair coin per pair, from the pair id and the search seed (SplitMix64).
    bool coin(std::int64_t pair, std::uint64_t seed)
    {
      std::uint64_t z = static_cast<std::uint64_t>(pair) ^ (seed * 0x9E3779B97F4A7C15ull);
      z += 0x9E3779B97F4A7C15ull;
      z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
      z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
      z ^= z >> 31;
      return (z >> 63) != 0;
    }

    enum class Entrapment { Real, Trap, Shared };
    Entrapment entrapmentClass(std::string_view group, const std::string& tag)
    {
      bool any_trap = false, any_real = false;
      std::size_t start = 0;
      while (start <= group.size())
      {
        const std::size_t end = std::min(group.find(';', start), group.size());
        const std::string_view member = group.substr(start, end - start);
        if (!member.empty()) { (member.substr(0, tag.size()) == tag ? any_trap : any_real) = true; }
        start = end + 1;
      }
      if (any_trap && any_real) { return Entrapment::Shared; }
      return any_trap ? Entrapment::Trap : Entrapment::Real;
    }
  }

  std::size_t PeakGroups::add(const SearchSet& set, std::size_t precursor, std::int64_t feature_id, const float* values,
                              float apex_rt_s, float rt_start_s, float rt_stop_s, float one_over_k0)
  {
    if (precursor >= set.size())
    { throw std::out_of_range("PeakGroups::add: precursor " + std::to_string(precursor) + " of " + std::to_string(set.size())); }
    const std::size_t r = scores.append(static_cast<std::int64_t>(precursor), set.pairOf(precursor),
                                        set.isDecoy(precursor), feature_id, values);
    apex_rt.push_back(apex_rt_s);
    rt_start.push_back(rt_start_s);
    rt_stop.push_back(rt_stop_s);
    im.push_back(one_over_k0);
    return r;
  }

  void PeakGroups::validate(const SearchSet& set) const
  {
    scores.validate();
    const std::size_t n = rows();
    if (apex_rt.size() != n || rt_start.size() != n || rt_stop.size() != n || im.size() != n)
    {
      throw std::invalid_argument("PeakGroups: " + std::to_string(n) + " score rows but " + std::to_string(apex_rt.size()) +
                                  " apex RTs, " + std::to_string(rt_start.size()) + "/" + std::to_string(rt_stop.size()) +
                                  " RT bounds and " + std::to_string(im.size()) + " 1/K0 values");
    }
    for (std::size_t r = 0; r < n; ++r)
    {
      const std::int64_t g = scores.group[r];
      if (g < 0 || static_cast<std::size_t>(g) >= set.size())
      { throw std::invalid_argument("PeakGroups: row " + std::to_string(r) + " names precursor " + std::to_string(g) + " of " + std::to_string(set.size())); }
      const std::size_t i = static_cast<std::size_t>(g);
      if (scores.pair[r] != set.pairOf(i) || (scores.is_decoy[r] != 0) != set.isDecoy(i))
      { throw std::invalid_argument("PeakGroups: row " + std::to_string(r) + " carries a pair id or label that precursor " + std::to_string(g) + " does not have"); }
    }
  }

  void PeakGroups::clear()
  {
    scores.clear();
    std::vector<float>().swap(apex_rt);
    std::vector<float>().swap(rt_start);
    std::vector<float>().swap(rt_stop);
    std::vector<float>().swap(im);
  }

  const std::vector<std::string>& rtImScoreNames()
  {
    static const std::vector<std::string> names = {"var_norm_rt_score", "var_im_delta_score", "var_im_ms1_delta_score"};
    return names;
  }

  void checkExtraction(const SearchSet& set, const PeakGroups& groups, const SearchParams& params)
  {
    (void)params;
    const std::size_t pairs = set.pairs();
    if (groups.rows() == 0)
    {
      throw SearchAbort("search: extraction found no peak group for any of the " + std::to_string(set.size()) +
                        " precursors (" + std::to_string(pairs) + " target-decoy pairs); nothing to score");
    }
    std::vector<char> seen(set.size(), 0);
    for (const std::int64_t g : groups.scores.group)
    {
      if (g >= 0 && static_cast<std::size_t>(g) < seen.size()) { seen[static_cast<std::size_t>(g)] = 1; }
    }
    const auto targets = static_cast<std::size_t>(std::count(seen.begin(), seen.begin() + static_cast<std::ptrdiff_t>(pairs), 1));
    const auto decoys = static_cast<std::size_t>(std::count(seen.begin() + static_cast<std::ptrdiff_t>(pairs), seen.end(), 1));
    const double ratio = targets > 0 ? static_cast<double>(decoys) / static_cast<double>(targets) : 0.0;
    if (targets == 0 || decoys == 0 || ratio < SearchParams::decoy_ratio_low || ratio > SearchParams::decoy_ratio_high)
    {
      throw SearchAbort("search: " + std::to_string(decoys) + " decoy and " + std::to_string(targets) +
                        " target precursors have peak groups (of " + std::to_string(pairs) + " pairs searched, " +
                        std::to_string(groups.rows()) + " peak groups); the decoy:target ratio " + fixed(ratio, 3) +
                        " is outside [" + fixed(SearchParams::decoy_ratio_low, 2) + ", " +
                        fixed(SearchParams::decoy_ratio_high, 2) + "], so pairing or extraction is broken and no "
                        "q-value would mean anything");
    }
  }

  ScoringOutcome scorePeakGroups(const SearchSet& set, const PeakGroups& groups, const SearchParams& params)
  {
    groups.validate(set);
    const odia::core::Options options = classifierOptions(params);
    ScoringOutcome out;
    out.scored = odia::core::scoreAndControl(groups.scores, options);
    const auto& result = out.scored.groups;
    const auto& diag = out.scored.diagnostics;

    // Entity keys. A decoy carries its target's sequence and protein group, so
    // both classes of one key form a structural pair at every level.
    std::vector<std::string> peptide(result.size()), protein(result.size());
    std::unordered_map<std::uint32_t, std::string> canonical;
    for (std::size_t g = 0; g < result.size(); ++g)
    {
      const auto i = static_cast<std::size_t>(result[g].group);
      const std::uint32_t handle = set.library.precursors().modified_sequence[i];
      auto it = canonical.find(handle);
      if (it == canonical.end()) { it = canonical.emplace(handle, canonicalModifiedSequence(set.modifiedSequence(i))).first; }
      peptide[g] = it->second;
      protein[g] = std::string(set.proteinGroup(i));
    }
    const auto peptides = odia::core::entityQValues(out.scored, peptide, false);
    const auto proteins = odia::core::entityQValues(out.scored, protein, true);
    out.peptide_q = peptides.group_qvalue;
    out.protein_q = proteins.group_qvalue;
    out.peptide_pairs = peptides.pairs;
    out.protein_pairs = proteins.pairs;
    out.peptides_at_q = odia::core::countAtQ(peptides.entities, SearchParams::identification_q);
    out.proteins_at_q = odia::core::countAtQ(proteins.entities, SearchParams::identification_q);

    if (!params.entrapment_tag.empty())
    {
      out.entrapment = true;
      std::size_t db_real = 0, db_trap = 0, reported = 0, trapped = 0;
      for (std::size_t k = 0; k < set.pairs(); ++k)
      {
        const auto c = entrapmentClass(set.proteinGroup(k), params.entrapment_tag);
        if (c == Entrapment::Real) { ++db_real; }
        else if (c == Entrapment::Trap) { ++db_trap; }
      }
      for (const auto& gr : result)
      {
        if (gr.is_decoy || !gr.winner || gr.qvalue > SearchParams::identification_q) { continue; }
        const auto c = entrapmentClass(set.proteinGroup(static_cast<std::size_t>(gr.group)), params.entrapment_tag);
        if (c == Entrapment::Shared) { ++out.entrapment_shared; continue; }
        ++reported;
        if (c == Entrapment::Trap) { ++trapped; }
      }
      out.entrapment_estimate = odia::core::entrapmentFdp(reported, trapped, db_real, db_trap);
      if (!out.entrapment_estimate.valid)
      { out.warnings.push_back("entrapment: no estimate (" + std::to_string(db_trap) + " entrapment and " + std::to_string(db_real) +
                               " real target precursors searched, " + std::to_string(reported) + " identified)"); }
      else if (out.entrapment_estimate.fdp > 0.015)
      { out.warnings.push_back("entrapment: combined FDP estimate " + percent(out.entrapment_estimate.fdp, 1.0) +
                               " at nominal 1% (" + std::to_string(trapped) + " entrapment of " + std::to_string(reported) + " identified)"); }
    }

    if (diag.pooled_vs_paired > SearchParams::pooled_vs_paired_warn || diag.pooled_vs_paired < 1.0 / SearchParams::pooled_vs_paired_warn)
    {
      out.warnings.push_back("pooled and paired precursor FDR disagree: " + std::to_string(diag.pooled_targets_at_q) + " vs " +
                             std::to_string(diag.targets_at_q) + " identifications at q <= 0.01; the pairs may not behave as "
                             "exchangeable target-decoy pairs");
    }
    if (diag.folds_unscaled > 0)
    { out.warnings.push_back(std::to_string(diag.folds_unscaled) + " cross-validation fold(s) could not be rescaled (too few precursors)"); }

    if (params.selftest)
    {
      out.selftest = true;
      out.selftest_limit = diag.targets_at_q / 100 + 1;
      // Label swap: decoys play targets. Real signal now sits on the "decoy"
      // side, so essentially nothing may pass.
      odia::core::ScoreTable swapped = groups.scores;
      for (auto& d : swapped.is_decoy) { d = d ? 0 : 1; }
      out.selftest_label_swap_ids = odia::core::scoreAndControl(swapped, options).diagnostics.targets_at_q;
      // Random labels: each pair's two labels exchanged on a fair coin. Targets
      // and decoys are then exchangeable by construction -- a pure null.
      odia::core::ScoreTable shuffled = groups.scores;
      for (std::size_t r = 0; r < shuffled.rows(); ++r)
      {
        if (coin(shuffled.pair[r], params.seed)) { shuffled.is_decoy[r] = shuffled.is_decoy[r] ? 0 : 1; }
      }
      out.selftest_random_label_ids = odia::core::scoreAndControl(shuffled, options).diagnostics.targets_at_q;
    }
    return out;
  }

  std::size_t identifications(const ScoringOutcome& outcome) { return outcome.scored.diagnostics.targets_at_q; }

  void checkGuards(const ScoringOutcome& outcome, const SearchParams& params)
  {
    const auto& d = outcome.scored.diagnostics;
    const std::string counts = std::to_string(d.target_groups) + " target and " + std::to_string(d.decoy_groups) +
                               " decoy precursors scored";
    if (d.iterations_trained == 0)
    {
      throw SearchAbort("search: the classifier learned no discriminant (0 iterations trained, " +
                        std::to_string(d.iterations_skipped) + " skipped for want of confident targets; " + counts + ", " +
                        std::to_string(d.features_used.size()) + " sub-scores); the scores would be the untrained initial score");
    }
    const std::size_t ids = d.targets_at_q;
    if (static_cast<double>(ids) > params.max_target_fraction * static_cast<double>(d.target_groups))
    {
      throw SearchAbort("search: " + std::to_string(ids) + " of " + std::to_string(d.target_groups) + " scored target precursors (" +
                        percent(static_cast<double>(ids), static_cast<double>(d.target_groups)) + ") pass q <= 0.01, more than "
                        "search:max_target_fraction " + fixed(params.max_target_fraction, 2) + " allows (" +
                        std::to_string(d.decoys_at_q) + " decoys pass); decoys this weak make the FDR estimate meaningless");
    }
    if (ids < params.min_ids)
    {
      throw SearchAbort("search: " + std::to_string(ids) + " target precursors pass q <= 0.01 (" + std::to_string(d.decoys_at_q) +
                        " decoys; " + counts + "), fewer than search:min_ids " + std::to_string(params.min_ids) +
                        "; too few for refinement or tuning");
    }
    if (outcome.selftest &&
        (outcome.selftest_label_swap_ids > outcome.selftest_limit || outcome.selftest_random_label_ids > outcome.selftest_limit))
    {
      throw SearchAbort("search: self-check failed: with labels swapped " + std::to_string(outcome.selftest_label_swap_ids) +
                        " and with random pair labels " + std::to_string(outcome.selftest_random_label_ids) +
                        " precursors pass q <= 0.01, where at most " + std::to_string(outcome.selftest_limit) +
                        " may (1% of the " + std::to_string(ids) + " identifications, plus one)");
    }
  }
}
