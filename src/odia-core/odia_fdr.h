// Copyright (c) 2026, Oliver Kohlbacher and the OpenDIAlyzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
//
// odia_fdr.h -- target-decoy q-values, p-values and PEP, at any level.
//
// Vendored from OpenDIAlyzer (src/odia_fdr.h, plus the q-value estimator that lived in
// src/odia_lda.h as lda_detail::assignQValues). Origin commit, file hashes and every change
// against the original are recorded in src/odia-core/MANIFEST.json. Dependency-free C++17:
// no OpenMS, no Eigen.
//
// One estimator, three levels. The precursor level ranks one best score per precursor, after
// CONCATENATED target-decoy competition within each pair (concatenatedCompetition); the peptide
// and protein levels roll those scores up to the entity and rank again. Pairs are always
// STRUCTURAL -- an integer pair id supplied by the caller -- never inferred from names.
//
//   1. CONTEXT FDR -- peptide- and protein-level q-values by rolling precursor scores up to
//      the entity and re-running target-decoy on that set (Rosenberger 2017, PyProphet).
//      A protein collects the best of many precursor draws, so its score is an extreme-value
//      statistic and its error rate is worse than the precursor one; a precursor q-value must
//      never be reported as a protein q-value.
//   2. PICKED competition (Savitski 2015; The 2022) -- each target competes with ITS OWN decoy
//      partner and enters the ranking once, as the winner. Removes the size bias of plain
//      target-decoy at the protein level: a pair's members have the same size by construction.
//   3. ENTRAPMENT FDP (Wen/Freestone/Keich/Noble 2025) -- a validation rig, not a model.

#ifndef ODIA_CORE_FDR_H
#define ODIA_CORE_FDR_H

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace odia::core
{

/// One item of a ranked list: a precursor's best peak group, or an entity.
struct RankedGroup
{
  std::size_t group_index = 0;   ///< caller's index of the item (also the final tie-break)
  std::size_t best_row = 0;      ///< caller's row behind the score (passed through)
  int label = 1;                 ///< 1 = target, 0 = decoy
  double score = 0.0;
  double qvalue = 1.0;
  double pvalue = 1.0;   ///< empirical p from the decoy null (TAIL probability at this score)
  double pep = 1.0;      ///< posterior error probability = LOCAL false-discovery rate at this score
};

/// Project the `pep` of a list sorted by score descending onto a non-decreasing sequence (PEP
/// must not fall as the score falls): weighted pool-adjacent-violators, the least-squares
/// isotonic fit. Items of equal score form one block and get one value.
///
/// ODIA swept a running MAX from the low-score end instead. That is an upper envelope of
/// everything BELOW an item, and the bottom of any list is at PEP 1, so every PEP came out 1.
inline void isotonicNonDecreasing(std::vector<RankedGroup>& ranked)
{
  struct Block { double sum; double weight; std::size_t end; };
  std::vector<Block> blocks;
  for (std::size_t begin = 0; begin < ranked.size();)
  {
    std::size_t end = begin + 1;
    while (end < ranked.size() && ranked[end].score == ranked[begin].score) { ++end; }
    double sum = 0.0;
    for (std::size_t i = begin; i < end; ++i) { sum += ranked[i].pep; }
    blocks.push_back({sum, static_cast<double>(end - begin), end});
    while (blocks.size() >= 2)
    {
      const Block& b = blocks[blocks.size() - 1];
      const Block& a = blocks[blocks.size() - 2];
      if (a.sum / a.weight <= b.sum / b.weight) { break; }
      const Block merged{a.sum + b.sum, a.weight + b.weight, b.end};
      blocks.pop_back();
      blocks.back() = merged;
    }
    begin = end;
  }
  std::size_t begin = 0;
  for (const Block& b : blocks)
  {
    const double value = std::min(1.0, std::max(0.0, b.sum / b.weight));
    for (std::size_t i = begin; i < b.end; ++i) { ranked[i].pep = value; }
    begin = b.end;
  }
}

/// How a ranked list turns counts into an FDR estimate at a threshold t, where T(t) and D(t) are
/// the targets and decoys at or above t and N_tar, N_dec those in the whole list.
enum class QEstimator
{
  /// pi0 * ((D + 1) / N_dec) / (T / N_tar): ODIA's estimator. On a list of pair WINNERS this is
  /// (D_win + 1) / N_dec_win divided by T_win / N_tar_win.
  Ratio,
  /// pi0 * (D + 1) / T: the classic target-decoy-competition count.
  Count
};

/// Target-decoy q-values, p-values and PEP for a list with one score per item.
///
/// Equal scores are one threshold, so the result does not depend on the order of ties. The list
/// is sorted in place (score descending, then group_index ascending).
///
///   q   : FDR(t) by @p estimator, monotonised from the bottom (a running minimum). The +1 is
///         Kall's finite-sample correction on the decoy count. A list with no decoys at all has
///         no null to scale by, so Ratio falls back to Count there: q = 1 / T, never 0.
///   p   : (decoys at or above + 1) / (N_dec + 1).
///   PEP : the local analogue of the q estimator in a sliding window of the ranked list, made
///         non-increasing in score by isotonic regression (isotonicNonDecreasing).
///
/// pi0 = 1 unless use_pi0 (Storey, lambda = 0.5). pi0 = 1 is the honest default.
inline void assignQValues(std::vector<RankedGroup>& ranked, bool use_pi0,
                          QEstimator estimator = QEstimator::Ratio)
{
  std::sort(ranked.begin(), ranked.end(), [](const RankedGroup& a, const RankedGroup& b) {
    if (a.score != b.score) { return a.score > b.score; }
    return a.group_index < b.group_index;
  });

  std::size_t Ntar = 0, Ndec = 0;
  for (const auto& r : ranked) { if (r.label == 1) { ++Ntar; } else { ++Ndec; } }
  // The ratio N_tar / N_dec by which the Ratio estimator scales decoy counts to expected null
  // targets; 1 for Count, and for a list without decoys (see above).
  const bool ratio = estimator == QEstimator::Ratio && Ndec > 0;
  const double scale = ratio ? static_cast<double>(Ntar) / static_cast<double>(Ndec) : 1.0;

  // Storey pi0: walking high->low, a target's empirical p-value from the decoy null is
  // decoys_seen / Ndec; null targets have ~uniform p, so the mass with p > lambda estimates pi0.
  double pi0 = 1.0;
  if (use_pi0 && Ntar > 0 && Ndec > 0)
  {
    const double lambda = 0.5;
    const std::size_t dthr = static_cast<std::size_t>(lambda * static_cast<double>(Ndec));
    std::size_t dec_seen = 0, tar_hi = 0;
    for (const auto& r : ranked)
    {
      if (r.label == 0) { ++dec_seen; }
      else if (dec_seen > dthr) { ++tar_hi; }
    }
    pi0 = static_cast<double>(tar_hi) / ((1.0 - lambda) * static_cast<double>(Ntar));
    if (!(pi0 > 0.0)) { pi0 = 1.0 / static_cast<double>(Ntar); }   // never 0
    if (pi0 > 1.0) { pi0 = 1.0; }                                   // never > 1
  }

  std::size_t targets = 0;
  std::size_t decoys = 0;
  for (std::size_t begin = 0; begin < ranked.size();)
  {
    std::size_t end = begin + 1;
    while (end < ranked.size() && ranked[end].score == ranked[begin].score) { ++end; }
    for (std::size_t i = begin; i < end; ++i)
    {
      if (ranked[i].label == 1) { ++targets; }
      else                      { ++decoys; }
    }
    // (D + 1) * scale / T is ((D + 1) / N_dec) / (T / N_tar) for Ratio and (D + 1) / T for Count.
    double fdr;
    if (targets == 0) { fdr = std::numeric_limits<double>::infinity(); }
    else
    {
      fdr = std::min(1.0, pi0 * (static_cast<double>(decoys) + 1.0) * scale /
                            static_cast<double>(targets));
    }
    for (std::size_t i = begin; i < end; ++i) { ranked[i].qvalue = fdr; }
    begin = end;
  }

  double running_min = 1.0;
  for (std::size_t end = ranked.size(); end > 0;)
  {
    std::size_t begin = end - 1;
    while (begin > 0 && ranked[begin - 1].score == ranked[end - 1].score) { --begin; }
    running_min = std::min(running_min, ranked[begin].qvalue);
    for (std::size_t i = begin; i < end; ++i) { ranked[i].qvalue = running_min; }
    end = begin;
  }

  // p-value: conservative empirical tail probability (Kall's +1 on both counts).
  {
    std::size_t dec_at_or_above = 0;
    for (std::size_t begin = 0; begin < ranked.size();)
    {
      std::size_t end = begin + 1;
      while (end < ranked.size() && ranked[end].score == ranked[begin].score) { ++end; }
      for (std::size_t i = begin; i < end; ++i) { if (ranked[i].label == 0) { ++dec_at_or_above; } }
      const double p = (Ndec > 0)
                         ? (static_cast<double>(dec_at_or_above) + 1.0) / (static_cast<double>(Ndec) + 1.0)
                         : 1.0;
      for (std::size_t i = begin; i < end; ++i) { ranked[i].pvalue = std::min(1.0, p); }
      begin = end;
    }
  }

  // PEP: in a score neighbourhood holding t targets and d decoys, the expected null-target count
  // is pi0 * d * scale (the same scale as q), so PEP ~ that over t. Sliding counts keep the pass
  // O(n); the window is ~0.5 % of the list and at least 101 wide.
  {
    const std::size_t n = ranked.size();
    const std::size_t half = std::max<std::size_t>(50, n / 200);
    std::size_t t = 0, d = 0, lo = 0, hi = 0;
    for (std::size_t i = 0; i < n; ++i)
    {
      const std::size_t new_lo = (i > half) ? i - half : 0;
      const std::size_t new_hi = std::min(n, i + half + 1);
      while (hi < new_hi) { if (ranked[hi].label == 1) { ++t; } else { ++d; } ++hi; }
      while (lo < new_lo) { if (ranked[lo].label == 1) { --t; } else { --d; } ++lo; }
      double pep = 1.0;
      if (t > 0) { pep = pi0 * static_cast<double>(d) * scale / static_cast<double>(t); }
      ranked[i].pep = std::min(1.0, std::max(0.0, pep));
    }
    isotonicNonDecreasing(ranked);
  }
}

// ------------------------------------------------------------------------------------------
// Structural pairing. A target and its decoy are a PAIR because the caller says so, by giving
// both the same integer pair id -- never because of how their names are spelled. ODIA paired a
// decoy with the target whose id it equals after stripping a 'DECOY_' prefix; DIALibGen's decoys
// are not named that way, and a naming convention is not a structure.
// ------------------------------------------------------------------------------------------

namespace detail
{

/// Winners of target-decoy pair competition over items (pair[i], label[i], score[i]).
/// An item whose pair id is negative, or whose pair has no member of the other label, has no
/// competitor and wins. Otherwise the target wins only if it scores STRICTLY higher: a tie
/// carries no evidence, and resolving it for the target would bias the estimate downward.
/// Throws std::invalid_argument if a pair id carries two targets or two decoys.
struct PairWinners
{
  std::vector<char> winner;
  std::size_t pairs_complete = 0;     ///< pairs with both members present
  std::size_t targets_unpaired = 0;   ///< targets with no decoy present
  std::size_t decoys_unpaired = 0;    ///< decoys with no target present
};

inline PairWinners pairWinners(const std::vector<std::int64_t>& pair, const std::vector<int>& label,
                               const std::vector<double>& score, const char* who)
{
  constexpr std::size_t none = std::numeric_limits<std::size_t>::max();
  PairWinners out;
  const std::size_t n = pair.size();
  out.winner.assign(n, 1);
  std::unordered_map<std::int64_t, std::pair<std::size_t, std::size_t>> members;   // target, decoy
  members.reserve(n);
  for (std::size_t i = 0; i < n; ++i)
  {
    if (pair[i] < 0) { (label[i] == 1 ? out.targets_unpaired : out.decoys_unpaired)++; continue; }
    auto& slot = members.emplace(pair[i], std::make_pair(none, none)).first->second;
    std::size_t& mine = label[i] == 1 ? slot.first : slot.second;
    if (mine != none)
    {
      throw std::invalid_argument(std::string(who) + ": pair " + std::to_string(pair[i]) +
                                  " has two " + (label[i] == 1 ? "targets" : "decoys") +
                                  " (items " + std::to_string(mine) + " and " +
                                  std::to_string(i) + ")");
    }
    mine = i;
  }
  for (const auto& kv : members)
  {
    const std::size_t t = kv.second.first, d = kv.second.second;
    if (t == none) { ++out.decoys_unpaired; continue; }
    if (d == none) { ++out.targets_unpaired; continue; }
    ++out.pairs_complete;
    const bool target_wins = score[t] > score[d];
    out.winner[t] = target_wins ? 1 : 0;
    out.winner[d] = target_wins ? 0 : 1;
  }
  return out;
}

} // namespace detail

/// Result of concatenated target-decoy competition, per input item.
struct Competition
{
  std::vector<char> winner;     ///< the item won its pair (or had no competitor)
  std::vector<double> qvalue;   ///< the winner's q; 1 for an item that lost its pair
  std::vector<double> pep;      ///< the winner's PEP; 1 for an item that lost its pair
  std::size_t pairs_complete = 0;
  std::size_t targets_unpaired = 0;
  std::size_t decoys_unpaired = 0;
  std::size_t target_winners = 0;
  std::size_t decoy_winners = 0;
};

/// CONCATENATED target-decoy competition: every pair collapses to its better-scoring member
/// (ties to the decoy), and q-values are computed over the winners alone by @p estimator,
/// monotonised. With the default, q at a threshold is (D_win + 1) / N_dec_win divided by
/// T_win / N_tar_win. pi0 is 1.
///
/// Why compete: a decoy shares its target's isolation window, retention time and much of its
/// fragment evidence, so a decoy lit up by its PRESENT target would otherwise enter the null
/// tail and distort it. Within the pair it simply loses.
///
/// @p pair, @p label (1 target, 0 decoy) and @p score are per item (e.g. per precursor, its best
/// d-score). Throws std::invalid_argument on size mismatch or a pair id with two targets or two
/// decoys; a negative pair id means "no partner".
inline Competition concatenatedCompetition(const std::vector<std::int64_t>& pair,
                                           const std::vector<int>& label,
                                           const std::vector<double>& score,
                                           QEstimator estimator = QEstimator::Ratio)
{
  const std::size_t n = pair.size();
  if (label.size() != n || score.size() != n)
  {
    throw std::invalid_argument("odia::core::concatenatedCompetition: " + std::to_string(n) +
                                " pair ids, " + std::to_string(label.size()) + " labels, " +
                                std::to_string(score.size()) + " scores");
  }
  const detail::PairWinners pw = detail::pairWinners(pair, label, score,
                                                     "odia::core::concatenatedCompetition");
  Competition out;
  out.winner = pw.winner;
  out.qvalue.assign(n, 1.0);
  out.pep.assign(n, 1.0);
  out.pairs_complete = pw.pairs_complete;
  out.targets_unpaired = pw.targets_unpaired;
  out.decoys_unpaired = pw.decoys_unpaired;

  std::vector<RankedGroup> ranked;
  ranked.reserve(n);
  for (std::size_t i = 0; i < n; ++i)
  {
    if (!pw.winner[i]) { continue; }
    RankedGroup r;
    r.group_index = i;
    r.best_row = i;
    r.label = label[i] == 1 ? 1 : 0;
    r.score = score[i];
    ranked.push_back(r);
    (r.label == 1 ? out.target_winners : out.decoy_winners)++;
  }
  assignQValues(ranked, false, estimator);
  for (const auto& r : ranked)
  {
    out.qvalue[r.group_index] = r.qvalue;
    out.pep[r.group_index] = r.pep;
  }
  return out;
}

/// ODIA's POOLED estimator, kept as a diagnostic: every item ranked, no competition, q by the
/// Ratio estimator. On a precursor list it is not the reported q; the ratio of its count at 1 %
/// to the concatenated count is a check (a large disagreement means pairs are not behaving as
/// exchangeable target-decoy pairs).
inline std::vector<double> pooledQValues(const std::vector<int>& label,
                                         const std::vector<double>& score, bool use_pi0 = false)
{
  if (label.size() != score.size())
  {
    throw std::invalid_argument("odia::core::pooledQValues: " + std::to_string(label.size()) +
                                " labels, " + std::to_string(score.size()) + " scores");
  }
  std::vector<RankedGroup> ranked(label.size());
  for (std::size_t i = 0; i < label.size(); ++i)
  {
    ranked[i].group_index = i;
    ranked[i].best_row = i;
    ranked[i].label = label[i] == 1 ? 1 : 0;
    ranked[i].score = score[i];
  }
  assignQValues(ranked, use_pi0, QEstimator::Ratio);
  std::vector<double> q(label.size(), 1.0);
  for (const auto& r : ranked) { q[r.group_index] = r.qvalue; }
  return q;
}

/// One item at whatever level is being controlled: a peptide, or a protein (group).
struct Entity
{
  std::string id;            ///< modified sequence, or protein accession / group key
  int label = 1;             ///< 1 = target, 0 = decoy
  double score = 0.0;        ///< best member score (see rollUp)
  std::int64_t pair = -1;    ///< shared by a target entity and its decoy; < 0 = no partner
  double qvalue = 1.0;
  double pvalue = 1.0;
  double pep = 1.0;
};

/// Roll member scores up to entities: each entity takes the score of its BEST member.
///
/// `member_key[i]` names the entity row i belongs to, `member_score[i]` is its score and
/// `member_label[i]` its class. Rows whose key is empty are skipped (unmapped precursors).
/// The maximum is the "best of n draws" that inflates the entity-level error rate relative to
/// the member one, which is why q-values must be recomputed at this level.
///
/// An entity is a (key, label) pair, so a decoy never merges into a target of the same key
/// (ODIA prefixed decoy keys to avoid exactly that). The target entity and the decoy entity of
/// one key form a structural pair: both get the key's index as their pair id. Use this when a
/// decoy carries its target's key (a protein group string, a target sequence); otherwise use
/// the overload with explicit pair ids. Entities come out in first-occurrence order.
inline std::vector<Entity> rollUp(const std::vector<std::string>& member_key,
                                  const std::vector<double>& member_score,
                                  const std::vector<int>& member_label)
{
  std::vector<Entity> out;
  const std::size_t n = member_key.size();
  if (member_score.size() != n || member_label.size() != n)
  {
    throw std::invalid_argument("odia::core::rollUp: " + std::to_string(n) + " keys, " +
                                std::to_string(member_score.size()) + " scores, " +
                                std::to_string(member_label.size()) + " labels");
  }
  std::unordered_map<std::string, std::int64_t> key_index;
  std::unordered_map<std::string, std::size_t> entity_index[2];   // by label: decoy, target
  key_index.reserve(n);
  for (std::size_t i = 0; i < n; ++i)
  {
    if (member_key[i].empty()) { continue; }
    const int label = member_label[i] == 1 ? 1 : 0;
    const auto key = key_index.emplace(member_key[i], static_cast<std::int64_t>(key_index.size()));
    const auto inserted = entity_index[label].emplace(member_key[i], out.size());
    if (inserted.second)
    {
      Entity e;
      e.id = member_key[i];
      e.label = label;
      e.score = member_score[i];
      e.pair = key.first->second;
      out.push_back(std::move(e));
    }
    else if (member_score[i] > out[inserted.first->second].score)
    {
      out[inserted.first->second].score = member_score[i];
    }
  }
  return out;
}

/// Roll up with EXPLICIT pair ids: `member_pair[i]` is the pair id of row i's entity. Every
/// member of one entity must carry the same pair id; a conflict throws std::invalid_argument.
inline std::vector<Entity> rollUp(const std::vector<std::string>& member_key,
                                  const std::vector<double>& member_score,
                                  const std::vector<int>& member_label,
                                  const std::vector<std::int64_t>& member_pair)
{
  if (member_pair.size() != member_key.size())
  {
    throw std::invalid_argument("odia::core::rollUp: " + std::to_string(member_key.size()) +
                                " keys, " + std::to_string(member_pair.size()) + " pair ids");
  }
  std::vector<Entity> out = rollUp(member_key, member_score, member_label);
  std::unordered_map<std::string, std::size_t> entity_index[2];
  for (std::size_t e = 0; e < out.size(); ++e) { entity_index[out[e].label].emplace(out[e].id, e); }
  std::vector<char> set(out.size(), 0);
  for (std::size_t i = 0; i < member_key.size(); ++i)
  {
    if (member_key[i].empty()) { continue; }
    const std::size_t e = entity_index[member_label[i] == 1 ? 1 : 0].at(member_key[i]);
    if (!set[e]) { out[e].pair = member_pair[i]; set[e] = 1; }
    else if (out[e].pair != member_pair[i])
    {
      throw std::invalid_argument("odia::core::rollUp: entity '" + member_key[i] +
                                  "' has members in pairs " + std::to_string(out[e].pair) +
                                  " and " + std::to_string(member_pair[i]));
    }
  }
  return out;
}

/// Target-decoy q-value / p-value / PEP over an entity set, identical estimator to the
/// precursor level. Entities are updated in place; their order is preserved.
inline void assignQValues(std::vector<Entity>& entities, bool use_pi0,
                          QEstimator estimator = QEstimator::Ratio)
{
  std::vector<RankedGroup> ranked;
  ranked.reserve(entities.size());
  for (std::size_t i = 0; i < entities.size(); ++i)
  {
    RankedGroup r;
    r.group_index = i;
    r.best_row = i;
    r.label = entities[i].label;
    r.score = entities[i].score;
    ranked.push_back(r);
  }
  assignQValues(ranked, use_pi0, estimator);
  for (const auto& r : ranked)
  {
    entities[r.group_index].qvalue = r.qvalue;
    entities[r.group_index].pvalue = r.pvalue;
    entities[r.group_index].pep = r.pep;
  }
}

/// Picked target-decoy competition (Savitski 2015; The 2022), paired STRUCTURALLY by
/// Entity::pair.
///
/// Each complete target/decoy pair collapses to ONE entry, the member with the better score
/// (ties to the DECOY); entities with a negative pair id or no partner pass through unchanged.
/// The output keeps input order (a pair appears where its winner stood). The caller then runs
/// assignQValues() on the result. A pair id carried by two targets or two decoys throws
/// std::invalid_argument. @p n_paired receives the number of complete pairs.
inline std::vector<Entity> pickedCompetition(const std::vector<Entity>& entities,
                                             std::size_t* n_paired = nullptr)
{
  std::vector<std::int64_t> pair(entities.size());
  std::vector<int> label(entities.size());
  std::vector<double> score(entities.size());
  for (std::size_t i = 0; i < entities.size(); ++i)
  {
    pair[i] = entities[i].pair;
    label[i] = entities[i].label;
    score[i] = entities[i].score;
  }
  const detail::PairWinners pw = detail::pairWinners(pair, label, score,
                                                     "odia::core::pickedCompetition");
  std::vector<Entity> out;
  out.reserve(entities.size());
  for (std::size_t i = 0; i < entities.size(); ++i)
  {
    if (pw.winner[i]) { out.push_back(entities[i]); }
  }
  if (n_paired) { *n_paired = pw.pairs_complete; }
  return out;
}

/// Result of the entrapment check. `valid` is false when the inputs cannot support an
/// estimate (no entrapment sequences in the database, or nothing reported).
struct EntrapmentEstimate
{
  std::size_t n_reported = 0;      ///< discoveries at the threshold (targets + entrapments)
  std::size_t n_entrapment = 0;    ///< of those, how many were entrapment sequences
  std::size_t db_target = 0;       ///< real target sequences searched
  std::size_t db_entrapment = 0;   ///< entrapment sequences searched
  double ratio = 0.0;              ///< db_entrapment / db_target
  double fdp = 1.0;                ///< estimated false discovery PROPORTION
  bool valid = false;
};

/// Combined entrapment estimator (Wen/Freestone/Keich/Noble, Nat Methods 2025):
///
///     FDP = (1 + 1/r) * n_entrapment / n_reported ,   r = db_entrapment / db_target
///
/// Entrapment sequences cannot be in the sample, so each entrapment discovery is false; they
/// sample the false-discovery process at rate r relative to the real targets, so n_entrapment/r
/// further false discoveries hide among the reported targets. A measurement of the reported
/// FDR, not a replacement for it. NOT the paired estimator, which needs an explicit
/// target<->entrapment pairing.
inline EntrapmentEstimate entrapmentFdp(std::size_t n_reported, std::size_t n_entrapment,
                                        std::size_t db_target, std::size_t db_entrapment)
{
  EntrapmentEstimate e;
  e.n_reported = n_reported;
  e.n_entrapment = n_entrapment;
  e.db_target = db_target;
  e.db_entrapment = db_entrapment;
  if (db_target == 0 || db_entrapment == 0 || n_reported == 0) { return e; }
  e.ratio = static_cast<double>(db_entrapment) / static_cast<double>(db_target);
  e.fdp = std::min(1.0, (1.0 + 1.0 / e.ratio) * static_cast<double>(n_entrapment) /
                          static_cast<double>(n_reported));
  e.valid = true;
  return e;
}

/// Targets at or below `q` (decoys never count as identifications).
inline std::size_t countAtQ(const std::vector<Entity>& entities, double q)
{
  return static_cast<std::size_t>(
    std::count_if(entities.begin(), entities.end(),
                  [&](const Entity& e) { return e.label == 1 && e.qvalue <= q; }));
}

} // namespace odia::core

#endif // ODIA_CORE_FDR_H
