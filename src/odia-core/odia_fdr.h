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
// One estimator, three levels. The precursor level ranks one best score per precursor; the
// peptide and protein levels roll those scores up to the entity and rank again:
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
#include <string>
#include <unordered_map>
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

/// Target-decoy q-values, p-values and PEP for a list with one score per item.
///
/// Equal scores are one threshold, so the result does not depend on the order of ties. The list
/// is sorted in place (score descending, then group_index ascending).
///
///   q   : FDR(t) = pi0 * ((D(t) + 1) / N_dec) / (T(t) / N_tar), monotonised from the bottom.
///         The +1 is Kall's finite-sample correction on the decoy count.
///   p   : (decoys at or above + 1) / (N_dec + 1).
///   PEP : the local analogue of the q estimator in a sliding window of the ranked list, made
///         non-increasing in score by an isotonic sweep.
///
/// pi0 = 1 unless use_pi0 (Storey, lambda = 0.5). pi0 = 1 is the honest default.
inline void assignQValues(std::vector<RankedGroup>& ranked, bool use_pi0)
{
  std::sort(ranked.begin(), ranked.end(), [](const RankedGroup& a, const RankedGroup& b) {
    if (a.score != b.score) { return a.score > b.score; }
    return a.group_index < b.group_index;
  });

  std::size_t Ntar = 0, Ndec = 0;
  for (const auto& r : ranked) { if (r.label == 1) { ++Ntar; } else { ++Ndec; } }

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
    double fdr;
    if (targets == 0) { fdr = std::numeric_limits<double>::infinity(); }
    else
    {
      const double dr = (Ndec > 0) ? (static_cast<double>(decoys) + 1.0) / static_cast<double>(Ndec) : 0.0;
      const double tr = static_cast<double>(targets) / static_cast<double>(Ntar);
      fdr = std::min(1.0, pi0 * dr / tr);
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
  // is pi0 * d * (Ntar / Ndec), so PEP ~ that over t. Sliding counts keep the pass O(n); the
  // window is ~0.5 % of the list and at least 101 wide. The raw local ratio is noisy, so it is
  // projected onto non-increasing-in-score by a running max from the low-score end.
  {
    const std::size_t n = ranked.size();
    const std::size_t half = std::max<std::size_t>(50, n / 200);
    const double scale = (Ndec > 0) ? (static_cast<double>(Ntar) / static_cast<double>(Ndec)) : 0.0;
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
    double running_max = 0.0;
    for (std::size_t end = n; end > 0; --end)
    {
      running_max = std::max(running_max, ranked[end - 1].pep);
      ranked[end - 1].pep = running_max;
    }
  }
}

/// One item at whatever level is being controlled: a peptide, or a protein (group).
struct Entity
{
  std::string id;            ///< modified sequence, or protein accession / group key
  int label = 1;             ///< 1 = target, 0 = decoy
  double score = 0.0;        ///< best member score (see rollUp)
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
inline std::vector<Entity> rollUp(const std::vector<std::string>& member_key,
                                  const std::vector<double>& member_score,
                                  const std::vector<int>& member_label)
{
  std::vector<Entity> out;
  const std::size_t n = member_key.size();
  if (member_score.size() != n || member_label.size() != n) { return out; }

  std::unordered_map<std::string, std::size_t> index;
  index.reserve(n);
  for (std::size_t i = 0; i < n; ++i)
  {
    if (member_key[i].empty()) { continue; }
    auto inserted = index.emplace(member_key[i], out.size());
    if (inserted.second)
    {
      Entity e;
      e.id = member_key[i];
      e.label = member_label[i] == 1 ? 1 : 0;
      e.score = member_score[i];
      out.push_back(std::move(e));
    }
    else if (member_score[i] > out[inserted.first->second].score)
    {
      out[inserted.first->second].score = member_score[i];
    }
  }
  return out;
}

/// Target-decoy q-value / p-value / PEP over an entity set, identical estimator to the
/// precursor level. Entities are updated in place; their order is preserved.
inline void assignQValues(std::vector<Entity>& entities, bool use_pi0)
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
  assignQValues(ranked, use_pi0);
  for (const auto& r : ranked)
  {
    entities[r.group_index].qvalue = r.qvalue;
    entities[r.group_index].pvalue = r.pvalue;
    entities[r.group_index].pep = r.pep;
  }
}

/// True if `id` is `decoy_tag` + something; writes the partner (target) id to `target_id`.
inline bool decoyPartnerId(const std::string& id, const std::string& decoy_tag,
                           std::string& target_id)
{
  if (decoy_tag.empty() || id.size() <= decoy_tag.size()) { return false; }
  if (id.compare(0, decoy_tag.size(), decoy_tag) != 0) { return false; }
  target_id = id.substr(decoy_tag.size());
  return true;
}

/// Picked target-decoy competition (Savitski 2015; The 2022).
///
/// Each target/decoy pair is collapsed to ONE entry carrying the better of the two scores and
/// the label of whichever won; entities with no partner pass through unchanged. The caller then
/// runs assignQValues() on the result. Ties go to the DECOY: a tie carries no evidence, and
/// resolving it for the target would bias the estimate downward.
///
/// Assumes the decoy id is `decoy_tag` + the target id.
inline std::vector<Entity> pickedCompetition(const std::vector<Entity>& entities,
                                             const std::string& decoy_tag,
                                             std::size_t* n_paired = nullptr)
{
  std::unordered_map<std::string, std::size_t> target_index;
  target_index.reserve(entities.size());
  for (std::size_t i = 0; i < entities.size(); ++i)
  {
    if (entities[i].label == 1) { target_index.emplace(entities[i].id, i); }
  }

  std::vector<Entity> out;
  out.reserve(entities.size());
  std::vector<char> consumed(entities.size(), 0);
  std::size_t paired = 0;

  for (std::size_t i = 0; i < entities.size(); ++i)
  {
    if (entities[i].label == 1) { continue; }
    std::string partner;
    if (!decoyPartnerId(entities[i].id, decoy_tag, partner)) { continue; }
    const auto it = target_index.find(partner);
    if (it == target_index.end()) { continue; }

    const Entity& decoy = entities[i];
    const Entity& target = entities[it->second];
    out.push_back(target.score > decoy.score ? target : decoy);
    consumed[i] = 1;
    consumed[it->second] = 1;
    ++paired;
  }
  for (std::size_t i = 0; i < entities.size(); ++i)
  {
    if (!consumed[i]) { out.push_back(entities[i]); }
  }
  if (n_paired) { *n_paired = paired; }
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
