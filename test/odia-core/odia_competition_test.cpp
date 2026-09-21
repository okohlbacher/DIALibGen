// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause
//
// Concatenated target-decoy competition at the precursor level, with STRUCTURAL pairs.
//
//   C1 ties go to the decoy; a loser never gets a q below 1.
//   C2 q-values on the winners, hand-computed, for both estimators.
//   C3 pairing comes from pair ids alone; a negative id or a missing partner means no competitor.
//   C4 a pair id carrying two targets (or two decoys) is refused.
//   C5 a winner list without decoys gets q = 1/T, not 0 (ODIA's Ratio estimator gave 0).
//   C6 the pooled estimator is kept, and differs from the concatenated one when decoys are lit
//      up by their own targets.

#include "synthetic.h"

#include <odia-core/odia_fdr.h>

#include <cmath>
#include <cstdio>
#include <random>
#include <stdexcept>
#include <vector>

using synth::check;

namespace
{
bool near(double a, double b) { return std::abs(a - b) < 1e-12; }
} // namespace

int main()
{
  using odia::core::QEstimator;

  // ---- C1: ties -------------------------------------------------------------------------
  {
    const auto c = odia::core::concatenatedCompetition({7, 7}, {1, 0}, {2.5, 2.5});
    check(!c.winner[0] && c.winner[1], "C1 an exact tie goes to the decoy");
    check(c.qvalue[0] == 1.0 && c.pep[0] == 1.0, "C1 the losing target has q = PEP = 1");
    check(c.pairs_complete == 1 && c.target_winners == 0 && c.decoy_winners == 1, "C1 counts");
    const auto d = odia::core::concatenatedCompetition({7, 7}, {1, 0}, {2.5000001, 2.5});
    check(d.winner[0] && !d.winner[1], "C1 a strictly better target wins");
  }

  // ---- C2: hand-computed q ---------------------------------------------------------------
  // Winners, best first: T x10 (20..11), D (10), T x5 (9..5), D x4 (4..1).
  // N_tar_win = 15, N_dec_win = 5.
  //   Ratio: top 10  -> ((0+1)/5) / (10/15) = 0.3   (monotonised: every one of them 0.3)
  //          next 5  -> ((1+1)/5) / (15/15) = 0.4   (the decoy between them also 0.4)
  //   Count: top 10  -> (0+1)/10 = 0.1; next 5 -> (1+1)/15 = 0.1333...
  {
    std::vector<std::int64_t> pair;
    std::vector<int> label;
    std::vector<double> score;
    std::int64_t p = 0;
    auto add = [&](int winner_label, double win, double lose) {
      pair.push_back(p); label.push_back(1); score.push_back(winner_label == 1 ? win : lose);
      pair.push_back(p); label.push_back(0); score.push_back(winner_label == 1 ? lose : win);
      ++p;
    };
    for (int i = 0; i < 10; ++i) { add(1, 20.0 - i, -100.0); }
    add(0, 10.0, 0.0);
    for (int i = 0; i < 5; ++i) { add(1, 9.0 - i, -100.0); }
    for (int i = 0; i < 4; ++i) { add(0, 4.0 - i, -50.0); }

    const auto r = odia::core::concatenatedCompetition(pair, label, score, QEstimator::Ratio);
    const auto c = odia::core::concatenatedCompetition(pair, label, score, QEstimator::Count);
    check(r.target_winners == 15 && r.decoy_winners == 5 && r.pairs_complete == 20, "C2 winner counts");
    bool ok_r = true, ok_c = true;
    for (std::size_t i = 0; i < pair.size(); ++i)
    {
      const std::int64_t pid = pair[i];
      if (!r.winner[i]) { ok_r = ok_r && r.qvalue[i] == 1.0; ok_c = ok_c && c.qvalue[i] == 1.0; continue; }
      if (pid < 10)       { ok_r = ok_r && near(r.qvalue[i], 0.3); ok_c = ok_c && near(c.qvalue[i], 0.1); }
      else if (pid <= 15) { ok_r = ok_r && near(r.qvalue[i], 0.4); ok_c = ok_c && near(c.qvalue[i], 2.0 / 15.0); }
    }
    check(ok_r, "C2 Ratio: (D_win+1)/N_dec_win / (T_win/N_tar_win), monotonised");
    check(ok_c, "C2 Count: (D_win+1)/T_win, monotonised");
  }

  // ---- C3: structure, not names ------------------------------------------------------------
  {
    // items: target of pair 3 (no decoy present), decoy of pair 4 (no target present),
    // two items with pair -1 (never compete), a complete pair 5.
    const auto c = odia::core::concatenatedCompetition({3, 4, -1, -1, 5, 5}, {1, 0, 1, 0, 1, 0},
                                                       {1.0, 2.0, 0.5, 9.0, 3.0, 1.0});
    check(c.winner == std::vector<char>({1, 1, 1, 1, 1, 0}), "C3 only a complete pair competes");
    check(c.pairs_complete == 1 && c.targets_unpaired == 2 && c.decoys_unpaired == 2, "C3 unpaired counts");
  }

  // ---- C4: a malformed pair is refused ----------------------------------------------------------
  {
    bool threw = false;
    try { odia::core::concatenatedCompetition({1, 1, 1}, {1, 0, 1}, {1.0, 2.0, 3.0}); }
    catch (const std::invalid_argument&) { threw = true; }
    check(threw, "C4 two targets in one pair throw");
    threw = false;
    try { odia::core::concatenatedCompetition({1, 1}, {0, 0}, {1.0, 2.0}); }
    catch (const std::invalid_argument&) { threw = true; }
    check(threw, "C4 two decoys in one pair throw");
  }

  // ---- C5: no decoy winners ---------------------------------------------------------------------
  {
    const auto c = odia::core::concatenatedCompetition({0, 0, 1, 1, 2, 2, 3, 3}, {1, 0, 1, 0, 1, 0, 1, 0},
                                                       {5, 1, 6, 1, 7, 1, 8, 1});
    bool ok = true;
    for (std::size_t i = 0; i < 8; i += 2) { ok = ok && near(c.qvalue[i], 0.25); }
    check(c.decoy_winners == 0 && ok, "C5 all-target winners: q = (0+1)/T = 0.25, not 0");
  }

  // ---- C6: pooled (diagnostic) vs concatenated -----------------------------------------------------
  // Decoys partly lit by their own target: a true target's decoy scores half the target's score
  // plus noise. Pooled, those decoys fill the null tail and push every q up; competing, each
  // loses to the target that lit it.
  {
    std::mt19937 rng(5);
    std::normal_distribution<double> nd(0.0, 1.0);
    std::vector<std::int64_t> pair;
    std::vector<int> label;
    std::vector<double> score;
    for (int p = 0; p < 4000; ++p)
    {
      const bool hit = p < 1500;
      const double t = (hit ? 5.0 : 0.0) + nd(rng);
      const double d = hit ? 0.5 * t + 0.5 * nd(rng) : nd(rng);
      pair.push_back(p); label.push_back(1); score.push_back(t);
      pair.push_back(p); label.push_back(0); score.push_back(d);
    }
    const auto c = odia::core::concatenatedCompetition(pair, label, score);
    const auto pooled = odia::core::pooledQValues(label, score);
    std::size_t ids_c = 0, ids_p = 0, false_c = 0;
    for (std::size_t i = 0; i < label.size(); ++i)
    {
      if (label[i] != 1) { continue; }
      if (c.winner[i] && c.qvalue[i] <= 0.01) { ++ids_c; if (pair[i] >= 1500) { ++false_c; } }
      if (pooled[i] <= 0.01) { ++ids_p; }
    }
    std::fprintf(stderr, "[C6] targets at q<=0.01: concatenated %zu (%zu null), pooled %zu\n",
                 ids_c, false_c, ids_p);
    check(ids_c > ids_p, "C6 competition removes decoys lit by their own targets from the null");
    check(ids_c > 0 && static_cast<double>(false_c) <= 0.02 * static_cast<double>(ids_c),
          "C6 and stays honest: at most 2 % of its IDs are null targets");
  }

  return synth::finish("odia_competition_test");
}
