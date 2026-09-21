// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause
//
// Label-blind fold rescaling: each fold's held-out scores are rescaled by the median and MAD of
// the best score of ALL its groups, never by statistics of its decoys.
//
// The decisive check is R3. A fold's model is trained on the OTHER folds only, so relabelling
// the pairs inside fold 0 (target <-> decoy) cannot change fold 0's model or its raw held-out
// scores. A label-blind rescale then leaves fold 0's final scores bit-identical. ODIA's rescale
// to the fold's decoy mean and sd changes them, because the set of decoys changed.

#include "synthetic.h"

#include <odia-core/odia_lda.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

using synth::check;

namespace
{
// Best score per group, grouped by fold.
std::vector<std::vector<double>> bestByFold(const synth::Data& d, const odia::core::LdaResult& r)
{
  std::vector<double> best(d.groups(), -1e300);
  std::vector<int> fold(d.groups(), -1);
  for (std::size_t i = 0; i < d.rows(); ++i)
  {
    const std::size_t g = static_cast<std::size_t>(d.group[i]);
    best[g] = std::max(best[g], r.dscore[i]);
    fold[g] = r.fold[i];
  }
  std::vector<std::vector<double>> out(static_cast<std::size_t>(r.n_folds));
  for (std::size_t g = 0; g < best.size(); ++g) { out[static_cast<std::size_t>(fold[g])].push_back(best[g]); }
  return out;
}
} // namespace

int main()
{
  // R1: the robust scale itself.
  {
    const auto a = odia::core::robustScale({1.0, 2.0, 3.0, 4.0, 100.0});
    check(a.valid && a.location == 3.0 && std::abs(a.scale - 1.4826) < 1e-12,
          "R1 median 3 and 1.4826 x MAD 1 for {1,2,3,4,100}");
    const auto b = odia::core::robustScale({4.0, 1.0, 3.0, 2.0});
    check(b.valid && b.location == 2.5 && std::abs(b.scale - 1.4826) < 1e-12,
          "R1 even count: the median averages the middle pair");
    check(!odia::core::robustScale({5.0, 5.0, 5.0}).valid, "R1 no spread -> not valid");
    check(!odia::core::robustScale({5.0}).valid, "R1 one value -> not valid");
  }

  const synth::Data d = synth::make({});
  const auto r = odia::core::scoreSemiSupervisedLDA(d.x, d.m, d.label, d.group, d.pair);

  // R2: after rescaling, each fold's best-per-group scores over ALL groups have median 0 and
  // normal-consistent MAD 1. (A decoy-only standardisation leaves the median well above 0 here:
  // 60 % of the targets carry signal.)
  {
    check(r.folds_unscaled == 0, "R2 every fold rescaled");
    const auto folds = bestByFold(d, r);
    for (std::size_t f = 0; f < folds.size(); ++f)
    {
      const auto rs = odia::core::robustScale(folds[f]);
      std::fprintf(stderr, "[R2] fold %zu: %zu groups, median %.3e, scale %.12f (was %.3f / %.3f)\n",
                   f, folds[f].size(), rs.location, rs.scale, r.fold_location[f], r.fold_scale[f]);
      check(std::abs(rs.location) < 1e-9, "R2 median of all held-out best scores is 0");
      check(std::abs(rs.scale - 1.0) < 1e-9, "R2 MAD scale of all held-out best scores is 1");
    }
  }

  // R3: label-blind. Swap target and decoy for every pair in fold 0; fold 0's scores must not move.
  {
    std::vector<int> flipped = d.label;
    std::size_t n_flipped = 0;
    for (std::size_t i = 0; i < d.rows(); ++i)
    {
      if (r.fold[i] == 0) { flipped[i] = 1 - flipped[i]; ++n_flipped; }
    }
    const auto s = odia::core::scoreSemiSupervisedLDA(d.x, d.m, flipped, d.group, d.pair);
    std::size_t moved = 0, others_moved = 0;
    for (std::size_t i = 0; i < d.rows(); ++i)
    {
      if (r.fold[i] == 0) { if (s.dscore[i] != r.dscore[i]) { ++moved; } }
      else if (s.dscore[i] != r.dscore[i]) { ++others_moved; }
    }
    std::fprintf(stderr, "[R3] relabelled %zu rows of fold 0: %zu of its scores moved; %zu rows of "
                         "other folds moved (their training data changed)\n",
                 n_flipped, moved, others_moved);
    check(n_flipped > 0, "R3 fold 0 is not empty");
    check(moved == 0, "R3 relabelling a held-out fold does not change its scores");
    check(others_moved > 0, "R3 sanity: the other folds did see the relabelled training data");
  }

  // R4: off means off.
  {
    odia::core::LdaParams p;
    p.normalize_folds = false;
    const auto s = odia::core::scoreSemiSupervisedLDA(d.x, d.m, d.label, d.group, d.pair, p);
    check(s.folds_unscaled == 0 && s.fold_scale == std::vector<double>(3, 1.0) &&
          s.fold_location == std::vector<double>(3, 0.0),
          "R4 normalize_folds=false reports identity transforms");
  }

  return synth::finish("odia_fold_rescale_test");
}
