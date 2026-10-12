// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause
//
// Pair-aware folds: both members of a target-decoy pair are held out together, and the fold of
// a pair depends on nothing but the pair id and the seed.
//
// ODIA shuffled target groups and decoy groups independently, so a target and its own decoy
// usually sat in different folds and the partner of a held-out precursor trained the model that
// scored it. That assignment fails F1 below for about two thirds of the pairs.

#include "synthetic.h"

#include <odia-core/odia_lda.h>

#include <algorithm>
#include <cstdio>
#include <map>
#include <stdexcept>
#include <vector>

using synth::check;

namespace
{
odia::core::LdaResult run(const synth::Data& d, const std::vector<int>& labels, unsigned seed = 42)
{
  odia::core::LdaParams p;
  p.seed = seed;
  return odia::core::scoreSemiSupervisedLDA(d.x, d.m, labels, d.group, d.pair, p);
}
} // namespace

int main()
{
  const synth::Data d = synth::make({});
  const auto r = run(d, d.label);

  // F1: every row of a pair -- target and decoy, all peak groups -- has ONE fold, and it is the
  // hash of the pair id.
  {
    std::map<std::int64_t, int> fold_of_pair;
    std::size_t split = 0, mismatch = 0;
    for (std::size_t i = 0; i < d.rows(); ++i)
    {
      const auto inserted = fold_of_pair.emplace(d.pair[i], r.fold[i]);
      if (!inserted.second && inserted.first->second != r.fold[i]) { ++split; }
      if (r.fold[i] != odia::core::foldOfPair(d.pair[i], 42, r.n_folds)) { ++mismatch; }
    }
    std::fprintf(stderr, "[F1] %zu rows in a different fold from their pair, %zu off the hash\n",
                 split, mismatch);
    check(r.n_folds == 3, "F1 three folds by default");
    check(split == 0, "F1 target and decoy of a pair share one fold");
    check(mismatch == 0, "F1 fold = foldOfPair(pair, seed, folds)");

    // Balanced in expectation: 2000 pairs over 3 folds.
    std::vector<int> size(3, 0);
    for (const auto& kv : fold_of_pair) { ++size[static_cast<std::size_t>(kv.second)]; }
    std::fprintf(stderr, "[F1] pairs per fold %d / %d / %d\n", size[0], size[1], size[2]);
    for (const int s : size) { check(s > 560 && s < 780, "F1 folds are roughly balanced"); }
  }

  // F2: label-blind. Swapping every label changes the model, never the folds.
  {
    std::vector<int> swapped(d.label.size());
    for (std::size_t i = 0; i < swapped.size(); ++i) { swapped[i] = 1 - d.label[i]; }
    const auto s = run(d, swapped);
    check(s.fold == r.fold, "F2 the fold assignment does not see labels");
  }

  // F3: subset-invariant. Drop every other pair (a smaller candidate set, another chunking):
  // the pairs that remain keep their folds.
  {
    synth::Data sub;
    sub.m = d.m;
    std::vector<int> kept_fold;
    for (std::size_t i = 0; i < d.rows(); ++i)
    {
      if (d.pair[i] % 2 != 0) { continue; }
      sub.x.insert(sub.x.end(), d.x.begin() + static_cast<std::ptrdiff_t>(i * d.m),
                   d.x.begin() + static_cast<std::ptrdiff_t>((i + 1) * d.m));
      sub.label.push_back(d.label[i]);
      sub.group.push_back(d.group[i]);
      sub.pair.push_back(d.pair[i]);
      kept_fold.push_back(r.fold[i]);
    }
    const auto s = run(sub, sub.label);
    check(s.fold == kept_fold, "F3 a pair keeps its fold when other pairs leave");
  }

  // F4: the seed is used, and the hash is the same number on every platform (pinned).
  {
    const auto s = run(d, d.label, 43);
    check(s.fold != r.fold, "F4 another seed gives another assignment");
    const int pinned[] = {odia::core::foldOfPair(0, 42, 3), odia::core::foldOfPair(1, 42, 3),
                          odia::core::foldOfPair(-7, 42, 3), odia::core::foldOfPair(123456789, 42, 3),
                          odia::core::foldOfPair(5, 7, 10)};
    std::fprintf(stderr, "[F4] pinned folds %d %d %d %d %d\n", pinned[0], pinned[1], pinned[2],
                 pinned[3], pinned[4]);
    check(pinned[0] == 2 && pinned[1] == 0 && pinned[2] == 2 && pinned[3] == 1 && pinned[4] == 0,
          "F4 foldOfPair values are pinned");
  }

  // F5: a group whose rows claim different pairs is refused, not silently split.
  {
    auto bad = d.pair;
    bad[1] = bad[0] + 1;   // second peak group of pair 0's target claims pair 1
    bool threw = false;
    try { odia::core::scoreSemiSupervisedLDA(d.x, d.m, d.label, d.group, bad); }
    catch (const std::invalid_argument&) { threw = true; }
    check(threw, "F5 a group with two pair ids throws");
  }

  // F6: the scoring still works on pair folds.
  check(r.n_iterations_trained > 0, "F6 the discriminant trains");

  return synth::finish("odia_pair_folds_test");
}
