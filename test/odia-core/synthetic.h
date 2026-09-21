// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause
//
// Synthetic score tables for the odia-core tests: target-decoy PAIRS of precursors, each with a
// few candidate peak groups. A fraction of the targets are true hits whose first peak group
// carries a shift on the informative sub-scores; every other row is N(0, 1) noise. Decoys are
// independent noise, so every decoy and every non-hit target is null.

#ifndef ODIA_CORE_TEST_SYNTHETIC_H
#define ODIA_CORE_TEST_SYNTHETIC_H

#include <cstdint>
#include <cstdio>
#include <random>
#include <vector>

namespace synth
{

struct Spec
{
  unsigned seed = 1;
  int pairs = 2000;
  int m = 8;                       // sub-scores
  int peaks = 4;                   // candidate peak groups per precursor
  double true_frac = 0.6;          // share of targets that are true hits (0 = pure null)
  double mu = 2.5;                 // shift on the informative sub-scores of a hit's peak group
  std::vector<int> informative = {0, 2, 5};
};

struct Data
{
  std::size_t m = 0;
  std::vector<double> x;               // row-major, rows x m
  std::vector<int> label;              // 1 = target, 0 = decoy
  std::vector<std::int64_t> group;     // target of pair p: 2p, its decoy: 2p + 1
  std::vector<std::int64_t> pair;      // p
  std::vector<std::int64_t> feature;   // unique within a group
  std::vector<char> group_true;        // by group id: a true hit
  int n_true = 0;

  std::size_t rows() const { return label.size(); }
  std::size_t groups() const { return group_true.size(); }
};

inline Data make(const Spec& spec)
{
  std::mt19937 rng(spec.seed);
  std::normal_distribution<double> noise(0.0, 1.0);
  std::uniform_real_distribution<double> u(0.0, 1.0);
  Data d;
  d.m = static_cast<std::size_t>(spec.m);
  d.group_true.assign(static_cast<std::size_t>(2 * spec.pairs), 0);
  for (int p = 0; p < spec.pairs; ++p)
  {
    const bool hit = spec.true_frac > 0.0 && u(rng) < spec.true_frac;
    d.n_true += hit ? 1 : 0;
    for (int member = 0; member < 2; ++member)
    {
      const bool target = member == 0;
      const std::int64_t g = 2 * static_cast<std::int64_t>(p) + member;
      d.group_true[static_cast<std::size_t>(g)] = (target && hit) ? 1 : 0;
      for (int k = 0; k < spec.peaks; ++k)
      {
        for (int j = 0; j < spec.m; ++j)
        {
          double v = noise(rng);
          if (target && hit && k == 0)
          {
            for (const int inf : spec.informative) { if (inf == j) { v += spec.mu; } }
          }
          d.x.push_back(v);
        }
        d.label.push_back(target ? 1 : 0);
        d.group.push_back(g);
        d.pair.push_back(p);
        d.feature.push_back(g * 1000 + 7 * k + 3);
      }
    }
  }
  return d;
}

inline int failures = 0;
inline void check(bool ok, const char* what)
{
  if (!ok) { std::fprintf(stderr, "FAIL: %s\n", what); ++failures; }
}
inline int finish(const char* name)
{
  if (failures) { std::fprintf(stderr, "%s FAILED (%d)\n", name, failures); return 1; }
  std::fprintf(stderr, "%s OK\n", name);
  return 0;
}

} // namespace synth

#endif // ODIA_CORE_TEST_SYNTHETIC_H
