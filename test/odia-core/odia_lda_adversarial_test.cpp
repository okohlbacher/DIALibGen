// Copyright (c) 2026, Oliver Kohlbacher and the OpenDIAlyzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
//
// odia_lda_adversarial_test.cpp -- the ways a semi-supervised target-decoy model cheats.
// Ported from OpenDIAlyzer src/odia_lda_adversarial_test.cpp (see src/odia-core/MANIFEST.json).
//
// Complements the oracle (odia_lda_test.cpp). Each property is a separate block; main()
// returns non-zero on any failure.

#include <odia-core/odia_lda.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <random>
#include <utility>
#include <vector>

namespace
{
struct Data
{
  std::vector<double> f;   // row-major, M columns
  int M = 0;
  std::vector<int> lab;
  std::vector<std::int64_t> g;
  std::vector<std::int64_t> pair;   // target i and decoy i form pair i
  std::vector<char> is_true;
  int n_true = 0;
};

// mu == 0 => pure null (no signal at all). Otherwise true_frac of target precursors carry
// signal (mu) on the informative dims, in exactly one of their peak groups.
Data make(unsigned seed, int M, int ntg, int ndec, int gpp, double true_frac, double mu,
          const std::vector<int>& info)
{
  std::mt19937 rng(seed);
  std::normal_distribution<double> nz(0.0, 1.0);
  std::uniform_real_distribution<double> u(0.0, 1.0);
  Data d;
  d.M = M;
  std::int64_t gid = 0;
  auto add = [&](bool tgt, bool tr) {
    const std::int64_t id = gid++;
    d.is_true.push_back(tr ? 1 : 0);
    for (int gg = 0; gg < gpp; ++gg)
    {
      std::vector<double> x(M);
      for (int j = 0; j < M; ++j) { x[j] = nz(rng); }
      if (tr && gg == 0) { for (int k : info) { x[k] += mu; } }
      d.f.insert(d.f.end(), x.begin(), x.end());
      d.lab.push_back(tgt ? 1 : 0);
      d.g.push_back(id);
      d.pair.push_back(tgt ? id : id - ntg);
    }
  };
  for (int i = 0; i < ntg; ++i) { const bool tr = (mu > 0.0) && (u(rng) < true_frac); d.n_true += tr; add(true, tr); }
  for (int i = 0; i < ndec; ++i) { add(false, false); }
  return d;
}

odia::core::LdaResult score(const Data& d)
{
  return odia::core::scoreSemiSupervisedLDA(d.f, static_cast<std::size_t>(d.M), d.lab, d.g, d.pair);
}

// best d-score row per precursor group -> its q, label, d
void bestPerGroup(const Data& d, const odia::core::LdaResult& s,
                  std::vector<double>& bq, std::vector<int>& bl, std::vector<double>& bd)
{
  std::int64_t ng = 0;
  for (auto g : d.g) { ng = std::max(ng, g + 1); }
  bq.assign(static_cast<std::size_t>(ng), 2.0);
  bl.assign(static_cast<std::size_t>(ng), -1);
  bd.assign(static_cast<std::size_t>(ng), -1e300);
  for (std::size_t i = 0; i < d.lab.size(); ++i)
  {
    const std::size_t g = static_cast<std::size_t>(d.g[i]);
    if (s.dscore[i] > bd[g]) { bd[g] = s.dscore[i]; bq[g] = s.qvalue[i]; bl[g] = d.lab[i]; }
  }
}
} // namespace

int main()
{
  int fails = 0;

  // P1 NULL CALIBRATION / no leakage: pure noise must yield almost no IDs at q < 0.01.
  {
    Data d = make(7, 8, 2000, 2000, 4, 0.0, 0.0, {});
    const auto s = score(d);
    std::vector<double> bq, bd; std::vector<int> bl; bestPerGroup(d, s, bq, bl, bd);
    int rt = 0, rd = 0;
    for (std::size_t id = 0; id < bq.size(); ++id) { if (bq[id] < 0.01) { (bl[id] == 1 ? rt : rd)++; } }
    const double frac = rt / 2000.0;
    std::fprintf(stderr, "[P1 null] targets@q<0.01=%d (%.3f of 2000), decoys=%d\n", rt, frac, rd);
    if (frac > 0.01) { std::fprintf(stderr, "  FAIL: leaks %d IDs on pure noise\n", rt); ++fails; }
  }

  // P2 DETERMINISM: same inputs + seed -> identical vectors.
  {
    Data d = make(11, 8, 1500, 1500, 3, 0.6, 2.2, {0, 2, 5});
    const auto a = score(d);
    const auto b = score(d);
    double md = 0.0;
    for (std::size_t i = 0; i < a.dscore.size(); ++i)
    {
      md = std::max(md, std::fabs(a.dscore[i] - b.dscore[i]));
      md = std::max(md, std::fabs(a.qvalue[i] - b.qvalue[i]));
    }
    std::fprintf(stderr, "[P2 determinism] max|delta|=%.2e\n", md);
    if (md > 1e-9) { std::fprintf(stderr, "  FAIL: nondeterministic\n"); ++fails; }
  }

  // P3 Q-MONOTONICITY: best-per-group sorted by d-score desc -> q non-decreasing.
  {
    Data d = make(13, 8, 1500, 1500, 3, 0.6, 2.2, {0, 2, 5});
    const auto s = score(d);
    std::vector<double> bq, bd; std::vector<int> bl; bestPerGroup(d, s, bq, bl, bd);
    std::vector<std::pair<double, double>> dq;
    for (std::size_t id = 0; id < bq.size(); ++id) { if (bl[id] >= 0) { dq.emplace_back(bd[id], bq[id]); } }
    std::sort(dq.begin(), dq.end(), [](const auto& x, const auto& y) { return x.first > y.first; });
    double running_max = -1.0; int viol = 0;
    for (const auto& p : dq) { if (p.second < running_max - 1e-9) { ++viol; } running_max = std::max(running_max, p.second); }
    std::fprintf(stderr, "[P3 monotonicity] violations=%d over %zu groups\n", viol, dq.size());
    if (viol > 0) { std::fprintf(stderr, "  FAIL: non-monotone q vs d-score\n"); ++fails; }
  }

  // P4 FDR CALIBRATION under known signal: controlled FDR + non-trivial recall.
  {
    Data d = make(17, 8, 2000, 2000, 4, 0.6, 2.5, {0, 2, 5});
    const auto s = score(d);
    std::vector<double> bq, bd; std::vector<int> bl; bestPerGroup(d, s, bq, bl, bd);
    int rt = 0, rd = 0, rtrue = 0;
    for (std::size_t id = 0; id < bq.size(); ++id)
    {
      if (bq[id] < 0.01) { if (bl[id] == 1) { ++rt; if (d.is_true[id]) { ++rtrue; } } else { ++rd; } }
    }
    const double fdr = rt ? static_cast<double>(rd) / rt : 1.0;
    const double recall = static_cast<double>(rtrue) / d.n_true;
    std::fprintf(stderr, "[P4 calibration] IDs=%d emp_FDR=%.3f recall=%.2f\n", rt, fdr, recall);
    if (fdr > 0.02)   { std::fprintf(stderr, "  FAIL: emp_FDR %.3f > 0.02\n", fdr); ++fails; }
    if (recall < 0.3) { std::fprintf(stderr, "  FAIL: recall %.2f < 0.3\n", recall); ++fails; }
  }

  if (fails) { std::fprintf(stderr, "adversarial test FAILED (%d property/properties)\n", fails); return 1; }
  std::fprintf(stderr, "adversarial test OK\n");
  return 0;
}
