// Copyright (c) 2026, Oliver Kohlbacher and the OpenDIAlyzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
//
// odia_lda_test.cpp -- oracle test for the semi-supervised LDA.
// Ported from OpenDIAlyzer src/odia_lda_test.cpp (see src/odia-core/MANIFEST.json).
//
// Synthetic target/decoy peak groups where a known subset of target precursors are true hits
// with elevated multivariate sub-scores; all other rows are null. A correct semi-supervised LDA
// must (a) rank true-hit rows above null rows, (b) recover a good fraction of true-hit
// precursors at q < 0.01, and (c) keep the empirical target-decoy FDR controlled.

#include <odia-core/odia_lda.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <random>
#include <utility>
#include <vector>

int main()
{
  std::mt19937 rng(123);
  std::normal_distribution<double> noise(0.0, 1.0);

  const int M = 8;                 // sub-scores
  const int n_target_prec = 2000;  // target precursors
  const int n_decoy_prec = 2000;   // decoy precursors
  const int groups_per_prec = 4;   // candidate peak groups per precursor
  const double true_frac = 0.6;    // fraction of target precursors that are real hits
  const double mu = 2.5;           // signal on the informative dims
  const int informative[] = {0, 2, 5};

  std::vector<double> feats;       // row-major, M columns
  std::vector<int> labels;
  std::vector<std::int64_t> group;
  std::vector<char> prec_is_true;
  std::int64_t gid = 0;

  auto add_prec = [&](bool is_target, bool is_true) {
    const std::int64_t id = gid++;
    prec_is_true.push_back(is_true ? 1 : 0);
    for (int g = 0; g < groups_per_prec; ++g)
    {
      std::vector<double> x(M);
      for (int j = 0; j < M; ++j) { x[j] = noise(rng); }
      if (is_true && g == 0) { for (int k : informative) { x[k] += mu; } }
      feats.insert(feats.end(), x.begin(), x.end());
      labels.push_back(is_target ? 1 : 0);
      group.push_back(id);
    }
  };

  std::uniform_real_distribution<double> u(0.0, 1.0);
  int n_true = 0;
  for (int i = 0; i < n_target_prec; ++i) { bool t = u(rng) < true_frac; n_true += t; add_prec(true, t); }
  for (int i = 0; i < n_decoy_prec; ++i)  { add_prec(false, false); }
  const std::size_t rows = labels.size();

  odia::core::LdaParams p;  // defaults
  const odia::core::LdaResult s = odia::core::scoreSemiSupervisedLDA(feats, M, labels, group, p);

  // ---- row order must not change the result ----
  // Fold assignment once used a group's first-occurrence POSITION, and row order is whatever the
  // parallel extraction produced, so identical input gave different ID counts across runs.
  {
    std::vector<std::size_t> perm(rows);
    for (std::size_t i = 0; i < perm.size(); ++i) { perm[i] = i; }
    std::mt19937 prng(999);
    std::shuffle(perm.begin(), perm.end(), prng);

    std::vector<double> f2(feats.size());
    std::vector<int> l2(rows);
    std::vector<std::int64_t> g2(rows);
    for (std::size_t i = 0; i < perm.size(); ++i)
    {
      std::copy(feats.begin() + static_cast<std::ptrdiff_t>(perm[i] * M),
                feats.begin() + static_cast<std::ptrdiff_t>((perm[i] + 1) * M),
                f2.begin() + static_cast<std::ptrdiff_t>(i * M));
      l2[i] = labels[perm[i]];
      g2[i] = group[perm[i]];
    }
    const odia::core::LdaResult s2 = odia::core::scoreSemiSupervisedLDA(f2, M, l2, g2, p);
    double worst = 0.0;
    for (std::size_t i = 0; i < perm.size(); ++i)
    {
      worst = std::max(worst, std::fabs(s2.dscore[i] - s.dscore[perm[i]]));
    }
    std::fprintf(stderr, "permuted-row d-score max |delta| = %.3e\n", worst);
    if (worst > 1e-9)
    {
      std::fprintf(stderr, "FAIL: row order changed the scores (max delta %.3e)\n", worst);
      return 1;
    }
  }

  // ---- basic contract ----
  if (s.dscore.size() != rows || s.qvalue.size() != rows)
  {
    std::fprintf(stderr, "FAIL: output size mismatch\n");
    return 1;
  }
  if (s.n_iterations_trained == 0)
  {
    std::fprintf(stderr, "FAIL: no iteration trained a discriminant\n");
    return 1;
  }

  // ---- separation: true-hit rows should outscore null rows on average ----
  double sum_true = 0, sum_null = 0; int n_t = 0, n_n = 0;
  for (std::size_t i = 0; i < rows; ++i)
  {
    const bool true_hit = (labels[i] == 1 && prec_is_true[static_cast<std::size_t>(group[i])] &&
                           (i % groups_per_prec) == 0);
    if (true_hit) { sum_true += s.dscore[i]; ++n_t; }
    else          { sum_null += s.dscore[i]; ++n_n; }
  }
  const double mean_true = sum_true / n_t, mean_null = sum_null / n_n;
  if (!(mean_true > mean_null))
  {
    std::fprintf(stderr, "FAIL: true-hit d-score %.3f not above null %.3f\n", mean_true, mean_null);
    return 1;
  }

  // ---- q < 0.01: recover a good fraction of true targets, controlled FDR ----
  const std::size_t ng = static_cast<std::size_t>(gid);
  std::vector<double> best_q(ng, 2.0);
  std::vector<int> best_lab(ng, -1);
  std::vector<double> best_d(ng, -1e300);
  for (std::size_t i = 0; i < rows; ++i)
  {
    const std::size_t g = static_cast<std::size_t>(group[i]);
    if (s.dscore[i] > best_d[g]) { best_d[g] = s.dscore[i]; best_q[g] = s.qvalue[i]; best_lab[g] = labels[i]; }
  }
  int rec_target = 0, rec_decoy = 0, rec_true = 0;
  for (std::size_t id = 0; id < ng; ++id)
  {
    if (best_q[id] < 0.01)
    {
      if (best_lab[id] == 1) { ++rec_target; if (prec_is_true[id]) { ++rec_true; } }
      else                   { ++rec_decoy; }
    }
  }
  const double emp_fdr = rec_target > 0 ? static_cast<double>(rec_decoy) / rec_target : 1.0;
  const double recall = static_cast<double>(rec_true) / n_true;
  std::fprintf(stderr, "recovered target-precursors@q<0.01=%d (true=%d/%d, recall=%.2f), decoys=%d, emp_FDR=%.3f\n",
               rec_target, rec_true, n_true, recall, rec_decoy, emp_fdr);

  if (recall < 0.35)   { std::fprintf(stderr, "FAIL: recall %.2f < 0.35\n", recall); return 1; }
  if (emp_fdr > 0.05)  { std::fprintf(stderr, "FAIL: empirical FDR %.3f > 0.05\n", emp_fdr); return 1; }

  // p-value, q-value and PEP must be three DISTINCT statistics:
  //   (a) they are not all identical, and PEP is not a constant;
  //   (b) PEP >= q at the same group -- local FDR cannot be below the average FDR above it;
  //   (c) PEP is non-increasing in the group's best d-score.
  // Checked per GROUP: every row carries its group's statistics, so a runner-up row with a low
  // d-score carries the PEP of its group's best row, and a row-level ordering proves nothing.
  // The ported ODIA check was row-level and passed only because every PEP was 1.
  {
    std::vector<std::size_t> best_row(ng, rows);
    for (std::size_t i = 0; i < rows; ++i)
    {
      const std::size_t g = static_cast<std::size_t>(group[i]);
      if (best_row[g] == rows || s.dscore[i] > s.dscore[best_row[g]]) { best_row[g] = i; }
    }
    std::size_t n_all_equal = 0, n_pep_below_q = 0, n_confident = 0;
    std::vector<std::pair<double, double>> ds;
    ds.reserve(ng);
    for (std::size_t g = 0; g < ng; ++g)
    {
      const std::size_t i = best_row[g];
      if (s.pvalue[i] == s.qvalue[i] && s.qvalue[i] == s.pep[i]) { ++n_all_equal; }
      if (s.pep[i] < s.qvalue[i] - 0.02) { ++n_pep_below_q; }   // both saturate near 1
      if (s.pep[i] < 0.05) { ++n_confident; }
      ds.emplace_back(s.dscore[i], s.pep[i]);
    }
    if (n_all_equal == ng)
    {
      std::fprintf(stderr, "FAIL: p-value, q-value and PEP are identical for all %zu groups\n", ng);
      return 1;
    }
    if (n_confident < static_cast<std::size_t>(rec_target) / 2)
    {
      std::fprintf(stderr, "FAIL: only %zu groups have PEP < 0.05 against %d targets at q < 0.01\n",
                   n_confident, rec_target);
      return 1;
    }
    if (n_pep_below_q > ng / 20)
    {
      std::fprintf(stderr, "FAIL: PEP < q-value for %zu/%zu groups\n", n_pep_below_q, ng);
      return 1;
    }
    std::sort(ds.begin(), ds.end(), [](const auto& a, const auto& b) { return a.first > b.first; });
    for (std::size_t i = 1; i < ds.size(); ++i)
    {
      if (ds[i].first < ds[i - 1].first && ds[i].second < ds[i - 1].second - 1e-9)
      {
        std::fprintf(stderr, "FAIL: PEP increases with d-score at rank %zu (%.4f -> %.4f)\n",
                     i, ds[i - 1].second, ds[i].second);
        return 1;
      }
    }
    std::fprintf(stderr, "p/q/PEP distinct OK (all-equal groups: %zu/%zu, PEP<0.05: %zu)\n",
                 n_all_equal, ng, n_confident);
  }

  std::fprintf(stderr, "odia_lda_test OK\n");
  return 0;
}
