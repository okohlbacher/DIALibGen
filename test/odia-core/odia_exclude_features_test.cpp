// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause
//
// Excluding named sub-scores from a fit.
//
// A discriminant that selects retention-time calibration anchors must not see the RT deviation
// score: anchors chosen by the agreement they are meant to correct bias the correction toward
// the uncorrected state. ODIA had a seed-only mask for its network classifier; here the column
// is removed from the whole fit, and the test proves it by rewriting the excluded column at
// random -- the result must not move by a single bit.
//
//   X1 an excluded column cannot influence the result
//   X2 without the exclusion the same rewrite does change it (the test has teeth)
//   X3 excluded, dropped and unmatched names are reported separately; matching is exact
//   X4 excluding every column is refused

#include "synthetic.h"

#include <odia-core/odia_core.h>

#include <cstdio>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

using synth::check;
using odia::core::ScoreTable;

namespace
{
// Eight sub-scores plus "var_norm_rt_score" in column 8, which carries a strong shift on the
// true hits' signal row -- the kind of score a calibration fit must not select anchors by.
ScoreTable withRtColumn(const synth::Data& d, unsigned rt_seed, double rt_shift)
{
  synth::Rng rng(rt_seed);
  ScoreTable t;
  std::vector<std::string> names;
  for (std::size_t j = 0; j < d.m; ++j) { names.push_back("var_s" + std::to_string(j)); }
  names.push_back("var_norm_rt_score");
  t.setColumns(names);
  for (std::size_t i = 0; i < d.rows(); ++i)
  {
    std::vector<double> row(d.x.begin() + static_cast<std::ptrdiff_t>(i * d.m),
                            d.x.begin() + static_cast<std::ptrdiff_t>((i + 1) * d.m));
    const bool signal_row = d.group_true[static_cast<std::size_t>(d.group[i])] && d.feature[i] % 1000 == 3;
    row.push_back(rng.normal() + (signal_row ? rt_shift : 0.0));
    t.append(d.group[i], d.pair[i], d.label[i] == 0, d.feature[i], row.data());
  }
  return t;
}

bool sameScores(const odia::core::ScoredResult& a, const odia::core::ScoredResult& b)
{
  if (a.dscore != b.dscore) { return false; }
  for (std::size_t g = 0; g < a.groups.size(); ++g)
  {
    if (a.groups[g].qvalue != b.groups[g].qvalue) { return false; }
  }
  return true;
}
} // namespace

int main()
{
  synth::Spec spec;
  spec.mu = 2.2;   // weaker signal on the ordinary sub-scores, so the RT column would dominate
  const synth::Data d = synth::make(spec);
  const ScoreTable a = withRtColumn(d, 1, 3.0);
  const ScoreTable b = withRtColumn(d, 2, -2.0);   // the same table, RT column rewritten

  odia::core::Options excl;
  excl.exclude_features = {"var_norm_rt_score", "var_im_delta_score"};

  // ---- X1 / X2 ------------------------------------------------------------------------------
  {
    const auto ra = odia::core::scoreAndControl(a, excl);
    const auto rb = odia::core::scoreAndControl(b, excl);
    check(sameScores(ra, rb), "X1 rewriting an excluded column changes nothing");

    const auto fa = odia::core::scoreAndControl(a);
    const auto fb = odia::core::scoreAndControl(b);
    check(!sameScores(fa, fb), "X2 without the exclusion the same rewrite changes the scores");
    std::fprintf(stderr, "[X1] IDs with the RT column %zu, excluded %zu\n",
                 fa.diagnostics.targets_at_q, ra.diagnostics.targets_at_q);
    check(fa.diagnostics.targets_at_q > ra.diagnostics.targets_at_q,
          "X2 sanity: the RT column did carry signal the excluded fit went without");
    check(ra.diagnostics.targets_at_q > 0 && ra.diagnostics.iterations_trained > 0,
          "X1 the excluded fit still trains and identifies on the other sub-scores");
  }

  // ---- X3 reporting --------------------------------------------------------------------------
  {
    ScoreTable t = a;
    // add a constant column so dropped and excluded can be told apart
    ScoreTable wide;
    std::vector<std::string> names = t.feature_names;
    names.push_back("var_constant");
    wide.setColumns(names);
    for (std::size_t i = 0; i < t.rows(); ++i)
    {
      std::vector<float> row(t.row(i), t.row(i) + t.width());
      row.push_back(1.0f);
      wide.append(t.group[i], t.pair[i], t.is_decoy[i] != 0, t.feature_id[i], row.data());
    }
    odia::core::Options o;
    o.exclude_features = {"var_norm_rt_score", "var_im_delta_score", "var_norm_rt"};
    const auto r = odia::core::scoreAndControl(wide, o);
    const auto& dg = r.diagnostics;
    check(dg.features_excluded == std::vector<std::string>({"var_norm_rt_score"}), "X3 excluded by exact name");
    check(dg.exclusions_unmatched == std::vector<std::string>({"var_im_delta_score", "var_norm_rt"}),
          "X3 names matching no column are reported (a prefix is not a match)");
    check(dg.features_dropped == std::vector<std::string>({"var_constant"}), "X3 dropped is not excluded");
    check(dg.features_used.size() == 8, "X3 the eight ordinary sub-scores are used");
  }

  // ---- X4 nothing left ----------------------------------------------------------------------------
  {
    odia::core::Options o;
    o.exclude_features = a.feature_names;
    bool threw = false;
    try { odia::core::scoreAndControl(a, o); }
    catch (const std::invalid_argument&) { threw = true; }
    check(threw, "X4 excluding every column is refused");
  }

  return synth::finish("odia_exclude_features_test");
}
