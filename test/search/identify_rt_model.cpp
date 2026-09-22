// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

// The RT calibration's choice between the line and LOWESS
// (Identifier::chooseRtModel), on synthetic calibration points with
// heavy-tailed scatter, as a predicted library's residuals are:
//
//   * a linear truth keeps the line;
//   * a bias confined to the top of the gradient -- the Astral run's case,
//     where the line ran a median 89 s late in the top library-RT decile --
//     switches to LOWESS, whose residuals there are centred again. (The rule
//     before the M3 review, LOWESS only when its 0.99-quantile window is 10 %
//     narrower, kept the line on the real Astral points, 321 s against 317 s:
//     the two windows came from different inlier sets. That comparison is
//     gone, so it is not re-enacted here.);
//   * too few points keep the line without fitting;
//   * the LOWESS outlier cut reports its fits, and never more than 10.

#include "synthetic_library.h"

#include <odia/search/Identifier.h>
#include <odia/search/RobustLine.h>

#include <nlohmann/json.hpp>

#include <OpenMS/DATASTRUCTURES/Param.h>

#include <cmath>
#include <cstdint>
#include <iostream>
#include <utility>
#include <vector>

using namespace ODIA::search;
using json = nlohmann::json;

namespace
{
  struct Rng
  {
    std::uint64_t state;
    double uniform()
    {
      std::uint64_t z = (state += 0x9E3779B97F4A7C15ull);
      z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
      z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
      z ^= z >> 31;
      return (static_cast<double>(z >> 11) + 0.5) / 9007199254740992.0;
    }
    double normal() { return std::sqrt(-2.0 * std::log(uniform())) * std::cos(6.283185307179586 * uniform()); }
  };

  /// (run seconds, assay RT) over a 3000 s run: assay = 0.03 x + bias(x) + noise,
  /// the noise 90 % N(0, 1) and 10 % N(0, 5); the bias rises from 0 at 2600 s
  /// to @p top_bias at 2800 s and stays there.
  std::vector<std::pair<double, double>> make(std::size_t n, double top_bias, std::uint64_t seed)
  {
    Rng r{seed};
    std::vector<std::pair<double, double>> out;
    for (std::size_t k = 0; k < n; ++k)
    {
      const double x = 3000.0 * (static_cast<double>(k) + 0.5) / static_cast<double>(n);
      const double bias = top_bias * std::min(1.0, std::max(0.0, (x - 2600.0) / 200.0));
      const double noise = (r.uniform() < 0.9 ? 1.0 : 5.0) * r.normal();
      out.emplace_back(x, 0.03 * x + bias + noise);
    }
    return out;
  }

  Identifier::RtModelChoice choose(const std::vector<std::pair<double, double>>& stock, json& record)
  {
    const RobustLine line = robustLine(stock, 3.0, 0.1);
    std::vector<std::pair<double, double>> points;
    for (std::size_t k = 0; k < stock.size(); ++k) { if (line.inlier[k]) { points.push_back(stock[k]); } }
    OpenMS::TransformationDescription t;
    t.setDataPoints(points);
    t.fitModel("linear", OpenMS::Param());
    Identifier::RtModelChoice c = Identifier::chooseRtModel(stock, points, t, 0.0, 3000.0);
    record = json::parse(c.record_json);
    return c;
  }
}

int main()
{
  json r;
  {
    const auto c = choose(make(1000, 0.0, 11), r);
    std::cout << "linear truth: " << r.dump() << "\n";
    CHECK(!c.lowess);
    CHECK(r["fits"].get<std::size_t>() >= 1 && r["fits"].get<std::size_t>() <= 10);
    CHECK(r.contains("converged") && r.contains("cv"));
  }
  {
    const auto c = choose(make(1000, 6.0, 12), r);
    std::cout << "top-decile bias: " << r.dump() << "\n";
    CHECK(c.lowess);
    CHECK(c.points.size() >= 200 && !c.fit.getDataPoints().empty());
    CHECK(r["cv"]["error_lowess"].get<double>() < 0.97 * r["cv"]["error_linear"].get<double>());
    // The chosen model's residuals in the top decile are centred again.
    double sum = 0.0;
    std::size_t n = 0;
    for (const auto& p : c.points)
    {
      if (p.first > 2700.0) { sum += p.second - c.fit.apply(p.first); ++n; }
    }
    CHECK(n > 50 && std::fabs(sum / static_cast<double>(n)) < 0.5);
  }
  {
    const auto c = choose(make(150, 8.0, 13), r);
    std::cout << "150 points: " << r.dump() << "\n";
    CHECK(!c.lowess && r["reason"] == "too few points");
  }
  if (synth::failures) { std::cerr << synth::failures << " check(s) failed\n"; return 1; }
  std::cout << "identify_rt_model: PASS\n";
  return 0;
}
