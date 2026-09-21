// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

/// A robust straight line through calibration points, for the RT calibration.
///
/// Why: OpenMS's iterative residual outlier removal (MRMRTNormalizer) stops as
/// soon as r^2 reaches search:calibration_min_rsq, so the points it returns
/// can still hold every wrong peak that did not push r^2 below that target --
/// on a synthetic run with a third of the seeds present, 40 of 192 kept
/// points were noise and the slope was 15 % off. This refit does not trust
/// the contaminated least-squares line: it starts from the repeated-median
/// line (Siegel; breakdown point 50 %), keeps the points within `cut` robust
/// standard deviations (1.4826 x MAD of the residuals), refits by least
/// squares on those, and repeats until the kept set is stable.
///
/// Standard library only, deterministic (no sampling), O(n^2) in the number
/// of points, which is at most a few hundred seeds.
#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <utility>
#include <vector>

namespace ODIA::search
{
  struct RobustLine
  {
    double intercept = std::numeric_limits<double>::quiet_NaN();
    double slope = std::numeric_limits<double>::quiet_NaN();
    double scale = std::numeric_limits<double>::quiet_NaN();   ///< robust SD of the inlier residuals (y units)
    std::vector<char> inlier;                                  ///< per input point
    std::size_t inliers = 0;
    int iterations = 0;
    bool valid() const { return std::isfinite(slope) && std::isfinite(intercept) && inliers >= 2; }
  };

  namespace detail
  {
    inline double median(std::vector<double> v)
    {
      if (v.empty()) { return std::numeric_limits<double>::quiet_NaN(); }
      const std::size_t h = v.size() / 2;
      std::nth_element(v.begin(), v.begin() + static_cast<std::ptrdiff_t>(h), v.end());
      const double upper = v[h];
      if (v.size() % 2 == 1) { return upper; }
      return 0.5 * (upper + *std::max_element(v.begin(), v.begin() + static_cast<std::ptrdiff_t>(h)));
    }
  }

  /// y = intercept + slope * x through @p points (x, y). @p cut is in robust
  /// standard deviations; @p min_scale is a floor on that SD (y units), so a
  /// near-perfect fit does not reject points for tiny residuals.
  inline RobustLine robustLine(const std::vector<std::pair<double, double>>& points, double cut = 3.0,
                               double min_scale = 0.0, int max_iterations = 20)
  {
    RobustLine out;
    const std::size_t n = points.size();
    out.inlier.assign(n, 0);
    if (n < 2) { return out; }

    // Repeated median: slope = med_i med_{j != i} (y_j - y_i) / (x_j - x_i).
    std::vector<double> per_point;
    per_point.reserve(n);
    std::vector<double> slopes;
    slopes.reserve(n);
    for (std::size_t i = 0; i < n; ++i)
    {
      slopes.clear();
      for (std::size_t j = 0; j < n; ++j)
      {
        const double dx = points[j].first - points[i].first;
        if (j != i && dx != 0.0) { slopes.push_back((points[j].second - points[i].second) / dx); }
      }
      if (!slopes.empty()) { per_point.push_back(detail::median(slopes)); }
    }
    if (per_point.empty()) { return out; }
    double slope = detail::median(per_point);
    std::vector<double> r(n);
    for (std::size_t i = 0; i < n; ++i) { r[i] = points[i].second - slope * points[i].first; }
    double intercept = detail::median(r);

    std::vector<char> keep(n, 1), previous;
    double scale = 0.0;
    for (int it = 0; it < max_iterations; ++it)
    {
      std::vector<double> abs_r;
      abs_r.reserve(n);
      for (std::size_t i = 0; i < n; ++i)
      {
        r[i] = points[i].second - (intercept + slope * points[i].first);
        if (it == 0 || keep[i]) { abs_r.push_back(std::fabs(r[i])); }
      }
      scale = std::max(min_scale, 1.4826 * detail::median(abs_r));
      previous = keep;
      std::size_t kept = 0;
      for (std::size_t i = 0; i < n; ++i)
      {
        keep[i] = (scale > 0.0 ? std::fabs(r[i]) <= cut * scale : r[i] == 0.0) ? 1 : 0;
        kept += keep[i];
      }
      out.iterations = it + 1;
      if (kept < 2) { break; }
      // Least squares on the kept points.
      double sx = 0, sy = 0, sxx = 0, sxy = 0, m = 0;
      for (std::size_t i = 0; i < n; ++i)
      {
        if (!keep[i]) { continue; }
        const double x = points[i].first, y = points[i].second;
        sx += x; sy += y; sxx += x * x; sxy += x * y; m += 1.0;
      }
      const double var = m * sxx - sx * sx;
      if (!(var > 0.0)) { break; }
      slope = (m * sxy - sx * sy) / var;
      intercept = (sy - slope * sx) / m;
      if (keep == previous && it > 0) { break; }
    }
    out.slope = slope;
    out.intercept = intercept;
    out.scale = scale;
    out.inlier = keep;
    out.inliers = static_cast<std::size_t>(std::count(keep.begin(), keep.end(), 1));
    return out;
  }
}
