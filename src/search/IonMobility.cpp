// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <odia/search/IonMobility.h>

#include <odia/search/SearchParams.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

namespace ODIA::search
{
  namespace
  {
    double quantileOf(std::vector<double> v, double q)
    {
      if (v.empty()) { return std::numeric_limits<double>::quiet_NaN(); }
      std::sort(v.begin(), v.end());
      const double pos = q * static_cast<double>(v.size() - 1);
      const auto a = static_cast<std::size_t>(std::floor(pos));
      const std::size_t b = std::min(v.size() - 1, a + 1);
      return v[a] + (pos - static_cast<double>(a)) * (v[b] - v[a]);
    }
  }

  MobilityApex mobilityApex(const std::vector<OpenSwath::SpectrumPtr>& spectra, const std::vector<double>& fragment_mz,
                            double ppm_half, double im_low, double im_high)
  {
    MobilityApex out;
    if (fragment_mz.empty() || !(im_high > im_low) || !(ppm_half > 0)) { return out; }
    const auto bins = static_cast<std::size_t>(std::ceil((im_high - im_low) / mobility_bin)) + 1;
    auto binOf = [&](double im) {
      return std::min(bins - 1, static_cast<std::size_t>(std::max(0.0, std::floor((im - im_low) / mobility_bin))));
    };
    // One raw mobilogram per fragment, and every matched peak (fragment, 1/K0, intensity).
    struct Hit { std::size_t fragment; double im, intensity; };
    std::vector<Hit> hits;
    std::vector<std::vector<double>> raw(fragment_mz.size(), std::vector<double>(bins, 0.0));
    for (const auto& s : spectra)
    {
      if (!s) { continue; }
      const auto mz_array = s->getMZArray();
      const auto intensity_array = s->getIntensityArray();
      const auto im_array = s->getDriftTimeArray();
      if (!mz_array || !intensity_array || !im_array) { continue; }
      const std::vector<double>& mz = mz_array->data;
      const std::vector<double>& in = intensity_array->data;
      const std::vector<double>& im = im_array->data;
      if (mz.size() != in.size() || mz.size() != im.size() || mz.empty()) { continue; }
      const bool sorted = std::is_sorted(mz.begin(), mz.end());
      for (std::size_t f = 0; f < fragment_mz.size(); ++f)
      {
        const double lo = fragment_mz[f] * (1.0 - ppm_half * 1e-6), hi = fragment_mz[f] * (1.0 + ppm_half * 1e-6);
        const std::size_t first = sorted ? static_cast<std::size_t>(std::lower_bound(mz.begin(), mz.end(), lo) - mz.begin()) : 0;
        for (std::size_t k = first; k < mz.size(); ++k)
        {
          if (mz[k] > hi) { if (sorted) { break; } continue; }
          if (mz[k] < lo || !(in[k] > 0) || !(im[k] >= im_low && im[k] <= im_high)) { continue; }
          raw[f][binOf(im[k])] += in[k];
          hits.push_back({f, im[k], in[k]});
        }
      }
    }
    if (hits.empty()) { return out; }

    // Smooth, scale each fragment to a maximum of 1, sum the votes.
    const int reach = static_cast<int>(std::ceil(3.0 * mobility_smoothing / mobility_bin));
    std::vector<double> kernel(2 * static_cast<std::size_t>(reach) + 1);
    for (int d = -reach; d <= reach; ++d)
    {
      const double x = d * mobility_bin / mobility_smoothing;
      kernel[static_cast<std::size_t>(d + reach)] = std::exp(-0.5 * x * x);
    }
    std::vector<double> votes(bins, 0.0), smooth(bins);
    for (const auto& m : raw)
    {
      std::fill(smooth.begin(), smooth.end(), 0.0);
      double top = 0.0;
      for (std::size_t b = 0; b < bins; ++b)
      {
        double v = 0.0;
        for (int d = -reach; d <= reach; ++d)
        {
          const long long j = static_cast<long long>(b) + d;
          if (j >= 0 && j < static_cast<long long>(bins)) { v += kernel[static_cast<std::size_t>(d + reach)] * m[static_cast<std::size_t>(j)]; }
        }
        smooth[b] = v;
        top = std::max(top, v);
      }
      if (!(top > 0)) { continue; }
      for (std::size_t b = 0; b < bins; ++b) { votes[b] += smooth[b] / top; }
    }
    std::size_t apex = 0;
    for (std::size_t b = 1; b < bins; ++b) { if (votes[b] > votes[apex]) { apex = b; } }
    const double centre = im_low + (static_cast<double>(apex) + 0.5) * mobility_bin;

    // Refine: the intensity-weighted mean of the peaks at the apex.
    double sum = 0.0, weighted = 0.0;
    std::vector<char> seen(fragment_mz.size(), 0);
    for (const Hit& h : hits)
    {
      if (std::fabs(h.im - centre) > mobility_support) { continue; }
      sum += h.intensity;
      weighted += h.intensity * h.im;
      seen[h.fragment] = 1;
    }
    if (!(sum > 0)) { return out; }
    out.im = weighted / sum;
    out.intensity = sum;
    out.fragments = static_cast<std::size_t>(std::count(seen.begin(), seen.end(), 1));
    return out;
  }

  MobilityApex mobilityAt(const std::vector<OpenSwath::SwathMap>& maps, double mz, double rt_s,
                          const std::vector<double>& fragment_mz, double ppm_half)
  {
    std::vector<OpenSwath::SpectrumPtr> spectra;
    double lo = std::numeric_limits<double>::infinity(), hi = -lo;
    std::size_t windows = 0;
    for (const auto& m : maps)
    {
      if (m.ms1 || !m.sptr || !(m.lower < mz && mz < m.upper)) { continue; }
      ++windows;
      lo = std::min(lo, m.imLower);
      hi = std::max(hi, m.imUpper);
      if (m.sptr->getNrSpectra() == 0) { continue; }
      for (auto& sp : m.sptr->getMultipleSpectra(rt_s, mobility_probe_spectra)) { spectra.push_back(std::move(sp)); }
    }
    MobilityApex out;
    if (!spectra.empty() && hi > lo) { out = mobilityApex(spectra, fragment_mz, ppm_half, lo, hi); }
    out.windows = windows;
    return out;
  }

  double libraryMobility(const Library& library, std::size_t i)
  {
    const auto& p = library.precursors();
    double k0 = i < p.im.size() ? static_cast<double>(p.im[i]) : std::numeric_limits<double>::quiet_NaN();
    if (!(std::isfinite(k0) && k0 > 0.0) && i < p.ccs.size() && std::isfinite(p.ccs[i]))
    { k0 = mobilityFromCCS(p.ccs[i], fromFixed(p.mz[i]), p.charge[i]); }
    return std::isfinite(k0) && k0 > 0.0 ? k0 : std::numeric_limits<double>::quiet_NaN();
  }

  double automaticImWindow(const std::vector<double>& abs_residuals, double& q99, double& robust_sd)
  {
    q99 = quantileOf(abs_residuals, 0.99);
    robust_sd = 1.4826 * quantileOf(abs_residuals, 0.5);
    double half = 0.0;
    if (std::isfinite(q99)) { half = std::max(half, q99); }
    if (std::isfinite(robust_sd)) { half = std::max(half, SearchParams::im_window_normal_quantile * robust_sd); }
    const double full = 2.0 * SearchParams::im_window_padding * half;
    return std::min(SearchParams::im_window_max, std::max(SearchParams::im_window_min, full));
  }

  double reportedMobility(double ms2, double ms1)
  {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    if (!(std::isfinite(ms2) && ms2 > 0.0)) { return nan; }
    if (std::isfinite(ms1) && ms1 > 0.0 && std::fabs(ms2 - ms1) > SearchParams::im_ms1_agreement) { return nan; }
    return ms2;
  }

  int windowOf(const std::vector<OpenSwath::SwathMap>& maps, double mz, double im)
  {
    int best = -1;
    for (std::size_t i = 0; i < maps.size(); ++i)
    {
      const auto& m = maps[i];
      if (m.ms1 || !(m.imLower < im && im < m.imUpper && m.lower < mz && mz < m.upper)) { continue; }
      if (best < 0) { best = static_cast<int>(i); continue; }
      const auto& b = maps[static_cast<std::size_t>(best)];
      if (std::fabs((b.imLower + b.imUpper) / 2 - im) > std::fabs((m.imLower + m.imUpper) / 2 - im)) { best = static_cast<int>(i); }
    }
    return best;
  }

  int assignWindow(const std::vector<OpenSwath::SwathMap>& maps, double mz, double im, double half_width, double& assign_im)
  {
    assign_im = im;
    const int held = windowOf(maps, mz, im);
    if (held >= 0 || !std::isfinite(im)) { return held; }
    int best = -1;
    double best_distance = std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < maps.size(); ++i)
    {
      const auto& m = maps[i];
      if (m.ms1 || !(m.lower < mz && mz < m.upper) || !(m.imUpper > m.imLower)) { continue; }
      const double distance = im <= m.imLower ? m.imLower - im : im - m.imUpper;
      if (distance < best_distance) { best_distance = distance; best = static_cast<int>(i); }
    }
    if (best < 0 || !(best_distance < half_width)) { return -1; }
    const auto& m = maps[static_cast<std::size_t>(best)];
    // Strictly inside, as OpenSWATH's pasef assignment requires, and by less
    // than any 1/K0 scan step, so no other window's range is entered.
    const double nudge = std::min(1e-6, (m.imUpper - m.imLower) / 4);
    assign_im = im <= m.imLower ? m.imLower + nudge : m.imUpper - nudge;
    return windowOf(maps, mz, assign_im);
  }
}
