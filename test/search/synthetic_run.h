// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

// A synthetic centroided DIA run with planted chromatographic peaks, and the
// library it was planted from. A fixture for the -run tests, not a physical
// model:
//
//   * library: tryptic-looking peptides (synth::peptides), charges 2 and 3,
//     up to six y ions with seeded relative intensities;
//   * signal: a share of the peptides elutes as Gaussian peaks at a known apex
//     -- a learnable function of the sequence (hydrophobicity and length) on a
//     gently curved gradient, of which the library RT is a noisy prediction --
//     all fragments co-eluting in their isolation window, and a precursor
//     isotope envelope in MS1;
//   * background, the null: every library precursor AND its shuffle decoy
//     (built exactly as the search builds them) gets one independent
//     interference trace per fragment and one at its precursor m/z, each at
//     its own random time near the precursor's expected RT. Absent targets
//     and decoys therefore have peak groups too -- badly co-eluting ones --
//     as on a real run, and the two classes see the same process;
//   * uniform random noise peaks in every spectrum.
//
// Everything is a function of the Spec and its seed: the random numbers come
// from SplitMix64, not from <random>'s distributions, whose output the
// standard leaves to the implementation.
#pragma once

#include "synthetic_library.h"

#include <odia/Library.h>
#include <odia/LibraryGenerator.h>

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/OPTIONS/PeakFileOptions.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/METADATA/Precursor.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <string>
#include <vector>

namespace synthrun
{
  struct Spec
  {
    std::size_t peptides = 800;       ///< library peptides; each at charge 2 and 3
    double planted_fraction = 0.375;  ///< of the peptides (both charges planted)
    std::uint64_t seed = 20260921;
    int windows = 20;                 ///< MS2 isolation windows, abutting
    double mz_low = 350.0, mz_high = 1000.0;
    double cycle_s = 0.5;             ///< one MS1 and every window per cycle
    double run_s = 300.0;             ///< acquisition from 0 to run_s
    double peak_sigma_s = 2.5;        ///< Gaussian elution profile
    double background_spread_s = 25.0;   ///< interference apex: expected RT +- this
    int noise_ms2 = 40, noise_ms1 = 150;   ///< uniform noise peaks per spectrum
    double mz_error_ppm = 2.0;        ///< SD of every peak's m/z error
  };

  /// One planted precursor: library index (targets only), apex and height.
  struct Planted
  {
    std::size_t index = 0;
    double apex_s = 0.0;
    double height = 0.0;
  };

  /// One Gaussian trace: m/z, apex, height, and the window it is acquired in
  /// (-1 = MS1).
  struct Trace
  {
    double mz, apex_s, height;
    int window;
  };

  struct Fixture
  {
    ODIA::Library library;            ///< targets only
    std::vector<Planted> planted;     ///< by library index
    std::vector<char> is_planted;     ///< per library precursor
    std::vector<Trace> traces;        ///< planted signal and background, in no particular order
    std::size_t background_traces = 0;
  };

  /// SplitMix64: a counter-based stream, identical on every platform.
  struct Rng
  {
    std::uint64_t state;
    explicit Rng(std::uint64_t seed) : state(seed) {}
    std::uint64_t next()
    {
      std::uint64_t z = (state += 0x9E3779B97F4A7C15ull);
      z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
      z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
      return z ^ (z >> 31);
    }
    double uniform() { return static_cast<double>(next() >> 11) * 0x1.0p-53; }   // [0, 1)
    double normal()
    {
      const double u1 = 1.0 - uniform(), u2 = uniform();
      return std::sqrt(-2.0 * std::log(u1)) * std::cos(6.283185307179586 * u2);
    }
  };

  inline std::uint64_t mixKey(std::uint64_t a, std::uint64_t b)
  {
    Rng r(a * 0x9E3779B97F4A7C15ull ^ b);
    return r.next();
  }

  /// RT (0 .. 100) -> apex in seconds, before noise: a gently curved gradient.
  inline double apexOf(const Spec& spec, double irt)
  {
    const double lo = 0.1 * spec.run_s, hi = 0.9 * spec.run_s;
    return lo + (hi - lo) * irt / 100.0 + 0.015 * spec.run_s * std::sin(3.141592653589793 * irt / 100.0);
  }

  /// The TRUE elution position of a peptide on 0 .. 100: mean Kyte-Doolittle
  /// hydrophobicity and length (the law test/refinement/synth_report.py uses),
  /// so a fine-tuned RT model can learn it. The library RT is this plus
  /// prediction noise.
  inline double trueRt(const std::string& sequence)
  {
    auto kd = [](char c) {
      switch (c)
      {
        case 'A': return 1.8; case 'R': return -4.5; case 'N': return -3.5; case 'D': return -3.5; case 'C': return 2.5;
        case 'Q': return -3.5; case 'E': return -3.5; case 'G': return -0.4; case 'H': return -3.2; case 'I': return 4.5;
        case 'L': return 3.8; case 'K': return -3.9; case 'M': return 1.9; case 'F': return 2.8; case 'P': return -1.6;
        case 'S': return -0.8; case 'T': return -0.7; case 'W': return -0.9; case 'Y': return -1.3; case 'V': return 4.2;
        default: return 0.0;
      }
    };
    double h = 0;
    for (const char c : sequence) { h += kd(c); }
    h /= static_cast<double>(std::max<std::size_t>(1, sequence.size()));
    const double minutes30 = 10.5 + 4.0 * h + 0.25 * static_cast<double>(sequence.size());
    return std::min(99.0, std::max(1.0, 100.0 * minutes30 / 30.0));
  }

  inline int windowOf(const Spec& spec, double mz)
  {
    if (!(mz > spec.mz_low && mz < spec.mz_high)) { return -2; }
    const double width = (spec.mz_high - spec.mz_low) / spec.windows;
    return std::min(spec.windows - 1, static_cast<int>((mz - spec.mz_low) / width));
  }

  /// The target library: peptides in synth::peptides order, charge 2 then 3.
  inline void buildLibrary(const Spec& spec, ODIA::Library& library)
  {
    auto& pre = library.precursors();
    auto& tr = library.transitions();
    const auto peptides = synth::peptides(spec.peptides, spec.seed);
    for (std::size_t n = 0; n < peptides.size(); ++n)
    {
      const auto& p = peptides[n];
      const OpenMS::AASequence seq = OpenMS::AASequence::fromString(p.sequence);
      pre.mz.push_back(ODIA::toFixed(seq.getMZ(p.charge)));
      // A prediction: the true position plus noise, the same for both charges.
      Rng predicted(mixKey(spec.seed + 4, n / 2));
      pre.irt.push_back(static_cast<float>(std::min(100.0, std::max(0.0, trueRt(p.sequence) + 1.5 * predicted.normal()))));
      pre.im.push_back(std::nanf(""));
      pre.ccs.push_back(std::nanf(""));
      pre.charge.push_back(static_cast<std::uint8_t>(p.charge));
      pre.decoy.push_back(0);
      pre.modified_sequence.push_back(library.strings().intern(p.sequence));
      pre.protein_group.push_back(library.strings().intern(p.protein));
      pre.transition_begin.push_back(static_cast<std::uint32_t>(tr.product_mz.size()));
      // y3 .. y(n-1), at most six, with seeded relative intensities (max 1).
      Rng r(mixKey(spec.seed + 1, n));
      std::vector<std::pair<double, int>> frags;
      for (std::size_t ord = 3; ord + 1 <= seq.size() && frags.size() < 6; ++ord)
      { frags.emplace_back(0.1 + 0.9 * r.uniform(), static_cast<int>(ord)); }
      double top = 0;
      for (const auto& f : frags) { top = std::max(top, f.first); }
      for (const auto& [intensity, ord] : frags)
      {
        tr.product_mz.push_back(ODIA::toFixed(seq.getSuffix(static_cast<OpenMS::Size>(ord)).getMZ(1, OpenMS::Residue::YIon)));
        tr.library_intensity.push_back(static_cast<float>(intensity / top));
        tr.type.push_back(ODIA::FragmentType::Y);
        tr.ordinal.push_back(static_cast<std::uint8_t>(ord));
        tr.charge.push_back(1);
        tr.loss.push_back(ODIA::LossType::None);
      }
      pre.transition_count.push_back(static_cast<std::uint32_t>(frags.size()));
    }
  }

  inline Fixture make(const Spec& spec)
  {
    Fixture fx;
    buildLibrary(spec, fx.library);
    const auto& pre = fx.library.precursors();
    const auto& tr = fx.library.transitions();
    const std::size_t n = fx.library.precursorCount();

    // Signal: per peptide (both charges share the apex).
    Rng noise(spec.seed ^ 0x5eedULL);
    std::vector<double> peptide_apex(spec.peptides, 0.0);
    std::vector<char> peptide_planted(spec.peptides, 0);
    for (std::size_t k = 0; k < spec.peptides; ++k)
    {
      peptide_planted[k] = static_cast<double>(mixKey(spec.seed, k) >> 11) * 0x1.0p-53 < spec.planted_fraction;
      const std::string sequence(fx.library.strings().get(pre.modified_sequence[2 * k]));
      peptide_apex[k] = apexOf(spec, trueRt(sequence)) + 1.0 * noise.normal();
    }
    fx.is_planted.assign(n, 0);
    for (std::size_t i = 0; i < n; ++i)
    {
      const std::size_t k = i / 2;
      const double mz = ODIA::fromFixed(pre.mz[i]);
      const double apex = peptide_apex[k];
      const int window = windowOf(spec, mz);
      if (!peptide_planted[k] || window < 0 || apex < 5 * spec.peak_sigma_s || apex > spec.run_s - 5 * spec.peak_sigma_s) { continue; }
      Rng r(mixKey(spec.seed + 2, i));
      const double height = std::exp(std::log(1e4) + (std::log(1e6) - std::log(1e4)) * r.uniform());
      fx.is_planted[i] = 1;
      fx.planted.push_back({i, apex, height});
      for (std::uint32_t f = 0; f < pre.transition_count[i]; ++f)
      {
        const std::uint32_t j = pre.transition_begin[i] + f;
        fx.traces.push_back({ODIA::fromFixed(tr.product_mz[j]), apex,
                             height * std::max(0.02, tr.library_intensity[j] * (1.0 + 0.15 * r.normal())), window});
      }
      const OpenMS::AASequence seq = OpenMS::AASequence::fromString(std::string(fx.library.strings().get(pre.modified_sequence[i])));
      const int z = pre.charge[i];
      const auto dist = seq.getFormula(OpenMS::Residue::Full, z).getIsotopeDistribution(OpenMS::CoarseIsotopePatternGenerator(3));
      for (std::size_t iso = 0; iso < dist.size(); ++iso)
      { fx.traces.push_back({mz + static_cast<double>(iso) * 1.0033548378 / z, apex, 3.0 * height * dist[iso].getIntensity(), -1}); }
    }

    // Background for targets and their shuffle decoys alike: one independent
    // trace per fragment and one at the precursor m/z, near the expected RT.
    ODIA::Library all;
    buildLibrary(spec, all);
    ODIA::LibraryGenerator::appendDecoys(all, ODIA::DecoyMethod::Shuffle, nullptr, 3, false);
    const auto& ap = all.precursors();
    const auto& at = all.transitions();
    for (std::size_t i = 0; i < all.precursorCount(); ++i)
    {
      const double mz = ODIA::fromFixed(ap.mz[i]);
      const int window = windowOf(spec, mz);
      if (window < 0) { continue; }
      Rng r(mixKey(spec.seed + 3, i));
      const double expected = apexOf(spec, ap.irt[i]);
      auto when = [&]() { return expected + spec.background_spread_s * (2.0 * r.uniform() - 1.0); };
      for (std::uint32_t f = 0; f < ap.transition_count[i]; ++f)
      {
        const std::uint32_t j = ap.transition_begin[i] + f;
        fx.traces.push_back({ODIA::fromFixed(at.product_mz[j]), when(), std::exp(std::log(1500.0) + 0.7 * r.normal()), window});
        ++fx.background_traces;
      }
      fx.traces.push_back({mz, when(), std::exp(std::log(3000.0) + 0.7 * r.normal()), -1});
      ++fx.background_traces;
    }
    return fx;
  }

  /// Write the run as centroided mzML (32-bit, zlib). MS2 windows abut on
  /// [mz_low, mz_high); each cycle is one MS1 spectrum and then every window.
  inline void writeMzML(const Spec& spec, const Fixture& fx, const std::string& path)
  {
    const double width = (spec.mz_high - spec.mz_low) / spec.windows;
    const int cycles = static_cast<int>(std::floor(spec.run_s / spec.cycle_s));
    const double reach = 4.0 * spec.peak_sigma_s;

    // Traces by acquisition channel (MS1 = 0, window w = w + 1), by apex.
    std::vector<std::vector<Trace>> channel(static_cast<std::size_t>(spec.windows) + 1);
    for (const auto& t : fx.traces) { channel[static_cast<std::size_t>(t.window + 1)].push_back(t); }
    for (auto& c : channel)
    { std::sort(c.begin(), c.end(), [](const Trace& a, const Trace& b) { return a.apex_s < b.apex_s || (a.apex_s == b.apex_s && a.mz < b.mz); }); }

    OpenMS::MSExperiment exp;
    Rng r(spec.seed ^ 0xD1AULL);
    auto addPeak = [&](OpenMS::MSSpectrum& s, double mz, double intensity) {
      if (!(intensity > 0)) { return; }
      OpenMS::Peak1D p;
      p.setMZ(mz * (1.0 + spec.mz_error_ppm * 1e-6 * r.normal()));
      p.setIntensity(static_cast<float>(intensity));
      s.push_back(p);
    };
    int native = 0;
    for (int c = 0; c < cycles; ++c)
    {
      for (int w = -1; w < spec.windows; ++w)
      {
        const double t = c * spec.cycle_s + (w + 1) * spec.cycle_s / (spec.windows + 1);
        OpenMS::MSSpectrum s;
        s.setRT(t);
        s.setMSLevel(w < 0 ? 1 : 2);
        s.setType(OpenMS::SpectrumSettings::SpectrumType::CENTROID);
        s.setNativeID("scan=" + std::to_string(++native));
        if (w >= 0)
        {
          OpenMS::Precursor prec;
          const double lo = spec.mz_low + w * width;
          prec.setMZ(lo + width / 2);
          prec.setIsolationWindowLowerOffset(width / 2);
          prec.setIsolationWindowUpperOffset(width / 2);
          s.setPrecursors({prec});
        }
        const auto& traces = channel[static_cast<std::size_t>(w + 1)];
        auto it = std::lower_bound(traces.begin(), traces.end(), t - reach,
                                   [](const Trace& a, double v) { return a.apex_s < v; });
        for (; it != traces.end() && it->apex_s <= t + reach; ++it)
        {
          const double elution = std::exp(-0.5 * std::pow((t - it->apex_s) / spec.peak_sigma_s, 2));
          addPeak(s, it->mz, it->height * elution * (1.0 + 0.05 * r.normal()));
        }
        const int n_noise = w < 0 ? spec.noise_ms1 : spec.noise_ms2;
        const double nlo = w < 0 ? spec.mz_low : 100.0, nhi = w < 0 ? spec.mz_high : 1500.0;
        for (int k = 0; k < n_noise; ++k)
        { addPeak(s, nlo + (nhi - nlo) * r.uniform(), std::exp(std::log(200.0) + r.normal())); }
        s.sortByPosition();
        exp.addSpectrum(s);
      }
    }
    OpenMS::MzMLFile file;
    file.getOptions().setMz32Bit(true);
    file.getOptions().setIntensity32Bit(true);
    file.getOptions().setCompression(true);
    file.store(path, exp);
  }

  /// Truth table: one line per planted precursor (library index, id, charge,
  /// precursor m/z, apex seconds, height).
  inline void writeTruth(const Fixture& fx, const std::string& path)
  {
    std::ofstream out(path);
    out << "index\tprecursor_id\tcharge\tprecursor_mz\tapex_s\theight\n";
    const auto& pre = fx.library.precursors();
    out.precision(10);
    for (const auto& p : fx.planted)
    {
      out << p.index << '\t' << fx.library.strings().get(pre.modified_sequence[p.index]) << static_cast<int>(pre.charge[p.index])
          << '\t' << static_cast<int>(pre.charge[p.index]) << '\t' << ODIA::fromFixed(pre.mz[p.index]) << '\t' << p.apex_s
          << '\t' << p.height << '\n';
    }
  }
}
