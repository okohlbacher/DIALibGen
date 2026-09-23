// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

// Ion mobility (M2) through stock OpenMS on a synthetic diaPASEF mzML
// (synthetic_run.h with im_bands 2: every isolation window acquired in two
// abutting 1/K0 bands, the library's 1/K0 MIS-calibrated by a known offset):
//
//   * loading: diaPASEF detected, the MS1 per-peak 1/K0 array seen;
//   * calibration: the library -> run 1/K0 line recovers the planted offset,
//     the automatic window lies in its clamp, and the calibrated 1/K0 puts
//     the seeds in their right band more often than the library's own;
//   * band assignment: planted precursors whose LIBRARY 1/K0 points at the
//     wrong band (or at none) are recovered, because extraction assigns by
//     the calibrated 1/K0 -- and with an identity map in its place they are not;
//   * the planted 1/K0 is recovered within 0.01 from the peak group's 1/K0;
//   * both members of every pair get the same drift time and precursor 1/K0,
//     hence the same band and the same 1/K0 range; decoys have peak groups
//     about as often as targets; the extraction does not depend on search:chunk;
//   * search:im_window -1 searches without ion mobility (no 1/K0, no IM
//     sub-scores), and a fixed width is used as given;
//   * a run WITHOUT ion mobility gives bitwise the same peak groups whatever
//     search:im_window says;
//   * mobilityApex, reportedMobility, windowOf and automaticImWindow on
//     hand-made input.

#include "synthetic_run.h"

#include <odia/search/AssayBuilder.h>
#include <odia/search/EvidencePrefilter.h>
#include <odia/search/Identifier.h>
#include <odia/search/IonMobility.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <string>
#include <vector>

using namespace ODIA::search;
namespace fs = std::filesystem;

namespace
{
  class Stages : public Identifier
  {
  public:
    explicit Stages(SearchParams p)
      : Identifier(std::move(p), [](const std::string& m) { std::cout << "  info: " << m << "\n"; },
                   [](const std::string& m) { std::cout << "  warn: " << m << "\n"; }) {}
    using Identifier::calibrate;
    using Identifier::extract;
    using Identifier::loadRun;
  };

  SearchParams params(double im_window = 0.0, std::size_t chunk = 20000)
  {
    SearchParams p;
    // The fixture plants the LIBRARY's fragments (see identify_extraction).
    p.intensities = Intensities::Library;
    p.threads = 2;
    p.chunk = chunk;
    p.im_window = im_window;
    return p;
  }

  /// Row indices in canonical (precursor, feature id) order.
  std::vector<std::size_t> canonical(const PeakGroups& g)
  {
    std::vector<std::size_t> order(g.rows());
    std::iota(order.begin(), order.end(), std::size_t{0});
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
      if (g.scores.group[a] != g.scores.group[b]) { return g.scores.group[a] < g.scores.group[b]; }
      return g.scores.feature_id[a] < g.scores.feature_id[b];
    });
    return order;
  }

  /// Bitwise equality of two peak-group tables, 1/K0 included.
  bool identical(const PeakGroups& a, const PeakGroups& b)
  {
    if (a.rows() != b.rows() || a.scores.feature_names != b.scores.feature_names)
    {
      std::cout << "  differ: " << a.rows() << " vs " << b.rows() << " rows\n";
      return false;
    }
    const auto oa = canonical(a), ob = canonical(b);
    const std::size_t w = a.scores.width();
    std::size_t differing = 0;
    for (std::size_t k = 0; k < oa.size(); ++k)
    {
      const std::size_t i = oa[k], j = ob[k];
      bool same = a.scores.group[i] == b.scores.group[j] && a.scores.feature_id[i] == b.scores.feature_id[j] &&
                  std::memcmp(a.scores.row(i), b.scores.row(j), w * sizeof(float)) == 0 &&
                  std::memcmp(&a.apex_rt[i], &b.apex_rt[j], sizeof(float)) == 0 &&
                  std::memcmp(&a.rt_start[i], &b.rt_start[j], sizeof(float)) == 0 &&
                  std::memcmp(&a.rt_stop[i], &b.rt_stop[j], sizeof(float)) == 0 &&
                  std::memcmp(&a.im[i], &b.im[j], sizeof(float)) == 0;
      if (!same && differing++ < 5) { std::cout << "  differ at canonical row " << k << "\n"; }
    }
    if (differing) { std::cout << "  " << differing << " of " << oa.size() << " rows differ\n"; }
    return differing == 0;
  }

  int bandOf(const synthrun::Spec& spec, double k0)
  {
    if (!(k0 > spec.im_low && k0 < spec.im_high)) { return -1; }
    return std::min(spec.im_bands - 1, static_cast<int>((k0 - spec.im_low) / ((spec.im_high - spec.im_low) / spec.im_bands)));
  }

  struct Recovery
  {
    std::vector<std::string> missed;                    ///< the first few planted targets not recovered, described
    std::size_t planted = 0, found = 0;                 ///< planted targets in the set, and with a peak group at the apex
    std::size_t wrong_band = 0, wrong_band_found = 0;   ///< ... whose LIBRARY 1/K0 is in another band or none
    std::vector<double> im_error;                       ///< |peak group 1/K0 - true 1/K0| of the found
    std::vector<double> apex_error;                     ///< |apex - planted apex| of the nearest peak group, every planted target with one
    std::size_t im_missing = 0;                         ///< found, but the peak group reports no 1/K0
  };

  /// Planted targets with a peak group within one cycle of their apex.
  Recovery recover(const synthrun::Spec& spec, const synthrun::Fixture& fx, const SearchSet& set, const PeakGroups& g,
                   const OpenMS::TransformationDescription* im_map = nullptr, const std::vector<OpenSwath::SwathMap>* maps = nullptr)
  {
    std::map<std::size_t, double> apex;
    for (const auto& p : fx.planted) { apex[p.index] = p.apex_s; }
    std::vector<int> best(set.size(), -1);
    std::vector<double> nearest(set.size(), std::numeric_limits<double>::infinity());
    for (std::size_t r = 0; r < g.rows(); ++r)
    {
      const auto i = static_cast<std::size_t>(g.scores.group[r]);
      if (set.isDecoy(i)) { continue; }
      const auto it = apex.find(set.source[i]);
      if (it == apex.end()) { continue; }
      const double d = std::fabs(g.apex_rt[r] - it->second);
      nearest[i] = std::min(nearest[i], d);
      if (d <= spec.cycle_s && (best[i] < 0 || d < std::fabs(g.apex_rt[static_cast<std::size_t>(best[i])] - it->second)))
      { best[i] = static_cast<int>(r); }
    }
    Recovery out;
    for (std::size_t k = 0; k < set.pairs(); ++k)
    {
      if (!apex.count(set.source[k])) { continue; }
      ++out.planted;
      if (std::isfinite(nearest[k])) { out.apex_error.push_back(nearest[k]); }
      const double truth = synthrun::trueMobility(spec, set.source[k]);
      const bool wrong = bandOf(spec, fx.library.precursors().im[set.source[k]]) != bandOf(spec, truth);
      out.wrong_band += wrong ? 1 : 0;
      if (best[k] < 0)
      {
        if (out.missed.size() < 12)
        {
          const double lib = fx.library.precursors().im[set.source[k]];
          const double mz = ODIA::fromFixed(set.library.precursors().mz[k]);
          std::string groups;
          for (std::size_t r = 0; r < g.rows(); ++r)
          {
            if (static_cast<std::size_t>(g.scores.group[r]) == k)
            { groups += " " + std::to_string(g.apex_rt[r]) + "@" + std::to_string(g.im[r]); }
          }
          out.missed.push_back(set.precursorId(k) + " mz " + std::to_string(mz) + " apex " + std::to_string(apex[set.source[k]]) +
                               " true " + std::to_string(truth) + " library " + std::to_string(lib) +
                               (im_map ? " calibrated " + std::to_string(im_map->apply(lib)) : std::string()) +
                               (maps && im_map ? " window " + std::to_string(windowOf(*maps, mz, im_map->apply(lib))) + " truth window " +
                                                   std::to_string(windowOf(*maps, mz, truth)) : std::string()) +
                               "; peak groups (apex@1/K0):" + (groups.empty() ? std::string(" none") : groups));
        }
        continue;
      }
      ++out.found;
      out.wrong_band_found += wrong ? 1 : 0;
      const float im = g.im[static_cast<std::size_t>(best[k])];
      if (std::isnan(im)) { ++out.im_missing; } else { out.im_error.push_back(std::fabs(im - truth)); }
    }
    std::sort(out.im_error.begin(), out.im_error.end());
    std::sort(out.apex_error.begin(), out.apex_error.end());
    return out;
  }

  double at(const std::vector<double>& sorted, double q)
  {
    return sorted.empty() ? std::nan("") : sorted[std::min(sorted.size() - 1, static_cast<std::size_t>(q * static_cast<double>(sorted.size())))];
  }
}

int main(int argc, char** argv)
{
  const fs::path dir = argc > 1 ? fs::path(argv[1]) : fs::temp_directory_path() / "dialibgen-identify-ion-mobility";
  fs::remove_all(dir);
  fs::create_directories(dir);

  // ---- 0. the pieces on hand-made input -------------------------------------------
  {
    // Six fragments at 1/K0 0.90 (three peaks each), one intense interference
    // at 1.10 on the first fragment's m/z, noise elsewhere.
    auto spectrum = std::make_shared<OpenSwath::Spectrum>();
    auto mz = std::make_shared<OpenSwath::BinaryDataArray>(), in = std::make_shared<OpenSwath::BinaryDataArray>(),
         k0 = std::make_shared<OpenSwath::BinaryDataArray>();
    k0->description = "Ion Mobility";
    std::vector<std::tuple<double, double, double>> peaks;
    const std::vector<double> fragments = {300.1, 400.2, 500.3, 600.4, 700.5, 800.6};
    for (const double f : fragments)
    {
      for (const double d : {-0.006, 0.0, 0.006}) { peaks.emplace_back(f, d == 0.0 ? 1000.0 : 600.0, 0.90 + d); }
    }
    peaks.emplace_back(300.1, 1e6, 1.10);
    peaks.emplace_back(450.0, 5e5, 0.70);
    std::sort(peaks.begin(), peaks.end());
    for (const auto& [m, i, k] : peaks) { mz->data.push_back(m); in->data.push_back(i); k0->data.push_back(k); }
    spectrum->setMZArray(mz);
    spectrum->setIntensityArray(in);
    spectrum->getDataArrays().push_back(k0);
    const MobilityApex a = mobilityApex({spectrum}, fragments, 10.0, 0.6, 1.4);
    std::cout << "mobilityApex: " << a.im << " from " << a.fragments << " fragments\n";
    CHECK(std::fabs(a.im - 0.90) < 0.002 && a.fragments == 6);
    CHECK(std::isnan(mobilityApex({spectrum}, {1234.5}, 10.0, 0.6, 1.4).im));
    CHECK(std::isnan(reportedMobility(-1.0, 0.9)) && std::isnan(reportedMobility(0.9, 0.95)));
    CHECK(reportedMobility(0.9, 0.91) == 0.9 && reportedMobility(0.9, -1.0) == 0.9 && reportedMobility(0.9, std::nan("")) == 0.9);
    double q99 = 0, sd = 0;
    CHECK(automaticImWindow(std::vector<double>(100, 0.0), q99, sd) == SearchParams::im_window_min);
    CHECK(automaticImWindow(std::vector<double>(100, 0.5), q99, sd) == SearchParams::im_window_max);
    const double w = automaticImWindow(std::vector<double>(100, 0.005), q99, sd);
    CHECK(std::fabs(w - 2.0 * SearchParams::im_window_padding * SearchParams::im_window_normal_quantile * 1.4826 * 0.005) < 1e-9);
    std::vector<OpenSwath::SwathMap> maps(3);
    maps[0].ms1 = true;
    maps[1].lower = 400; maps[1].upper = 425; maps[1].imLower = 0.6; maps[1].imUpper = 1.0;
    maps[2].lower = 400; maps[2].upper = 425; maps[2].imLower = 0.9; maps[2].imUpper = 1.3;
    CHECK(windowOf(maps, 410, 0.7) == 1 && windowOf(maps, 410, 1.2) == 2 && windowOf(maps, 410, 1.4) == -1);
    CHECK(windowOf(maps, 410, 0.94) == 1 && windowOf(maps, 410, 0.96) == 2 && windowOf(maps, 430, 0.8) == -1);
  }

  // ---- 1. the diaPASEF fixture -----------------------------------------------------
  synthrun::Spec spec;
  spec.im_bands = 2;
  spec.im_library_offset = 0.08;   // the library's 1/K0 reads 0.08 high: 11 % of the precursors in the wrong band
  const synthrun::Fixture fx = synthrun::make(spec);
  const std::string mzml = (dir / "pasef.mzML").string();
  synthrun::writeMzML(spec, fx, mzml);
  std::cout << "fixture: " << fx.library.precursorCount() << " precursors, " << fx.planted.size() << " planted, "
            << fs::file_size(mzml) / 1000 << " kB mzML, " << spec.windows << " windows x " << spec.im_bands << " bands\n";

  Stages stages(params());
  RunData run = stages.loadRun(mzml);
  std::size_t ms2 = 0;
  for (const auto& m : run.maps) { ms2 += m.ms1 ? 0 : 1; }
  CHECK(ms2 == static_cast<std::size_t>(spec.windows * spec.im_bands));
  CHECK(run.ion_mobility && run.ms1_ion_mobility);
  CHECK(run.provenance_json.find("\"im_extraction\":true") != std::string::npos);
  const std::vector<IsolationWindow> windows = Identifier::isolationWindows(run);
  const SearchSet set = EvidencePrefilter::select(fx.library, stages.params(), windows, run.maps,
                                                  [](const std::string& m) { std::cout << "  info: " << m << "\n"; });

  // ---- 2. calibration --------------------------------------------------------------
  const Calibration cal = stages.calibrate(fx.library, set, run);
  std::cout << "calibration: " << cal.points << " RT points, 1/K0 window " << cal.im_window << "\n  " << cal.provenance_json << "\n";
  CHECK(!cal.bootstrap && cal.points >= 50);
  CHECK(cal.im_window >= SearchParams::im_window_min && cal.im_window <= SearchParams::im_window_max);
  CHECK(cal.provenance_json.find("\"ion_mobility\"") != std::string::npos);
  {
    // The planted mis-calibration, undone: library 1.08 -> run 1.00.
    const double mapped = cal.im.apply(1.08), slope = cal.im.apply(1.2) - cal.im.apply(1.0);
    std::cout << "1/K0 map: 1.08 -> " << mapped << ", slope " << slope / 0.2 << "\n";
    CHECK(std::fabs(mapped - 1.00) < 0.005);
    CHECK(std::fabs(slope / 0.2 - 1.0) < 0.03);
    const auto j = cal.provenance_json.find("\"window_assignment\"");
    CHECK(j != std::string::npos);
  }

  // ---- 3. extraction: band assignment, 1/K0 recovery, pair symmetry ----------------
  PeakGroups g;
  stages.extract(set, run, cal, g);
  g.validate(set);
  {
    const Recovery r = recover(spec, fx, set, g, &cal.im, &run.maps);
    for (const auto& m : r.missed) { std::cout << "  missed: " << m << "\n"; }
    std::cout << "extraction: " << g.rows() << " peak groups, " << g.scores.width() << " sub-scores; planted " << r.planted
              << ", recovered " << r.found << "; library 1/K0 in the wrong band: " << r.wrong_band << ", recovered "
              << r.wrong_band_found << "; 1/K0 error median " << at(r.im_error, 0.5) << ", p95 " << at(r.im_error, 0.95)
              << ", none reported " << r.im_missing << "; apex error of the nearest peak group median " << at(r.apex_error, 0.5)
              << " s, p90 " << at(r.apex_error, 0.9) << " s, p95 " << at(r.apex_error, 0.95) << " s\n";
    CHECK(r.planted >= 300);
    CHECK(static_cast<double>(r.found) >= 0.95 * static_cast<double>(r.planted));
    CHECK(r.wrong_band >= 20);
    CHECK(static_cast<double>(r.wrong_band_found) >= 0.9 * static_cast<double>(r.wrong_band));
    CHECK(at(r.im_error, 0.5) < 0.004);
    CHECK(at(r.im_error, 0.95) < 0.01);
    CHECK(static_cast<double>(r.im_missing) <= 0.05 * static_cast<double>(r.found));
    for (const char* name : {"var_im_xcorr_shape", "var_im_xcorr_coelution", "var_im_delta_score", "var_im_ms1_delta_score"})
    { CHECK(std::find(g.scores.feature_names.begin(), g.scores.feature_names.end(), name) != g.scores.feature_names.end()); }
    std::size_t targets_with = 0, decoys_with = 0;
    std::vector<char> seen(set.size(), 0);
    for (const auto k : g.scores.group) { seen[static_cast<std::size_t>(k)] = 1; }
    for (std::size_t i = 0; i < set.size(); ++i) { (set.isDecoy(i) ? decoys_with : targets_with) += seen[i]; }
    std::cout << "precursors with peak groups: " << targets_with << " targets, " << decoys_with << " decoys\n";
    CHECK(decoys_with > 0 && static_cast<double>(decoys_with) >= 0.8 * static_cast<double>(targets_with) &&
          static_cast<double>(decoys_with) <= 1.25 * static_cast<double>(targets_with));
  }
  {
    // Both members of every pair: the same drift time and precursor 1/K0, so
    // the same window and the same 1/K0 range.
    AssayOptions o;
    o.ion_mobility = true;
    o.im_map = [&cal](double k0) { return cal.im.apply(k0); };
    std::vector<std::size_t> all(set.size());
    std::iota(all.begin(), all.end(), std::size_t{0});
    const AssayChunk a = AssayBuilder::build(set, all, o);
    std::map<std::string, double> drift;
    for (const auto& c : a.experiment.compounds) { drift[c.id] = c.drift_time; }
    std::size_t pairs_equal = 0, transitions_equal = 0, transitions = 0;
    for (std::size_t k = 0; k < set.pairs(); ++k)
    {
      const double t = drift.at(set.precursorId(k)), d = drift.at(set.precursorId(set.partnerOf(k)));
      pairs_equal += (t == d && t > 0 && windowOf(run.maps, ODIA::fromFixed(set.library.precursors().mz[k]), t) ==
                                          windowOf(run.maps, ODIA::fromFixed(set.library.precursors().mz[set.partnerOf(k)]), d)) ? 1 : 0;
    }
    for (const auto& tr : a.experiment.transitions)
    {
      ++transitions;
      transitions_equal += tr.precursor_im == drift.at(tr.peptide_ref) ? 1 : 0;
    }
    std::cout << "pair symmetry: " << pairs_equal << " of " << set.pairs() << " pairs share drift time and window; "
              << transitions_equal << " of " << transitions << " transitions carry their compound's 1/K0; " << a.missing_im
              << " without 1/K0\n";
    CHECK(pairs_equal == set.pairs() && transitions_equal == transitions && a.missing_im == 0);
  }
  {
    // The assay's 1/K0 is what assigns the band. With the library's own 1/K0
    // (an identity map) and a window wide enough to hold its 0.08 error, the
    // precursors whose library 1/K0 lies in their true band are found and
    // those whose library 1/K0 points at the other band (or none) are lost:
    // they are extracted where their signal is not. With the calibrated 1/K0
    // (above) both kinds are found.
    Calibration uncalibrated = cal;
    uncalibrated.im = OpenMS::TransformationDescription();
    uncalibrated.im_window = 0.2;
    PeakGroups u;
    stages.extract(set, run, uncalibrated, u);
    const Recovery r = recover(spec, fx, set, u);
    const std::size_t right = r.planted - r.wrong_band, right_found = r.found - r.wrong_band_found;
    std::cout << "with the library's 1/K0 and a 0.2 window: recovered " << right_found << " of " << right
              << " in their band, " << r.wrong_band_found << " of " << r.wrong_band << " in the wrong band\n";
    CHECK(static_cast<double>(right_found) >= 0.9 * static_cast<double>(right));
    CHECK(static_cast<double>(r.wrong_band_found) <= 0.2 * static_cast<double>(r.wrong_band));
  }
  {
    // Chunking changes nothing.
    Stages chunked(params(0.0, 300));
    PeakGroups c;
    chunked.extract(set, run, cal, c);
    std::cout << "chunk 300: " << c.rows() << " peak groups\n";
    CHECK(identical(g, c));
  }

  // ---- 4. search:im_window: fixed, and off -----------------------------------------------
  {
    Stages fixed(params(0.07));
    const Calibration f = fixed.calibrate(fx.library, set, run);
    CHECK(f.im_window == 0.07);
    Stages off(params(-1.0));
    RunData r = off.loadRun(mzml);
    CHECK(r.ion_mobility && r.provenance_json.find("\"im_extraction\":false") != std::string::npos);
    const Calibration c = off.calibrate(fx.library, set, r);
    CHECK(c.im_window == -1.0 && c.provenance_json.find("\"ion_mobility\"") == std::string::npos);
    PeakGroups o;
    off.extract(set, r, c, o);
    std::size_t with_im = 0;
    for (const float v : o.im) { with_im += std::isnan(v) ? 0 : 1; }
    const Recovery ro = recover(spec, fx, set, o);
    std::cout << "search:im_window -1: " << o.rows() << " peak groups, " << with_im << " with a 1/K0; planted recovered "
              << ro.found << " of " << ro.planted << ", apex error median " << at(ro.apex_error, 0.5) << " s, p90 "
              << at(ro.apex_error, 0.9) << " s\n";
    CHECK(o.rows() > 0 && with_im == 0);
    CHECK(std::none_of(o.scores.feature_names.begin(), o.scores.feature_names.end(),
                       [](const std::string& n) { return n.rfind("var_im_", 0) == 0; }));
    r.maps.clear();
  }
  run.maps.clear();

  // ---- 5. a run without ion mobility: search:im_window changes nothing --------------------
  {
    const synthrun::Spec plain;
    const synthrun::Fixture pf = synthrun::make(plain);
    const std::string file = (dir / "plain.mzML").string();
    synthrun::writeMzML(plain, pf, file);
    std::vector<PeakGroups> tables;
    for (const double setting : {0.0, 0.05, -1.0})
    {
      Stages s(params(setting));
      RunData r = s.loadRun(file);
      CHECK(!r.ion_mobility && r.provenance_json.find("\"im_extraction\":false") != std::string::npos);
      CHECK(r.provenance_json.find("im_range") == std::string::npos);
      const SearchSet ps = EvidencePrefilter::select(pf.library, s.params(), Identifier::isolationWindows(r), r.maps,
                                                     [](const std::string&) {});
      const Calibration c = s.calibrate(pf.library, ps, r);
      CHECK(c.im_window == -1.0 && c.provenance_json.find("ion_mobility") == std::string::npos);
      tables.emplace_back();
      s.extract(ps, r, c, tables.back());
      r.maps.clear();
    }
    std::cout << "no ion mobility: " << tables[0].rows() << " peak groups at every search:im_window\n";
    CHECK(identical(tables[0], tables[1]) && identical(tables[0], tables[2]));
  }

  fs::remove_all(dir);
  if (synth::failures) { std::cerr << synth::failures << " check(s) failed\n"; return 1; }
  std::cout << "identify_ion_mobility: PASS\n";
  return 0;
}
