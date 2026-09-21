// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

// The run-facing stages on a synthetic mzML with planted peaks
// (synthetic_run.h), through stock OpenMS: loadRun, calibrate, extract.
//
//   * calibration finds the planted RT map (not the bootstrap), with sane
//     windows, and refuses a run without signal unless allow_bootstrap;
//   * extraction recovers >= 95 % of the planted precursors with a peak
//     group within one cycle of the planted apex, and gives decoys peak
//     groups about as often as targets;
//   * the peak-group table does not depend on search:chunk;
//   * cache mode loads into a scratch directory under search:cache_dir and
//     recovers the same planted precursors;
//   * search:candidates evidence: the prefilter's seeds calibrate the run,
//     their prefilter RT agrees with the calibration points, and extraction
//     of the evidence set recovers the planted precursors it holds, with the
//     same peak groups at two chunk sizes.

#include "synthetic_run.h"

#include <odia/search/EvidencePrefilter.h>
#include <odia/search/Identifier.h>
#include <odia/search/RobustLine.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <iostream>
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

  SearchParams params(std::size_t chunk = 20000)
  {
    SearchParams p;
    p.threads = 2;
    p.chunk = chunk;
    return p;
  }

  /// Planted targets of @p set with a peak group within @p tolerance seconds of the apex.
  std::pair<std::size_t, std::size_t> recovered(const synthrun::Fixture& fx, const SearchSet& set, const PeakGroups& groups,
                                                double tolerance)
  {
    std::map<std::size_t, double> apex;   // library index -> planted apex
    for (const auto& p : fx.planted) { apex[p.index] = p.apex_s; }
    std::vector<char> hit(set.size(), 0);
    for (std::size_t r = 0; r < groups.rows(); ++r)
    {
      const auto i = static_cast<std::size_t>(groups.scores.group[r]);
      if (set.isDecoy(i)) { continue; }
      const auto it = apex.find(set.source[i]);
      if (it != apex.end() && std::fabs(groups.apex_rt[r] - it->second) <= tolerance) { hit[i] = 1; }
    }
    std::size_t planted = 0, found = 0;
    for (std::size_t k = 0; k < set.pairs(); ++k)
    {
      if (apex.count(set.source[k])) { ++planted; found += hit[k]; }
    }
    return {found, planted};
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

  /// Bitwise equality in canonical order; the first difference is printed.
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
      std::string what;
      if (a.scores.group[i] != b.scores.group[j] || a.scores.feature_id[i] != b.scores.feature_id[j]) { what = "ids"; }
      for (std::size_t c = 0; c < w && what.empty(); ++c)
      {
        if (std::memcmp(a.scores.row(i) + c, b.scores.row(j) + c, sizeof(float)) != 0)
        { what = a.scores.feature_names[c] + " " + std::to_string(a.scores.row(i)[c]) + " vs " + std::to_string(b.scores.row(j)[c]); }
      }
      if (what.empty() && (std::memcmp(&a.apex_rt[i], &b.apex_rt[j], sizeof(float)) != 0 ||
                           std::memcmp(&a.rt_start[i], &b.rt_start[j], sizeof(float)) != 0 ||
                           std::memcmp(&a.rt_stop[i], &b.rt_stop[j], sizeof(float)) != 0)) { what = "RT"; }
      if (!what.empty() && differing++ < 5)
      {
        std::cout << "  differ at row " << k << " (precursor " << a.scores.group[i] << ", feature " << a.scores.feature_id[i]
                  << ", apex " << a.apex_rt[i] << "/" << b.apex_rt[j] << "): " << what << "\n";
      }
    }
    if (differing) { std::cout << "  " << differing << " of " << oa.size() << " rows differ\n"; }
    return differing == 0;
  }
}

int main(int argc, char** argv)
{
  const fs::path dir = argc > 1 ? fs::path(argv[1]) : fs::temp_directory_path() / "dialibgen-identify-extraction";
  fs::remove_all(dir);
  fs::create_directories(dir);

  // ---- 0. the robust calibration line: 40 % of the points are noise -------------
  {
    synthrun::Rng r(3);
    std::vector<std::pair<double, double>> points;
    std::vector<char> truth;
    for (int i = 0; i < 300; ++i)
    {
      const double x = 3000.0 * r.uniform();
      const bool real = i % 5 >= 2;
      points.emplace_back(x, real ? 5.0 + 0.03 * x + 0.5 * r.normal() : 100.0 * r.uniform());
      truth.push_back(real ? 1 : 0);
    }
    const RobustLine line = robustLine(points, 3.0, 0.1);
    std::size_t real_in = 0, noise_in = 0;
    for (std::size_t k = 0; k < points.size(); ++k) { (truth[k] ? real_in : noise_in) += line.inlier[k]; }
    std::cout << "robust line: slope " << line.slope << " intercept " << line.intercept << " scale " << line.scale << ", inliers "
              << real_in << " real + " << noise_in << " noise, " << line.iterations << " iterations\n";
    CHECK(line.valid());
    CHECK(std::fabs(line.slope - 0.03) < 0.001 && std::fabs(line.intercept - 5.0) < 1.0);
    CHECK(real_in >= 170 && noise_in <= 10);
    CHECK(!robustLine({{1.0, 2.0}}).valid());
  }

  // The task's example scale: 1,600 library precursors, about 600 planted.
  const synthrun::Spec spec;
  const synthrun::Fixture fx = synthrun::make(spec);
  const std::string mzml = (dir / "run.mzML").string();
  synthrun::writeMzML(spec, fx, mzml);
  std::cout << "fixture: " << fx.library.precursorCount() << " precursors, " << fx.planted.size() << " planted, "
            << fx.background_traces << " background traces, " << fs::file_size(mzml) / 1000 << " kB mzML\n";
  CHECK(fx.planted.size() >= 450 && fx.planted.size() <= 750);

  // ---- 1. load, calibrate, extract ---------------------------------------------
  Stages stages(params());
  const SearchSet set = CandidateSelector::select(fx.library, stages.params());
  RunData run = stages.loadRun(mzml);
  std::size_t ms2 = 0;
  for (const auto& m : run.maps) { ms2 += m.ms1 ? 0 : 1; }
  CHECK(ms2 == static_cast<std::size_t>(spec.windows));
  CHECK(run.maps.size() == ms2 + 1);
  CHECK(run.read_mode == "normal" && run.scratch.empty() && !run.ion_mobility);

  const Calibration cal = stages.calibrate(fx.library, set, run);
  std::cout << "calibration: " << cal.points << " points from " << cal.seeds << " seeds, r^2 " << cal.rsq << ", RT window "
            << cal.rt_window << " s, m/z " << cal.mz_ppm << " / " << cal.ms1_mz_ppm << " ppm\n  " << cal.provenance_json << "\n";
  CHECK(!cal.bootstrap);
  CHECK(cal.points >= 20);
  CHECK(cal.rsq > 0.95);
  CHECK(cal.rt_window >= 30.0 && cal.rt_window <= 90.0);
  CHECK(cal.mz_ppm >= 4.0 && cal.mz_ppm <= 20.0);
  CHECK(cal.im_window == -1.0);
  {
    // The map must put planted precursors where they are.
    OpenMS::TransformationDescription inverse = cal.rt;
    inverse.invert();
    std::vector<double> error;
    for (const auto& p : fx.planted)
    { error.push_back(std::fabs(inverse.apply(set.rt_scale.toAssay(fx.library.precursors().irt[p.index])) - p.apex_s)); }
    std::sort(error.begin(), error.end());
    std::cout << "planted apex vs calibrated RT: median " << error[error.size() / 2] << " s, p95 "
              << error[error.size() * 95 / 100] << " s\n";
    CHECK(error[error.size() / 2] < 4.0);
    CHECK(error[error.size() * 95 / 100] < cal.rt_window / 2);
  }

  PeakGroups a;
  stages.extract(set, run, cal, a);
  a.validate(set);
  const auto [found, planted] = recovered(fx, set, a, spec.cycle_s);
  std::size_t targets_with = 0, decoys_with = 0;
  {
    std::vector<char> seen(set.size(), 0);
    for (const auto g : a.scores.group) { seen[static_cast<std::size_t>(g)] = 1; }
    for (std::size_t i = 0; i < set.size(); ++i) { (set.isDecoy(i) ? decoys_with : targets_with) += seen[i]; }
  }
  std::cout << "extraction: " << a.rows() << " peak groups, " << a.scores.width() << " sub-scores; planted recovered within "
            << spec.cycle_s << " s: " << found << " of " << planted << "; precursors with peak groups: " << targets_with
            << " targets, " << decoys_with << " decoys\n";
  CHECK(planted >= 400);
  CHECK(static_cast<double>(found) >= 0.95 * static_cast<double>(planted));
  CHECK(decoys_with > 0 && static_cast<double>(decoys_with) >= 0.8 * static_cast<double>(targets_with) &&
        static_cast<double>(decoys_with) <= 1.25 * static_cast<double>(targets_with));
  for (const float v : a.im) { CHECK(std::isnan(v)); }
  CHECK(std::find(a.scores.feature_names.begin(), a.scores.feature_names.end(), "var_norm_rt_score") != a.scores.feature_names.end());
  // Every configured sub-score is actually FILLED on a run with MS1: a column
  // stock OpenSWATH never computes for this setup would be dropped silently
  // and its evidence would be missing (the MS1-MS2 co-elution scores were).
  {
    const std::size_t width = a.scores.width();
    std::size_t ms1_columns = 0;
    for (std::size_t c = 0; c < width; ++c)
    {
      std::size_t finite = 0;
      for (std::size_t r = 0; r < a.rows(); ++r) { finite += std::isfinite(a.scores.values[r * width + c]) ? 1 : 0; }
      const std::string& name = a.scores.feature_names[c];
      ms1_columns += name.rfind("var_ms1_", 0) == 0 ? 1 : 0;
      if (finite == 0) { std::cout << "  never filled: " << name << "\n"; }
      CHECK(finite > 0);
    }
    CHECK(ms1_columns == 6);
  }

  // ---- 2. neither repetition nor chunking changes a single value ---------------------
  {
    PeakGroups again;
    stages.extract(set, run, cal, again);
    std::cout << "repeat: " << again.rows() << " peak groups\n";
    CHECK(identical(a, again));
    Stages chunked(params(250));
    PeakGroups b;
    chunked.extract(set, run, cal, b);
    std::cout << "chunk 250: " << b.rows() << " peak groups\n";
    CHECK(identical(a, b));
  }
  run.maps.clear();

  // ---- 3. cache mode: scratch directory under search:cache_dir ------------------------
  {
    SearchParams p = params(1000);
    p.readoptions = ReadMode::Cache;
    p.cache_dir = (dir / "cache").string();
    Stages cached(p);
    RunData c = cached.loadRun(mzml);
    CHECK(c.read_mode == "cache");
    CHECK(!c.scratch.empty() && fs::exists(c.scratch) && c.scratch.parent_path() == fs::path(p.cache_dir));
    PeakGroups g;
    cached.extract(set, c, cal, g);
    const auto [f, n] = recovered(fx, set, g, spec.cycle_s);
    std::cout << "cache mode: " << g.rows() << " peak groups, planted recovered " << f << " of " << n << "\n";
    CHECK(n == planted && static_cast<double>(f) >= 0.95 * static_cast<double>(n));
    c.maps.clear();
    fs::remove_all(c.scratch);
  }

  // ---- 4. a run without signal: calibration refuses, or bootstraps when told to ------
  {
    synthrun::Spec empty = spec;
    empty.planted_fraction = 0.0;
    const synthrun::Fixture none = synthrun::make(empty);
    const std::string blank = (dir / "blank.mzML").string();
    synthrun::writeMzML(empty, none, blank);
    Stages strict(params());
    const SearchSet s = CandidateSelector::select(none.library, strict.params());
    RunData r = strict.loadRun(blank);
    std::string message;
    try { (void)strict.calibrate(none.library, s, r); }
    catch (const SearchAbort& e) { message = e.what(); }
    std::cout << "no signal: " << message << "\n";
    CHECK(message.find("RT calibration failed") != std::string::npos);
    SearchParams p = params();
    p.allow_bootstrap = true;
    Stages lenient(p);
    const Calibration b = lenient.calibrate(none.library, s, r);
    CHECK(b.bootstrap);
    CHECK(b.rt_window > 0 && b.mz_ppm == 20.0);
  }

  // ---- 5. a file that is not a run ----------------------------------------------------
  {
    Stages stray(params());
    bool threw = false;
    try { (void)stray.loadRun((dir / "absent.mzML").string()); } catch (const std::exception&) { threw = true; }
    CHECK(threw);
  }

  // ---- 6. search:candidates evidence: seeds, calibration, extraction ------------------
  {
    SearchParams p = params();
    CHECK(p.candidates == "evidence");   // the default
    Stages ev(p);
    RunData r = ev.loadRun(mzml);
    const std::vector<IsolationWindow> windows = Identifier::isolationWindows(r);
    const SearchSet es = EvidencePrefilter::select(fx.library, ev.params(), windows, r.maps,
                                                   [](const std::string& m) { std::cout << "  info: " << m << "\n"; });
    const Calibration c = ev.calibrate(fx.library, es, r);
    std::cout << "evidence calibration: " << c.points << " points from " << c.seeds << " seeds, r^2 " << c.rsq << ", RT window "
              << c.rt_window << " s\n  " << c.provenance_json << "\n";
    CHECK(es.seeds.size() >= 100);
    CHECK(c.provenance_json.find("\"seed_rule\":\"evidence\"") != std::string::npos);
    CHECK(c.provenance_json.find("\"prefilter_agreement\"") != std::string::npos);
    CHECK(c.provenance_json.find("\"nonlinear\"") != std::string::npos);
    CHECK(!c.bootstrap && c.points >= 50 && c.rsq > 0.95);
    CHECK(c.rt_window >= 30.0 && c.rt_window <= 90.0);
    {
      OpenMS::TransformationDescription inverse = c.rt;
      inverse.invert();
      std::vector<double> error;
      for (const auto& pl : fx.planted)
      { error.push_back(std::fabs(inverse.apply(es.rt_scale.toAssay(fx.library.precursors().irt[pl.index])) - pl.apex_s)); }
      std::sort(error.begin(), error.end());
      std::cout << "planted apex vs evidence-calibrated RT: median " << error[error.size() / 2] << " s, p95 "
                << error[error.size() * 95 / 100] << " s\n";
      CHECK(error[error.size() / 2] < 4.0);
      CHECK(error[error.size() * 95 / 100] < c.rt_window / 2);
    }
    PeakGroups g1, g2;
    ev.extract(es, r, c, g1);
    g1.validate(es);
    Stages small(params(300));
    small.extract(es, r, c, g2);
    const auto [f, n] = recovered(fx, es, g1, spec.cycle_s);
    std::cout << "evidence set: " << es.pairs() << " pairs, planted among them " << n << ", recovered " << f << "; "
              << g1.rows() << " peak groups (chunk 300: " << g2.rows() << ")\n";
    CHECK(n >= 400 && static_cast<double>(f) >= 0.95 * static_cast<double>(n));
    CHECK(identical(g1, g2));
    r.maps.clear();
  }

  fs::remove_all(dir);
  if (synth::failures) { std::cerr << synth::failures << " check(s) failed\n"; return 1; }
  std::cout << "identify_extraction: PASS\n";
  return 0;
}
