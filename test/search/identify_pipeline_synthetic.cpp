// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

// Identifier end to end on a SYNTHETIC run: the three run-facing seams are
// replaced by a generator of peak groups (signal on a known share of targets,
// noise everywhere else), so candidate selection, decoys, scoring, the
// entity roll-ups, every run-level guard, the report and the provenance are
// exercised without OpenMS extraction.

#include "synthetic_library.h"

#include <odia/LibraryRefiner.h>
#include <odia/search/Identifier.h>

#include <nlohmann/json.hpp>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <random>
#include <string>
#include <vector>

using namespace ODIA::search;
namespace fs = std::filesystem;
using json = nlohmann::json;

namespace
{
  std::uint64_t mix(std::uint64_t z)
  {
    z += 0x9E3779B97F4A7C15ull;
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
    return z ^ (z >> 31);
  }

  struct Plan
  {
    int present_percent = 40;   ///< targets whose first peak group carries signal
    float signal = 2.0f;
    bool decoys = true;         ///< false: decoys get no peak groups (a broken extraction)
  };

  class SyntheticRun : public Identifier
  {
  public:
    SyntheticRun(SearchParams p, Plan plan)
      : Identifier(std::move(p), [](const std::string&) {}, [](const std::string&) {}), plan_(plan) {}

    fs::path scratch;   ///< what loadRun created
    std::size_t present = 0;

    bool isPresent(const SearchSet& set, std::size_t i) const
    {
      return !set.isDecoy(i) && static_cast<int>(mix(static_cast<std::uint64_t>(set.pairOf(i)) + 17) % 100) < plan_.present_percent;
    }

  protected:
    RunData loadRun(const std::string& path) override
    {
      RunData r;
      r.path = path;
      r.name = "synthetic";
      r.read_mode = "synthetic";
      static int counter = 0;
      r.scratch = fs::temp_directory_path() / ("dialibgen-synthetic-scratch-" + std::to_string(++counter));
      fs::create_directories(r.scratch);
      std::ofstream(r.scratch / "window_0.cached") << "cache\n";
      scratch = r.scratch;
      return r;
    }

    Calibration calibrate(const ODIA::Library&, const SearchSet& set, RunData&) override
    {
      Calibration c;
      c.seeds = set.pairs();
      c.points = set.pairs();
      c.rsq = 0.99;
      c.coverage = 1.0;
      c.rt_window = 300.0;
      c.mz_ppm = 15.0;
      c.ms1_mz_ppm = 10.0;
      return c;
    }

    void extract(const SearchSet& set, RunData&, const Calibration&, PeakGroups& out) override
    {
      out.setColumns({"var_xcorr_shape_weighted", "var_xcorr_coelution_weighted", "var_library_corr", "var_library_sangle",
                      "var_log_sn_score", "var_massdev_score_weighted", "var_isotope_correlation_score", "var_norm_rt_score"});
      present = 0;
      // Backwards, to show the row order extraction produces does not matter.
      for (std::size_t n = set.size(); n-- > 0;)
      {
        if (set.isDecoy(n) && !plan_.decoys) { continue; }
        const bool signal = isPresent(set, n);
        present += signal ? 1 : 0;
        const double apex_min = 5.0 + 0.3 * set.library.precursors().irt[n];
        for (int f = 0; f < 3; ++f)
        {
          std::mt19937_64 rng(mix(n * 8 + static_cast<std::size_t>(f)));
          std::normal_distribution<float> noise(0.0f, 1.0f);
          float v[8];
          for (float& x : v) { x = noise(rng); }
          if (signal && f == 0) { for (int c = 0; c < 6; ++c) { v[c] += plan_.signal; } }
          const float apex = static_cast<float>(apex_min * 60.0) + 20.0f * static_cast<float>(f);
          out.add(set, n, f, v, apex, apex - 6.0f, apex + 6.0f, std::nanf(""));
        }
      }
    }

  private:
    Plan plan_;
  };

  SearchParams params(int threads = 1)
  {
    SearchParams p;
    p.subset = 1000;
    p.max_pairs = 0;
    p.min_ids = 100;
    p.threads = threads;
    return p;
  }

  std::string bytes(const fs::path& f)
  {
    std::ifstream in(f, std::ios::binary);
    return std::string(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
  }

  /// identify() must throw SearchAbort mentioning @p expected, write no report
  /// and leave no scratch directory behind.
  void expectAbort(const ODIA::Library& lib, SearchParams p, Plan plan, const fs::path& out, const std::string& expected)
  {
    SyntheticRun run(std::move(p), plan);
    std::string message;
    try { (void)run.identify(lib, "synthetic.mzML", out.string(), "test"); }
    catch (const SearchAbort& e) { message = e.what(); }
    catch (const std::exception& e) { message = std::string("NOT A GUARD: ") + e.what(); }
    std::cout << "abort (" << expected << "): " << message << "\n";
    CHECK(message.find(expected) != std::string::npos);
    CHECK(!fs::exists(out));
    CHECK(!run.scratch.empty() && !fs::exists(run.scratch));
  }
}

int main(int argc, char** argv)
{
  const fs::path dir = argc > 1 ? fs::path(argv[1]) : fs::temp_directory_path() / "dialibgen-identify-synthetic";
  fs::remove_all(dir);
  fs::create_directories(dir);
  const ODIA::Library lib = synth::library(synth::peptides(1500, 21));

  // ---- 1. a good run ----------------------------------------------------------
  const fs::path a = dir / "a.ids.parquet";
  SyntheticRun run(params(1), Plan());
  const IdentificationResult result = run.identify(lib, "synthetic.mzML", a.string(), "test");
  CHECK(fs::exists(a));
  CHECK(!fs::exists(run.scratch));
  const json prov = json::parse(result.provenance_json);
  std::cout << "identified " << result.identified << " (" << run.present << " targets with signal); report "
            << result.report_rows << " rows, " << result.report_decoys << " decoys\n";
  CHECK(result.identified >= 100);
  CHECK(result.identified <= run.present + run.present / 20 + 5);
  CHECK(result.run_name == "synthetic");
  CHECK(prov["experimental"] == true);
  CHECK(prov["candidates"]["pairs"] == prov["scoring"]["target_precursors"]);
  CHECK(prov["identifications"]["precursors"] == result.identified);
  CHECK(prov["report"]["rows"] == result.report_rows);
  CHECK(prov["settings"]["decoys"] == "shuffle");
  CHECK(prov["scoring"]["iterations_trained"].get<int>() > 0);
  CHECK(prov["selftest"].is_null() && prov["entrapment"].is_null());
  CHECK(prov["run"]["read_mode"] == "synthetic");
  CHECK(prov["calibration"]["rsq"] == 0.99);

  // The report through refine's reader: false discoveries stay near 1 %.
  {
    ODIA::RefineParams p;
    ODIA::RefineStats st;
    const auto obs = ODIA::LibraryRefiner::readObservations(a.string(), p, st);
    CHECK(st.run == "synthetic");
    CHECK(st.ids_decoy == result.report_decoys);
    std::size_t false_ids = 0, ids = 0;
    // Rebuild which keys carried signal from the SearchSet the run searched.
    const SearchSet set = CandidateSelector::select(lib, params(1));
    for (std::size_t k = 0; k < set.pairs(); ++k)
    {
      const auto it = obs.find(ODIA::LibraryRefiner::key(set.modifiedSequence(k), set.charge(k)));
      if (it == obs.end() || it->second.q > 0.01) { continue; }
      ++ids;
      false_ids += run.isPresent(set, k) ? 0 : 1;
      const double expected_min = 5.0 + 0.3 * set.library.precursors().irt[k];
      CHECK(std::fabs(it->second.rt - expected_min) < 1.0);   // minutes, one of the three peak groups
    }
    std::cout << "refine reader: " << ids << " precursors pass all three gates, " << false_ids << " without signal\n";
    CHECK(ids > 0);
    CHECK(static_cast<double>(false_ids) <= 0.03 * static_cast<double>(ids) + 2.0);
  }

  // ---- 2. thread count does not change the report ------------------------------
  {
    const fs::path b = dir / "b.ids.parquet";
    SyntheticRun again(params(3), Plan());
    const auto r = again.identify(lib, "synthetic.mzML", b.string(), "test");
    CHECK(r.identified == result.identified);
    CHECK(bytes(a) == bytes(b));
  }

  // ---- 3. self-checks and the RT/IM ablation -------------------------------------
  {
    SearchParams p = params(2);
    p.selftest = true;
    p.rt_im_scores = false;
    SyntheticRun checked(p, Plan());
    const auto r = checked.identify(lib, "synthetic.mzML", (dir / "c.ids.parquet").string(), "test");
    const json j = json::parse(r.provenance_json);
    std::cout << "self-check: " << j["selftest"].dump() << "\n";
    CHECK(!j["selftest"].is_null());
    CHECK(j["selftest"]["label_swap_ids"].get<std::size_t>() <= j["selftest"]["limit"].get<std::size_t>());
    CHECK(j["selftest"]["random_label_ids"].get<std::size_t>() <= j["selftest"]["limit"].get<std::size_t>());
    CHECK(j["scoring"]["features_excluded"] == json::array({"var_norm_rt_score"}));
  }

  // ---- 4. every guard aborts with counts, writes nothing, cleans up ----------------
  {
    SearchParams p = params(1);
    p.min_ids = 100000;
    expectAbort(lib, p, Plan(), dir / "min_ids.parquet", "fewer than search:min_ids");
  }
  {
    Plan broken;
    broken.decoys = false;
    expectAbort(lib, params(1), broken, dir / "ratio.parquet", "decoy:target ratio");
  }
  {
    Plan easy;
    easy.present_percent = 90;
    easy.signal = 4.0f;
    expectAbort(lib, params(1), easy, dir / "fraction.parquet", "search:max_target_fraction");
  }
  {
    Plan null;
    null.present_percent = 0;
    expectAbort(lib, params(1), null, dir / "null.parquet", "search: ");
  }
  {
    ScoringOutcome untrained;
    untrained.scored.diagnostics.iterations_trained = 0;
    untrained.scored.diagnostics.target_groups = 1000;
    std::string message;
    try { checkGuards(untrained, params(1)); } catch (const SearchAbort& e) { message = e.what(); }
    CHECK(message.find("learned no discriminant") != std::string::npos);
  }
  // An existing report is refused before the run is touched.
  {
    SyntheticRun refused(params(1), Plan());
    bool threw = false;
    try { (void)refused.identify(lib, "synthetic.mzML", a.string(), "test"); }
    catch (const std::runtime_error& e) { threw = std::string(e.what()).find("refusing to overwrite") != std::string::npos; }
    CHECK(threw);
    CHECK(refused.scratch.empty());
  }

  fs::remove_all(dir);
  if (synth::failures) { std::cerr << synth::failures << " check(s) failed\n"; return 1; }
  std::cout << "identify_pipeline_synthetic: PASS\n";
  return 0;
}
