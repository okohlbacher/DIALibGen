// Copyright (c) 2026, Oliver Kohlbacher and the DIALibRefine authors.
// SPDX-License-Identifier: BSD-3-Clause

// Fragment-intensity replacement, tested without any fixture: the library and
// the observations are built here, so every expected value is visible beside
// the assertion that checks it.
//
// The test that matters most is mz_mismatch_throws. A fragment numbering that
// is off by one still agrees on IDENTITY for most fragments, every counter
// still looks healthy, and every intensity lands on the wrong fragment. Only
// the m/z cross-check notices.

#include <odia/LibraryRefiner.h>

#include <arrow/api.h>
#include <arrow/config.h>
#include <arrow/compute/api.h>
#include <arrow/io/file.h>
#include <arrow/util/config.h>
#include <parquet/arrow/reader.h>
#include <parquet/arrow/writer.h>

#include <chrono>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <locale>
#include <string>
#include <vector>

using namespace ODIA;

namespace
{
  int failures = 0;
  void check(bool ok, const std::string& what)
  { if (!ok) { ++failures; std::cerr << "FAIL: " << what << "\n"; } }

  struct Frag { FragmentType type; int ordinal; int charge; double mz; float intensity; };

  /// One target precursor "PEPTIDEK"/2 with the given transitions.
  Library makeLibrary(const std::vector<Frag>& frags)
  {
    Library lib;
    auto& p = lib.precursors();
    auto& t = lib.transitions();
    p.mz.push_back(toFixed(500.25)); p.irt.push_back(0.5f); p.im.push_back(1.0f); p.ccs.push_back(400.0f);
    p.charge.push_back(2); p.decoy.push_back(0);
    p.modified_sequence.push_back(lib.strings().intern("PEPTIDEK"));
    p.protein_group.push_back(lib.strings().intern("sp|P1|TEST"));
    p.transition_begin.push_back(0);
    p.transition_count.push_back(static_cast<std::uint32_t>(frags.size()));
    for (const Frag& f : frags)
    {
      t.product_mz.push_back(toFixed(f.mz)); t.library_intensity.push_back(f.intensity);
      t.type.push_back(f.type); t.ordinal.push_back(static_cast<std::uint8_t>(f.ordinal));
      t.charge.push_back(static_cast<std::int8_t>(f.charge)); t.loss.push_back(LossType::None);
    }
    return lib;
  }

  LibraryRefiner::ObsMap observe(const std::string& info, const std::string& quant, const std::string& corr)
  {
    Observation o;
    o.rt = 42.0f; o.im = 0.9f; o.q = 0.001f;
    check(parseFragmentInfo(info, quant, corr, o.frags), "fixture triple parses");
    LibraryRefiner::ObsMap m;
    m.emplace(LibraryRefiner::key("PEPTIDEK", 2), o);
    return m;
  }

  const std::vector<Frag> LIB = {
    {FragmentType::Y, 5, 1, 600.30, 1.00f}, {FragmentType::Y, 4, 1, 500.25, 0.80f},
    {FragmentType::B, 3, 1, 300.15, 0.60f}, {FragmentType::Y, 6, 1, 700.35, 0.40f},
  };

  void parser()
  {
    std::vector<ObservedFragment> f;
    std::size_t bad = 0;
    check(parseFragmentInfo("b9^1/668.37;y7^2/330.17;", "10;20;", "0.9;-0.1;", f, &bad) && f.size() == 2 && bad == 0,
          "a well-formed triple parses, and the trailing ';' makes no phantom fragment");
    check(f.size() == 2 && f[1].type == FragmentType::Y && f[1].ordinal == 7 && f[1].charge == 2 && f[1].correlation < 0,
          "type, ordinal, charge and a NEGATIVE correlation survive");
    f.assign(1, ObservedFragment{});
    check(!parseFragmentInfo("b9^1/668.37;y7^1/659.35;", "10;", "0.9;0.8;", f) && f.size() == 1,
          "ragged lists are a ROW failure and write nothing");
    bad = 0;
    check(parseFragmentInfo("by9^1/668.37;b0^1/1.0;b3^9/1.0;b3^1/x;y7^1/659.35;", "1;1;1;1;5;", "1;1;1;1;1;", f, &bad)
          && f.size() == 1 && bad == 4, "an unknown series, ordinal 0, charge 9 and a non-numeric m/z are bad tokens, not b ions");
    check(parseFragmentInfo("y7^1/659.35;", "5;", "", f) && f.size() == 1 && std::isnan(f[0].correlation),
          "an absent correlation column gives NaN, not zero");
    for (const std::string token : {"0x1p3", " 5", "5 ", "+-5", "++5", "--5", "-+5", "5,1", "nan", "inf"})
    {
      for (int column = 0; column < 3; ++column)
      {
        bad = 0;
        check(parseFragmentInfo(column == 0 ? "y7^1/" + token : "y7^1/659.35",
                                column == 1 ? token : "5", column == 2 ? token : "0.9", f, &bad) &&
              f.empty() && bad == 1, "non-decimal fragment token is rejected in column " + std::to_string(column) + ": " + token);
      }
    }
    struct CommaDecimal : std::numpunct<char> { char do_decimal_point() const override { return ','; } };
    const auto original = std::locale::global(std::locale(std::locale::classic(), new CommaDecimal));
    check(parseFragmentInfo("y7^1/+6.5935e2", "+5.0e1", "-.1", f) && f.size() == 1 &&
          std::abs(f[0].mz - 659.35) < 1e-9 && f[0].quant == 50 && std::abs(f[0].correlation + 0.1f) < 1e-7f,
          "decimal scientific notation uses the classic locale under a comma-decimal global locale");
    std::locale::global(original);
  }

  void replaces_and_reranks()
  {
    Library lib = makeLibrary(LIB);
    // y4 is the observed base peak; b3 anti-correlates; y6 has no area.
    const auto obs = observe("y5^1/600.30;y4^1/500.25;b3^1/300.15;y6^1/700.35;", "500;1000;800;0;", "0.9;0.95;-0.2;0.9;");
    RefineParams p; p.write_rt = false; p.write_intensity = true; p.intensity_min_fragments = 2;
    RefineStats st;
    LibraryRefiner::refine(lib, obs, p, st);
    const auto& t = lib.transitions();
    check(st.intensity_replaced_precursors == 1 && lib.precursors().transition_count[0] == 2, "two trusted fragments survive");
    check(st.intensity_gated_correlation == 1 && st.intensity_gated_zero_quant == 1, "each gate reports its own rejection");
    check(t.type[0] == FragmentType::Y && t.ordinal[0] == 4 && t.ordinal[1] == 5, "re-ranked by OBSERVED intensity: y4 now leads");
    // library_max: the observed maximum takes the value the precursor's maximum already had (1.0).
    check(std::abs(t.library_intensity[0] - 1.0f) < 1e-6f && std::abs(t.library_intensity[1] - 0.5f) < 1e-6f,
          "library_max scaling: 1000 -> 1.0 and 500 -> 0.5");
    check(std::abs(fromFixed(t.product_mz[0]) - 500.25) < 1e-3, "a transition keeps its OWN m/z when it moves");
    check(st.intensity_matched_transitions + st.intensity_mz_mismatch + st.intensity_unmatched_in_library +
          st.intensity_loss_bearing == st.intensity_candidate_transitions, "the four fates of a library transition are exhaustive");
    check(st.intensity_rank_agreement == 0.0, "rank agreement records that the predicted base peak (y5) was wrong here");
  }

  void no_restrict_preserves_counts()
  {
    Library lib = makeLibrary(LIB);
    const auto obs = observe("y5^1/600.30;y4^1/500.25;b3^1/300.15;y6^1/700.35;", "500;1000;800;0;", "0.9;0.95;-0.2;0.9;");
    RefineParams p; p.write_rt = false; p.write_intensity = true; p.intensity_restrict = false;
    RefineStats st;
    bool threw = false;
    // Not every transition is trusted, so under no-restrict this precursor keeps its
    // predictions -- and then NOTHING was replaced, which is an error by design.
    try { LibraryRefiner::refine(lib, obs, p, st); } catch (const std::runtime_error&) { threw = true; }
    check(threw && lib.precursors().transition_count[0] == 4 && lib.transitions().library_intensity[0] == 1.00f,
          "no-restrict never changes a transition count, and a run that replaces nothing is refused");
  }

  void decoy_intensities_are_symmetric()
  {
    Library lib = makeLibrary(LIB).subsetByIndex({0, 0});
    lib.precursors().decoy[1] = 1;
    auto& transitions = lib.transitions();
    for (std::size_t j = LIB.size(); j < transitions.product_mz.size(); ++j)
    { transitions.product_mz[j] += toFixed(37.0); }
    // A decoy has its own m/z and may arrive in a different transition order.
    auto reverse_decoy = [&](auto& values) { std::reverse(values.begin() + LIB.size(), values.end()); };
    reverse_decoy(transitions.product_mz); reverse_decoy(transitions.library_intensity);
    reverse_decoy(transitions.type); reverse_decoy(transitions.ordinal);
    reverse_decoy(transitions.charge); reverse_decoy(transitions.loss);
    const auto obs = observe("y5^1/600.30;y4^1/500.25;b3^1/300.15;y6^1/700.35;",
                             "500;1000;500;0;", "0.9;0.9;0.9;0.9;");
    RefineParams p; p.write_rt = false; p.write_intensity = true;
    RefineStats st;
    LibraryRefiner::refine(lib, obs, p, st);
    check(st.intensity_replaced_precursors == 1 && st.intensity_replaced_decoys == 1 &&
          lib.precursors().transition_count == std::vector<std::uint32_t>({3, 3}),
          "target and shifted decoy both keep the three trusted fragment identities");
    const auto& t = lib.transitions();
    check(t.ordinal == std::vector<std::uint8_t>({4, 3, 5, 4, 3, 5}) &&
          t.library_intensity == std::vector<float>({1.0f, 0.5f, 0.5f, 1.0f, 0.5f, 0.5f}),
          "target and decoy use identical observed values and identity tie-breaking order");
    for (std::size_t j = 0; j < 3 && t.product_mz.size() == 6; ++j)
    {
      check(t.type[j] == t.type[j + 3] && t.ordinal[j] == t.ordinal[j + 3] &&
            t.charge[j] == t.charge[j + 3] && t.loss[j] == t.loss[j + 3],
            "paired transitions retain identical fragment identities");
      check(t.product_mz[j + 3] - t.product_mz[j] == toFixed(37.0),
            "decoy transitions retain their own shifted m/z after reordering");
    }

    const auto all_trusted = observe("y5^1/600.30;y4^1/500.25;b3^1/300.15;y6^1/700.35;",
                                     "500;1000;500;250;", "0.9;0.9;0.9;0.9;");
    for (int scenario = 0; scenario < 4; ++scenario)
    {
      lib = makeLibrary(LIB).subsetByIndex({0, 0});
      lib.precursors().decoy[1] = 1;
      auto& missing = lib.transitions();
      if (scenario < 2)
      {
        --lib.precursors().transition_count[1];
        missing.product_mz.pop_back(); missing.library_intensity.pop_back();
        missing.type.pop_back(); missing.ordinal.pop_back(); missing.charge.pop_back(); missing.loss.pop_back();
      }
      else { missing.ordinal.back() = 5; } // Duplicate y5 must not hide the absent y6 by preserving the count.
      const Library before = lib.subsetByIndex({0, 1});
      p.intensity_restrict = scenario % 2; st = RefineStats{};
      std::string error;
      try { LibraryRefiner::refine(lib, all_trusted, p, st); }
      catch (const std::runtime_error& e) { error = e.what(); }
      check(error.find("replaced no precursor") != std::string::npos && st.intensity_decoy_asymmetry == 1 &&
            st.intensity_replaced_precursors == 0 && st.intensity_replaced_decoys == 0,
            "a decoy missing a kept target fragment reverts both spectra, also under no-restrict");
      check(lib.precursors().transition_begin == before.precursors().transition_begin &&
            lib.precursors().transition_count == before.precursors().transition_count &&
            lib.transitions().product_mz == before.transitions().product_mz &&
            lib.transitions().library_intensity == before.transitions().library_intensity &&
            lib.transitions().type == before.transitions().type && lib.transitions().ordinal == before.transitions().ordinal &&
            lib.transitions().charge == before.transitions().charge && lib.transitions().loss == before.transitions().loss,
            "refused asymmetric replacement leaves both rows' complete transition arrays unchanged");
    }
  }

  void mixed_rt_units_refused()
  {
    Library lib = makeLibrary(LIB).subsetByIndex({0, 0});
    lib.precursors().modified_sequence[1] = lib.strings().intern("UNSEENPEPK");
    const auto obs = observe("y5^1/600.30;", "500;", "0.9;");
    RefineParams p; p.filter = false;
    RefineStats st;
    std::string error;
    try { LibraryRefiner::refine(lib, obs, p, st); }
    catch (const std::runtime_error& e) { error = e.what(); }
    check(error.find("mix reference-run minutes") != std::string::npos &&
          lib.precursors().irt == std::vector<float>({0.5f, 0.5f}),
          "shared refinement API refuses mixed RT units before modifying either precursor");
    lib.precursors().irt = {15.0f, 20.0f};
    p.library_rt_in_minutes = true;
    LibraryRefiner::refine(lib, obs, p, st);
    check(lib.precursors().irt == std::vector<float>({42.0f, 20.0f}),
          "known whole-library minutes permit observed RT while retaining unmatched predictions");
  }

  void mz_mismatch_throws()
  {
    Library lib = makeLibrary(LIB);
    // The same fragments, numbered one higher: identities still collide with the
    // library's (y5, y6 exist on both sides) but the m/z belong to other fragments.
    const auto obs = observe("y6^1/600.30;y5^1/500.25;b4^1/300.15;y7^1/700.35;", "500;1000;800;300;", "0.9;0.95;0.9;0.9;");
    RefineParams p; p.write_rt = false; p.write_intensity = true; p.intensity_min_fragments = 1;
    RefineStats st;
    std::string msg;
    try { LibraryRefiner::refine(lib, obs, p, st); } catch (const std::runtime_error& e) { msg = e.what(); }
    check(msg.find("ppm") != std::string::npos, "a shifted fragment numbering is REFUSED, naming the ppm check: " + msg);
    check(lib.transitions().library_intensity[0] == 1.00f, "and nothing was written before refusing");
  }

  void mixed_provenance_refused()
  {
    Library lib = makeLibrary(LIB);
    const auto obs = observe("y5^1/600.30;", "500;", "0.9;");
    RefineParams p; p.write_rt = false; p.write_intensity = true; p.filter = false;
    RefineStats st;
    bool threw = false;
    try { LibraryRefiner::refine(lib, obs, p, st); } catch (const std::runtime_error&) { threw = true; }
    check(threw, "-write_intensity with the filter off is refused without -allow_mixed_intensity");
  }

  void duplicate_keys_count_once()
  {
    Library lib = makeLibrary(LIB).subsetByIndex({0, 0});
    auto obs = observe("y5^1/600.30;", "500;", "0.9;");
    obs.emplace(LibraryRefiner::key("ABSENTPEPK", 2), Observation{});
    RefineParams p; p.write_rt = false;
    RefineStats st;
    LibraryRefiner::refine(lib, obs, p, st);
    check(st.matched == 2 && st.ids_unmatched == 1 && st.match_fraction == 0.5,
          "matching two target rows on one key covers only one of two reference keys");
    p.min_match_fraction = 0.75;
    st = RefineStats{};
    bool refused = false;
    try { LibraryRefiner::refine(lib, obs, p, st); }
    catch (const std::runtime_error&) { refused = true; }
    check(refused, "duplicate targets cannot bypass minimum reference coverage");
  }

  void loss_bearing_rank_and_preservation()
  {
    Library lib = makeLibrary(LIB);
    lib.transitions().loss[0] = LossType::Water;
    const auto obs = observe("y4^1/500.25;b3^1/300.15;y6^1/700.35;", "1000;800;300;", "0.9;0.9;0.9;");
    RefineParams p; p.write_rt = false; p.write_intensity = true;
    RefineStats st;
    LibraryRefiner::refine(lib, obs, p, st);
    check(st.intensity_loss_bearing == 1 && st.intensity_rank_agreement == 1.0,
          "rank agreement compares no-loss transitions when the original maximum carries a loss");
    check(lib.transitions().library_intensity.front() == 1.0f,
          "library_max normalization still uses the original spectrum maximum");
    lib = makeLibrary(LIB); lib.transitions().loss[0] = LossType::Water;
    p.intensity_restrict = false; st = RefineStats{};
    std::string message;
    try { LibraryRefiner::refine(lib, obs, p, st); }
    catch (const std::runtime_error& e) { message = e.what(); }
    check(message.find("neutral-loss") != std::string::npos && lib.transitionCount() == LIB.size(),
          "no-restrict explains why a neutral-loss precursor keeps its full predicted spectrum");
  }

  void report_layouts()
  {
    const auto path = std::filesystem::temp_directory_path() /
      ("dialibgen-fragments-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) + ".parquet");
    auto status = [](const arrow::Status& s) { if (!s.ok()) { throw std::runtime_error(s.ToString()); } };
    auto write = [&](bool numbered, const std::string& omit = "", const std::string& first_id = "y5^1/600.30",
                     double quantity = 500.0)
    {
      std::vector<std::shared_ptr<arrow::Field>> fields;
      std::vector<std::shared_ptr<arrow::Array>> arrays;
      auto text = [&](const std::string& name, const std::vector<std::string>& values)
      {
        if (name == omit) { return; }
        arrow::StringBuilder b; status(b.AppendValues(values));
        std::shared_ptr<arrow::Array> a; status(b.Finish(&a));
        fields.push_back(arrow::field(name, arrow::utf8())); arrays.push_back(a);
      };
      auto numbers = [&](const std::string& name, const std::vector<double>& values)
      {
        if (name == omit) { return; }
        arrow::DoubleBuilder b; status(b.AppendValues(values));
        std::shared_ptr<arrow::Array> a; status(b.Finish(&a));
        fields.push_back(arrow::field(name, arrow::float64())); arrays.push_back(a);
      };
      text("Modified.Sequence", {"PEPTIDEK", "ENTRAPK", "BADQPEPK", "DECOYPEPK", "NANQPEPK"});
      text("Protein.Group", {"P1", "ENTRAP_P1", "P2", "DECOY_P2", "P3"});
      numbers("Precursor.Charge", {2, 2, 2, 2, 2});
      numbers("Q.Value", {0.001, 0.001, 0.02, 0.001, std::numeric_limits<double>::quiet_NaN()});
      numbers("Global.Q.Value", {0.001, 0.001, 0.001, 0.001, 0.001});
      numbers("PG.Q.Value", {0.001, 0.001, 0.001, 0.001, 0.001});
      numbers("Decoy", {0, 0, 0, 1, 0});
      numbers("RT", {42, 42, 42, 42, 42});
      if (numbered)
      {
        // Sparse slot numbers are valid; identities, not slot positions, match the library.
        const std::vector<std::string> prefix = {"Fr.0", "Fr.2", "Fr.4", "Fr.9"};
        const std::vector<std::string> ids = {first_id, "y4^1/500.25", "b3^1/300.15", "y6^1/700.35"};
        const std::vector<double> quant = {quantity, 1000, 800, 0}, corr = {0.9, 0.95, -0.2, 0.9};
        for (std::size_t i = 0; i < prefix.size(); ++i)
        {
          text(prefix[i] + ".Id", std::vector<std::string>(5, ids[i]));
          numbers(prefix[i] + ".Quantity", std::vector<double>(5, quant[i]));
          numbers(prefix[i] + ".Score", std::vector<double>(5, corr[i]));
        }
        text("Fr.12.Id", std::vector<std::string>(5, ""));
        numbers("Fr.12.Quantity", std::vector<double>(5, 0));
        numbers("Fr.12.Score", std::vector<double>(5, 0));
      }
      else
      {
        text("Fragment.Info", std::vector<std::string>(5, "y5^1/600.30;y4^1/500.25;b3^1/300.15;y6^1/700.35;"));
        text("Fragment.Quant.Raw", std::vector<std::string>(5, "500;1000;800;0;"));
        text("Fragment.Correlations", std::vector<std::string>(5, "0.9;0.95;-0.2;0.9;"));
      }
      auto output = arrow::io::FileOutputStream::Open(path.string());
      if (!output.ok()) { throw std::runtime_error(output.status().ToString()); }
      status(parquet::arrow::WriteTable(*arrow::Table::Make(arrow::schema(fields), arrays),
        arrow::default_memory_pool(), *output, 5));
      status((*output)->Close());
    };
    RefineParams p; p.write_rt = false; p.write_intensity = true; p.intensity_min_fragments = 2;
    RefineStats old_stats, new_stats;
    std::cerr << "report layouts: writing DIA-NN 1.9 fixture\n";
    write(false);
    {
      std::cerr << "report layouts: opening Arrow input\n";
      auto input = arrow::io::ReadableFile::Open(path.string());
      if (!input.ok()) { throw std::runtime_error(input.status().ToString()); }
      std::cerr << "report layouts: opening Parquet reader\n";
      auto reader = parquet::arrow::OpenFile(*input, arrow::default_memory_pool());
      if (!reader.ok()) { throw std::runtime_error(reader.status().ToString()); }
      std::shared_ptr<arrow::Table> table;
      std::cerr << "report layouts: reading Parquet table\n";
      status((*reader)->ReadTable(&table));
      std::cerr << "report layouts: casting numeric ChunkedArray\n";
      auto numbers = arrow::compute::Cast(table->GetColumnByName("Precursor.Charge"), arrow::float64());
      if (!numbers.ok()) { throw std::runtime_error(numbers.status().ToString()); }
      check(numbers->chunked_array()->length() == 5, "numeric cast preserves all rows");
      std::cerr << "report layouts: casting string ChunkedArray\n";
      auto strings = arrow::compute::Cast(table->GetColumnByName("Modified.Sequence"), arrow::utf8());
      if (!strings.ok()) { throw std::runtime_error(strings.status().ToString()); }
      check(strings->chunked_array()->length() == 5, "string cast preserves all rows");
      std::cerr << "report layouts: Arrow input and casts passed\n";
    }
    std::cerr << "report layouts: reading DIA-NN 1.9 fixture\n";
    const auto old_obs = LibraryRefiner::readObservations(path.string(), p, old_stats);
    std::cerr << "report layouts: writing DIA-NN 2 fixture\n";
    write(true);
    std::cerr << "report layouts: reading DIA-NN 2 fixture\n";
    const auto new_obs = LibraryRefiner::readObservations(path.string(), p, new_stats);
    check(new_obs.size() == 2 && new_obs.count(LibraryRefiner::key("ENTRAPK", 2)) == 1,
      "DIA-NN 2 applies the same gates to target and ENTRAP proteins; rejects q, NaN and decoy rows");
    check(new_stats.ids_q_above == 1 && new_stats.ids_q_invalid == 1 && new_stats.ids_decoy == 1,
      "DIA-NN 2 reports each row rejection");
    Library old_lib = makeLibrary(LIB), new_lib = makeLibrary(LIB);
    LibraryRefiner::refine(old_lib, old_obs, p, old_stats);
    LibraryRefiner::refine(new_lib, new_obs, p, new_stats);
    check(old_lib.transitions().library_intensity == new_lib.transitions().library_intensity &&
          old_lib.transitions().product_mz == new_lib.transitions().product_mz &&
          old_stats.intensity_matched_transitions == new_stats.intensity_matched_transitions &&
          new_stats.intensity_bad_tokens == 0 && new_stats.intensity_replaced_precursors == 1,
      "DIA-NN 1.9 lists and 2.x numbered columns yield the same intensities, m/z and counters");

    auto refused = [&](const std::string& missing)
    {
      write(true, missing); RefineStats s;
      try { LibraryRefiner::readObservations(path.string(), p, s); } catch (const std::runtime_error&) { return true; }
      return false;
    };
    check(refused("Fr.2.Quantity"), "DIA-NN 2 refuses an incomplete fragment triplet");
    check(refused("Fr.2.Id"), "DIA-NN 2 refuses quantity/score without identity");
    check(refused("Fr.2.Score"), "DIA-NN 2 refuses a missing Score when the quality gate is on");
    p.intensity_min_correlation = -1;
    RefineStats ungated;
    check(LibraryRefiner::readObservations(path.string(), p, ungated).at(LibraryRefiner::key("PEPTIDEK", 2)).frags.size() == 4,
      "a missing Score is accepted only when its quality gate is explicitly disabled");
    p.intensity_min_correlation = 0;

    p.write_rt = true;
    check(refused("RT"), "requested RT replacement refuses a report without an RT column");
    p.write_rt = false;
    RefineStats no_rt;
    check(LibraryRefiner::readObservations(path.string(), p, no_rt).size() == 2,
          "filtering without RT replacement permits an absent RT column");
    p.write_im = true;
    check(refused("IM"), "requested mobility replacement refuses a report without an IM column");
    p.write_im = false;

    write(true, "", "y5^1/601.30"); RefineStats shifted; Library shifted_lib = makeLibrary(LIB);
    const auto shifted_obs = LibraryRefiner::readObservations(path.string(), p, shifted);
    bool mz_refused = false;
    try { LibraryRefiner::refine(shifted_lib, shifted_obs, p, shifted); }
    catch (const std::runtime_error& e) { mz_refused = std::string(e.what()).find("ppm") != std::string::npos; }
    check(mz_refused && shifted_lib.transitions().library_intensity[0] == 1.0f,
      "DIA-NN 2 identities retain the mandatory m/z check before writing");
    write(true, "", "by5^1/600.30"); RefineStats malformed;
    const auto malformed_obs = LibraryRefiner::readObservations(path.string(), p, malformed);
    check(malformed.intensity_bad_tokens == 2 && malformed_obs.at(LibraryRefiner::key("PEPTIDEK", 2)).frags.size() == 3,
      "DIA-NN 2 rejects malformed identities for each passing row");
    write(true, "", "y5^1/600.30", std::numeric_limits<double>::infinity()); RefineStats nonfinite;
    LibraryRefiner::readObservations(path.string(), p, nonfinite);
    check(nonfinite.intensity_bad_tokens == 2, "DIA-NN 2 rejects nonfinite fragment quantities");

    // Real fragment fields from diann-bench-3.0/out/D0_r2.parquet, row 0.
    // DIA-NN 2.0 Academia, compiled Jan 28 2025 05:36:10;
    // --report-lib-info --export-quant. Report SHA-256:
    // 634897a30468ca1a6f9aa802fb874aeb6b12216164d339ca533c0b08cd19fb50
    // Only fragment fields are copied; the joining key and gates below are synthetic.
    struct RealFragment { const char* id; FragmentType type; int ordinal; double mz; float quant, score; };
    const RealFragment captured[] = {
      {"y16^1/1364.691650", FragmentType::Y, 16, 1364.691650, 318.015380859375f, 0.5786356925964355f},
      {"b8^1/569.304199", FragmentType::B, 8, 569.304199, 432.0181579589844f, 0.10540767014026642f},
      {"b7^1/498.267059", FragmentType::B, 7, 498.267059, 341.01593017578125f, 0.2240317314863205f},
      {"b9^1/668.372620", FragmentType::B, 9, 668.372620, 390.0226745605469f, 0.48129919171333313f},
      {"b6^1/427.229950", FragmentType::B, 6, 427.229950, 467.0227966308594f, 0.32008859515190125f},
      {"y17^1/1463.760132", FragmentType::Y, 17, 1463.760132, 0.0f, 0.0f},
      {"b5^1/356.192841", FragmentType::B, 5, 356.192841, 172.00811767578125f, -0.10410183668136597f},
      {"y13^1/1109.569824", FragmentType::Y, 13, 1109.569824, 60.00340270996094f, 0.022481275722384453f},
      {"y12^1/1052.548340", FragmentType::Y, 12, 1052.548340, 53.00226974487305f, 0.17406094074249268f},
      {"b4^1/285.155731", FragmentType::B, 4, 285.155731, 274.0181579589844f, 0.3053339421749115f},
      {"y18^1/1534.797241", FragmentType::Y, 18, 1534.797241, 0.0f, 0.0f},
      {"y7^1/659.347107", FragmentType::Y, 7, 659.347107, 123.00689697265625f, 0.16552139818668365f},
    };
    std::vector<std::shared_ptr<arrow::Field>> fields;
    std::vector<std::shared_ptr<arrow::Array>> arrays;
    auto add = [&](const std::string& name, auto& builder, auto value)
    {
      status(builder.Append(value));
      std::shared_ptr<arrow::Array> array; status(builder.Finish(&array));
      fields.push_back(arrow::field(name, array->type())); arrays.push_back(array);
    };
    arrow::StringBuilder text;
    arrow::Int64Builder integer;
    arrow::FloatBuilder number;
    add("Modified.Sequence", text, "PEPTIDEK");
    add("Precursor.Charge", integer, 2);
    add("Decoy", integer, 0);
    add("Q.Value", number, 0.001f);
    for (std::size_t i = 0; i < std::size(captured); ++i)
    {
      const std::string prefix = "Fr." + std::to_string(i);
      add(prefix + ".Index", integer, static_cast<std::int64_t>(i));
      add(prefix + ".Id", text, captured[i].id);
      add(prefix + ".Quantity", number, captured[i].quant);
      add(prefix + ".Score", number, captured[i].score);
    }
    auto output = arrow::io::FileOutputStream::Open(path.string());
    if (!output.ok()) { throw std::runtime_error(output.status().ToString()); }
    status(parquet::arrow::WriteTable(*arrow::Table::Make(arrow::schema(fields), arrays),
      arrow::default_memory_pool(), *output, 1));
    status((*output)->Close());
    RefineParams real_params; real_params.write_rt = false; real_params.write_intensity = true;
    real_params.q_global = real_params.q_protein = 1;
    RefineStats real_stats;
    const auto real_obs = LibraryRefiner::readObservations(path.string(), real_params, real_stats);
    check(real_obs.size() == 1 && real_stats.intensity_bad_tokens == 0,
          "real DIA-NN 2 fragment row loads through Parquet without bad tokens");
    const auto& fragments = real_obs.at(LibraryRefiner::key("PEPTIDEK", 2)).frags;
    check(fragments.size() == std::size(captured), "all 12 real fragment slots survive loading");
    for (const auto& expected : captured)
    {
      const auto found = std::find_if(fragments.begin(), fragments.end(), [&](const auto& f) {
        return f.type == expected.type && f.ordinal == expected.ordinal && f.charge == 1;
      });
      check(found != fragments.end() && std::abs(found->mz - expected.mz) < 1e-9 &&
            std::abs(found->quant - expected.quant) <= 1e-6f * std::abs(expected.quant) &&
            std::abs(found->correlation - expected.score) <= 1e-6f * std::abs(expected.score),
            std::string("real fragment preserves identity, m/z, quantity and score: ") + expected.id);
    }
    std::filesystem::remove(path);
  }
}

int main()
{
  if (arrow::GetBuildInfo().version_string != ARROW_VERSION_STRING)
  {
    std::cerr << "Arrow header/runtime mismatch: " << ARROW_VERSION_STRING << " / "
              << arrow::GetBuildInfo().version_string << '\n';
    return 1;
  }
  std::cerr << "intensity_match: parser\n";
  parser();
  std::cerr << "intensity_match: replacement and reranking\n";
  replaces_and_reranks();
  decoy_intensities_are_symmetric();
  std::cerr << "intensity_match: preservation and refusal guards\n";
  no_restrict_preserves_counts();
  mz_mismatch_throws();
  mixed_provenance_refused();
  mixed_rt_units_refused();
  duplicate_keys_count_once();
  loss_bearing_rank_and_preservation();
  std::cerr << "intensity_match: report layouts\n";
  report_layouts();
  if (failures) { std::cerr << failures << " failure(s)\n"; return 1; }
  std::cout << "intensity_match: ok\n";
  return 0;
}
