// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

// The evidence prefilter (search:candidates evidence) and the chunking:
//
//   1. the decoy fragments the prefilter indexes are exactly those
//      appendSearchDecoys builds, and those the searched set carries;
//   2. on a synthetic run with planted peaks, planted targets reach the depth
//      and are seen at their apex; null decoys do not;
//   3. label symmetry: indexing every target fragment as its decoy's and vice
//      versa swaps the evidence exactly and keeps the identical pair set;
//   4. pair-union: a pair is kept when either member passes, never otherwise;
//   5. the cap: exactly max_pairs pairs, pairs never split, and on a null
//      (members exchangeable) it changes the target:decoy ratio of passing
//      members by < 1 %; the ratio guard fires on a target-only selection;
//   6. determinism: the same evidence and the same set at 1 and 4 threads;
//   7. chunks never split a pair and each spreads over many isolation windows.

#include "synthetic_run.h"

#include <odia/search/AssayBuilder.h>
#include <odia/search/EvidencePrefilter.h>
#include <odia/search/Identifier.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <algorithm>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

using namespace ODIA::search;
namespace fs = std::filesystem;

namespace
{
  class Loader : public Identifier
  {
  public:
    explicit Loader(SearchParams p) : Identifier(std::move(p), [](const std::string&) {}, [](const std::string&) {}) {}
    using Identifier::loadRun;
  };

  void threads(int n)
  {
#ifdef _OPENMP
    omp_set_num_threads(n);
#else
    (void)n;
#endif
  }

  bool same(const MemberEvidence& a, const MemberEvidence& b)
  {
    return a.depth == b.depth && a.spectra == b.spectra && std::memcmp(&a.rt, &b.rt, sizeof(float)) == 0 &&
           std::memcmp(&a.intensity, &b.intensity, sizeof(float)) == 0;
  }

  bool sameEvidence(const std::vector<MemberEvidence>& a, const std::vector<MemberEvidence>& b)
  {
    if (a.size() != b.size()) { return false; }
    for (std::size_t m = 0; m < a.size(); ++m) { if (!same(a[m], b[m])) { return false; } }
    return true;
  }

  /// The m/z (as the prefilter stores them) of precursor @p i's top F
  /// transitions by library intensity, ties to the earlier transition.
  std::vector<float> topFragments(const ODIA::Library& lib, std::size_t i)
  {
    const auto& p = lib.precursors();
    const auto& t = lib.transitions();
    std::vector<std::uint32_t> order(p.transition_count[i]);
    for (std::uint32_t k = 0; k < order.size(); ++k) { order[k] = p.transition_begin[i] + k; }
    std::stable_sort(order.begin(), order.end(), [&](std::uint32_t a, std::uint32_t b) { return t.library_intensity[a] > t.library_intensity[b]; });
    std::vector<float> out;
    for (std::size_t k = 0; k < order.size() && k < SearchParams::prefilter_fragments; ++k)
    { out.push_back(static_cast<float>(ODIA::fromFixed(t.product_mz[order[k]]))); }
    return out;
  }

  std::vector<float> stored(const PrefilterPairs& u, std::size_t member)
  {
    const float* f = u.mz(member);
    return std::vector<float>(f, f + u.fragments[member / 2]);
  }
}

int main(int argc, char** argv)
{
  const fs::path dir = argc > 1 ? fs::path(argv[1]) : fs::temp_directory_path() / "dialibgen-identify-prefilter";
  fs::remove_all(dir);
  fs::create_directories(dir);

  // ---- 1. the indexed decoy fragments are the searched decoy fragments --------------
  {
    auto precursors = synth::peptides(600, 13);
    for (const char* s : {"EEEEEEK", "LLLLLLLLLLR", "ALILILK", "GIILLIR"})   // no decoy: unshufflable, isobaric
    {
      synth::Precursor p;
      p.sequence = s;
      p.protein = "ODD";
      p.rt = 40.0f;
      for (int z : {2, 3}) { p.charge = z; precursors.push_back(p); }
    }
    ODIA::Library lib = synth::library(precursors);
    // A slot no decoy can reproduce (an unknown loss): dropped from both assays.
    {
      auto& t = lib.transitions();
      t.loss[lib.precursors().transition_begin[0] + 1] = ODIA::LossType::Other;
    }
    for (const char* method : {"shuffle", "pseudo_reverse"})
    {
      SearchParams params;
      params.decoys = parseSearchDecoyMethod(method);
      params.max_pairs = 0;
      const DecoyRules rules = CandidateSelector::decoyRules(lib, params);
      // appendSearchDecoys on a copy holding the same targets.
      ODIA::Library copy = lib.subsetByIndex([&] {
        std::vector<std::size_t> all(lib.precursorCount());
        for (std::size_t i = 0; i < all.size(); ++i) { all[i] = i; }
        return all;
      }());
      const DecoyBuild built = appendSearchDecoys(copy, rules);
      std::size_t made = 0, compared = 0, mismatches = 0;
      for (std::size_t i = 0; i < lib.precursorCount(); ++i)
      {
        const DecoyAssay a = searchDecoy(lib, i, rules);
        CHECK(a.outcome == built.outcome[i]);
        if (a.outcome != DecoyOutcome::Made) { CHECK(a.slots.empty() && a.mz.empty()); continue; }
        CHECK(a.redrawn == built.redrawn[i]);
        const std::size_t d = lib.precursorCount() + made++;
        const auto& cp = copy.precursors();
        const auto& ct = copy.transitions();
        CHECK(cp.decoy[d] == 1 && cp.transition_count[d] == a.slots.size() && cp.transition_count[i] == a.slots.size());
        for (std::size_t k = 0; k < a.slots.size(); ++k)
        {
          ++compared;
          mismatches += ct.product_mz[cp.transition_begin[d] + k] != a.mz[k];                                   // decoy
          mismatches += ct.product_mz[cp.transition_begin[i] + k] != lib.transitions().product_mz[a.slots[k]];  // target
          mismatches += ct.charge[cp.transition_begin[d] + k] != a.charge[k];
        }
      }
      std::cout << method << ": " << made << " decoys, " << compared << " fragments compared with appendSearchDecoys, "
                << mismatches << " differ\n";
      CHECK(made == built.made && made > 1000 && mismatches == 0);
      CHECK(built.slots_dropped >= 1);

      // The prefilter's index, and the searched set built from the same targets.
      const PrefilterPairs universe = EvidencePrefilter::pairs(lib, params, {});
      CHECK(universe.size() == made);
      CHECK(universe.stats.no_decoy == lib.precursorCount() - made);
      std::vector<std::pair<std::uint64_t, std::size_t>> chosen(universe.targets);
      const SearchSet set = CandidateSelector::fromTargets(lib, params, chosen, universe.stats);
      CHECK(set.pairs() == made);
      std::map<std::size_t, std::size_t> pair_of_source;
      for (std::size_t k = 0; k < universe.size(); ++k) { pair_of_source[universe.targets[k].second] = k; }
      std::size_t checked = 0, differ = 0;
      for (std::size_t k = 0; k < set.pairs(); ++k)
      {
        const std::size_t u = pair_of_source.at(set.source[k]);
        differ += topFragments(set.library, k) != stored(universe, 2 * u);                       // target
        differ += topFragments(set.library, set.partnerOf(k)) != stored(universe, 2 * u + 1);    // decoy
        ++checked;
      }
      std::cout << "  prefilter index vs searched set: " << checked << " pairs, " << differ << " members differ\n";
      CHECK(differ == 0);
    }
  }

  // ---- the synthetic run ---------------------------------------------------------------
  synthrun::Spec spec;
  spec.peptides = 600;
  const synthrun::Fixture fx = synthrun::make(spec);
  const std::string mzml = (dir / "run.mzML").string();
  synthrun::writeMzML(spec, fx, mzml);
  SearchParams params;
  params.threads = 4;
  Loader loader(params);
  RunData run = loader.loadRun(mzml);
  const std::vector<IsolationWindow> windows = Identifier::isolationWindows(run);
  std::cout << "fixture: " << fx.library.precursorCount() << " precursors, " << fx.planted.size() << " planted, "
            << windows.size() << " windows\n";

  threads(4);
  const PrefilterPairs universe = EvidencePrefilter::pairs(fx.library, params, windows);
  SweepStats sw;
  const std::vector<MemberEvidence> evidence = EvidencePrefilter::sweep(universe, run.maps, params, &sw);
  std::cout << "sweep: " << sw.spectra << " spectra, " << sw.maps << " maps, " << sw.entries << " index entries, "
            << sw.peaks << " peaks\n";
  CHECK(sw.maps == static_cast<std::size_t>(spec.windows) && sw.windows == sw.maps);
  CHECK(evidence.size() == 2 * universe.size());

  // ---- 2. planted targets pass and are seen at their apex; null decoys do not --------
  {
    std::map<std::size_t, double> apex;
    for (const auto& p : fx.planted) { apex[p.index] = p.apex_s; }
    std::size_t planted = 0, passed = 0, at_apex = 0, null_decoys = 0, null_decoys_passing = 0, null_targets_passing = 0;
    for (std::size_t k = 0; k < universe.size(); ++k)
    {
      const auto it = apex.find(universe.targets[k].second);
      const MemberEvidence& t = evidence[2 * k];
      const MemberEvidence& d = evidence[2 * k + 1];
      if (it != apex.end())
      {
        ++planted;
        passed += t.depth >= std::min<int>(params.prefilter_depth, universe.fragments[k]) ? 1 : 0;
        at_apex += std::fabs(t.rt - it->second) <= 2.0 * spec.peak_sigma_s ? 1 : 0;
      }
      else
      {
        ++null_decoys;
        null_decoys_passing += d.depth >= params.prefilter_depth ? 1 : 0;
        null_targets_passing += t.depth >= params.prefilter_depth ? 1 : 0;
      }
    }
    std::cout << "planted: " << passed << " of " << planted << " reach the depth, " << at_apex << " seen within "
              << 2.0 * spec.peak_sigma_s << " s of the apex; unplanted pairs: " << null_targets_passing << " targets and "
              << null_decoys_passing << " decoys of " << null_decoys << " pass\n";
    CHECK(planted >= 300);
    CHECK(static_cast<double>(passed) >= 0.95 * static_cast<double>(planted));
    CHECK(static_cast<double>(at_apex) >= 0.95 * static_cast<double>(planted));
    // The fixture's background gives every member its own interference traces
    // near its expected RT (synthetic_run.h), so null members of both classes
    // reach the depth by chance -- equally often, since the two classes see
    // the same process.
    const double spread = 4.0 * std::sqrt(static_cast<double>(null_targets_passing + null_decoys_passing)) + 5.0;
    CHECK(std::fabs(static_cast<double>(null_targets_passing) - static_cast<double>(null_decoys_passing)) <= spread);
    // Seeds: targets that pass themselves, one per peptide, best first: the
    // strongest evidence is planted signal.
    const RtScale scale = CandidateSelector::rtScale(fx.library);
    const auto seeds = EvidencePrefilter::seeds(fx.library, universe, evidence, params, scale);
    std::set<std::string> peptides;
    std::size_t planted_top = 0;
    for (std::size_t k = 0; k < seeds.size(); ++k)
    {
      const auto& s = seeds[k];
      CHECK(s.depth >= params.prefilter_depth);
      CHECK(peptides.insert(std::string(fx.library.strings().get(fx.library.precursors().modified_sequence[s.index]))).second);
      if (k < 100) { planted_top += apex.count(s.index); }
      if (k > 0) { CHECK(s.depth <= seeds[k - 1].depth); }
    }
    std::cout << "seeds: " << seeds.size() << ", " << planted_top << " of the best 100 planted\n";
    CHECK(seeds.size() >= 150 && planted_top >= 95);
  }

  // ---- 3. label symmetry ---------------------------------------------------------------
  {
    const std::vector<MemberEvidence> swapped = EvidencePrefilter::sweep(universe, run.maps, params, nullptr, true);
    std::size_t differ = 0, asymmetric = 0;
    for (std::size_t k = 0; k < universe.size(); ++k)
    {
      differ += !same(swapped[2 * k], evidence[2 * k + 1]) || !same(swapped[2 * k + 1], evidence[2 * k]);
      asymmetric += evidence[2 * k].depth != evidence[2 * k + 1].depth;
    }
    std::cout << "label swap: " << differ << " pairs differ from the exact swap (" << asymmetric
              << " pairs whose members differ in depth)\n";
    CHECK(differ == 0);
    CHECK(asymmetric > 200);   // the swap is a real test: most planted pairs are lopsided
    for (const std::size_t cap : {std::size_t{0}, std::size_t{150}})
    {
      SearchParams p = params;
      p.max_pairs = cap;
      const PrefilterSelection a = EvidencePrefilter::choose(universe, evidence, p, windows);
      const PrefilterSelection b = EvidencePrefilter::choose(universe, swapped, p, windows);
      std::cout << "  cap " << cap << ": " << a.kept.size() << " pairs kept, swapped labels " << b.kept.size() << "\n";
      CHECK(a.kept == b.kept);
      CHECK(a.targets_passing == b.decoys_passing && a.decoys_passing == b.targets_passing);
      CHECK(!a.kept.empty());
    }
  }

  // ---- 4. pair-union --------------------------------------------------------------------
  {
    PrefilterPairs u;
    std::vector<MemberEvidence> e;
    auto add = [&](int target_depth, int decoy_depth) {
      const std::size_t k = u.targets.size();
      u.targets.emplace_back(1000 + k, k);
      u.fragments.push_back(6);
      u.fragment_mz.insert(u.fragment_mz.end(), 2 * SearchParams::prefilter_fragments, 500.0f);
      u.precursor_mz.push_back(500.0);
      u.library_rt.push_back(static_cast<float>(k));
      u.charge.push_back(2);
      MemberEvidence t, d;
      t.depth = static_cast<std::uint8_t>(target_depth);
      d.depth = static_cast<std::uint8_t>(decoy_depth);
      t.spectra = target_depth >= 5 ? 1 : 0;
      d.spectra = decoy_depth >= 5 ? 1 : 0;
      e.push_back(t);
      e.push_back(d);
    };
    add(6, 0);   // target only
    add(0, 6);   // decoy only
    add(2, 2);   // neither
    add(5, 5);   // both
    add(4, 3);   // neither
    SearchParams p;
    p.prefilter_depth = 5;
    const PrefilterSelection s = EvidencePrefilter::choose(u, e, p, {});
    CHECK((s.kept == std::vector<std::size_t>{0, 1, 3}));
    CHECK(s.targets_passing == 2 && s.decoys_passing == 2 && s.both_passing == 1 && s.union_pairs == 3);
    CHECK(s.depth_targets[6] == 1 && s.depth_decoys[6] == 1 && s.depth_targets[0] == 1);
  }

  // ---- 5. the cap, on a null of exchangeable members; the ratio guard ------------------------
  {
    PrefilterPairs u;
    std::vector<MemberEvidence> e;
    synthrun::Rng r(77);
    const std::size_t n = 2000000;
    std::vector<IsolationWindow> w;
    for (int k = 0; k < 40; ++k) { w.push_back({400.0 + 15.0 * k, 415.0 + 15.0 * k}); }
    for (std::size_t k = 0; k < n; ++k)
    {
      u.targets.emplace_back(r.next(), k);
      u.fragments.push_back(6);
      u.precursor_mz.push_back(400.5 + 599.0 * r.uniform());
      u.library_rt.push_back(static_cast<float>(100.0 * r.uniform()));
      u.charge.push_back(static_cast<std::uint8_t>(2 + r.next() % 3));
      for (int m = 0; m < 2; ++m)
      {
        MemberEvidence x;
        const double v = r.uniform();
        x.depth = static_cast<std::uint8_t>(v < 0.55 ? 2 : v < 0.75 ? 4 : v < 0.9 ? 5 : 6);
        x.spectra = x.depth >= 5 ? static_cast<std::uint32_t>(1 + r.next() % 12) : 0;
        e.push_back(x);
      }
    }
    u.fragment_mz.assign(2 * n * SearchParams::prefilter_fragments, 500.0f);
    SearchParams p;
    p.max_pairs = 300000;
    const PrefilterSelection s = EvidencePrefilter::choose(u, e, p, w);
    const double before = static_cast<double>(s.targets_passing) / static_cast<double>(s.decoys_passing);
    const double after = static_cast<double>(s.kept_targets_passing) / static_cast<double>(s.kept_decoys_passing);
    std::cout << "null cap: " << s.union_pairs << " union pairs -> " << s.kept.size() << " over " << s.strata
              << " strata; passing target:decoy " << before << " before, " << after << " after ("
              << 100.0 * std::fabs(after / before - 1.0) << " % change)\n";
    CHECK(s.kept.size() == p.max_pairs && s.capped == s.union_pairs - p.max_pairs);
    CHECK(std::fabs(after / before - 1.0) < 0.01);
    CHECK(std::is_sorted(s.kept.begin(), s.kept.end()));
    // Every library-RT decile keeps its share of the cap.
    std::vector<std::size_t> decile(10, 0);
    for (const std::size_t k : s.kept) { ++decile[std::min<std::size_t>(9, static_cast<std::size_t>(u.library_rt[k] / 10.0f))]; }
    for (const std::size_t c : decile) { CHECK(c > 0.09 * p.max_pairs && c < 0.11 * p.max_pairs); }
    // The ratio guard: whole pairs pass, a target-only (or decoy-heavy) selection does not.
    bool fired = false;
    EvidencePrefilter::checkRatio(s.kept.size(), s.kept.size());
    try { EvidencePrefilter::checkRatio(s.kept.size(), 0); } catch (const SearchAbort&) { fired = true; }
    CHECK(fired);
    fired = false;
    try { EvidencePrefilter::checkRatio(1000, 1300); } catch (const SearchAbort&) { fired = true; }
    CHECK(fired);
  }

  // ---- 6. determinism across threads ---------------------------------------------------
  {
    threads(1);
    const PrefilterPairs one = EvidencePrefilter::pairs(fx.library, params, windows);
    const std::vector<MemberEvidence> e1 = EvidencePrefilter::sweep(one, run.maps, params);
    const SearchSet s1 = EvidencePrefilter::select(fx.library, params, windows, run.maps, nullptr);
    threads(4);
    const SearchSet s4 = EvidencePrefilter::select(fx.library, params, windows, run.maps, nullptr);
    CHECK(one.targets == universe.targets && one.fragment_mz.size() == universe.fragment_mz.size() &&
          std::memcmp(one.fragment_mz.data(), universe.fragment_mz.data(), one.fragment_mz.size() * sizeof(float)) == 0);
    CHECK(sameEvidence(e1, evidence));
    CHECK(s1.source == s4.source && s1.draw == s4.draw && s1.prefilter_json == s4.prefilter_json);
    CHECK(s1.library.transitions().product_mz == s4.library.transitions().product_mz);
    CHECK(s1.seeds.size() == s4.seeds.size());
    for (std::size_t k = 0; k < std::min(s1.seeds.size(), s4.seeds.size()); ++k)
    { CHECK(s1.seeds[k].index == s4.seeds[k].index && s1.seeds[k].rt_s == s4.seeds[k].rt_s); }
    std::cout << "threads 1 vs 4: " << s1.pairs() << " pairs, " << s1.seeds.size() << " seeds, identical\n";
    CHECK(s1.pairs() > 0);
    // The entrapment estimate's database ratio comes from the pairs BEFORE the
    // prefilter, which keeps present targets preferentially.
    CHECK(!s1.entrapment_db_universe);
    SearchParams tagged = params;
    tagged.entrapment_tag = "SYNPROT_1";
    const SearchSet st = EvidencePrefilter::select(fx.library, tagged, windows, run.maps, nullptr);
    std::size_t trap = 0;
    for (const auto& t : universe.targets)
    { trap += std::string(fx.library.strings().get(fx.library.precursors().protein_group[t.second])).rfind("SYNPROT_1", 0) == 0; }
    std::cout << "entrapment base: " << st.entrapment_db_trap << " tagged and " << st.entrapment_db_real << " untagged of "
              << universe.size() << " pairs before the prefilter (" << st.pairs() << " searched)\n";
    CHECK(st.entrapment_db_universe && st.entrapment_db_trap == trap && trap > 0);
    CHECK(st.entrapment_db_real + st.entrapment_db_trap == universe.size());
    CHECK(st.source == s1.source);   // the tag changes nothing that is searched
    // Every searched pair is whole, and inside a window.
    for (std::size_t k = 0; k < s1.pairs(); ++k)
    {
      CHECK(s1.library.precursors().decoy[k] == 0 && s1.library.precursors().decoy[s1.partnerOf(k)] == 1);
      CHECK(CandidateSelector::inWindow(ODIA::fromFixed(s1.library.precursors().mz[k]), windows));
    }
  }
  run.maps.clear();

  // ---- 7. chunks: whole pairs, spread over the windows ------------------------------------
  {
    const ODIA::Library lib = synth::library(synth::peptides(3000, 5));
    SearchParams p;
    p.candidates = "random";
    p.max_pairs = 0;
    std::vector<IsolationWindow> w;
    for (int k = 0; k < 20; ++k) { w.push_back({350.0 + 32.5 * k, 350.0 + 32.5 * (k + 1)}); }
    const SearchSet set = CandidateSelector::select(lib, p, w);
    auto windowOf = [&](std::size_t i) {
      const double mz = ODIA::fromFixed(set.library.precursors().mz[i]);
      for (std::size_t k = 0; k < w.size(); ++k) { if (w[k].lower < mz && mz < w[k].upper) { return k; } }
      return w.size();
    };
    const std::size_t occupied = [&] {
      std::set<std::size_t> s;
      for (std::size_t k = 0; k < set.pairs(); ++k) { s.insert(windowOf(k)); }
      return s.size();
    }();
    for (const auto& [chunk, batch] : std::vector<std::pair<std::size_t, std::size_t>>{{1000, 100}, {1000, 500}, {101, 500}, {20000, 500}})
    {
      const auto parts = AssayBuilder::chunks(set, chunk, batch);
      std::vector<char> seen(set.size(), 0);
      std::size_t min_spread = w.size();
      for (const auto& part : parts)
      {
        CHECK(!part.empty() && part.size() <= chunk && part.size() % 2 == 0);
        const std::size_t half = part.size() / 2;
        std::set<std::size_t> spread;
        for (std::size_t j = 0; j < half; ++j)
        {
          CHECK(!set.isDecoy(part[j]) && part[half + j] == set.partnerOf(part[j]));   // targets, then their decoys
          CHECK(j == 0 || part[j - 1] < part[j]);
          spread.insert(windowOf(part[j]));
        }
        for (const std::size_t i : part) { CHECK(seen[i] == 0); seen[i] = 1; }
        min_spread = std::min(min_spread, spread.size());
      }
      CHECK(std::count(seen.begin(), seen.end(), 1) == static_cast<std::ptrdiff_t>(set.size()));
      std::cout << "chunks of " << chunk << " (batch " << batch << "): " << parts.size() << " chunks, each over at least "
                << min_spread << " of " << occupied << " occupied windows\n";
      if (chunk == 1000 && batch == 100) { CHECK(min_spread >= 8); }   // 10 pieces of 50 pairs, dealt across m/z
      if (chunk == 20000) { CHECK(parts.size() == 1 && min_spread == occupied); }
    }
  }

  fs::remove_all(dir);
  if (synth::failures) { std::cerr << synth::failures << " check(s) failed\n"; return 1; }
  std::cout << "identify_prefilter: PASS\n";
  return 0;
}
