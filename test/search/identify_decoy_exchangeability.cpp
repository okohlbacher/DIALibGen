// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

// Search decoys must be exchangeable with null (absent) targets: a known-null
// pair must be won by its target and by its decoy equally often. This is the
// test of search:intensities predicted against search:intensities library.
//
//   * the rule itself: PredictedAssays::rank orders a member's fragments by
//     its own prediction; pairCount is a symmetric function of the pair;
//   * the prefilter: every member's indexed fragments are its OWN top
//     fragments by the model (library: the decoy re-uses its target's slots);
//   * exchangeability on a synthetic run whose spectra hold the predicted top
//     fragments of "real" peptides (and noise), searched with a library of
//     random null peptides predicted by the same toy model (toy_fragment_model.h:
//     Pro-directed y ions are the most intense). The target of a null pair
//     then carries fragments chosen because they are the ones real peptides
//     fragment into; with library intensities its decoy puts other
//     compositions into those slots and loses: the null pairs' prefilter
//     evidence is target-favoured (z > 3 here, as on the acceptance runs:
//     81 and 79). With predicted intensities the decoy picks its own and the
//     sign test is balanced (|z| < 3);
//   * both members predicted even where the library already holds the model's
//     own assays, and a library RICHER than its own prediction (a build with
//     its own fragment cap) giving its targets no advantage: section 3b, the
//     library-assay shortcut this branch removed;
//   * the entrapment tag rule on UniProt-style ids (sp|ENTRAP_...|...).

#include "synthetic_library.h"
#include "toy_fragment_model.h"

#include <odia/search/EvidencePrefilter.h>
#include <odia/search/PredictedAssays.h>
#include <odia/search/Scoring.h>
#include <odia/search/SearchDecoys.h>

#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMS.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

using namespace ODIA::search;

namespace
{
  struct Rng
  {
    std::uint64_t state;
    std::uint64_t next()
    {
      std::uint64_t z = (state += 0x9E3779B97F4A7C15ull);
      z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
      z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
      return z ^ (z >> 31);
    }
    double uniform() { return static_cast<double>(next() >> 11) / 9007199254740992.0; }
  };

  /// A tryptic-looking peptide; P three times as frequent as the other
  /// residues, so that many peptides carry a short Pro-directed y ion
  /// (y2 = PK, PR) whose m/z recurs across the "proteome". Residues are
  /// drawn independently (no rule against repeats): a shuffle of such a
  /// peptide is then distributed exactly as another one, which is what a
  /// null target is, so any asymmetry the test sees comes from the search.
  std::string peptide(Rng& rng)
  {
    static const std::string aa = "ADEFGHILMNPPPQSTVWY";
    const std::size_t len = 8 + rng.next() % 9;
    std::string s;
    while (s.size() + 1 < len) { s.push_back(aa[rng.next() % aa.size()]); }
    s.push_back(rng.next() % 2 ? 'K' : 'R');
    return s;
  }

  /// A member's whole fragment universe by @p model, in the search's rank
  /// order; @p above_floor receives how many of them are above the floor.
  std::vector<PredictedFragment> ranked(FragmentModel& model, const std::string& sequence, int charge, std::size_t& above_floor,
                                        double lo = 200.0, double hi = 1800.0)
  {
    const OpenMS::AASequence seq = OpenMS::AASequence::fromString(sequence);
    const auto spectra = model.predict({seq}, {charge});
    std::vector<PredictedFragment> out;
    above_floor = 0;
    PredictedAssays::rank(seq, charge, spectra[0], lo, hi, out, above_floor);
    return out;
  }

  /// A member's top @p n fragments by the model, the rule of the search.
  std::vector<PredictedFragment> top(FragmentModel& model, const std::string& sequence, int charge, std::size_t n,
                                     double lo = 200.0, double hi = 1800.0)
  {
    std::size_t above = 0;
    std::vector<PredictedFragment> out = ranked(model, sequence, charge, above, lo, hi);
    if (out.size() > n) { out.resize(n); }
    return out;
  }

  /// A library of @p sequences at charge 2, each with its top @p fragments
  /// fragments by the model: what a generator would write with that model and
  /// that fragment cap.
  ODIA::Library library(FragmentModel& model, const std::vector<std::string>& sequences, const std::string& protein_prefix,
                        std::size_t fragments = 6)
  {
    ODIA::Library lib;
    auto& pre = lib.precursors();
    auto& tr = lib.transitions();
    for (std::size_t k = 0; k < sequences.size(); ++k)
    {
      const auto assay = top(model, sequences[k], 2, fragments);
      pre.mz.push_back(ODIA::toFixed(OpenMS::AASequence::fromString(sequences[k]).getMZ(2)));
      pre.irt.push_back(static_cast<float>(k % 100));
      pre.im.push_back(std::nanf(""));
      pre.ccs.push_back(std::nanf(""));
      pre.charge.push_back(2);
      pre.decoy.push_back(0);
      pre.modified_sequence.push_back(lib.strings().intern(sequences[k]));
      pre.protein_group.push_back(lib.strings().intern(protein_prefix + std::to_string(k / 3)));
      pre.transition_begin.push_back(static_cast<std::uint32_t>(tr.product_mz.size()));
      for (const auto& f : assay)
      {
        tr.product_mz.push_back(ODIA::toFixed(f.mz));
        tr.library_intensity.push_back(f.intensity);
        tr.type.push_back(f.type);
        tr.ordinal.push_back(f.ordinal);
        tr.charge.push_back(f.charge);
        tr.loss.push_back(ODIA::LossType::None);
      }
      pre.transition_count.push_back(static_cast<std::uint32_t>(assay.size()));
    }
    return lib;
  }

  /// A synthetic run: 300 MS2 spectra, each holding the indexed fragments of a
  /// dozen "real" peptides (the ones the library's null targets are NOT) and
  /// 300 noise peaks.
  std::vector<OpenSwath::SwathMap> syntheticRun(const std::vector<std::vector<PredictedFragment>>& real_top, Rng& rng)
  {
    auto experiment = std::make_shared<OpenMS::PeakMap>();
    for (int s = 0; s < 300; ++s)
    {
      OpenMS::MSSpectrum spectrum;
      spectrum.setMSLevel(2);
      spectrum.setRT(static_cast<double>(s));
      std::vector<std::pair<double, double>> peaks;
      for (int r = 0; r < 12; ++r)
      {
        for (const auto& f : real_top[rng.next() % real_top.size()]) { peaks.emplace_back(f.mz, 1000.0); }
      }
      for (int q = 0; q < 300; ++q) { peaks.emplace_back(200.0 + 1300.0 * rng.uniform(), 10.0 + 490.0 * rng.uniform()); }
      std::sort(peaks.begin(), peaks.end());
      for (const auto& [mz, intensity] : peaks)
      {
        OpenMS::Peak1D p;
        p.setMZ(mz);
        p.setIntensity(static_cast<float>(intensity));
        spectrum.push_back(p);
      }
      experiment->addSpectrum(spectrum);
    }
    OpenSwath::SwathMap map(250.0, 2000.0, 1125.0, false);
    map.sptr = std::make_shared<OpenMS::SpectrumAccessOpenMS>(experiment);
    return {map};
  }

  struct Sign
  {
    std::size_t target = 0, decoy = 0;
    double z() const { return (target + decoy) > 0 ? (static_cast<double>(target) - static_cast<double>(decoy)) / std::sqrt(static_cast<double>(target + decoy)) : 0.0; }
  };
}

int main()
{
  toy::ProlineModel model;

  // ---- 1. the rule on one member, and the pair's count ----------------------------------
  {
    // AGPLDEFK: the y ion starting at P (y6, PLDEFK) leads the ranking.
    const auto ranked = top(model, "AGPLDEFK", 2, 100);
    CHECK(!ranked.empty() && ranked[0].type == ODIA::FragmentType::Y && ranked[0].ordinal == 6 && ranked[0].charge == 1);
    for (std::size_t k = 1; k < ranked.size(); ++k) { CHECK(ranked[k - 1].intensity >= ranked[k].intensity); }
    for (const auto& f : ranked) { CHECK(f.mz >= 200.0 && f.mz <= 1800.0); }
    CHECK(PredictedAssays::pairCount(12, 10, 8, 20, 20, 3) == 8);
    CHECK(PredictedAssays::pairCount(6, 10, 10, 20, 20, 3) == 6);
    CHECK(PredictedAssays::pairCount(12, 1, 10, 20, 20, 3) == 3);
    CHECK(PredictedAssays::pairCount(12, 1, 1, 2, 20, 3) == 0);
    for (std::size_t a : {0u, 2u, 5u, 9u})
    {
      for (std::size_t b : {0u, 3u, 7u, 12u})
      { CHECK(PredictedAssays::pairCount(9, a, b, 11, 14, 3) == PredictedAssays::pairCount(9, b, a, 14, 11, 3)); }
    }
  }

  // ---- 2. a null library, a synthetic run of other ("real") peptides -------------------
  Rng rng{20260922};
  std::vector<std::string> nulls, reals;
  for (std::size_t k = 0; k < 12000; ++k) { nulls.push_back(peptide(rng)); }
  for (std::size_t k = 0; k < 400; ++k) { reals.push_back(peptide(rng)); }
  const ODIA::Library lib = library(model, nulls, "sp|ENTRAP_P");
  std::vector<std::vector<PredictedFragment>> real_top;
  for (const auto& s : reals) { real_top.push_back(top(model, s, 2, 6)); }

  const std::vector<OpenSwath::SwathMap> maps = syntheticRun(real_top, rng);
  const std::vector<IsolationWindow> windows = {{250.0, 2000.0}};

  // ---- 3. library vs predicted: the null pairs' prefilter evidence ----------------------
  Sign by_mode[2];
  for (const Intensities mode : {Intensities::Library, Intensities::Predicted})
  {
    SearchParams p;
    p.intensities = mode;
    p.max_pairs = 0;
    p.prefilter_depth = 1;   // spectra = spectra with any indexed fragment matched
    p.threads = 2;
    const bool predicted = mode == Intensities::Predicted;
    const std::size_t seen_before = model.peptides_seen;
    const PrefilterPairs universe = EvidencePrefilter::pairs(lib, p, windows, predicted ? &model : nullptr);
    if (predicted)
    {
      // Both members of every pair are predicted, on this library too: it
      // carries the toy model's own top six, which changes nothing about the
      // rule and so must not change what the search predicts.
      const std::size_t predicted_now = model.peptides_seen - seen_before;
      std::cout << "predicted " << predicted_now << " peptides for " << universe.size() << " pairs\n";
      CHECK(predicted_now >= 2 * universe.size() && predicted_now <= 2 * universe.processed);
    }
    const std::vector<MemberEvidence> evidence = EvidencePrefilter::sweep(universe, maps, p);
    Sign& sign = by_mode[predicted ? 1 : 0];
    auto strength = [](const MemberEvidence& e) { return (static_cast<std::uint64_t>(e.depth) << 32) | e.spectra; };
    std::uint64_t spectra_t = 0, spectra_d = 0;
    for (std::size_t k = 0; k < universe.size(); ++k)
    {
      const auto t = strength(evidence[2 * k]), d = strength(evidence[2 * k + 1]);
      if (t > d) { ++sign.target; }
      else if (d > t) { ++sign.decoy; }
      spectra_t += evidence[2 * k].spectra;
      spectra_d += evidence[2 * k + 1].spectra;
    }
    std::cout << toString(mode) << ": " << universe.size() << " null pairs; stronger target " << sign.target << ", stronger decoy "
              << sign.decoy << " (z " << sign.z() << "); spectra matched " << spectra_t << " by targets, " << spectra_d << " by decoys\n";

    // Every member's indexed fragments are its own top fragments by the
    // model (predicted), or the target's slots (library).
    const DecoyRules rules = CandidateSelector::decoyRules(lib, p);
    const double lo = ODIA::fromFixed(rules.fragment_min), hi = ODIA::fromFixed(rules.fragment_max);
    std::size_t own = 0, checked = 0, same_slots = 0;
    for (std::size_t k = 0; k < universe.size(); k += 97)
    {
      const std::size_t i = universe.targets[k].second;
      const DecoyAssay a = searchDecoy(lib, i, rules);
      CHECK(a.outcome == DecoyOutcome::Made && !a.sequence.empty());
      const auto t = top(model, nulls[i], 2, SearchParams::prefilter_fragments, lo, hi);
      const auto d = top(model, a.sequence, 2, SearchParams::prefilter_fragments, lo, hi);
      const float* tm = universe.mz(2 * k);
      const float* dm = universe.mz(2 * k + 1);
      bool target_own = true, decoy_own = true;
      for (std::size_t j = 0; j < universe.fragments[k]; ++j)
      {
        target_own = target_own && std::fabs(tm[j] - static_cast<float>(t[j].mz)) < 1e-3f;
        decoy_own = decoy_own && std::fabs(dm[j] - static_cast<float>(d[j].mz)) < 1e-3f;
      }
      CHECK(target_own);   // the library holds the toy model's top six: both modes index them
      own += decoy_own ? 1 : 0;
      ++checked;
      // Library: the decoy's j-th indexed fragment is the target's j-th slot recomputed.
      same_slots += (!decoy_own) ? 1 : 0;
    }
    if (predicted) { CHECK(own == checked); }
    else { CHECK(same_slots > checked / 2); }
  }
  // The "before" state fails the sign test, the fix passes it.
  CHECK(by_mode[0].z() > 3.0);
  CHECK(std::fabs(by_mode[1].z()) < 3.0);

  // ---- 3b. a library richer than the model's own prediction ------------------------------
  //
  // The regression (docs/design/built-in-identification.md, "The
  // library-assay shortcut, removed"): a library written by another build
  // holds up to its own fragment cap per target -- MORE transitions than the
  // model itself
  // puts above the floor for that target. Letting such targets keep their
  // library assays, and predicting only the decoys, reads the TARGET's
  // library transition count where the pair's count rule must read the
  // target's own above-floor count, so the pair comes out with more fragments
  // than the rule allows, and the extra ones are the target's library's. On a
  // timsTOF library whose sampled targets ALL carried exactly the model's
  // prediction it changed 21 % of the pairs and the entrapment FDP went from
  // 0.74 % to 1.30 %. So: both members' own predictions, the rule's count.
  {
    toy::SteepProlineModel steep;
    // Twelve fragments per target, a generator's fragment cap; the steep
    // model puts fewer than that above the floor for most of these peptides.
    const ODIA::Library rich = library(steep, nulls, "sp|ENTRAP_P", 12);
    std::vector<std::vector<PredictedFragment>> steep_real;
    for (const auto& s : reals) { steep_real.push_back(top(steep, s, 2, SearchParams::prefilter_fragments)); }
    Rng run_rng{20260923};
    const std::vector<OpenSwath::SwathMap> steep_maps = syntheticRun(steep_real, run_rng);

    SearchParams p;
    p.intensities = Intensities::Predicted;
    p.max_pairs = 0;
    p.prefilter_depth = 1;
    p.threads = 2;
    const PrefilterPairs universe = EvidencePrefilter::pairs(rich, p, windows, &steep);
    const DecoyRules rules = CandidateSelector::decoyRules(rich, p);
    const double lo = ODIA::fromFixed(rules.fragment_min), hi = ODIA::fromFixed(rules.fragment_max);
    const auto& pre = rich.precursors();

    std::size_t checked = 0, by_rule = 0, own_assays = 0, richer = 0, kept_differs = 0;
    for (std::size_t k = 0; k < universe.size(); k += 7)
    {
      const std::size_t i = universe.targets[k].second;
      const DecoyAssay a = searchDecoy(rich, i, rules);
      CHECK(a.outcome == DecoyOutcome::Made && !a.sequence.empty());
      std::size_t above_t = 0, above_d = 0;
      const auto t = ranked(steep, nulls[i], 2, above_t, lo, hi);
      const auto d = ranked(steep, a.sequence, 2, above_d, lo, hi);
      const std::size_t lib_count = pre.transition_count[i];
      // The rule, read from the two members' OWN predictions.
      const std::size_t want =
        PredictedAssays::pairCount(lib_count, above_t, above_d, t.size(), d.size(), SearchParams::min_assay_fragments);
      // What keeping the target's library assay gives instead: the whole
      // library assay counts as the target's above-floor fragments.
      const std::size_t kept =
        PredictedAssays::pairCount(lib_count, lib_count, above_d, lib_count, d.size(), SearchParams::min_assay_fragments);
      richer += above_t < lib_count ? 1 : 0;
      kept_differs += kept != want ? 1 : 0;
      by_rule += universe.assay_fragments[k] == want ? 1 : 0;
      // Both members' indexed fragments are the top of their own prediction.
      const float* tm = universe.mz(2 * k);
      const float* dm = universe.mz(2 * k + 1);
      bool own = universe.fragments[k] == std::min<std::size_t>(SearchParams::prefilter_fragments, want);
      for (std::size_t j = 0; own && j < universe.fragments[k]; ++j)
      {
        own = std::fabs(tm[j] - static_cast<float>(t[j].mz)) < 1e-3f && std::fabs(dm[j] - static_cast<float>(d[j].mz)) < 1e-3f;
      }
      own_assays += own ? 1 : 0;
      ++checked;
    }
    std::cout << "rich library: " << checked << " of " << universe.size() << " pairs checked; " << richer
              << " whose target's library assay holds more transitions than its own prediction has above the floor; "
              << kept_differs << " whose count keeping it would change; " << by_rule << " with the rule's count, " << own_assays
              << " with both members' own assays\n";
    CHECK(checked > 200);
    CHECK(richer > checked / 5);          // the library really is the longer list
    CHECK(kept_differs > checked / 10);   // and keeping it really would change the count
    CHECK(by_rule == checked);
    CHECK(own_assays == checked);

    // ... and the null pairs' prefilter evidence stays balanced.
    const std::vector<MemberEvidence> evidence = EvidencePrefilter::sweep(universe, steep_maps, p);
    auto strength = [](const MemberEvidence& e) { return (static_cast<std::uint64_t>(e.depth) << 32) | e.spectra; };
    Sign sign;
    for (std::size_t k = 0; k < universe.size(); ++k)
    {
      const auto t = strength(evidence[2 * k]), d = strength(evidence[2 * k + 1]);
      if (t > d) { ++sign.target; }
      else if (d > t) { ++sign.decoy; }
    }
    std::cout << "rich library: " << universe.size() << " null pairs; stronger target " << sign.target << ", stronger decoy "
              << sign.decoy << " (z " << sign.z() << ")\n";
    CHECK(std::fabs(sign.z()) < 3.0);
  }

  // ---- 4. the searched set's assays are the predicted ones -------------------------------
  {
    SearchParams p;
    p.intensities = Intensities::Predicted;
    p.max_pairs = 300;
    p.candidates = "random";
    p.threads = 2;
    const SearchSet set = EvidencePrefilter::selectRandom(lib, p, windows, model);
    CHECK(set.pairs() == 300);
    CHECK(set.stats.drawn == set.stats.no_decoy + set.stats.capped + set.stats.pairs);
    const auto& pre = set.library.precursors();
    const auto& tr = set.library.transitions();
    const DecoyRules rules = CandidateSelector::decoyRules(lib, p);
    const double lo = ODIA::fromFixed(rules.fragment_min), hi = ODIA::fromFixed(rules.fragment_max);
    std::size_t matched = 0;
    for (std::size_t k = 0; k < set.pairs(); ++k)
    {
      const std::size_t decoy = set.partnerOf(k);
      CHECK(pre.transition_count[k] == pre.transition_count[decoy]);
      const DecoyAssay a = searchDecoy(lib, set.source[k], rules);
      const auto t = top(model, std::string(set.modifiedSequence(k)), 2, pre.transition_count[k], lo, hi);
      const auto d = top(model, a.sequence, 2, pre.transition_count[decoy], lo, hi);
      bool same = t.size() == pre.transition_count[k] && d.size() == pre.transition_count[decoy];
      for (std::uint32_t j = 0; same && j < pre.transition_count[k]; ++j)
      {
        same = std::fabs(ODIA::fromFixed(tr.product_mz[pre.transition_begin[k] + j]) - t[j].mz) < 1e-4 &&
               std::fabs(ODIA::fromFixed(tr.product_mz[pre.transition_begin[decoy] + j]) - d[j].mz) < 1e-4 &&
               tr.library_intensity[pre.transition_begin[decoy] + j] == d[j].intensity;
      }
      matched += same ? 1 : 0;
    }
    std::cout << "searched set: " << matched << " of " << set.pairs() << " pairs carry both members' own predicted assays\n";
    CHECK(matched == set.pairs());
    // Deterministic: the same set again.
    const SearchSet again = EvidencePrefilter::selectRandom(lib, p, windows, model);
    CHECK(again.library.transitions().product_mz == set.library.transitions().product_mz);
    CHECK(again.library.transitions().library_intensity == set.library.transitions().library_intensity);
  }

  // ---- 5. the entrapment tag: at the start of an id or after a '|' ----------------------
  {
    const std::string tag = "ENTRAP_";
    CHECK(entrapmentClass("ENTRAP_P12345", tag) == Entrapment::Trap);
    CHECK(entrapmentClass("sp|ENTRAP_Q96QE4|ENTRAP_LR37B_HUMAN", tag) == Entrapment::Trap);
    CHECK(entrapmentClass("tr|A0A024|ENTRAP_X_HUMAN", tag) == Entrapment::Trap);
    CHECK(entrapmentClass("sp|ENTRAP_Q1|A;sp|ENTRAP_Q2|B", tag) == Entrapment::Trap);
    CHECK(entrapmentClass("sp|ENTRAP_Q1|A;sp|P12345|B", tag) == Entrapment::Shared);
    CHECK(entrapmentClass("sp|P12345|ALBU_HUMAN", tag) == Entrapment::Real);
    CHECK(entrapmentClass("sp|P12345|XENTRAP_HUMAN", tag) == Entrapment::Real);   // not after a '|'
    CHECK(entrapmentClass("P12345", tag) == Entrapment::Real);
    CHECK(entrapmentClass("sp|ENTRAP_Q1|A", "") == Entrapment::Real);
  }

  if (synth::failures) { std::cerr << synth::failures << " check(s) failed\n"; return 1; }
  std::cout << "identify_decoy_exchangeability: PASS\n";
  return 0;
}
