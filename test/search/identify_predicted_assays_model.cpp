// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

// search:intensities predicted with the real PeptDeep MS2 model. BOTH members
// of every pair are predicted from their own sequences, whatever the library
// holds, and the pair's count is read from the two members' own predictions:
//
//   * every member's assay is its own sequence's prediction by the same model,
//     ranked by the same rule -- the decoy's from searchDecoy's sequence, and
//     the target's too;
//   * on a library DIALibGen generated with that very model
//     (LibraryGenerator::predictFragmentIntensities, same instrument and NCE)
//     the two agree fragment for fragment: the search's rule mirrors the
//     generator's;
//   * both members carry the same number of fragments, and that number is
//     PredictedAssays::pairCount of the two members' OWN above-floor counts.
//     Checked on a second library written PAST the search's intensity floor,
//     as another build's is: its targets hold more transitions than the model
//     puts above the floor for them, so a target that keeps its library assay
//     brings the library's count into the pair's rule and the pair comes out
//     with more fragments than the rule allows. That is the regression --
//     docs/design/built-in-identification.md, "The library-assay shortcut,
//     removed"; on the timsTOF library it moved 21 % of the pairs and the
//     entrapment FDP from 0.74 % to 1.30 %.
//
//   identify_predicted_assays_model <peptdeep_ms2_dynamic.onnx>

#include "synthetic_library.h"

#include <odia/LibraryGenerator.h>
#include <odia/search/EvidencePrefilter.h>
#include <odia/search/PredictedAssays.h>

#include <OpenMS/CHEMISTRY/AASequence.h>

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

using namespace ODIA::search;

namespace
{
  /// A copy of @p in whose every target carries the model's top @p fragments
  /// of its OWN ranked universe -- floor or no floor. A library another build
  /// wrote looks like this: it kept a fixed number of transitions per
  /// precursor, where this search counts only those above 1e-4 of the base
  /// peak (PredictedAssays::predicted_floor). DIALibGen's own generator cuts
  /// at that same floor, so its libraries never show the difference.
  ODIA::Library pastTheFloor(const ODIA::Library& in, FragmentModel& model, const DecoyRules& rules, std::size_t fragments)
  {
    const auto& ip = in.precursors();
    const double lo = ODIA::fromFixed(rules.fragment_min), hi = ODIA::fromFixed(rules.fragment_max);
    std::vector<OpenMS::AASequence> peptides;
    std::vector<int> charges;
    for (std::size_t i = 0; i < in.precursorCount(); ++i)
    {
      peptides.push_back(OpenMS::AASequence::fromString(std::string(in.strings().get(ip.modified_sequence[i]))));
      charges.push_back(static_cast<int>(ip.charge[i]));
    }
    const auto spectra = model.predict(peptides, charges);

    ODIA::Library out;
    auto& pre = out.precursors();
    auto& tr = out.transitions();
    std::vector<PredictedFragment> ranked;
    for (std::size_t i = 0; i < in.precursorCount(); ++i)
    {
      std::size_t above = 0;
      PredictedAssays::rank(peptides[i], charges[i], spectra[i], lo, hi, ranked, above);
      const std::size_t n = std::min(fragments, ranked.size());
      pre.mz.push_back(ip.mz[i]);
      pre.irt.push_back(ip.irt[i]);
      pre.im.push_back(ip.im[i]);
      pre.ccs.push_back(ip.ccs[i]);
      pre.charge.push_back(ip.charge[i]);
      pre.decoy.push_back(ip.decoy[i]);
      pre.modified_sequence.push_back(out.strings().intern(std::string(in.strings().get(ip.modified_sequence[i]))));
      pre.protein_group.push_back(out.strings().intern(std::string(in.strings().get(ip.protein_group[i]))));
      pre.transition_begin.push_back(static_cast<std::uint32_t>(tr.product_mz.size()));
      for (std::size_t j = 0; j < n; ++j)
      {
        tr.product_mz.push_back(ODIA::toFixed(ranked[j].mz));
        tr.library_intensity.push_back(ranked[j].intensity);
        tr.type.push_back(ranked[j].type);
        tr.ordinal.push_back(ranked[j].ordinal);
        tr.charge.push_back(static_cast<std::uint8_t>(ranked[j].charge));
        tr.loss.push_back(ODIA::LossType::None);
      }
      pre.transition_count.push_back(static_cast<std::uint32_t>(n));
    }
    return out;
  }

  struct Tally
  {
    std::size_t pairs = 0, targets_as_library = 0, targets_own = 0, decoys_own = 0, fragments_differing = 0;
    std::size_t counts_by_rule = 0, longer_than_own = 0, kept_would_differ = 0;
  };

  /// Search @p lib and check every pair of the searched set: the target's
  /// assay against its library assay, both members' against their own
  /// predictions, and the count against the pair's rule.
  Tally verify(const std::string& name, const ODIA::Library& lib, FragmentModel& model, const SearchParams& p)
  {
    const SearchSet set = EvidencePrefilter::selectRandom(lib, p, {}, model);
    std::cout << name << ": searched " << set.pairs() << " pairs of " << lib.precursorCount() << " library targets; no decoy "
              << set.stats.no_decoy << " (copy " << set.stats.decoy_copy << ", too few " << set.stats.decoy_too_few_fragments
              << ", unpredictable " << set.stats.decoy_unpredictable << ")\n";

    const DecoyRules rules = CandidateSelector::decoyRules(lib, p);
    const double lo = ODIA::fromFixed(rules.fragment_min), hi = ODIA::fromFixed(rules.fragment_max);
    const auto& lp = lib.precursors();
    const auto& lt = lib.transitions();
    const auto& sp = set.library.precursors();
    const auto& st = set.library.transitions();
    Tally t;
    t.pairs = set.pairs();
    for (std::size_t k = 0; k < set.pairs(); ++k)
    {
      const std::size_t decoy = set.partnerOf(k), i = set.source[k];
      const std::uint32_t n = sp.transition_count[k];
      CHECK(n == sp.transition_count[decoy] && n >= SearchParams::min_assay_fragments && n <= lp.transition_count[i]);
      // The target: its library assay, in library order, cut to n -- the same
      // model wrote both lists, so they agree.
      bool same = true;
      for (std::uint32_t j = 0; j < n; ++j)
      {
        const std::uint32_t a = sp.transition_begin[k] + j, b = lp.transition_begin[i] + j;
        const bool ion = st.type[a] == lt.type[b] && st.ordinal[a] == lt.ordinal[b] && st.charge[a] == lt.charge[b];
        const bool value =
          std::fabs(st.library_intensity[a] - lt.library_intensity[b]) <= 1e-4f * std::max(1.0f, lt.library_intensity[b]);
        same = same && ion && value;
        t.fragments_differing += (ion && value) ? 0 : 1;
      }
      t.targets_as_library += same ? 1 : 0;
      // Both members through the model, from their own sequences, one batch.
      const DecoyAssay a = searchDecoy(lib, i, rules);
      const OpenMS::AASequence target_seq = OpenMS::AASequence::fromString(std::string(lib.strings().get(lp.modified_sequence[i])));
      const OpenMS::AASequence decoy_seq = OpenMS::AASequence::fromString(a.sequence);
      const int z = static_cast<int>(sp.charge[k]);
      const auto spectra = model.predict({target_seq, decoy_seq}, {z, z});
      std::vector<PredictedFragment> target_ranked, decoy_ranked;
      std::size_t above_t = 0, above_d = 0;
      PredictedAssays::rank(target_seq, z, spectra[0], lo, hi, target_ranked, above_t);
      PredictedAssays::rank(decoy_seq, z, spectra[1], lo, hi, decoy_ranked, above_d);
      auto assayIsOwn = [&](std::size_t member, const std::vector<PredictedFragment>& ranked) {
        if (ranked.size() < n) { return false; }
        for (std::uint32_t j = 0; j < n; ++j)
        {
          const std::uint32_t d = sp.transition_begin[member] + j;
          if (st.type[d] != ranked[j].type || st.ordinal[d] != ranked[j].ordinal || st.charge[d] != ranked[j].charge ||
              std::fabs(ODIA::fromFixed(st.product_mz[d]) - ranked[j].mz) >= 1e-4)
          { return false; }
        }
        return true;
      };
      t.targets_own += assayIsOwn(k, target_ranked) ? 1 : 0;
      t.decoys_own += assayIsOwn(decoy, decoy_ranked) ? 1 : 0;
      // The count, from the two members' own above-floor counts -- and what
      // keeping the target's library assay would make of it, which puts the
      // library's transition count where the target's own count belongs.
      const std::size_t library_count = lp.transition_count[i];
      const std::size_t want = PredictedAssays::pairCount(library_count, above_t, above_d, target_ranked.size(),
                                                          decoy_ranked.size(), SearchParams::min_assay_fragments);
      const std::size_t kept = PredictedAssays::pairCount(library_count, library_count, above_d, library_count,
                                                          decoy_ranked.size(), SearchParams::min_assay_fragments);
      t.counts_by_rule += n == want ? 1 : 0;
      t.longer_than_own += above_t < library_count ? 1 : 0;
      t.kept_would_differ += kept != want ? 1 : 0;
    }
    std::cout << name << ": targets as in the library " << t.targets_as_library << " of " << t.pairs << " ("
              << t.fragments_differing << " fragments differ); by their own prediction: targets " << t.targets_own
              << ", decoys " << t.decoys_own << "; counts by the pair's rule " << t.counts_by_rule << " ("
              << t.longer_than_own << " targets whose library assay outlasts their own above-floor list, "
              << t.kept_would_differ << " pairs whose count keeping it would change)\n";
    // Batches of a different composition may move a prediction in its last
    // bits, which can swap two nearly tied fragments: allow a handful.
    const double all = 0.98 * static_cast<double>(t.pairs);
    CHECK(static_cast<double>(t.targets_as_library) >= all);
    CHECK(static_cast<double>(t.targets_own) >= all);
    CHECK(static_cast<double>(t.decoys_own) >= all);
    CHECK(static_cast<double>(t.counts_by_rule) >= all);
    return t;
  }
}

int main(int argc, char** argv)
{
  if (argc < 2) { std::cerr << "usage: identify_predicted_assays_model <peptdeep_ms2_dynamic.onnx>\n"; return 2; }
  const std::string model_path = argv[1];

  ODIA::Library lib = synth::library(synth::peptides(250, 7));
  ODIA::DigestParams digest;
  digest.fragment_mz_min = 200.0;
  digest.fragment_mz_max = 1800.0;
  digest.max_fragment_charge = 2;
  digest.min_fragments = 3;
  digest.max_fragments = 12;
  const std::size_t unpredicted = ODIA::LibraryGenerator::predictFragmentIntensities(lib, model_path, digest, 30.0f, "QE", false, 2);
  CHECK(unpredicted == 0);

  PeptDeepFragmentModel model(model_path, "QE", 30.0, 2);
  CHECK(model.describe() == "PeptDeep MS2, instrument QE, NCE 30");
  SearchParams p;
  p.intensities = Intensities::Predicted;
  p.candidates = "random";
  p.max_pairs = 0;
  p.threads = 2;

  // 1. DIALibGen's own library: the search's rule mirrors the generator's.
  const Tally own = verify("generated library", lib, model, p);
  CHECK(own.pairs > 400);

  // 2. A library written past the floor, as another build's is: the targets
  //    hold more transitions than the model puts above the floor for them.
  //    Both members are still predicted and the count is still the rule's --
  //    which is what the library-assay shortcut broke.
  const ODIA::Library rich = pastTheFloor(lib, model, CandidateSelector::decoyRules(lib, p), 12);
  const Tally past = verify("library past the floor", rich, model, p);
  CHECK(past.pairs > 400);
  CHECK(past.longer_than_own > past.pairs / 5);      // the library really is the longer list
  CHECK(past.kept_would_differ > past.pairs / 20);   // and keeping it really would change the count

  if (synth::failures) { std::cerr << synth::failures << " check(s) failed\n"; return 1; }
  std::cout << "identify_predicted_assays_model: PASS\n";
  return 0;
}
