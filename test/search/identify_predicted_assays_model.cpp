// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

// search:intensities predicted with the real PeptDeep MS2 model:
//
//   * the search's rule mirrors the generator's: on a library whose targets
//     LibraryGenerator::predictFragmentIntensities predicted (same model,
//     instrument and NCE), every searched target's assay is its library
//     assay, cut to the pair's count;
//   * every decoy's assay is its own sequence's prediction: the same model
//     on the decoy peptide (searchDecoy's sequence), ranked by the same rule;
//   * both members of a pair carry the same number of fragments.
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
  const SearchSet set = EvidencePrefilter::selectRandom(lib, p, {}, model);
  std::cout << "searched " << set.pairs() << " pairs of " << lib.precursorCount() << " library targets; no decoy "
            << set.stats.no_decoy << " (copy " << set.stats.decoy_copy << ", too few " << set.stats.decoy_too_few_fragments
            << ", unpredictable " << set.stats.decoy_unpredictable << ")\n";
  CHECK(set.pairs() > 400);

  const DecoyRules rules = CandidateSelector::decoyRules(lib, p);
  const double lo = ODIA::fromFixed(rules.fragment_min), hi = ODIA::fromFixed(rules.fragment_max);
  const auto& lp = lib.precursors();
  const auto& lt = lib.transitions();
  const auto& sp = set.library.precursors();
  const auto& st = set.library.transitions();
  std::size_t targets_as_library = 0, decoys_own = 0, fragments_differing = 0;
  for (std::size_t k = 0; k < set.pairs(); ++k)
  {
    const std::size_t decoy = set.partnerOf(k), i = set.source[k];
    const std::uint32_t n = sp.transition_count[k];
    CHECK(n == sp.transition_count[decoy] && n >= SearchParams::min_assay_fragments && n <= lp.transition_count[i]);
    // The target: its library assay, in library order, cut to n.
    bool same = true;
    for (std::uint32_t j = 0; j < n; ++j)
    {
      const std::uint32_t a = sp.transition_begin[k] + j, b = lp.transition_begin[i] + j;
      const bool ion = st.type[a] == lt.type[b] && st.ordinal[a] == lt.ordinal[b] && st.charge[a] == lt.charge[b];
      const bool value = std::fabs(st.library_intensity[a] - lt.library_intensity[b]) <= 1e-4f * std::max(1.0f, lt.library_intensity[b]);
      same = same && ion && value;
      fragments_differing += (ion && value) ? 0 : 1;
    }
    targets_as_library += same ? 1 : 0;
    // The decoy: the model on its own sequence, the same rule.
    const DecoyAssay a = searchDecoy(lib, i, rules);
    const OpenMS::AASequence seq = OpenMS::AASequence::fromString(a.sequence);
    const auto spectra = model.predict({seq}, {static_cast<int>(sp.charge[decoy])});
    std::vector<PredictedFragment> ranked;
    std::size_t above = 0;
    PredictedAssays::rank(seq, sp.charge[decoy], spectra[0], lo, hi, ranked, above);
    bool own = ranked.size() >= n;
    for (std::uint32_t j = 0; own && j < n; ++j)
    {
      const std::uint32_t d = sp.transition_begin[decoy] + j;
      own = st.type[d] == ranked[j].type && st.ordinal[d] == ranked[j].ordinal && st.charge[d] == ranked[j].charge &&
            std::fabs(ODIA::fromFixed(st.product_mz[d]) - ranked[j].mz) < 1e-4;
    }
    decoys_own += own ? 1 : 0;
  }
  std::cout << "targets as in the library: " << targets_as_library << " of " << set.pairs() << " (" << fragments_differing
            << " fragments differ); decoys by their own prediction: " << decoys_own << "\n";
  // Batches of a different composition may move a prediction in its last
  // bits, which can swap two nearly tied fragments: allow a handful.
  CHECK(static_cast<double>(targets_as_library) >= 0.98 * static_cast<double>(set.pairs()));
  CHECK(static_cast<double>(decoys_own) >= 0.98 * static_cast<double>(set.pairs()));

  if (synth::failures) { std::cerr << synth::failures << " check(s) failed\n"; return 1; }
  std::cout << "identify_predicted_assays_model: PASS\n";
  return 0;
}
