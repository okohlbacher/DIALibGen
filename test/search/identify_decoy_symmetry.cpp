// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

// The search decoys differ from their targets in fragment m/z only, and in
// nothing a target could never show. Each check here is a way a decoy set can
// be weaker (or stronger) than null targets by construction, found on a real
// Astral run:
//   1. a target whose interior cannot be rearranged (EEEEEEK) must not get a
//      substituted residue instead -- it has no decoy and leaves the search;
//   2. every decoy fragment lies in the m/z range of the library's target
//      fragments, where the generator put every target fragment;
//   3. no decoy is a fragment-level copy of its target (I = L swaps);
//   4. only methods that keep both termini are offered: reverse moves the
//      C-terminal K/R and is refused like mutate;
//   5. each terminus keeps TWO residues, so the short ions a proteome shares
//      (y1, y2, b1, b2 and their complements) are shared by the pair -- with
//      one kept, null pairs went to whichever member drew the more common
//      short composition (SearchDecoys.h).

#include "synthetic_library.h"

#include <odia/search/CandidateSelector.h>
#include <odia/search/SearchDecoys.h>

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/Residue.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

using namespace ODIA::search;

namespace
{
  std::vector<double> fragments(const SearchSet& set, std::size_t i)
  {
    const auto& p = set.library.precursors();
    const auto& t = set.library.transitions();
    std::vector<double> mz;
    for (std::uint32_t k = 0; k < p.transition_count[i]; ++k)
    { mz.push_back(ODIA::fromFixed(t.product_mz[p.transition_begin[i] + k])); }
    return mz;
  }

  /// Every decoy fragment within @p ppm of some target fragment.
  bool isCopy(const std::vector<double>& target, const std::vector<double>& decoy, double ppm)
  {
    for (const double d : decoy)
    {
      bool near = false;
      for (const double t : target) { near = near || std::fabs(d - t) <= d * ppm * 1e-6; }
      if (!near) { return false; }
    }
    return true;
  }

  SearchParams params(const char* method)
  {
    SearchParams p;
    p.subset = 0;
    p.max_pairs = 0;
    p.decoys = parseSearchDecoyMethod(method);
    return p;
  }
}

int main()
{
  // A background of ordinary peptides, which sets a wide fragment range.
  auto precursors = synth::peptides(400, 11);
  auto add = [&precursors](const std::string& sequence, const std::string& protein) {
    synth::Precursor p;
    p.sequence = sequence;
    p.protein = protein;
    p.rt = 50.0f;
    for (int z : {2, 3}) { p.charge = z; precursors.push_back(p); }
  };
  // 1. Interiors of one residue: no rearrangement exists.
  const std::vector<std::string> unshufflable = {"EEEEEEK", "LLLLLLLLLLR", "GSSSSSR", "LGGGGGR"};
  for (const auto& s : unshufflable) { add(s, "RUNS"); }
  // 3. Interiors of I and L only: every rearrangement has the target's masses.
  const std::vector<std::string> isobaric = {"ALILILK", "GIILLIR", "SLLIIIK"};
  for (const auto& s : isobaric) { add(s, "ISOBARIC"); }

  for (const char* method : {"shuffle", "pseudo_reverse"})
  {
    const ODIA::Library lib = synth::library(precursors);
    const SearchSet set = CandidateSelector::select(lib, params(method));
    std::cout << method << ": " << set.pairs() << " pairs of " << set.stats.drawn << " drawn, " << set.stats.no_decoy
              << " without a decoy\n";
    CHECK(set.pairs() > 700);

    // The input library's target fragment range.
    double lo = 1e300, hi = 0;
    {
      const auto& p = lib.precursors();
      const auto& t = lib.transitions();
      for (std::size_t i = 0; i < lib.precursorCount(); ++i)
      {
        for (std::uint32_t k = 0; k < p.transition_count[i]; ++k)
        {
          const double mz = ODIA::fromFixed(t.product_mz[p.transition_begin[i] + k]);
          lo = std::min(lo, mz);
          hi = std::max(hi, mz);
        }
      }
    }

    std::size_t outside = 0, copies = 0;
    std::set<std::string> searched;
    for (std::size_t k = 0; k < set.pairs(); ++k)
    {
      const std::size_t d = set.partnerOf(k);
      searched.insert(std::string(set.modifiedSequence(k)));
      const auto target = fragments(set, k), decoy = fragments(set, d);
      CHECK(target.size() == decoy.size());
      for (const double mz : decoy) { outside += (mz < lo || mz > hi) ? 1 : 0; }
      copies += isCopy(target, decoy, 10.0) ? 1 : 0;
    }
    std::cout << "  decoy fragments outside the target range [" << lo << ", " << hi << "]: " << outside
              << "; fragment-level copies: " << copies << "\n";
    CHECK(outside == 0);   // 2.
    CHECK(copies == 0);    // 3.
    for (const auto& s : unshufflable) { CHECK(searched.count(s) == 0); }   // 1.
    for (const auto& s : isobaric) { CHECK(searched.count(s) == 0); }       // 3.
  }

  // 2, with pressure: a narrow library whose W-rich peptides can only be
  // rearranged into range by keeping the W residues away from the C-terminus.
  {
    std::vector<synth::Precursor> narrow;
    const std::vector<std::string> sequences = {
      "WWWWWWGGGGGGGGR", "WWWWWGGGGGGGGGK", "WWWWGGGGGGGGGGR", "WWWWWWGGGGGGGK", "WWWWWGGGGGGGGR",
      "WWWWWWWGGGGGGGR", "WWWWWWAGGGGGGGK", "WWWWWSGGGGGGGGR", "WWWWWWGSGGGGGGK", "WWWWWWGGAGGGGGR",
      "AGGGGSGGAGGSGAK", "SGGAGGSGGAGGSGR", "GGSGGAGGSGGAGGK"};
    for (std::size_t n = 0; n < sequences.size(); ++n)
    {
      synth::Precursor p;
      p.sequence = sequences[n];
      p.protein = "NARROW_" + std::to_string(n);
      p.rt = static_cast<float>(n);
      for (int z : {2, 3}) { p.charge = z; narrow.push_back(p); }
    }
    const ODIA::Library lib = synth::library(narrow);
    double lo = 1e300, hi = 0;
    for (std::size_t i = 0; i < lib.transitionCount(); ++i)
    {
      lo = std::min(lo, ODIA::fromFixed(lib.transitions().product_mz[i]));
      hi = std::max(hi, ODIA::fromFixed(lib.transitions().product_mz[i]));
    }
    for (const char* method : {"shuffle", "pseudo_reverse"})
    {
      const SearchSet set = CandidateSelector::select(lib, params(method));
      std::size_t outside = 0;
      for (std::size_t k = 0; k < set.pairs(); ++k)
      {
        for (const double mz : fragments(set, set.partnerOf(k))) { outside += (mz < lo || mz > hi) ? 1 : 0; }
      }
      std::cout << "narrow library, " << method << ": " << set.pairs() << " pairs, " << set.stats.no_decoy
                << " without a decoy, " << outside << " decoy fragments outside [" << lo << ", " << hi << "]\n";
      CHECK(outside == 0);
    }
  }

  // 5. Both termini: TWO residues each, so the short ions a whole proteome
  //    shares (y1, y2, b1, b2, and with them y(n-1), y(n-2), b(n-1), b(n-2))
  //    are the pair's, not one member's. With one residue kept, a decoy
  //    traded its target's y2 for another composition, and a null pair was
  //    decided by which of the two is the more common one (a 4 % excess of
  //    decoys at the prefilter depth on an entrapment library; targets whose
  //    second-to-last residue is P won 1.73 : 1).
  {
    const ODIA::Library lib = synth::library(precursors);
    const DecoyRules rules = CandidateSelector::decoyRules(lib, params("shuffle"));
    std::size_t made = 0, interior_differs = 0;
    for (std::size_t i = 0; i < lib.precursorCount(); ++i)
    {
      const DecoyAssay a = searchDecoy(lib, i, rules);
      if (a.outcome != DecoyOutcome::Made) { continue; }
      ++made;
      const auto target = OpenMS::AASequence::fromString(std::string(lib.strings().get(lib.precursors().modified_sequence[i])));
      const auto decoy = OpenMS::AASequence::fromString(a.sequence);
      CHECK(decoy.size() == target.size());
      CHECK(decoy.getPrefix(SEARCH_DECOY_KEEP_NTERM).toString() == target.getPrefix(SEARCH_DECOY_KEEP_NTERM).toString());
      CHECK(decoy.getSuffix(SEARCH_DECOY_KEEP_CTERM).toString() == target.getSuffix(SEARCH_DECOY_KEEP_CTERM).toString());
      // y1, y2, b1, b2 and their complements are then the same ions in both.
      for (const std::size_t o : {std::size_t{1}, std::size_t{2}})
      {
        for (int z : {1, 2})
        {
          CHECK(std::fabs(decoy.getSuffix(o).getMZ(z, OpenMS::Residue::YIon) - target.getSuffix(o).getMZ(z, OpenMS::Residue::YIon)) < 1e-9);
          CHECK(std::fabs(decoy.getPrefix(o).getMZ(z, OpenMS::Residue::BIon) - target.getPrefix(o).getMZ(z, OpenMS::Residue::BIon)) < 1e-9);
          const std::size_t big = target.size() - o;
          CHECK(std::fabs(decoy.getPrefix(big).getMZ(z, OpenMS::Residue::BIon) - target.getPrefix(big).getMZ(z, OpenMS::Residue::BIon)) < 1e-9);
          CHECK(std::fabs(decoy.getSuffix(big).getMZ(z, OpenMS::Residue::YIon) - target.getSuffix(big).getMZ(z, OpenMS::Residue::YIon)) < 1e-9);
        }
      }
      interior_differs += decoy.toString() != target.toString() ? 1 : 0;
    }
    std::cout << "terminal residues kept: " << made << " decoys, " << interior_differs << " differ from their target\n";
    CHECK(made > 700);
    CHECK(interior_differs == made);
  }

  // 4. Only methods that keep both termini.
  for (const char* refused : {"reverse", "mutate", "none", "anything"})
  {
    bool threw = false;
    try { (void)parseSearchDecoyMethod(refused); } catch (const std::invalid_argument&) { threw = true; }
    CHECK(threw);
  }
  {
    SearchParams p;
    p.decoys = ODIA::DecoyMethod::Reverse;
    bool threw = false;
    try { p.validate(); } catch (const std::invalid_argument&) { threw = true; }
    CHECK(threw);
  }

  if (synth::failures) { std::cerr << synth::failures << " check(s) failed\n"; return 1; }
  std::cout << "identify_decoy_symmetry: PASS\n";
  return 0;
}
