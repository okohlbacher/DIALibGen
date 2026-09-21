// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

// The paired random subset, the in-memory decoys and the sequence-blind assays:
// every selected target has exactly one decoy, the draw is deterministic per
// seed and blind to labels and row order, and assays carry no sequence.

#include "synthetic_library.h"

#include <odia/search/AssayBuilder.h>
#include <odia/search/CandidateSelector.h>

#include <algorithm>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

using namespace ODIA::search;

namespace
{
  std::vector<std::string> ids(const SearchSet& s)
  {
    std::vector<std::string> v;
    for (std::size_t i = 0; i < s.size(); ++i) { v.push_back(s.precursorId(i)); }
    return v;
  }

  SearchParams params(std::size_t subset, std::size_t max_pairs, std::uint64_t seed = 42)
  {
    SearchParams p;
    p.subset = subset;
    p.max_pairs = max_pairs;
    p.seed = seed;
    return p;
  }

  bool sameTransitions(const SearchSet& a, const SearchSet& b)
  {
    const auto& x = a.library.transitions();
    const auto& y = b.library.transitions();
    return x.product_mz == y.product_mz && x.library_intensity == y.library_intensity && x.ordinal == y.ordinal;
  }
}

int main()
{
  const auto precursors = synth::peptides(3000, 7);   // 6000 targets
  const ODIA::Library lib = synth::library(precursors);

  // ---- 1. pairing ---------------------------------------------------------
  const SearchSet set = CandidateSelector::select(lib, params(1000, 0));
  const auto& st = set.stats;
  std::cout << "selected " << st.pairs << " pairs of " << st.drawn << " drawn (" << st.no_decoy << " without decoy) from "
            << st.eligible << " eligible\n";
  CHECK(st.targets == precursors.size());
  CHECK(st.eligible == precursors.size());
  CHECK(st.drawn == 1000);
  CHECK(st.drawn == st.no_decoy + st.capped + st.pairs);
  CHECK(st.pairs >= 990);
  CHECK(set.size() == 2 * set.pairs());
  CHECK(set.library.precursorCount() == set.size());
  {
    const auto& p = set.library.precursors();
    std::map<std::string, int> target_keys, decoy_keys;
    for (std::size_t k = 0; k < set.pairs(); ++k)
    {
      const std::size_t t = k, d = set.pairs() + k;
      CHECK(p.decoy[t] == 0 && p.decoy[d] == 1);
      CHECK(!set.isDecoy(t) && set.isDecoy(d));
      CHECK(set.pairOf(t) == static_cast<std::int64_t>(k) && set.pairOf(d) == set.pairOf(t));
      CHECK(set.partnerOf(t) == d && set.partnerOf(d) == t);
      CHECK(set.modifiedSequence(t) == set.modifiedSequence(d));
      CHECK(set.charge(t) == set.charge(d));
      CHECK(set.proteinGroup(t) == set.proteinGroup(d));
      CHECK(p.mz[t] == p.mz[d]);
      CHECK(p.irt[t] == p.irt[d]);
      CHECK(set.source[t] == set.source[d]);
      CHECK(set.precursorId(d) == set.precursorId(t) + "_decoy");
      // The draw key is a function of what target and decoy share.
      CHECK(set.draw[k] == CandidateSelector::drawKey(set.modifiedSequence(d), set.charge(d), 42));
      CHECK(set.draw[k] == CandidateSelector::drawKey(set.modifiedSequence(t), set.charge(t), 42));
      // The decoy is a different assay: its fragment m/z moved.
      bool moved = false;
      for (std::uint32_t j = 0; j < std::min(p.transition_count[t], p.transition_count[d]); ++j)
      {
        moved = moved || set.library.transitions().product_mz[p.transition_begin[t] + j] !=
                         set.library.transitions().product_mz[p.transition_begin[d] + j];
      }
      CHECK(moved);
      CHECK(p.transition_count[d] >= SearchParams::min_assay_fragments);
      ++target_keys[set.precursorId(t)];
      ++decoy_keys[set.precursorId(t)];
      // Source points at the matching target of the input library.
      const auto& src = lib.precursors();
      CHECK(lib.strings().get(src.modified_sequence[set.source[t]]) == set.modifiedSequence(t));
      CHECK(src.charge[set.source[t]] == p.charge[t] && src.decoy[set.source[t]] == 0);
    }
    // Exactly one decoy per target: every key once among targets, once among decoys.
    CHECK(target_keys.size() == set.pairs() && decoy_keys.size() == set.pairs());
    for (const auto& kv : target_keys) { CHECK(kv.second == 1); }
    // Targets are ordered by precursor m/z.
    for (std::size_t k = 1; k < set.pairs(); ++k) { CHECK(p.mz[k - 1] <= p.mz[k]); }
  }
  // The RT scale is the INPUT library's range, not the subset's.
  {
    float lo = 1e9f, hi = -1e9f;
    for (const auto& pr : precursors) { lo = std::min(lo, pr.rt); hi = std::max(hi, pr.rt); }
    CHECK(set.rt_scale.min == lo && set.rt_scale.max == hi);
  }

  // ---- 2. deterministic per seed, different across seeds ------------------
  const SearchSet again = CandidateSelector::select(lib, params(1000, 0));
  CHECK(ids(again) == ids(set));
  CHECK(again.source == set.source && again.draw == set.draw);
  CHECK(sameTransitions(again, set));
  const SearchSet other = CandidateSelector::select(lib, params(1000, 0, 43));
  {
    const auto a = ids(set), b = ids(other);
    const std::set<std::string> sa(a.begin(), a.end());
    std::size_t shared = 0;
    for (const auto& id : b) { shared += sa.count(id); }
    std::cout << "seed 42 vs 43: " << shared / 2 << " of " << set.pairs() << " pairs shared\n";
    CHECK(other.pairs() >= 990);
    CHECK(shared < a.size() / 2);   // a 1/6 subset: expect ~1/6 overlap
  }

  // ---- 3. blind to row order and to decoys already in the file -------------
  {
    auto shuffled = precursors;
    std::mt19937_64 rng(99);
    for (std::size_t k = shuffled.size(); k > 1; --k) { std::swap(shuffled[k - 1], shuffled[rng() % k]); }
    auto with_decoys = shuffled;
    for (const auto& pr : precursors)
    {
      auto d = pr;
      d.decoy = true;
      with_decoys.push_back(d);
    }
    const ODIA::Library lib2 = synth::library(with_decoys);
    const SearchSet s2 = CandidateSelector::select(lib2, params(1000, 0));
    CHECK(ids(s2) == ids(set));
    CHECK(s2.draw == set.draw);
    CHECK(sameTransitions(s2, set));
    CHECK(s2.stats.library_decoys_ignored == precursors.size());
    CHECK(s2.stats.targets == precursors.size());
  }

  // ---- 4. the cap keeps the lowest draw keys --------------------------------
  {
    const SearchSet capped = CandidateSelector::select(lib, params(1000, 300));
    CHECK(capped.pairs() == 300 && capped.stats.pairs == 300);
    CHECK(capped.stats.drawn == 1000);
    CHECK(capped.stats.drawn == capped.stats.no_decoy + capped.stats.capped + capped.stats.pairs);
    CHECK(capped.stats.capped >= 1000 - 300 - set.stats.no_decoy);
    std::vector<std::uint64_t> keys(set.draw);
    std::sort(keys.begin(), keys.end());
    std::vector<std::uint64_t> kept(capped.draw);
    std::sort(kept.begin(), kept.end());
    CHECK(std::equal(kept.begin(), kept.end(), keys.begin()));
    // subset 0 = every eligible target
    const SearchSet all = CandidateSelector::select(lib, params(0, 0));
    CHECK(all.stats.drawn == precursors.size());
  }
  // ... also when many drawn targets get no decoy: palindromes have no reversed
  // decoy, so the cap-aware construction must extend its draw and still keep
  // exactly the lowest-key pairs.
  {
    auto mixed = synth::peptides(200, 5);
    std::mt19937_64 rng(8);
    static const std::string aa = "ADEFGHILMNPQSTVWY";
    for (std::size_t k = 0; k < 200; ++k)
    {
      std::string half;
      while (half.size() < 4 + k % 4) { half.push_back(aa[rng() % aa.size()]); }
      synth::Precursor p;
      p.sequence = "K" + half + std::string(half.rbegin(), half.rend()) + "K";
      p.protein = "PALINDROME_" + std::to_string(k);
      p.rt = static_cast<float>(k) / 2.0f;
      for (int z : {2, 3}) { p.charge = z; mixed.push_back(p); }
    }
    const ODIA::Library pal = synth::library(mixed);
    SearchParams p = params(0, 300);
    p.decoys = parseSearchDecoyMethod("reverse");
    const SearchSet capped = CandidateSelector::select(pal, p);
    p.max_pairs = 0;
    const SearchSet full = CandidateSelector::select(pal, p);
    std::cout << "reverse decoys: " << full.stats.no_decoy << " of " << full.stats.drawn << " drawn targets without a decoy\n";
    CHECK(full.stats.no_decoy >= 350);
    CHECK(capped.pairs() == 300);
    CHECK(capped.stats.no_decoy > 0);
    CHECK(capped.stats.drawn == capped.stats.no_decoy + capped.stats.capped + capped.stats.pairs);
    std::vector<std::uint64_t> keys(full.draw), kept(capped.draw);
    std::sort(keys.begin(), keys.end());
    std::sort(kept.begin(), kept.end());
    CHECK(std::equal(kept.begin(), kept.end(), keys.begin()));
    for (std::size_t k = 0; k < capped.pairs(); ++k)
    {
      const std::string s(capped.modifiedSequence(k));
      CHECK(s != std::string(s.rbegin(), s.rend()));
    }
  }

  // ---- 5. exclusions ----------------------------------------------------------
  {
    auto odd = synth::peptides(50, 3);
    auto few = odd[0]; few.sequence = "PEPTIDEAFGHK"; few.fragments = 2;          // too few transitions
    auto nan_rt = odd[1]; nan_rt.sequence = "PEPTIDEAFGMK"; nan_rt.rt = std::nanf("");
    auto dup1 = odd[2]; dup1.sequence = "DAPTIDEAFGMR";
    auto dup2 = dup1; dup2.protein = "OTHER";
    odd.push_back(few); odd.push_back(nan_rt); odd.push_back(dup1); odd.push_back(dup2);
    const SearchSet s = CandidateSelector::select(synth::library(odd), params(0, 0));
    CHECK(s.stats.ineligible_fragments == 1);
    CHECK(s.stats.ineligible_rt == 1);
    CHECK(s.stats.duplicate_key == 2);
    CHECK(s.stats.eligible == 100);
  }

  // ---- 6. decoy methods ---------------------------------------------------------
  {
    bool refused = false;
    try { (void)parseSearchDecoyMethod("mutate"); } catch (const std::invalid_argument&) { refused = true; }
    CHECK(refused);
    refused = false;
    SearchParams p = params(100, 0);
    p.decoys = ODIA::DecoyMethod::Mutate;
    try { (void)CandidateSelector::select(lib, p); } catch (const std::invalid_argument&) { refused = true; }
    CHECK(refused);
    for (const char* m : {"pseudo_reverse", "reverse"})
    {
      p.decoys = parseSearchDecoyMethod(m);
      const SearchSet s = CandidateSelector::select(lib, p);
      CHECK(s.pairs() >= 90);
      CHECK(ids(s).size() == 2 * s.pairs());
    }
  }

  // ---- 7. assays: sequence-blind, verbatim protein groups, [0, 100] RT --------
  {
    std::vector<std::size_t> all(set.size());
    for (std::size_t i = 0; i < all.size(); ++i) { all[i] = i; }
    const AssayChunk chunk = AssayBuilder::build(set, all);
    const auto& exp = chunk.experiment;
    CHECK(exp.compounds.size() == set.size());
    CHECK(exp.transitions.size() == set.library.transitionCount());
    std::set<std::string> compound_ids;
    std::size_t k = 0;
    for (std::size_t c = 0; c < exp.compounds.size(); ++c)
    {
      const auto& comp = exp.compounds[c];
      const std::size_t i = chunk.precursor[c];
      CHECK(comp.sequence.empty());
      CHECK(comp.id == set.precursorId(i));
      compound_ids.insert(comp.id);
      CHECK(comp.protein_refs.size() == 1 && comp.protein_refs[0] == std::string(set.proteinGroup(i)));
      CHECK(comp.charge == set.charge(i));
      CHECK(comp.drift_time == -1.0);
      CHECK(comp.rt >= 0.0 && comp.rt <= 100.0);
      CHECK(std::fabs(comp.rt - set.rt_scale.toAssay(set.library.precursors().irt[i])) < 1e-9);
      CHECK(comp.peptide_group_label == set.precursorId(set.isDecoy(i) ? set.partnerOf(i) : i));
      for (std::uint32_t j = 0; j < set.library.precursors().transition_count[i]; ++j, ++k)
      {
        const auto& t = exp.transitions[k];
        CHECK(t.peptide_ref == comp.id);
        CHECK(t.decoy == set.isDecoy(i));
        CHECK(t.detecting_transition && t.quantifying_transition && !t.identifying_transition);
        CHECK(t.precursor_im == -1.0);
      }
    }
    CHECK(compound_ids.size() == set.size());
    // Ion mobility on: a library without 1/K0 or CCS is counted, not guessed.
    AssayOptions im;
    im.ion_mobility = true;
    CHECK(AssayBuilder::build(set, {0, 1, 2}, im).missing_im == 3);

    const auto parts = AssayBuilder::chunks(set, 101);
    std::size_t covered = 0;
    for (const auto& part : parts)
    {
      CHECK(part.size() <= 101 && part.size() % 2 == 0);
      const std::set<std::size_t> in(part.begin(), part.end());
      for (const std::size_t i : part) { CHECK(in.count(set.partnerOf(i)) == 1); }
      covered += part.size();
    }
    CHECK(covered == set.size());
  }

  if (synth::failures) { std::cerr << synth::failures << " check(s) failed\n"; return 1; }
  std::cout << "identify_subset_pairing: PASS\n";
  return 0;
}
