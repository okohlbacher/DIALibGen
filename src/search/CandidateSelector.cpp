// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <odia/search/CandidateSelector.h>

#include <odia/search/SearchDecoys.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>

namespace ODIA::search
{
  namespace
  {
    /// SplitMix64's finaliser: FNV-1a alone mixes its last bytes (the charge)
    /// poorly into the high bits, which are what a ranking compares.
    std::uint64_t finalise(std::uint64_t z)
    {
      z += 0x9E3779B97F4A7C15ull;
      z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
      z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
      return z ^ (z >> 31);
    }

    /// Identity of a precursor within one library: interned sequence handle and charge.
    std::uint64_t identity(const Library& library, std::size_t i)
    {
      return (static_cast<std::uint64_t>(library.precursors().modified_sequence[i]) << 8) |
             library.precursors().charge[i];
    }

    constexpr std::size_t none = std::numeric_limits<std::size_t>::max();
  }

  std::string_view SearchSet::modifiedSequence(std::size_t i) const
  { return library.strings().get(library.precursors().modified_sequence[i]); }

  std::string_view SearchSet::proteinGroup(std::size_t i) const
  { return library.strings().get(library.precursors().protein_group[i]); }

  int SearchSet::charge(std::size_t i) const { return library.precursors().charge[i]; }

  std::string SearchSet::precursorId(std::size_t i) const
  {
    std::string id(modifiedSequence(i));
    id += std::to_string(charge(i));
    if (isDecoy(i)) { id += "_decoy"; }
    return id;
  }

  std::uint64_t CandidateSelector::drawKey(std::string_view modified_sequence, int charge, std::uint64_t seed)
  {
    std::uint64_t h = 14695981039346656037ull;
    auto mix = [&h](unsigned char c) { h ^= c; h *= 1099511628211ull; };
    for (const char c : modified_sequence) { mix(static_cast<unsigned char>(c)); }
    mix('/');
    for (const char c : std::to_string(charge)) { mix(static_cast<unsigned char>(c)); }
    return finalise(h ^ finalise(seed));
  }

  RtScale CandidateSelector::rtScale(const Library& library)
  {
    const auto& p = library.precursors();
    RtScale s;
    bool any = false;
    for (std::size_t i = 0; i < library.precursorCount(); ++i)
    {
      if (p.decoy[i] || !std::isfinite(p.irt[i])) { continue; }
      const double v = p.irt[i];
      if (!any) { s.min = s.max = v; any = true; }
      else { s.min = std::min(s.min, v); s.max = std::max(s.max, v); }
    }
    if (!any) { throw std::invalid_argument("search: no target precursor of the library has a finite RT"); }
    if (!(s.max > s.min))
    { throw std::invalid_argument("search: every library RT is " + std::to_string(s.min) + "; an RT scale needs a range"); }
    return s;
  }

  bool CandidateSelector::inWindow(double mz, const std::vector<IsolationWindow>& windows)
  {
    for (const auto& w : windows)
    {
      if (w.lower < mz && mz < w.upper) { return true; }
    }
    return false;
  }

  SearchSet CandidateSelector::select(const Library& library, const SearchParams& params,
                                      const std::vector<IsolationWindow>& windows)
  {
    params.validate();
    SearchSet out;
    SelectionStats& st = out.stats;
    const auto& pre = library.precursors();
    const std::size_t library_size = library.precursorCount();
    st.library_precursors = library_size;
    out.rt_scale = rtScale(library);
    // Decoy fragments must stay where target fragments can be: the input
    // library's target fragment range, read from targets only.
    DecoyRules rules;
    rules.method = params.decoys;
    rules.min_fragments = SearchParams::min_assay_fragments;
    std::tie(rules.fragment_min, rules.fragment_max) = targetFragmentRange(library);
    st.fragment_mz_min = fromFixed(rules.fragment_min);
    st.fragment_mz_max = fromFixed(rules.fragment_max);

    // 1. Eligible targets. Every exclusion is a property of the target alone,
    //    which its decoy would share, so it removes whole pairs. The window
    //    check reads the precursor m/z, which a decoy inherits: a pair outside
    //    every isolation window could never be extracted, and would only use
    //    up search:max_pairs.
    st.windows = windows.size();
    std::vector<std::pair<std::uint64_t, std::size_t>> by_identity;
    for (std::size_t i = 0; i < library_size; ++i)
    {
      if (pre.decoy[i]) { ++st.library_decoys_ignored; continue; }
      ++st.targets;
      if (pre.mz[i] == MZ_INVALID) { ++st.ineligible_mz; continue; }
      if (pre.charge[i] == 0) { ++st.ineligible_charge; continue; }
      if (!std::isfinite(pre.irt[i])) { ++st.ineligible_rt; continue; }
      if (pre.transition_count[i] < SearchParams::min_assay_fragments) { ++st.ineligible_fragments; continue; }
      if (!windows.empty() && !inWindow(fromFixed(pre.mz[i]), windows)) { ++st.ineligible_window; continue; }
      by_identity.emplace_back(identity(library, i), i);
    }
    // A (sequence, charge) held by two targets cannot be paired structurally:
    // both decoys would carry the same key. Excluding ALL of them keeps the
    // rule independent of row order.
    std::sort(by_identity.begin(), by_identity.end());
    std::vector<std::pair<std::uint64_t, std::size_t>> draws;   // (draw key, library index)
    draws.reserve(by_identity.size());
    for (std::size_t a = 0; a < by_identity.size();)
    {
      std::size_t b = a + 1;
      while (b < by_identity.size() && by_identity[b].first == by_identity[a].first) { ++b; }
      if (b - a > 1) { st.duplicate_key += b - a; }
      else
      {
        const std::size_t i = by_identity[a].second;
        draws.emplace_back(drawKey(library.strings().get(pre.modified_sequence[i]), pre.charge[i], params.seed), i);
      }
      a = b;
    }
    std::vector<std::pair<std::uint64_t, std::size_t>>().swap(by_identity);
    st.eligible = draws.size();

    // 2. The draw: the `subset` lowest keys. Ties (never observed; 64-bit keys)
    //    fall back to the key itself, which is unique after step 1.
    auto less = [&](const std::pair<std::uint64_t, std::size_t>& a, const std::pair<std::uint64_t, std::size_t>& b) {
      if (a.first != b.first) { return a.first < b.first; }
      const auto sa = library.strings().get(pre.modified_sequence[a.second]);
      const auto sb = library.strings().get(pre.modified_sequence[b.second]);
      if (sa != sb) { return sa < sb; }
      return pre.charge[a.second] < pre.charge[b.second];
    };
    if (params.subset != 0 && params.subset < draws.size())
    {
      std::nth_element(draws.begin(), draws.begin() + static_cast<std::ptrdiff_t>(params.subset), draws.end(), less);
      draws.resize(params.subset);
    }
    std::sort(draws.begin(), draws.end(), less);
    st.drawn = draws.size();

    // 3. Decoys, in draw order, for as many drawn targets as search:max_pairs
    //    needs -- the result is the same as building every decoy and keeping
    //    the lowest max_pairs pairs, without building the ones the cap drops.
    //    Built on a library holding exactly those targets. Its strings are
    //    re-interned, so equal sequences now share one handle even if the
    //    input's did not; a duplicate the input's handles hid is caught here
    //    and excluded like the others.
    const std::size_t cap = params.max_pairs != 0 ? params.max_pairs : std::numeric_limits<std::size_t>::max();
    std::size_t n = std::min(draws.size(), cap == std::numeric_limits<std::size_t>::max() ? draws.size() : cap + cap / 50 + 16);
    Library sub;
    std::vector<std::size_t> decoy_of;
    std::size_t made = 0;
    DecoyBuild built;
    for (;;)
    {
      std::vector<std::size_t> drawn;
      drawn.reserve(n);
      for (std::size_t t = 0; t < n; ++t) { drawn.push_back(draws[t].second); }
      sub = library.subsetByIndex(drawn);
      {
        std::unordered_map<std::uint64_t, std::size_t> seen;
        std::vector<char> duplicate(n, 0);
        for (std::size_t t = 0; t < n; ++t)
        {
          const auto ins = seen.emplace(identity(sub, t), t);
          if (!ins.second) { duplicate[t] = 1; duplicate[ins.first->second] = 1; }
        }
        if (std::find(duplicate.begin(), duplicate.end(), 1) != duplicate.end())
        {
          std::vector<std::pair<std::uint64_t, std::size_t>> kept;
          for (std::size_t t = 0; t < draws.size(); ++t)
          {
            if (t < n && duplicate[t]) { ++st.duplicate_key; --st.drawn; continue; }
            kept.push_back(draws[t]);
          }
          draws.swap(kept);
          n = std::min(n, draws.size());
          continue;
        }
      }
      built = appendSearchDecoys(sub, rules);
      made = built.made;
      const auto skipped = static_cast<std::size_t>(
        std::count_if(built.outcome.begin(), built.outcome.end(), [](DecoyOutcome o) { return o != DecoyOutcome::Made; }));
      if (made + skipped != n || built.outcome.size() != n || sub.precursorCount() != n + made)
      {
        throw std::logic_error("search: decoy construction accounted for " + std::to_string(made) + " decoys and " +
                               std::to_string(skipped) + " failures among " + std::to_string(n) + " targets");
      }
      if (made >= cap || n == draws.size()) { break; }
      n = std::min(draws.size(), n + 2 * (cap - made) + 16);   // too many failures: extend and rebuild
    }

    // Structural pairing: a decoy stores its target's sequence handle and charge.
    std::unordered_map<std::uint64_t, std::size_t> target_of;
    target_of.reserve(n);
    for (std::size_t t = 0; t < n; ++t) { target_of.emplace(identity(sub, t), t); }
    decoy_of.assign(n, none);
    for (std::size_t d = n; d < sub.precursorCount(); ++d)
    {
      const auto it = target_of.find(identity(sub, d));
      if (it == target_of.end() || decoy_of[it->second] != none)
      { throw std::logic_error("search: decoy " + std::to_string(d) + " has no unique target"); }
      decoy_of[it->second] = d;
    }

    // 4. Whole pairs only, in draw order, up to the cap. Drawn targets past the
    //    last kept pair are `capped`, so drawn = no_decoy + capped + pairs
    //    whatever n the loop above settled on.
    std::vector<std::size_t> pair_target;
    pair_target.reserve(std::min(made, cap));
    std::size_t walked = 0;
    for (; walked < n && pair_target.size() < cap; ++walked)
    {
      switch (built.outcome[walked])
      {
        case DecoyOutcome::Made: break;
        case DecoyOutcome::Unparsable: ++st.decoy_unparsable; break;
        case DecoyOutcome::Unshufflable: ++st.decoy_unshufflable; break;
        case DecoyOutcome::OutOfRange: ++st.decoy_out_of_range; break;
        case DecoyOutcome::Copy: ++st.decoy_copy; break;
        case DecoyOutcome::TooFewFragments: ++st.decoy_too_few_fragments; break;
      }
      if (decoy_of[walked] == none) { ++st.no_decoy; }
      else
      {
        pair_target.push_back(walked);
        if (built.redrawn[walked] > 0) { ++st.decoy_redrawn; }
        st.fragment_slots_dropped += pre.transition_count[draws[walked].second] - sub.precursors().transition_count[walked];
      }
    }
    if (st.no_decoy != st.decoy_unparsable + st.decoy_unshufflable + st.decoy_out_of_range + st.decoy_copy +
                       st.decoy_too_few_fragments)
    { throw std::logic_error("search: decoy failures do not add up to the targets without a decoy"); }
    st.capped = draws.size() - walked;
    st.pairs = pair_target.size();

    // 5. Targets by precursor m/z (ties keep draw order), decoys in the same order.
    std::stable_sort(pair_target.begin(), pair_target.end(), [&](std::size_t a, std::size_t b) {
      return sub.precursors().mz[a] < sub.precursors().mz[b];
    });
    std::vector<std::size_t> order(pair_target);
    for (const std::size_t t : pair_target) { order.push_back(decoy_of[t]); }
    out.library = sub.subsetByIndex(order);
    out.source.resize(order.size());
    out.draw.resize(pair_target.size());
    for (std::size_t k = 0; k < pair_target.size(); ++k)
    {
      out.source[k] = out.source[k + pair_target.size()] = draws[pair_target[k]].second;
      out.draw[k] = draws[pair_target[k]].first;
    }
    return out;
  }
}
