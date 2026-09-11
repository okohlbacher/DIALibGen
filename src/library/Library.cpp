// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <unordered_map>
#include <odia/Library.h>

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <cstring>
#include <numeric>

namespace ODIA
{

  FragmentType parseFragmentType(std::string_view s)
  {
    if (s.empty()) { return FragmentType::Unknown; }
    switch (s.front())
    {
      case 'a': case 'A': return FragmentType::A;
      case 'b': case 'B': return FragmentType::B;
      case 'c': case 'C': return FragmentType::C;
      case 'x': case 'X': return FragmentType::X;
      case 'y': case 'Y': return FragmentType::Y;
      case 'z': case 'Z': return FragmentType::Z;
      case 'p': case 'P': return FragmentType::Precursor;
      default: return FragmentType::Unknown;
    }
  }

  LossType parseLossType(std::string_view s)
  {
    if (s.empty() || s == "noloss" || s == "none" || s == "-") { return LossType::None; }
    if (s == "H2O" || s == "h2o") { return LossType::Water; }
    if (s == "NH3" || s == "nh3") { return LossType::Ammonia; }
    // HPO3 and H3PO4 differ by exactly one water (79.96633 vs 97.97690 Da).
    // Aliasing them was harmless while the label was only round-tripped; it
    // stopped being harmless once the label started feeding a mass.
    if (s == "H3PO4") { return LossType::Phospho; }
    if (s == "HPO3") { return LossType::Metaphosphate; }
    if (s == "CO" || s == "co") { return LossType::CO; }
    return LossType::Other;
  }

  std::string_view toString(FragmentType t)
  {
    switch (t)
    {
      case FragmentType::A: return "a";
      case FragmentType::B: return "b";
      case FragmentType::C: return "c";
      case FragmentType::X: return "x";
      case FragmentType::Y: return "y";
      case FragmentType::Z: return "z";
      case FragmentType::Precursor: return "p";
      case FragmentType::Unknown: break;
    }
    return "?";
  }

  std::string_view toString(LossType l)
  {
    switch (l)
    {
      case LossType::None: return "noloss";
      case LossType::Water: return "H2O";
      case LossType::Ammonia: return "NH3";
      case LossType::Phospho: return "H3PO4";
      case LossType::Metaphosphate: return "HPO3";
      case LossType::CO: return "CO";
      case LossType::Other: break;
    }
    return "other";
  }

  // -------------------------------------------------------------- StringArena

  void StringArena::reserve(std::size_t entries, std::size_t chars)
  {
    entries_.reserve(entries);
    lookup_.reserve(entries);
    if (chars > BLOCK && blocks_.empty())
    {
      blocks_.emplace_back(chars);
      block_used_ = 0;
    }
  }

  char* StringArena::allocate(std::size_t n)
  {
    if (blocks_.empty() || block_used_ + n > blocks_.back().size())
    {
      blocks_.emplace_back(std::max(n, BLOCK));
      block_used_ = 0;
    }
    char* p = blocks_.back().data() + block_used_;
    block_used_ += n;
    return p;
  }

  std::uint32_t StringArena::intern(std::string_view s)
  {
    if (auto it = lookup_.find(s); it != lookup_.end()) { return it->second; }

    char* p = allocate(s.size());
    // memcpy with a null source is UB even for a zero count, and an absent
    // column yields an empty string_view whose data() is null.
    if (!s.empty()) { std::memcpy(p, s.data(), s.size()); }
    bytes_ += s.size();

    entries_.push_back({p, static_cast<std::uint32_t>(s.size())});
    const auto handle = static_cast<std::uint32_t>(entries_.size() - 1);
    lookup_.emplace(std::string_view(p, s.size()), handle);
    return handle;
  }

  std::string_view StringArena::get(std::uint32_t handle) const
  {
    if (handle >= entries_.size()) { return {}; }
    return std::string_view(entries_[handle].data, entries_[handle].length);
  }

  void StringArena::releaseLookup()
  {
    // swap-with-empty, because clear() on an unordered_map keeps the bucket
    // array -- which is the larger half at this size.
    std::unordered_map<std::string_view, std::uint32_t>().swap(lookup_);
  }

  std::size_t StringArena::footprintBytes() const
  {
    std::size_t blocks = 0;
    for (const auto& b : blocks_) { blocks += b.capacity(); }
    // One hash node per distinct string, plus the bucket array.
    //
    // The model was key view + value + a next pointer, with a comment claiming
    // it was "exact enough to stop the figure being wrong by a factor". It was
    // wrong by a factor: releasing this map dropped peak RSS 1.28 -> 1.13 GiB
    // (154 MB) on a 4.26 M-precursor library while the accounting predicted
    // only 47 MiB. libstdc++ caches the hash in the node and pads to alignment,
    // and the allocator rounds on top of that. 48 B is measured against that
    // delta rather than derived from sizeof.
    const std::size_t node = 48;
    const std::size_t map = lookup_.size() * node + lookup_.bucket_count() * sizeof(void*);
    return blocks + entries_.capacity() * sizeof(Entry) + map;
  }

  // ------------------------------------------------------------------ Library

  std::size_t Library::decoyCount() const
  {
    return static_cast<std::size_t>(
      std::count(precursors_.decoy.begin(), precursors_.decoy.end(), std::uint8_t{1}));
  }

  std::size_t Library::invalidMzCount() const
  {
    return static_cast<std::size_t>(
      std::count(precursors_.mz.begin(), precursors_.mz.end(), MZ_INVALID));
  }

  std::size_t Library::invalidMzTransitionCount() const
  {
    return static_cast<std::size_t>(
      std::count(transitions_.product_mz.begin(), transitions_.product_mz.end(), MZ_INVALID));
  }

  void Library::reserve(std::size_t precursors, std::size_t transitions)
  {
    auto& p = precursors_;
    p.mz.reserve(precursors);
    p.irt.reserve(precursors);
    p.im.reserve(precursors);
    p.ccs.reserve(precursors);
    p.charge.reserve(precursors);
    p.decoy.reserve(precursors);
    p.modified_sequence.reserve(precursors);
    p.protein_group.reserve(precursors);
    p.transition_begin.reserve(precursors);
    p.transition_count.reserve(precursors);

    auto& t = transitions_;
    t.product_mz.reserve(transitions);
    t.library_intensity.reserve(transitions);
    t.type.reserve(transitions);
    t.ordinal.reserve(transitions);
    t.charge.reserve(transitions);
    t.loss.reserve(transitions);
  }

  void Library::shrinkToFit()
  {
    auto& p = precursors_;
    p.mz.shrink_to_fit(); p.irt.shrink_to_fit(); p.im.shrink_to_fit();
    p.ccs.shrink_to_fit();
    p.charge.shrink_to_fit(); p.decoy.shrink_to_fit();
    p.modified_sequence.shrink_to_fit(); p.protein_group.shrink_to_fit();
    p.transition_begin.shrink_to_fit(); p.transition_count.shrink_to_fit();

    auto& t = transitions_;
    t.product_mz.shrink_to_fit(); t.library_intensity.shrink_to_fit();
    t.type.shrink_to_fit(); t.ordinal.shrink_to_fit();
    t.charge.shrink_to_fit(); t.loss.shrink_to_fit();
  }

  std::size_t Library::dropDecoys()
  {
    const std::size_t n = precursorCount();
    std::size_t removed = 0;
    for (std::size_t i = 0; i < n; ++i) { if (precursors_.decoy[i]) { ++removed; } }
    if (removed == 0) { return 0; }

    // Rebuild rather than erase in place: transitions are a CSR layout, so
    // removing a precursor means recomputing every later transition_begin.
    // This mirrors sortByPrecursorMz's rebuild, including copying every
    // per-transition array -- omitting one there segfaulted the sort, and the
    // same omission here would do the same.
    PrecursorArrays np;
    TransitionArrays nt;
    const std::size_t keep = n - removed;
    np.mz.reserve(keep); np.irt.reserve(keep); np.im.reserve(keep);
    np.ccs.reserve(keep); np.charge.reserve(keep); np.decoy.reserve(keep);
    np.modified_sequence.reserve(keep); np.protein_group.reserve(keep);
    np.transition_begin.reserve(keep); np.transition_count.reserve(keep);

    for (std::size_t i = 0; i < n; ++i)
    {
      if (precursors_.decoy[i]) { continue; }
      np.mz.push_back(precursors_.mz[i]);
      np.irt.push_back(precursors_.irt[i]);
      np.im.push_back(precursors_.im[i]);
      np.ccs.push_back(precursors_.ccs.empty()
                         ? std::numeric_limits<float>::quiet_NaN()
                         : precursors_.ccs[i]);
      np.charge.push_back(precursors_.charge[i]);
      np.decoy.push_back(0);
      np.modified_sequence.push_back(precursors_.modified_sequence[i]);
      np.protein_group.push_back(precursors_.protein_group[i]);

      const std::uint32_t begin = precursors_.transition_begin[i];
      const std::uint32_t count = precursors_.transition_count[i];
      np.transition_begin.push_back(static_cast<std::uint32_t>(nt.product_mz.size()));
      np.transition_count.push_back(count);
      for (std::uint32_t k = 0; k < count; ++k)
      {
        const std::uint32_t s = begin + k;
        nt.product_mz.push_back(transitions_.product_mz[s]);
        nt.library_intensity.push_back(transitions_.library_intensity[s]);
        nt.type.push_back(transitions_.type[s]);
        nt.ordinal.push_back(transitions_.ordinal[s]);
        nt.charge.push_back(transitions_.charge[s]);
        nt.loss.push_back(transitions_.loss[s]);
      }
    }
    precursors_ = std::move(np);
    transitions_ = std::move(nt);
    return removed;
  }

  Library Library::subsetByIndex(const std::vector<std::size_t>& keep) const
  {
    Library out;
    // Handles are RE-INTERNED, not copied. `StringArena::Entry` holds a raw
    // pointer into the arena's own blocks, so copying the arena produces
    // entries that resolve into THIS library's storage -- fine while the parent
    // outlives the subset, a segfault the moment it does not. Re-interning also
    // means the subset carries only the strings it actually uses, which for the
    // -min_library_fragments filter is the point.
    const std::size_t n = precursorCount();
    std::unordered_map<std::uint32_t, std::uint32_t> remap;
    const auto rehome = [&](std::uint32_t h) {
      if (h == StringArena::npos) { return h; }
      const auto it = remap.find(h);
      if (it != remap.end()) { return it->second; }
      const std::uint32_t fresh = out.strings_.intern(strings_.get(h));
      remap.emplace(h, fresh);
      return fresh;
    };
    out.precursors_.mz.reserve(keep.size());
    out.precursors_.irt.reserve(keep.size());
    out.precursors_.im.reserve(keep.size());
    out.precursors_.ccs.reserve(keep.size());
    out.precursors_.charge.reserve(keep.size());
    out.precursors_.decoy.reserve(keep.size());
    out.precursors_.modified_sequence.reserve(keep.size());
    out.precursors_.protein_group.reserve(keep.size());
    out.precursors_.transition_begin.reserve(keep.size());
    out.precursors_.transition_count.reserve(keep.size());

    for (const std::size_t i : keep)
    {
      // An out-of-range index would read past the arrays and produce a library
      // that looks valid. Refuse instead -- the caller built this list from a
      // set of its own and a mismatch is a bug there, not bad input.
      if (i >= n)
      {
        throw std::out_of_range("Library::subsetByIndex: precursor " +
                                std::to_string(i) + " of " + std::to_string(n));
      }
      out.precursors_.mz.push_back(precursors_.mz[i]);
      out.precursors_.irt.push_back(precursors_.irt[i]);
      out.precursors_.im.push_back(precursors_.im[i]);
      out.precursors_.ccs.push_back(precursors_.ccs.empty()
                                      ? std::numeric_limits<float>::quiet_NaN()
                                      : precursors_.ccs[i]);
      out.precursors_.charge.push_back(precursors_.charge[i]);
      out.precursors_.decoy.push_back(precursors_.decoy[i]);
      out.precursors_.modified_sequence.push_back(rehome(precursors_.modified_sequence[i]));
      out.precursors_.protein_group.push_back(rehome(precursors_.protein_group[i]));

      const std::uint32_t begin = precursors_.transition_begin[i];
      const std::uint32_t count = precursors_.transition_count[i];
      out.precursors_.transition_begin.push_back(
        static_cast<std::uint32_t>(out.transitions_.product_mz.size()));
      out.precursors_.transition_count.push_back(count);
      for (std::uint32_t k = 0; k < count; ++k)
      {
        const std::uint32_t s = begin + k;
        out.transitions_.product_mz.push_back(transitions_.product_mz[s]);
        out.transitions_.library_intensity.push_back(transitions_.library_intensity[s]);
        out.transitions_.type.push_back(transitions_.type[s]);
        out.transitions_.ordinal.push_back(transitions_.ordinal[s]);
        out.transitions_.charge.push_back(transitions_.charge[s]);
        out.transitions_.loss.push_back(transitions_.loss[s]);
      }
    }
    // Gathering destroys m/z order even when `keep` was ascending, because the
    // source order is by m/z only if the source was sorted. Say so rather than
    // inherit a stale flag; every window slice depends on it.
    out.sorted_by_mz_ = false;
    return out;
  }

  void Library::sortByPrecursorMz()
  {
    const std::size_t n = precursorCount();
    std::vector<std::uint32_t> order(n);
    std::iota(order.begin(), order.end(), 0u);
    std::sort(order.begin(), order.end(), [this](std::uint32_t a, std::uint32_t b) {
      if (precursors_.mz[a] != precursors_.mz[b]) { return precursors_.mz[a] < precursors_.mz[b]; }
      // Ties broken deterministically: run-to-run reproducibility is a
      // precondition for every later measurement, not a nicety.
      if (precursors_.charge[a] != precursors_.charge[b]) { return precursors_.charge[a] < precursors_.charge[b]; }
      if (precursors_.decoy[a] != precursors_.decoy[b]) { return precursors_.decoy[a] < precursors_.decoy[b]; }
      // On the sequence TEXT, not the arena handle. The handle is
      // first-appearance order, so tie-breaking on it made the sort a function
      // of how the input happened to be ordered: re-sorting our own output
      // reordered 4,532 rows of the fixture.
      return strings_.get(precursors_.modified_sequence[a])
           < strings_.get(precursors_.modified_sequence[b]);
    });

    PrecursorArrays np;
    TransitionArrays nt;
    np.mz.resize(n); np.irt.resize(n); np.im.resize(n); np.ccs.resize(n);
    np.charge.resize(n); np.decoy.resize(n);
    np.modified_sequence.resize(n); np.protein_group.resize(n);
    np.transition_begin.resize(n); np.transition_count.resize(n);
    nt.product_mz.reserve(transitionCount());
    nt.library_intensity.reserve(transitionCount());
    nt.type.reserve(transitionCount());
    nt.ordinal.reserve(transitionCount());
    nt.charge.reserve(transitionCount());
    nt.loss.reserve(transitionCount());

    for (std::size_t i = 0; i < n; ++i)
    {
      const std::uint32_t src = order[i];
      np.mz[i] = precursors_.mz[src];
      np.irt[i] = precursors_.irt[src];
      np.im[i] = precursors_.im[src];
      np.ccs[i] = precursors_.ccs.empty() ? std::numeric_limits<float>::quiet_NaN()
                                          : precursors_.ccs[src];
      np.charge[i] = precursors_.charge[src];
      np.decoy[i] = precursors_.decoy[src];
      np.modified_sequence[i] = precursors_.modified_sequence[src];
      np.protein_group[i] = precursors_.protein_group[src];

      const std::uint32_t begin = precursors_.transition_begin[src];
      const std::uint32_t count = precursors_.transition_count[src];
      np.transition_begin[i] = static_cast<std::uint32_t>(nt.product_mz.size());
      np.transition_count[i] = count;
      for (std::uint32_t k = 0; k < count; ++k)
      {
        const std::uint32_t s = begin + k;
        nt.product_mz.push_back(transitions_.product_mz[s]);
        nt.library_intensity.push_back(transitions_.library_intensity[s]);
        nt.type.push_back(transitions_.type[s]);
        nt.ordinal.push_back(transitions_.ordinal[s]);
        nt.charge.push_back(transitions_.charge[s]);
        nt.loss.push_back(transitions_.loss[s]);
      }
    }

    precursors_ = std::move(np);
    transitions_ = std::move(nt);
    sorted_by_mz_ = true;
  }

  std::size_t Library::lowerBound(double mz_low) const
  {
    if (!sorted_by_mz_)
    {
      throw std::logic_error("Library::lowerBound called before sortByPrecursorMz");
    }
    const MzFixed key = toFixed(mz_low);
    const auto it = std::lower_bound(precursors_.mz.begin(), precursors_.mz.end(), key);
    return static_cast<std::size_t>(it - precursors_.mz.begin());
  }

  std::size_t Library::footprintBytes() const
  {
    auto vec = [](const auto& v) { return v.capacity() * sizeof(typename std::decay_t<decltype(v)>::value_type); };
    const auto& p = precursors_;
    const auto& t = transitions_;
    return vec(p.mz) + vec(p.irt) + vec(p.im) + vec(p.ccs) + vec(p.charge) + vec(p.decoy)
         + vec(p.modified_sequence) + vec(p.protein_group)
         + vec(p.transition_begin) + vec(p.transition_count)
         + vec(t.product_mz) + vec(t.library_intensity) + vec(t.type)
         + vec(t.ordinal) + vec(t.charge) + vec(t.loss)
         + strings_.footprintBytes();
  }

} // namespace ODIA
