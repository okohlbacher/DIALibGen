// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <odia/search/SearchDecoys.h>

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>

namespace ODIA::search
{
  namespace
  {
    using OpenMS::AASequence;
    using OpenMS::Residue;

    /// Neutral-loss mass in Da; NaN for a loss no decoy can reproduce.
    double lossMass(LossType loss)
    {
      static const double water = OpenMS::EmpiricalFormula("H2O").getMonoWeight();
      static const double ammonia = OpenMS::EmpiricalFormula("NH3").getMonoWeight();
      static const double phospho = OpenMS::EmpiricalFormula("H3PO4").getMonoWeight();
      static const double metaphosphate = OpenMS::EmpiricalFormula("HPO3").getMonoWeight();
      static const double carbon_monoxide = OpenMS::EmpiricalFormula("CO").getMonoWeight();
      switch (loss)
      {
        case LossType::None: return 0.0;
        case LossType::Water: return water;
        case LossType::Ammonia: return ammonia;
        case LossType::Phospho: return phospho;
        case LossType::Metaphosphate: return metaphosphate;
        case LossType::CO: return carbon_monoxide;
        case LossType::Other: break;
      }
      return std::numeric_limits<double>::quiet_NaN();
    }

    /// One fragment slot of a target: where it sits in the transition arrays
    /// and how to recompute it for a rearranged sequence.
    struct Slot
    {
      std::uint32_t index;
      FragmentType type;
      std::size_t ordinal;
      int charge;
      double loss;
    };

    /// m/z of @p slot on @p s; NaN when it cannot be computed.
    double fragmentMz(const AASequence& s, const Slot& slot)
    {
      try
      {
        double mz = 0.0;
        switch (slot.type)
        {
          case FragmentType::A: mz = s.getPrefix(slot.ordinal).getMZ(slot.charge, Residue::AIon); break;
          case FragmentType::B: mz = s.getPrefix(slot.ordinal).getMZ(slot.charge, Residue::BIon); break;
          case FragmentType::C: mz = s.getPrefix(slot.ordinal).getMZ(slot.charge, Residue::CIon); break;
          case FragmentType::X: mz = s.getSuffix(slot.ordinal).getMZ(slot.charge, Residue::XIon); break;
          case FragmentType::Y: mz = s.getSuffix(slot.ordinal).getMZ(slot.charge, Residue::YIon); break;
          case FragmentType::Z: mz = s.getSuffix(slot.ordinal).getMZ(slot.charge, Residue::ZIon); break;
          default: return std::numeric_limits<double>::quiet_NaN();
        }
        return mz - slot.loss / slot.charge;
      }
      catch (const std::exception&) { return std::numeric_limits<double>::quiet_NaN(); }
    }

    bool reproducible(FragmentType t)
    {
      return t == FragmentType::A || t == FragmentType::B || t == FragmentType::C || t == FragmentType::X ||
             t == FragmentType::Y || t == FragmentType::Z;
    }

    struct Token
    {
      std::string text;
      bool modified;
    };

    enum class Try { Ok, Same, Range, Copy, Invalid };
  }

  const char* toString(DecoyOutcome o)
  {
    switch (o)
    {
      case DecoyOutcome::Made: return "made";
      case DecoyOutcome::Unparsable: return "unparsable";
      case DecoyOutcome::Unshufflable: return "unshufflable";
      case DecoyOutcome::OutOfRange: return "out_of_range";
      case DecoyOutcome::Copy: return "copy";
      case DecoyOutcome::TooFewFragments: return "too_few_fragments";
      case DecoyOutcome::Unpredictable: return "unpredictable";
    }
    return "?";
  }

  std::pair<MzFixed, MzFixed> targetFragmentRange(const Library& library)
  {
    const auto& p = library.precursors();
    const auto& t = library.transitions();
    MzFixed lo = std::numeric_limits<MzFixed>::max(), hi = 0;
    for (std::size_t i = 0; i < library.precursorCount(); ++i)
    {
      if (p.decoy[i]) { continue; }
      for (std::uint32_t k = 0; k < p.transition_count[i]; ++k)
      {
        const MzFixed v = t.product_mz[p.transition_begin[i] + k];
        if (v == MZ_INVALID) { continue; }
        lo = std::min(lo, v);
        hi = std::max(hi, v);
      }
    }
    if (hi == 0) { throw std::invalid_argument("search: the library's targets carry no representable fragment m/z"); }
    return {lo, hi};
  }

  DecoyAssay searchDecoy(const Library& library, std::size_t i, const DecoyRules& rules)
  {
    if (rules.method != DecoyMethod::Shuffle && rules.method != DecoyMethod::PseudoReverse)
    { throw std::invalid_argument("search decoys: only shuffle and pseudo_reverse keep both termini of the target"); }
    if (i >= library.precursorCount())
    { throw std::out_of_range("search decoys: precursor " + std::to_string(i) + " of " + std::to_string(library.precursorCount())); }
    if (library.precursors().decoy[i]) { throw std::invalid_argument("search decoys: precursor " + std::to_string(i) + " is a decoy"); }
    DecoyAssay out;
    const std::size_t keep_n = DECOY_KEEP_NTERM, keep_c = DECOY_KEEP_CTERM;
    {
      const std::string sequence(library.strings().get(library.precursors().modified_sequence[i]));
      AASequence target;
      std::vector<Token> tokens;
      try
      {
        target = AASequence::fromString(sequence);
        // Residue by residue, modifications attached, so a rearrangement moves
        // a modified residue whole. Terminal modifications stay terminal.
        for (OpenMS::Size k = 0; k < target.size(); ++k)
        {
          const auto& residue = target[k];
          auto single = AASequence::fromString(residue.getOneLetterCode());
          if (residue.isModified()) { single.setModification(0, residue.getModification()->getId()); }
          tokens.push_back({single.toString(), residue.isModified()});
        }
      }
      catch (const std::exception&) { out.outcome = DecoyOutcome::Unparsable; return out; }
      if (tokens.size() <= keep_n + keep_c + 1) { out.outcome = DecoyOutcome::Unparsable; return out; }

      // The slots a decoy can reproduce; the others go from both assays.
      const auto& pre = library.precursors();
      const auto& tr = library.transitions();
      const std::uint32_t begin = pre.transition_begin[i], count = pre.transition_count[i];
      std::vector<Slot> slots;
      std::vector<double> target_mz;
      for (std::uint32_t k = 0; k < count; ++k)
      {
        const std::uint32_t s = begin + k;
        const std::size_t ordinal = tr.ordinal[s];
        const double loss = lossMass(tr.loss[s]);
        if (!reproducible(tr.type[s]) || ordinal == 0 || ordinal >= tokens.size() || !std::isfinite(loss) ||
            tr.product_mz[s] == MZ_INVALID)
        { continue; }
        slots.push_back({s, tr.type[s], ordinal, tr.charge[s] == 0 ? 1 : static_cast<int>(tr.charge[s]), loss});
        target_mz.push_back(fromFixed(tr.product_mz[s]));
      }
      if (slots.size() < std::max<std::size_t>(1, rules.min_fragments))
      { out.outcome = DecoyOutcome::TooFewFragments; return out; }
      std::vector<double> sorted_target(target_mz);
      std::sort(sorted_target.begin(), sorted_target.end());

      std::string original;
      for (std::size_t k = keep_n; k + keep_c < tokens.size(); ++k) { original += tokens[k].text; }
      std::vector<MzFixed> decoy_mz(slots.size());
      std::string decoy_sequence;

      auto attempt = [&](const std::vector<Token>& mid) {
        std::string middle;
        for (const auto& tok : mid) { middle += tok.text; }
        if (middle == original) { return Try::Same; }
        std::string rearranged;
        for (std::size_t k = 0; k < keep_n; ++k) { rearranged += tokens[k].text; }
        rearranged += middle;
        for (std::size_t k = tokens.size() - keep_c; k < tokens.size(); ++k) { rearranged += tokens[k].text; }
        AASequence decoy;
        try
        {
          decoy = AASequence::fromString(rearranged);
          if (target.hasNTerminalModification()) { decoy.setNTerminalModification(target.getNTerminalModification()->getId()); }
          if (target.hasCTerminalModification()) { decoy.setCTerminalModification(target.getCTerminalModification()->getId()); }
        }
        catch (const std::exception&) { return Try::Invalid; }
        if (decoy == target) { return Try::Same; }
        bool copy = true;
        for (std::size_t k = 0; k < slots.size(); ++k)
        {
          const double mz = fragmentMz(decoy, slots[k]);
          const MzFixed fixed = std::isfinite(mz) ? toFixed(mz) : MZ_INVALID;
          if (fixed == MZ_INVALID) { return Try::Invalid; }
          if (fixed < rules.fragment_min || fixed > rules.fragment_max) { return Try::Range; }
          decoy_mz[k] = fixed;
          if (copy)
          {
            // Does any target fragment lie within copy_ppm of this one?
            const double tol = mz * rules.copy_ppm * 1e-6;
            const auto it = std::lower_bound(sorted_target.begin(), sorted_target.end(), mz - tol);
            copy = it != sorted_target.end() && *it <= mz + tol;
          }
        }
        if (!copy) { decoy_sequence = decoy.toString(); }
        return copy ? Try::Copy : Try::Ok;
      };

      std::vector<Token> mid(tokens.begin() + static_cast<std::ptrdiff_t>(keep_n),
                             tokens.end() - static_cast<std::ptrdiff_t>(keep_c));
      Try result = Try::Same;
      std::size_t differing = 0, range = 0, copies = 0, invalid = 0;
      if (rules.method == DecoyMethod::PseudoReverse)
      {
        std::reverse(mid.begin(), mid.end());
        result = attempt(mid);
        if (result != Try::Same && result != Try::Ok) { ++differing; }
        range += result == Try::Range;
        copies += result == Try::Copy;
        invalid += result == Try::Invalid;
      }
      else
      {
        // Deterministic in the sequence (FNV-1a), so the same library gives
        // the same decoys on every machine: the first arrangement is the one
        // LibraryGenerator's shuffle draws. Fisher-Yates with rejection
        // sampling, since std::shuffle is not specified across standard libraries.
        std::uint64_t seed = 14695981039346656037ull;
        for (const char ch : sequence) { seed = (seed ^ static_cast<std::uint8_t>(ch)) * 1099511628211ull; }
        std::mt19937_64 rng(seed);
        for (int a = 0; a < rules.arrangements; ++a)
        {
          for (std::size_t k = mid.size(); k > 1; --k)
          {
            const std::uint64_t bound = k;
            const std::uint64_t threshold = (std::uint64_t{0} - bound) % bound;
            std::uint64_t draw;
            do { draw = rng(); } while (draw < threshold);
            std::swap(mid[k - 1], mid[draw % bound]);
          }
          result = attempt(mid);
          if (result == Try::Ok) { break; }
          if (result == Try::Same) { continue; }
          ++differing;
          range += result == Try::Range;
          copies += result == Try::Copy;
          invalid += result == Try::Invalid;
        }
      }
      if (result != Try::Ok)
      {
        if (differing == 0) { out.outcome = DecoyOutcome::Unshufflable; }
        else if (invalid == differing) { out.outcome = DecoyOutcome::Unparsable; }
        else { out.outcome = range >= copies ? DecoyOutcome::OutOfRange : DecoyOutcome::Copy; }
        return out;
      }
      out.redrawn = static_cast<std::uint16_t>(std::min<std::size_t>(differing, 65535));
      out.slots.reserve(slots.size());
      out.charge.reserve(slots.size());
      for (const Slot& slot : slots)
      {
        out.slots.push_back(slot.index);
        out.charge.push_back(static_cast<std::int8_t>(slot.charge));
      }
      out.mz = std::move(decoy_mz);
      out.sequence = std::move(decoy_sequence);
    }
    return out;
  }

  DecoyBuild appendSearchDecoys(Library& library, const DecoyRules& rules)
  {
    if (rules.method != DecoyMethod::Shuffle && rules.method != DecoyMethod::PseudoReverse)
    { throw std::invalid_argument("search decoys: only shuffle and pseudo_reverse keep both termini of the target"); }
    const std::size_t n = library.precursorCount();
    for (std::size_t i = 0; i < n; ++i)
    {
      if (library.precursors().decoy[i]) { throw std::invalid_argument("search decoys: the library must hold targets only"); }
    }
    DecoyBuild out;
    out.outcome.assign(n, DecoyOutcome::Made);
    out.redrawn.assign(n, 0);

    // Each decoy is a function of its own target alone: computed in parallel,
    // applied in library order, so the result does not depend on threads.
    std::vector<DecoyAssay> assays(n);
    {
      const Library& read = library;
      const auto count = static_cast<std::ptrdiff_t>(n);
      std::string failure;
#pragma omp parallel for schedule(dynamic, 64)
      for (std::ptrdiff_t k = 0; k < count; ++k)
      {
        try { assays[static_cast<std::size_t>(k)] = searchDecoy(read, static_cast<std::size_t>(k), rules); }
        catch (const std::exception& e)
        {
#pragma omp critical(odia_search_decoys)
          { if (failure.empty()) { failure = e.what(); } }
        }
      }
      if (!failure.empty()) { throw std::runtime_error("search decoys: " + failure); }
    }

    auto& p = library.precursors();
    auto& t = library.transitions();
    for (std::size_t i = 0; i < n; ++i)
    {
      DecoyAssay& a = assays[i];
      out.outcome[i] = a.outcome;
      if (a.outcome != DecoyOutcome::Made) { continue; }
      out.redrawn[i] = a.redrawn;
      const std::uint32_t begin = p.transition_begin[i], count = p.transition_count[i];
      const std::size_t kept = a.slots.size();
      // Drop unreproducible slots from the target first, so both assays carry
      // the same slots in the same order: slot k then sits at begin + k.
      if (kept < count)
      {
        for (std::size_t k = 0; k < kept; ++k)
        {
          const std::uint32_t from = a.slots[k], to = begin + static_cast<std::uint32_t>(k);
          if (from == to) { continue; }
          t.product_mz[to] = t.product_mz[from];
          t.library_intensity[to] = t.library_intensity[from];
          t.type[to] = t.type[from];
          t.ordinal[to] = t.ordinal[from];
          t.charge[to] = t.charge[from];
          t.loss[to] = t.loss[from];
        }
        out.slots_dropped += count - kept;
        p.transition_count[i] = static_cast<std::uint32_t>(kept);
      }
      Library::checkTransitionCapacity(t.product_mz.size(), kept);
      const auto new_begin = static_cast<std::uint32_t>(t.product_mz.size());
      for (std::size_t k = 0; k < kept; ++k)
      {
        const std::uint32_t s = begin + static_cast<std::uint32_t>(k);
        t.product_mz.push_back(a.mz[k]);
        t.library_intensity.push_back(t.library_intensity[s]);
        t.type.push_back(t.type[s]);
        t.ordinal.push_back(t.ordinal[s]);
        t.charge.push_back(a.charge[k]);
        t.loss.push_back(t.loss[s]);
      }
      // Everything else is the target's, the precursor m/z included (the field
      // convention; a shuffle keeps the composition, so it is also correct).
      p.mz.push_back(p.mz[i]);
      p.irt.push_back(p.irt[i]);
      p.im.push_back(i < p.im.size() ? p.im[i] : std::numeric_limits<float>::quiet_NaN());
      p.ccs.push_back(i < p.ccs.size() ? p.ccs[i] : std::numeric_limits<float>::quiet_NaN());
      p.charge.push_back(p.charge[i]);
      p.decoy.push_back(1);
      p.modified_sequence.push_back(p.modified_sequence[i]);
      p.protein_group.push_back(p.protein_group[i]);
      p.transition_begin.push_back(new_begin);
      p.transition_count.push_back(static_cast<std::uint32_t>(kept));
      ++out.made;
      a = DecoyAssay();   // release as we go
    }
    if (out.made) { library.markUnsorted(); }
    return out;
  }
}
