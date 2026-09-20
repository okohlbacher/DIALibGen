// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <cmath>
#include <cstdint>
#include <limits>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace ODIA
{

  /// m/z stored as fixed-point in units of 1e-5 Th.
  ///
  /// Four bytes, like float, but the error is distributed better for this use.
  /// Tolerances in DIA are relative (ppm), and fixed-point error is constant in
  /// absolute terms, so the scale is chosen against the *bottom* of the m/z
  /// range: at 1e-5 Th the half-quantum is 0.033 ppm at m/z 150 and 0.0025 ppm
  /// at m/z 2000, against a flat ~0.03-0.05 ppm for float32. Range is
  /// 0 - 42,949 Th.
  using MzFixed = std::uint32_t;

  inline constexpr double MZ_QUANTUM = 1e-5;

  /// Largest representable m/z: 2^32 quanta of 1e-5 Th.
  inline constexpr double MZ_MAX = 42949.67295;

  /// Sentinel for "no valid m/z". Distinguishable from a real value, which is
  /// never 0 in practice.
  inline constexpr MzFixed MZ_INVALID = 0;

  /// Convert to fixed-point, or MZ_INVALID if the value cannot be represented.
  ///
  /// The domain check is not defensive padding: a NaN or out-of-range double
  /// converted by a bare cast is undefined behaviour, and on this toolchain it
  /// silently produced m/z 0, 42,448 and 2,764 Th from an empty field, "-500.1"
  /// and "1e9" respectively. Those precursors then sort to the ends of the
  /// library and never match a window, so a truncated or mis-columned input
  /// loads without complaint and quietly loses peptides.
  inline MzFixed toFixed(double mz)
  {
    if (!(mz > 0.0) || mz > MZ_MAX) { return MZ_INVALID; }   // false for NaN
    return static_cast<MzFixed>(mz / MZ_QUANTUM + 0.5);
  }

  inline constexpr double fromFixed(MzFixed v)
  {
    return static_cast<double>(v) * MZ_QUANTUM;
  }

  /// Nitrogen, the drift gas every timsTOF in reach runs, and the drift-tube
  /// temperature its calibration assumes.
  inline constexpr double DRIFT_GAS_MASS = 28.0134;
  inline constexpr double DRIFT_GAS_TEMPERATURE_K = 305.0;

  /// Mason-Schamp, reduced: 1/K0 = CCS * sqrt(mu * T) / (18509 * z), with
  /// mu the reduced mass of the ion and the drift gas.
  ///
  /// The two quantities are NOT interchangeable and this is not a unit
  /// conversion: CCS is a property of the ion, 1/K0 is what one instrument
  /// measures for it under its own calibration. ODIA used to emit angstroms
  /// and stop for exactly that reason. It now emits both, because the
  /// alternative measured worse: shipping no mobility at all costs a diaPASEF
  /// consumer the entire mobility dimension, whereas the derived value is good
  /// to 2.9% against 37,193 measured 1/K0 values on S08 -- and every consumer
  /// recalibrates mobility against its own run regardless.
  ///
  /// Fitted against those measurements, the coefficient comes out 1039.07
  /// where this formula gives 18509/sqrt(305) = 1059.82 -- 2.0% apart, which is
  /// inside the run-to-run calibration spread and is why the textbook constant
  /// is used rather than a number fitted to one file.
  inline double mobilityFromCCS(double ccs, double mz, int charge)
  {
    if (charge <= 0) { return std::numeric_limits<double>::quiet_NaN(); }
    const double m_ion = mz * charge;
    const double mu = m_ion * DRIFT_GAS_MASS / (m_ion + DRIFT_GAS_MASS);
    return ccs * std::sqrt(mu * DRIFT_GAS_TEMPERATURE_K) / (18509.0 * charge);
  }

  /// The inverse, for a library that carries 1/K0 but no CCS.
  inline double ccsFromMobility(double one_over_k0, double mz, int charge)
  {
    if (charge <= 0) { return std::numeric_limits<double>::quiet_NaN(); }
    const double m_ion = mz * charge;
    const double mu = m_ion * DRIFT_GAS_MASS / (m_ion + DRIFT_GAS_MASS);
    return one_over_k0 * 18509.0 * charge / std::sqrt(mu * DRIFT_GAS_TEMPERATURE_K);
  }

  enum class FragmentType : std::uint8_t
  {
    Unknown = 0, A, B, C, X, Y, Z, Precursor
  };

  enum class LossType : std::uint8_t
  {
    None = 0, Water, Ammonia, Phospho, Metaphosphate, CO, Other
  };

  FragmentType parseFragmentType(std::string_view s);
  LossType parseLossType(std::string_view s);
  std::string_view toString(FragmentType t);
  std::string_view toString(LossType l);

  /// Interned strings addressed by a 4-byte handle.
  ///
  /// The library's string fields repeat massively -- a peptide reference has
  /// ~7.1 M distinct values across ~78.6 M transitions, a protein group far
  /// fewer. Storing std::string per row is what produced ~471 M allocations and
  /// a 32.65 GB library load in the predecessor project; the strings themselves
  /// were a fifth of that, the rest was arena fragmentation caused by the
  /// allocation count. One arena plus 4-byte handles removes both.
  class StringArena
  {
  public:
    static constexpr std::uint32_t npos = 0xFFFFFFFFu;

    StringArena() = default;
    StringArena(StringArena&&) = default;
    StringArena& operator=(StringArena&&) = default;

    /// NOT copyable, and the compiler has to say so.
    ///
    /// `Entry::data` is a raw pointer into `blocks_`. A default copy allocates
    /// fresh blocks and then copies `entries_` VERBATIM, so every handle in the
    /// copy resolves into the SOURCE arena's storage. The copy is correct for
    /// exactly as long as the source outlives it and wrong the instant it does
    /// not -- silently, because the memory is still mapped until the source's
    /// blocks are freed.
    ///
    /// `Library::subsetByIndex` copied the arena this way. It survived because
    /// its only caller held the parent alive for the subset's whole lifetime.
    /// The first caller to write `library = library.subsetByIndex(...)` -- the
    /// `-min_library_fragments` filter -- destroyed the parent through the
    /// move-assign and took the arena with it, and the next `strings().get()`
    /// segfaulted inside memcpy. Re-intern into a fresh arena instead; that is
    /// also what makes a subset's arena hold only the strings the subset uses.
    StringArena(const StringArena&) = delete;
    StringArena& operator=(const StringArena&) = delete;

    /// Returns a handle for @p s, storing it only if it is new.
    std::uint32_t intern(std::string_view s);

    /// Drop the interning index, keeping the strings.
    ///
    /// `lookup_` exists only to answer "have I seen this string before" while a
    /// file is being read. Nothing consults it afterwards -- handles resolve
    /// through `entries_` -- but it is a member, so it lived for the library's
    /// whole lifetime.
    ///
    /// At 4.26 M distinct sequences that is one 28-byte node each plus a bucket
    /// array, ~150 MB of a ~1.03 GiB library retained to answer a question
    /// nobody asks again. Measured, not estimated: `footprintBytes()` already
    /// counts the map explicitly.
    ///
    /// After this, `intern()` would re-intern rather than deduplicate, so it
    /// must only be called when loading is finished. Calling it twice is safe.
    void releaseLookup();

    std::string_view get(std::uint32_t handle) const;

    std::size_t size() const { return entries_.size(); }
    std::size_t bytes() const { return bytes_; }

    /// Everything the arena holds, including the interning map and unused block
    /// capacity -- not just the characters. Reporting only bytes() understated
    /// the arena by 2.6x, which matters because D3's whole point is that the
    /// memory claim is measured rather than asserted.
    std::size_t footprintBytes() const;

    void reserve(std::size_t entries, std::size_t chars);

  private:
    /// Characters live in fixed blocks that are never reallocated.
    ///
    /// The lookup map is keyed by views into this storage, so the storage must
    /// have stable addresses. Backing it with a growing std::string silently
    /// breaks interning: on reallocation every key in the map dangles, lookups
    /// stop matching, and the same string is stored again -- which shows up as
    /// an inflated distinct-string count rather than as a crash.
    static constexpr std::size_t BLOCK = 1u << 20;

    struct Entry
    {
      const char* data;
      std::uint32_t length;
    };

    char* allocate(std::size_t n);

    std::vector<std::vector<char>> blocks_;
    std::size_t block_used_ = 0;
    std::size_t bytes_ = 0;
    std::vector<Entry> entries_;
    std::unordered_map<std::string_view, std::uint32_t> lookup_;
  };

  /// An assay library as parallel arrays.
  ///
  /// Layout is struct-of-arrays with a CSR relation: precursors own a
  /// contiguous run of transitions via (transition_begin, transition_count).
  /// Sorting precursors by m/z therefore makes each isolation window a
  /// contiguous slice rather than a gather, which is what lets the library be
  /// paged rather than held resident.
  ///
  /// Deliberately absent: a per-transition name. It is derived from the peptide
  /// reference and the fragment annotation, and the annotation is itself
  /// reconstructible from (type, ordinal, charge, loss), which cost 4 bytes
  /// here. Storing it would be ~78.5 M distinct strings.
  class Library
  {
  public:
    struct PrecursorArrays
    {
      std::vector<MzFixed> mz;
      std::vector<float> irt;              ///< library prediction, not an observation
      std::vector<float> im;               ///< 1/K0 as read from a file; NaN when absent

      /// Predicted collision cross-section, in square angstroms. NaN when
      /// absent.
      ///
      /// A separate field from `im` on purpose. CCS and the 1/K0 an instrument
      /// reports are different quantities, related by the Mason-Schamp equation
      /// through the drift gas and the instrument's calibration. ODIA does not
      /// convert -- that is done downstream, where the instrument is known --
      /// so writing angstroms into a field consumers read as 1/K0 would be
      /// exactly the silent unit error the separation avoids.
      std::vector<float> ccs;
      std::vector<std::uint8_t> charge;
      std::vector<std::uint8_t> decoy;
      std::vector<std::uint32_t> modified_sequence;  ///< StringArena handle
      std::vector<std::uint32_t> protein_group;      ///< StringArena handle
      std::vector<std::uint32_t> transition_begin;
      std::vector<std::uint32_t> transition_count;
    };

    struct TransitionArrays
    {
      std::vector<MzFixed> product_mz;
      std::vector<float> library_intensity;   ///< relative; f32 is more than it carries
      std::vector<FragmentType> type;
      std::vector<std::uint8_t> ordinal;
      std::vector<std::int8_t> charge;
      std::vector<LossType> loss;
    };

    std::size_t precursorCount() const { return precursors_.mz.size(); }
    std::size_t transitionCount() const { return transitions_.product_mz.size(); }
    std::size_t decoyCount() const;

    /// Precursors whose m/z could not be represented (see MZ_INVALID). They can
    /// never match an isolation window, so they are silently lost work unless
    /// reported.
    std::size_t invalidMzCount() const;

    /// Transitions whose m/z could not be represented. Reachable in more ways
    /// than the precursor case -- subtracting a neutral loss can drive a small
    /// fragment to or below zero.
    std::size_t invalidMzTransitionCount() const;

    const PrecursorArrays& precursors() const { return precursors_; }
    const TransitionArrays& transitions() const { return transitions_; }
    const StringArena& strings() const { return strings_; }

    // Clearing sorted_by_mz_ here would be wrong: the non-const overload is what
    // any read-only `auto& p = lib.precursors()` binds to on a non-const
    // Library, so the ordinary reading idiom would silently disarm lowerBound.
    // A caller that actually reorders must say so with markUnsorted().
    PrecursorArrays& precursors() { return precursors_; }
    TransitionArrays& transitions() { return transitions_; }
    StringArena& strings() { return strings_; }

    /// Reorder precursors by ascending m/z, rebuilding the CSR index.
    ///
    /// After this, the precursors of an isolation window are a contiguous range
    /// locatable by binary search.
    void sortByPrecursorMz();

    /// Drop every decoy precursor and its transitions, keeping targets in
    /// order. Returns how many were removed.
    ///
    /// Exists so a CACHED library can be re-decoyed with a different method
    /// without repeating the ~19 minutes of retention-time, fragment-intensity
    /// and CCS inference that produced its targets. None of that inference
    /// depends on how decoys are made.
    std::size_t dropDecoys();

    /// A new library holding only @p keep, in the order given, with their
    /// transitions gathered and the string arena carried over whole.
    ///
    /// Exists for the CiRT seed: a
    /// few hundred precursors can be searched blind over the entire gradient,
    /// which is exactly what is unaffordable for 10^7. Searching them means
    /// running the real extractor and picker over them, and those take a
    /// Library -- so the subset has to BE one, not a mask the callers of
    /// extraction would each have to honour.
    ///
    /// The arena is copied rather than re-interned: handles stay valid, the
    /// cost is one arena (~32 MB on our libraries) against the correctness of
    /// every string handle in the gathered precursors, and a seed library is
    /// short-lived.
    Library subsetByIndex(const std::vector<std::size_t>& keep) const;

    /// First precursor with m/z >= @p mz_low, for window slicing.
    ///
    /// Throws if the library is not sorted. Returning 0 instead would be
    /// indistinguishable from a legitimate answer, and every window would
    /// silently slice from the start of the library.
    std::size_t lowerBound(double mz_low) const;

    /// Declare that the arrays have been reordered behind our back.
    void markUnsorted() { sorted_by_mz_ = false; }

    bool isSortedByMz() const { return sorted_by_mz_; }

    /// Bytes held by the arrays and the arena, for measuring memory use.
    std::size_t footprintBytes() const;

    /// CSR offsets/counts are 32-bit: at most 4,294,967,295 transitions per library.
    /// Check before allocation or append; subtraction avoids overflow in current + additional.
    static void checkTransitionCapacity(std::size_t current, std::size_t additional = 0);

    void reserve(std::size_t precursors, std::size_t transitions);

    /// Release the growth slack in every array.
    ///
    /// Called after loading so footprintBytes() reports what is held rather
    /// than what std::vector happened to allocate on the way -- otherwise the
    /// figure varies with the input path rather than the data.
    void shrinkToFit();

  private:
    PrecursorArrays precursors_;
    TransitionArrays transitions_;
    StringArena strings_;
    bool sorted_by_mz_ = false;
  };

} // namespace ODIA
