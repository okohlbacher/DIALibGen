// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

/// The search's in-memory decoys: one per selected target, built so that a
/// decoy differs from its target in fragment m/z ONLY, and in nothing a target
/// could never show.
///
/// The library generator's LibraryGenerator::appendDecoys is not used, for
/// three reasons measured on a full Orbitrap Astral run:
///   * its shuffle falls back to substituting a residue from a fixed table
///     when the interior cannot be rearranged (EEEEEEK, LLLLLLLLLLR); that
///     table is not ours to use, and the substituted decoy no longer matches
///     the precursor m/z it keeps. Here such a target has no decoy and leaves
///     the search with it (no_decoy, "unshufflable");
///   * it never range-checks the recomputed fragments, while every target
///     fragment lies in the generator's fragment window: 0.9 % of decoy
///     fragments fell where no target fragment can be, some outside the
///     acquired MS2 scan range, which makes decoys weaker than null targets
///     (anti-conservative). Here every decoy fragment must lie within the
///     m/z range of the library's TARGET fragments;
///   * it accepts any rearrangement whose residue string differs, so an I<->L
///     swap, or a permutation that leaves every selected fragment's
///     composition unchanged, gives a decoy with its target's fragment masses
///     (a copy, which ties with its target and loses true identifications).
///     Here such an arrangement is rejected.
/// A rejected arrangement is re-drawn (shuffle) up to a fixed number of times;
/// pseudo_reverse has one arrangement only. Every rule reads the target alone,
/// so a failure removes a whole pair, label-blind.
///
/// A fragment slot no decoy can reproduce (an unknown ion type or neutral loss)
/// is dropped from the target AND the decoy, so both assays keep the same slots.
#pragma once

#include <odia/Library.h>
#include <odia/LibraryGenerator.h>

#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

namespace ODIA::search
{
  /// Why a target got no search decoy.
  enum class DecoyOutcome : std::uint8_t
  {
    Made = 0,
    Unparsable,       ///< the sequence does not parse, or is too short to keep both termini and rearrange
    Unshufflable,     ///< no arrangement of the interior differs from the target's (EEEEEEK)
    OutOfRange,       ///< every arrangement tried puts a fragment outside the targets' fragment m/z range
    Copy,             ///< every arrangement tried reproduces the target's fragment masses (I = L)
    TooFewFragments   ///< fewer reproducible fragment slots than the assay minimum
  };

  const char* toString(DecoyOutcome o);

  struct DecoyRules
  {
    DecoyMethod method = DecoyMethod::Shuffle;   ///< Shuffle or PseudoReverse
    /// The m/z range every decoy fragment must lie in: that of the library's
    /// target fragments (targetFragmentRange), inclusive, in fixed-point units.
    MzFixed fragment_min = 0;
    MzFixed fragment_max = 0;
    /// A decoy whose every fragment lies within this of one of its target's
    /// fragments is a copy.
    double copy_ppm = 10.0;
    std::size_t min_fragments = 3;
    /// Shuffle: arrangements drawn before the target is given up.
    int arrangements = 50;
  };

  struct DecoyBuild
  {
    std::size_t made = 0;
    std::vector<DecoyOutcome> outcome;   ///< per target, in library order
    std::vector<std::uint16_t> redrawn;  ///< per target: arrangements rejected for range or copy before one was accepted
    std::size_t slots_dropped = 0;       ///< target fragment slots dropped from both classes (not reproducible)
  };

  /// The inclusive m/z range of @p library's target fragments (representable
  /// values only). Throws std::invalid_argument when there is none.
  std::pair<MzFixed, MzFixed> targetFragmentRange(const Library& library);

  /// Append one decoy per target of @p library, which must hold targets only,
  /// after them and in their order (decoy k belongs to the k-th target that
  /// got one). A decoy copies its target's modified sequence, charge, protein
  /// group, precursor m/z, RT, 1/K0, CCS, fragment slots and intensities; its
  /// fragment m/z are recomputed from the rearranged sequence. Deterministic
  /// in the sequence. Throws std::invalid_argument on a decoy in @p library or
  /// on a method other than Shuffle or PseudoReverse.
  DecoyBuild appendSearchDecoys(Library& library, const DecoyRules& rules);
}
