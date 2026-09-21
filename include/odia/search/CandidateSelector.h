// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

/// Which precursors one search looks at, and their in-memory decoys.
///
/// M1 selection is a deterministic, PAIRED, label-blind random subset: every
/// eligible target is ranked by a hash of the key its decoy will share
/// ("<modified sequence>/<charge>"), salted by search:seed, and the lowest
/// search:subset keys are drawn. Decoys are then built for the drawn targets
/// only, in memory, with LibraryGenerator::appendDecoys on a subset library;
/// decoys already present in the library file are never searched. A target
/// whose decoy cannot be built leaves the search together with its would-be
/// decoy, so what is searched is always whole pairs.
///
/// Nothing here reads the run, and nothing reads a label to decide what is
/// kept: the draw key is identical for a target and its decoy.
#pragma once

#include <odia/Library.h>
#include <odia/search/SearchParams.h>

#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

namespace ODIA::search
{
  /// Library RT -> the assay scale [0, 100], by the INPUT library's own RT
  /// range (targets, finite values). One scale for the whole search, whatever
  /// the subset: it removes the 0..1 / iRT / minutes ambiguity of library RT
  /// columns without looking at the run.
  struct RtScale
  {
    double min = 0.0;
    double max = 0.0;
    double toAssay(double library_rt) const { return 100.0 * (library_rt - min) / (max - min); }
    double toLibrary(double assay_rt) const { return min + assay_rt * (max - min) / 100.0; }
  };

  /// What selection did, count by count. drawn = no_decoy + capped + pairs.
  struct SelectionStats
  {
    std::size_t library_precursors = 0;
    std::size_t library_decoys_ignored = 0;   ///< decoys already in the file: never searched
    std::size_t targets = 0;                  ///< target precursors in the library
    std::size_t ineligible_mz = 0;            ///< precursor m/z not representable
    std::size_t ineligible_charge = 0;        ///< charge 0
    std::size_t ineligible_rt = 0;            ///< library RT not finite
    std::size_t ineligible_fragments = 0;     ///< fewer than SearchParams::min_assay_fragments transitions
    std::size_t duplicate_key = 0;            ///< targets sharing (sequence, charge) with another: all excluded
    std::size_t eligible = 0;
    std::size_t drawn = 0;                    ///< after search:subset
    std::size_t no_decoy = 0;                 ///< drawn targets for which no decoy could be built
    std::size_t capped = 0;                   ///< pairs removed by search:max_pairs
    std::size_t pairs = 0;                    ///< searched target-decoy pairs
  };

  /// The precursors of one search: targets at [0, pairs()), their decoys at
  /// [pairs(), 2 * pairs()); target k and decoy pairs() + k form pair k. Both
  /// blocks are in the same order, by target precursor m/z and then draw key,
  /// so a chunk of pairs is an m/z-contiguous slice of the run's windows.
  ///
  /// Decoys carry their target's modified sequence, charge, protein group,
  /// precursor m/z, RT and 1/K0 (LibraryGenerator::appendDecoys); only their
  /// fragment m/z differ.
  struct SearchSet
  {
    Library library;                   ///< the searched precursors (a new library; the input is never modified)
    std::vector<std::size_t> source;   ///< per precursor: index of the target (for a decoy: of its target) in the input library
    std::vector<std::uint64_t> draw;   ///< per pair: its draw key
    RtScale rt_scale;
    SelectionStats stats;

    std::size_t size() const { return source.size(); }
    std::size_t pairs() const { return source.size() / 2; }
    bool isDecoy(std::size_t i) const { return i >= pairs(); }
    /// Pair id of precursor @p i, shared by a target and its decoy: 0 .. pairs() - 1.
    std::int64_t pairOf(std::size_t i) const { return static_cast<std::int64_t>(isDecoy(i) ? i - pairs() : i); }
    std::size_t partnerOf(std::size_t i) const { return isDecoy(i) ? i - pairs() : i + pairs(); }
    std::string_view modifiedSequence(std::size_t i) const;
    std::string_view proteinGroup(std::size_t i) const;
    int charge(std::size_t i) const;
    /// "<modified sequence><charge>", plus "_decoy" for a decoy: DIALibGen's
    /// Precursor.Id convention, unique within the set.
    std::string precursorId(std::size_t i) const;
  };

  class CandidateSelector
  {
  public:
    /// FNV-1a64 of "<modified sequence>/<charge>", seeded with @p seed.
    /// The same for a target and for its decoy, whatever the method.
    static std::uint64_t drawKey(std::string_view modified_sequence, int charge, std::uint64_t seed);

    /// The input library's RT range over finite target values. Throws
    /// std::invalid_argument when it is empty or a single point.
    static RtScale rtScale(const Library& library);

    /// Select pairs and build their decoys. Deterministic in (library content,
    /// params); independent of the library's row order. Throws
    /// std::invalid_argument on settings the selector cannot honour.
    static SearchSet select(const Library& library, const SearchParams& params);
  };
}
