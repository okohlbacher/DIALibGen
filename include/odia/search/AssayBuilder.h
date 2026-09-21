// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

/// SearchSet slices -> OpenSWATH assays (OpenSwath::LightTargetedExperiment).
///
/// Built per extraction chunk straight from the library arrays, never for the
/// whole library. Peptide sequences are left EMPTY for targets and decoys
/// alike: OpenSWATH then computes no sequence-derived sub-score, which keeps
/// the scoring symmetric although a DIALibGen decoy stores its target's
/// sequence. Protein references are the library's Protein.Group strings
/// verbatim (the tuner's cohort key).
#pragma once

#include <odia/search/CandidateSelector.h>

#include <OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>

#include <cstddef>
#include <functional>
#include <string>
#include <vector>

namespace ODIA::search
{
  struct AssayOptions
  {
    /// Give compounds and transitions a drift time (diaPASEF). Off: -1, which
    /// is OpenSWATH's "no ion mobility".
    bool ion_mobility = false;
    /// Library 1/K0 -> assay drift time, e.g. the inverted IM calibration.
    /// Identity when empty.
    std::function<double(double)> im_map;
  };

  struct AssayChunk
  {
    OpenSwath::LightTargetedExperiment experiment;
    std::vector<std::size_t> precursor;   ///< per compound: its SearchSet index
    /// Compounds left without a drift time although ion_mobility was on (no
    /// 1/K0 and no CCS in the library). On diaPASEF data OpenSWATH silently
    /// extracts nothing for them, so the caller must refuse or count them.
    std::size_t missing_im = 0;
  };

  class AssayBuilder
  {
  public:
    /// Compound id of precursor @p i: SearchSet::precursorId.
    static std::string compoundId(const SearchSet& set, std::size_t i) { return set.precursorId(i); }

    /// Assays for the precursors @p which (SearchSet indices), in that order.
    /// Compound: id, charge, rt = rt_scale.toAssay(library RT) on [0, 100],
    /// drift_time, sequence "", protein_refs {Protein.Group}, peptide_group_label
    /// = the pair's target id. One transition per library transition:
    /// "<id>_<k>", detecting and quantifying, not identifying. One LightProtein
    /// per distinct protein group of the chunk.
    static AssayChunk build(const SearchSet& set, const std::vector<std::size_t>& which,
                            const AssayOptions& options = AssayOptions());

    /// Split the set into chunks of at most @p precursors precursors (at least
    /// one pair each), every pair in the same chunk as its partner: chunk c
    /// holds the targets of pairs [a, b) followed by their decoys.
    static std::vector<std::vector<std::size_t>> chunks(const SearchSet& set, std::size_t precursors);
  };
}
