// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <odia/search/AssayBuilder.h>

#include <algorithm>
#include <cmath>
#include <set>
#include <stdexcept>
#include <string>

namespace ODIA::search
{
  namespace
  {
    /// One compound and its transitions for library precursor @p i. Both
    /// build() and buildTargets() go through here, so the seed assays the
    /// calibration searches have exactly the layout the extraction uses.
    void appendAssay(OpenSwath::LightTargetedExperiment& exp, const Library& library, std::size_t i,
                     const std::string& id, const std::string& group_label, const std::string& protein_group,
                     double assay_rt, double drift_time, bool decoy)
    {
      const auto& p = library.precursors();
      const auto& t = library.transitions();
      OpenSwath::LightCompound c;
      c.id = id;
      // Sequence-blind on purpose, for BOTH classes: a DIALibGen decoy carries
      // its target's sequence, so any sequence-derived sub-score would score
      // the decoy as its target.
      c.sequence.clear();
      c.charge = p.charge[i];
      c.rt = assay_rt;
      c.drift_time = drift_time;
      c.protein_refs = {protein_group};
      c.peptide_group_label = group_label;

      const double precursor_mz = fromFixed(p.mz[i]);
      const std::uint32_t begin = p.transition_begin[i];
      for (std::uint32_t k = 0; k < p.transition_count[i]; ++k)
      {
        const std::uint32_t j = begin + k;
        OpenSwath::LightTransition tr;
        tr.transition_name = c.id + "_" + std::to_string(k);
        tr.peptide_ref = c.id;
        tr.library_intensity = t.library_intensity[j];
        tr.product_mz = fromFixed(t.product_mz[j]);
        tr.precursor_mz = precursor_mz;
        // Window assignment reads the TRANSITION's 1/K0, extraction the
        // compound's; they must agree.
        tr.precursor_im = c.drift_time;
        tr.fragment_charge = t.charge[j] == 0 ? 1 : t.charge[j];
        tr.decoy = decoy;
        tr.detecting_transition = true;
        tr.quantifying_transition = true;
        tr.identifying_transition = false;
        exp.transitions.push_back(std::move(tr));
      }
      exp.compounds.push_back(std::move(c));
    }
  }

  AssayChunk AssayBuilder::build(const SearchSet& set, const std::vector<std::size_t>& which, const AssayOptions& options)
  {
    const auto& p = set.library.precursors();
    AssayChunk out;
    out.precursor = which;
    auto& exp = out.experiment;
    exp.compounds.reserve(which.size());
    std::size_t transitions = 0;
    for (const std::size_t i : which)
    {
      if (i >= set.size()) { throw std::out_of_range("AssayBuilder: precursor " + std::to_string(i) + " of " + std::to_string(set.size())); }
      transitions += p.transition_count[i];
    }
    exp.transitions.reserve(transitions);
    std::set<std::string> proteins;

    for (const std::size_t i : which)
    {
      double drift = -1.0;
      if (options.ion_mobility)
      {
        double k0 = (i < p.im.size()) ? static_cast<double>(p.im[i]) : std::nan("");
        if (!std::isfinite(k0) && i < p.ccs.size() && std::isfinite(p.ccs[i]))
        { k0 = mobilityFromCCS(p.ccs[i], fromFixed(p.mz[i]), p.charge[i]); }
        if (std::isfinite(k0) && k0 > 0.0)
        {
          if (options.im_map) { k0 = options.im_map(k0); }
          drift = k0;
        }
        else { ++out.missing_im; }
      }
      const std::string pg(set.proteinGroup(i));
      proteins.insert(pg);
      appendAssay(exp, set.library, i, set.precursorId(i), set.precursorId(set.isDecoy(i) ? set.partnerOf(i) : i), pg,
                  set.rt_scale.toAssay(p.irt[i]), drift, set.isDecoy(i));
    }
    exp.proteins.reserve(proteins.size());
    for (const auto& pg : proteins) { exp.proteins.push_back(OpenSwath::LightProtein{pg, ""}); }
    return out;
  }

  OpenSwath::LightTargetedExperiment AssayBuilder::buildTargets(const Library& library, const std::vector<std::size_t>& which,
                                                                const RtScale& scale)
  {
    const auto& p = library.precursors();
    OpenSwath::LightTargetedExperiment exp;
    exp.compounds.reserve(which.size());
    std::set<std::string> proteins;
    for (const std::size_t i : which)
    {
      if (i >= library.precursorCount())
      { throw std::out_of_range("AssayBuilder: precursor " + std::to_string(i) + " of " + std::to_string(library.precursorCount())); }
      if (p.decoy[i]) { throw std::invalid_argument("AssayBuilder: calibration seeds are targets; precursor " + std::to_string(i) + " is a decoy"); }
      std::string id(library.strings().get(p.modified_sequence[i]));
      id += std::to_string(p.charge[i]);
      const std::string pg(library.strings().get(p.protein_group[i]));
      proteins.insert(pg);
      appendAssay(exp, library, i, id, id, pg, scale.toAssay(p.irt[i]), -1.0, false);
    }
    exp.proteins.reserve(proteins.size());
    for (const auto& pg : proteins) { exp.proteins.push_back(OpenSwath::LightProtein{pg, ""}); }
    return exp;
  }

  std::vector<std::vector<std::size_t>> AssayBuilder::chunks(const SearchSet& set, std::size_t precursors)
  {
    const std::size_t per_chunk = std::max<std::size_t>(1, precursors / 2);   // pairs per chunk
    std::vector<std::vector<std::size_t>> out;
    for (std::size_t a = 0; a < set.pairs(); a += per_chunk)
    {
      const std::size_t b = std::min(set.pairs(), a + per_chunk);
      std::vector<std::size_t> chunk;
      chunk.reserve(2 * (b - a));
      for (std::size_t k = a; k < b; ++k) { chunk.push_back(k); }
      for (std::size_t k = a; k < b; ++k) { chunk.push_back(set.pairs() + k); }
      out.push_back(std::move(chunk));
    }
    return out;
  }
}
