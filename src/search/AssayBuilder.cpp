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
  AssayChunk AssayBuilder::build(const SearchSet& set, const std::vector<std::size_t>& which, const AssayOptions& options)
  {
    const auto& p = set.library.precursors();
    const auto& t = set.library.transitions();
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
      OpenSwath::LightCompound c;
      c.id = set.precursorId(i);
      // Sequence-blind on purpose, for BOTH classes: a DIALibGen decoy carries
      // its target's sequence, so any sequence-derived sub-score would score
      // the decoy as its target.
      c.sequence.clear();
      c.charge = set.charge(i);
      c.rt = set.rt_scale.toAssay(p.irt[i]);
      c.drift_time = -1.0;
      if (options.ion_mobility)
      {
        double k0 = (i < p.im.size()) ? static_cast<double>(p.im[i]) : std::nan("");
        if (!std::isfinite(k0) && i < p.ccs.size() && std::isfinite(p.ccs[i]))
        { k0 = mobilityFromCCS(p.ccs[i], fromFixed(p.mz[i]), p.charge[i]); }
        if (std::isfinite(k0) && k0 > 0.0)
        {
          if (options.im_map) { k0 = options.im_map(k0); }
          c.drift_time = k0;
        }
        else { ++out.missing_im; }
      }
      const std::string pg(set.proteinGroup(i));
      c.protein_refs = {pg};
      proteins.insert(pg);
      c.peptide_group_label = set.precursorId(set.isDecoy(i) ? set.partnerOf(i) : i);

      const double precursor_mz = fromFixed(p.mz[i]);
      const bool decoy = set.isDecoy(i);
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
    exp.proteins.reserve(proteins.size());
    for (const auto& pg : proteins) { exp.proteins.push_back(OpenSwath::LightProtein{pg, ""}); }
    return out;
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
