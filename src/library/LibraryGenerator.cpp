// Copyright (c) 2026, Oliver Kohlbacher and the DIALibraryGenerator authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <odia/LibraryGenerator.h>
#include <odia/PeptDeepPredictor.h>

#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/ModifiedPeptideGenerator.h>
#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/FORMAT/FASTAFile.h>

#include <iostream>
#include <algorithm>
#include <iomanip>
#include <locale>
#include <sstream>
#include <thread>
#include <fstream>
#include <string_view>
#include <cmath>
#include <limits>
#include <map>
#include <random>
#include <set>
#include <stdexcept>
#include <unordered_map>

using namespace OpenMS;

namespace ODIA
{
  namespace
  {
    /// DIA-NN's substitution table: residue -> replacement.
    ///
    /// Fixed, not random, and reproduced here so the decoys match what the
    /// reference engine would produce.
    constexpr std::string_view MUTATE_FROM = "GAVLIFMPWSCTYHKRQEND";
    constexpr std::string_view MUTATE_TO   = "LLLVVLLLLTSSSSLLNDQE";
    static_assert(MUTATE_FROM.size() == MUTATE_TO.size(),
                  "the substitution tables must correspond position by position");

    /// True if @p list, a ';'-delimited set, already holds @p accession.
    bool containsAccession(std::string_view list, std::string_view accession)
    {
      std::size_t pos = 0;
      while (pos <= list.size())
      {
        const std::size_t sep = list.find(';', pos);
        const std::size_t end = sep == std::string_view::npos ? list.size() : sep;
        if (list.substr(pos, end - pos) == accession) { return true; }
        if (sep == std::string_view::npos) { break; }
        pos = sep + 1;
      }
      return false;
    }

    /// Neutral-loss mass for a label, in Da.
    ///
    /// Recomputing a decoy fragment from the bare ion series and then copying
    /// the target's loss label onto it leaves the m/z ~18 Da (H2O), ~17 (NH3)
    /// or ~98 (H3PO4) too high. That is the same anti-conservative asymmetry the
    /// shift-based version had -- target m/z come from the file and include the
    /// loss, decoy m/z would not -- restricted to loss-bearing transitions. The
    /// b/y-only fixture cannot show it; a phospho library is wrong by ~98 Th.
    double lossMass(LossType loss)
    {
      static const double water = EmpiricalFormula("H2O").getMonoWeight();
      static const double ammonia = EmpiricalFormula("NH3").getMonoWeight();
      static const double phospho = EmpiricalFormula("H3PO4").getMonoWeight();
      static const double metaphosphate = EmpiricalFormula("HPO3").getMonoWeight();
      static const double carbon_monoxide = EmpiricalFormula("CO").getMonoWeight();
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
      return std::nan("");   // unknown loss: caller drops the transition
    }

    char mutateResidue(char aa)
    {
      const auto pos = MUTATE_FROM.find(aa);
      return pos == std::string_view::npos ? aa : MUTATE_TO[pos];
    }

    /// One fragment before it is written into the library arrays.
    struct Fragment
    {
      double mz;
      FragmentType type;
      std::uint8_t ordinal;
      std::int8_t charge;
    };

    void enumerateAllFragments(const AASequence& peptide, const DigestParams& p,
                               int precursor_charge, std::vector<Fragment>& out)
    {
      out.clear();
      const Size n = peptide.size();
      if (n < 2) { return; }

      // Up to the precursor's own charge, not one below it.
      //
      // Requiring the complementary fragment to keep a charge is the textbook
      // assumption, and it is what this did. Measured against DIA-NN's library
      // on the human proteome, it is also wrong for this purpose: 21.6% of the
      // fragments DIA-NN keeps for a 2+ precursor are themselves 2+, and 35.6%
      // for a 3+ precursor. Those were not ions we ranked poorly -- they were
      // ions we never enumerated, and they accounted for 44.7% of the assays
      // DIA-NN carried that we did not. The PeptDeep MS2 model predicts the
      // z2 channels for every precursor charge, so the restriction was ours
      // alone.
      const int max_z = std::min(p.max_fragment_charge, std::max(1, precursor_charge));
      for (int z = 1; z <= max_z; ++z)
      {
        for (Size i = 1; i < n; ++i)
        {
          const double b = peptide.getPrefix(i).getMZ(z, Residue::BIon);
          if (b >= p.fragment_mz_min && b <= p.fragment_mz_max)
          {
            out.push_back({b, FragmentType::B, static_cast<std::uint8_t>(i),
                           static_cast<std::int8_t>(z)});
          }
          const double y = peptide.getSuffix(i).getMZ(z, Residue::YIon);
          if (y >= p.fragment_mz_min && y <= p.fragment_mz_max)
          {
            out.push_back({y, FragmentType::Y, static_cast<std::uint8_t>(i),
                           static_cast<std::int8_t>(z)});
          }
        }
      }

    }

    void enumerateFragments(const AASequence& peptide, const DigestParams& p,
                            int precursor_charge, std::vector<Fragment>& out)
    {
      enumerateAllFragments(peptide, p, precursor_charge, out);

      // Without predicted intensities there is no basis for ranking, so the cap
      // is applied by descending m/z: high-m/z fragments carry more sequence and
      // sit in a less crowded part of the spectrum. predictFragmentIntensities
      // replaces this, and must re-enumerate to do so -- the top twelve by m/z
      // are not the top twelve by intensity.
      std::sort(out.begin(), out.end(), [](const Fragment& a, const Fragment& b) {
        if (a.mz != b.mz) { return a.mz > b.mz; }
        if (a.type != b.type) { return a.type < b.type; }
        if (a.ordinal != b.ordinal) { return a.ordinal < b.ordinal; }
        return a.charge < b.charge;
      });
      if (out.size() > p.max_fragments) { out.resize(p.max_fragments); }
    }
  } // namespace

  DecoyMethod parseDecoyMethod(const std::string& s)
  {
    if (s == "none") { return DecoyMethod::None; }
    if (s == "pseudo_reverse") { return DecoyMethod::PseudoReverse; }
    if (s == "reverse") { return DecoyMethod::Reverse; }
    if (s == "shuffle") { return DecoyMethod::Shuffle; }
    return DecoyMethod::Mutate;
  }

  std::string_view toString(DecoyMethod m)
  {
    switch (m)
    {
      case DecoyMethod::None: return "none";
      case DecoyMethod::Mutate: return "mutate";
      case DecoyMethod::PseudoReverse: return "pseudo_reverse";
      case DecoyMethod::Reverse: return "reverse";
      case DecoyMethod::Shuffle: return "shuffle";
    }
    return "mutate";
  }

  LibraryGenerator::Stats LibraryGenerator::generate(const std::string& fasta_file,
                                                     const DigestParams& params,
                                                     Library& library)
  {
    Stats stats;

    std::vector<FASTAFile::FASTAEntry> entries;
    FASTAFile().load(fasta_file, entries);
    stats.proteins = entries.size();

    ProteaseDigestion digestion;
    digestion.setEnzyme(params.enzyme);
    digestion.setMissedCleavages(params.missed_cleavages);

    // Peptides are deduplicated across proteins; a peptide seen in several
    // proteins accumulates them, which is what protein inference later needs.
    std::map<std::string, std::string> peptide_to_proteins;
    std::vector<AASequence> peptides;
    for (const auto& entry : entries)
    {
      AASequence protein;
      try { protein = AASequence::fromString(entry.sequence); }
      catch (const std::exception&) { continue; }  // non-standard residues

      peptides.clear();
      digestion.digest(protein, peptides, params.min_length, params.max_length);

      // The initiator methionine is cleaved from most proteins before they are
      // ever seen, so the N-terminal peptide observed in a spectrum usually
      // starts one residue in. Digesting the truncated sequence as well is what
      // produces it. Only the peptides that actually differ are added: every
      // peptide not spanning position 1 is identical in both digests, and
      // adding it twice would inflate nothing but the work.
      if (params.n_terminal_methionine_excision && !entry.sequence.empty() &&
          entry.sequence.front() == 'M')
      {
        try
        {
          const auto excised = AASequence::fromString(entry.sequence.substr(1));
          std::vector<AASequence> more;
          digestion.digest(excised, more, params.min_length, params.max_length);
          // The first peptide of the excised digest is the only one that can be
          // new; the rest repeat the untruncated digest exactly.
          if (!more.empty()) { peptides.push_back(more.front()); }
        }
        catch (const std::exception&) {}   // non-standard residues, as above
      }

      for (const auto& pep : peptides)
      {
        auto& proteins = peptide_to_proteins[pep.toUnmodifiedString()];
        // Compare on ';' boundaries, not by substring: a plain find() drops P1
        // when P12 is already present, which is routine with isoform accessions
        // (P12345 / P12345-2) and with the sp|X|Y form we emit.
        if (proteins.empty()) { proteins = entry.identifier; }
        else if (!containsAccession(proteins, entry.identifier))
        {
          proteins += ';';
          proteins += entry.identifier;
        }
      }
    }
    stats.peptides = peptide_to_proteins.size();

    // StringList is std::vector<OpenMS::String>, NOT std::vector<std::string>:
    // the two are distinct types and no implicit conversion exists between the
    // vectors even though String converts to and from std::string element-wise.
    // Built explicitly so this compiles against a stock OpenMS release and not
    // only against a tree where the alias happens to be spelled differently.
    auto toStringList = [](const std::vector<std::string>& v) {
      return StringList(v.begin(), v.end());
    };
    ModifiedPeptideGenerator::MapToResidueType fixed_map =
      ModifiedPeptideGenerator::getModifications(toStringList(params.fixed_modifications));
    ModifiedPeptideGenerator::MapToResidueType variable_map =
      ModifiedPeptideGenerator::getModifications(toStringList(params.variable_modifications));

    library.reserve(peptide_to_proteins.size() * params.charges.size(),
                    peptide_to_proteins.size() * params.charges.size() * params.max_fragments);

    std::vector<Fragment> fragments;
    std::vector<AASequence> forms;
    for (const auto& [sequence, proteins] : peptide_to_proteins)
    {
      AASequence base;
      try { base = AASequence::fromString(sequence); }
      catch (const std::exception&) { continue; }

      ModifiedPeptideGenerator::applyFixedModifications(fixed_map, base);

      forms.clear();
      forms.push_back(base);
      if (!params.variable_modifications.empty())
      {
        ModifiedPeptideGenerator::applyVariableModifications(
          variable_map, base, params.max_variable_modifications, forms, true);
      }

      const std::uint32_t protein_handle = library.strings().intern(proteins);

      for (const auto& form : forms)
      {
        const std::uint32_t sequence_handle =
          library.strings().intern(form.toString());

        for (int z : params.charges)
        {
          const double precursor_mz = form.getMZ(z);
          if (precursor_mz < params.precursor_mz_min || precursor_mz > params.precursor_mz_max)
          {
            ++stats.dropped_precursor_mz;
            continue;
          }

          enumerateFragments(form, params, z, fragments);
          if (fragments.size() < params.min_fragments)
          {
            ++stats.dropped_too_few_fragments;
            continue;
          }

          auto& p = library.precursors();
          auto& t = library.transitions();
          p.mz.push_back(toFixed(precursor_mz));
          p.irt.push_back(std::nanf(""));   // filled in by prediction
          p.im.push_back(std::nanf(""));
          p.ccs.push_back(std::nanf(""));   // likewise
          p.charge.push_back(static_cast<std::uint8_t>(z));
          p.decoy.push_back(0);
          p.modified_sequence.push_back(sequence_handle);
          p.protein_group.push_back(protein_handle);
          p.transition_begin.push_back(static_cast<std::uint32_t>(t.product_mz.size()));
          p.transition_count.push_back(static_cast<std::uint32_t>(fragments.size()));

          for (const auto& f : fragments)
          {
            t.product_mz.push_back(toFixed(f.mz));
            t.library_intensity.push_back(1.0f);  // placeholder until prediction
            t.type.push_back(f.type);
            t.ordinal.push_back(f.ordinal);
            t.charge.push_back(f.charge);
            t.loss.push_back(LossType::None);
          }
        }
      }
    }

    stats.precursors = library.precursorCount();
    stats.transitions = library.transitionCount();
    return stats;
  }


  namespace
  {
    /// Sessions to run inference across, from what the caller asked for.
    ///
    /// Capped rather than taken as given: each session holds its own copy of
    /// the weights plus an ONNX arena sized for a MAX_BATCH_ROWS batch, so
    /// memory grows about half a gigabyte per session while the throughput
    /// curve flattens. Measured on the MS2 stage of a 60k-precursor build:
    ///
    ///     sessions    MS2      peak RSS
    ///            1   409.7 s   0.86 GiB
    ///            8    50.2 s   5.02 GiB   8.2x
    ///           16    34.6 s   8.82 GiB  11.8x
    ///           32    32.5 s  16.11 GiB  12.6x
    ///
    /// 16 is the knee. Going on to 32 buys 6% more speed for 83% more memory,
    /// which on a shared node is how a library build starts failing in the
    /// allocator instead of finishing slightly later.
    void writeDeviceLine(const char* stage, const PeptDeepPredictor& p)
    {
      std::cerr << "[" << stage << "] provider "
                << (p.provider() == PeptDeepPredictor::Provider::CUDA ? "CUDA" : "CPU")
                << " device " << p.device()
                << " sessions " << p.sessionCount() << "\n";
    }

    int inferenceSessions(unsigned requested)
    {
      const unsigned hardware = std::max(1u, std::thread::hardware_concurrency());
      const unsigned want = requested == 0 ? hardware : requested;
      return static_cast<int>(std::min(want, 16u));
    }
  } // namespace

  float LibraryGenerator::freeCysteineRtOffset(std::size_t free_cysteines)
  {
    // Pooled medians of (observed - predicted)/slope over the two benchmarks.
    // Nothing beyond three: 3+ carried 58 peptides on Astral and none on S08,
    // and the effect has plainly saturated by then.
    switch (free_cysteines)
    {
      case 0:  return 0.0f;
      case 1:  return 0.0343f;
      case 2:  return 0.0645f;
      default: return 0.0580f;
    }
  }

  namespace
  {
    /// Cysteines carrying no modification -- the ones the RT model cannot see.
    std::size_t freeCysteines(const AASequence& peptide)
    {
      std::size_t n = 0;
      for (std::size_t i = 0; i < peptide.size(); ++i)
      {
        if (peptide[i].getOneLetterCode() == "C" && !peptide[i].isModified()) { ++n; }
      }
      return n;
    }
  } // namespace

  std::size_t LibraryGenerator::predictRetentionTimes(Library& library,
                                                      const std::string& rt_model_path,
                                                      bool prefer_gpu, unsigned sessions,
                                                      bool free_cysteine_rt_correction)
  {
    // One prediction per distinct sequence. The RT model has no charge input,
    // so predicting per precursor would repeat identical work for every charge
    // state of the same peptide.
    std::map<std::uint32_t, std::size_t> handle_to_slot;
    std::vector<AASequence> unique_peptides;
    std::vector<std::uint32_t> slot_handle;

    const auto& p = library.precursors();
    for (std::size_t i = 0; i < library.precursorCount(); ++i)
    {
      const auto handle = p.modified_sequence[i];
      if (handle_to_slot.count(handle)) { continue; }
      AASequence peptide;
      try
      {
        peptide = AASequence::fromString(std::string(library.strings().get(handle)));
      }
      catch (const std::exception&)
      {
        continue;   // left NaN, and counted below
      }
      handle_to_slot.emplace(handle, unique_peptides.size());
      unique_peptides.push_back(std::move(peptide));
      slot_handle.push_back(handle);
    }
    if (unique_peptides.empty()) { return library.precursorCount(); }

    PeptDeepPredictor predictor(rt_model_path, prefer_gpu, 1, inferenceSessions(sessions));
    // State the device: a silent CPU fallback turns a 4-minute stage
    // into a 25-minute one with nothing in the log to explain it.
    writeDeviceLine("RT", predictor);
    std::vector<PeptDeepPredictor::Failure> failures;
    const auto predicted = predictor.predictRT(unique_peptides, &failures);

    std::map<std::uint32_t, float> by_handle;
    for (std::size_t slot = 0; slot < predicted.size(); ++slot)
    {
      float value = predicted[slot];
      // Applied here, in the model's own units, so that any downstream affine
      // rescale to an iRT gauge carries it through unchanged.
      if (free_cysteine_rt_correction && !std::isnan(value))
      {
        value += freeCysteineRtOffset(freeCysteines(unique_peptides[slot]));
      }
      by_handle.emplace(slot_handle[slot], value);
    }

    std::size_t unpredicted = 0;
    auto& precursors = library.precursors();
    for (std::size_t i = 0; i < library.precursorCount(); ++i)
    {
      const auto it = by_handle.find(precursors.modified_sequence[i]);
      const float value = it == by_handle.end() ? std::nanf("") : it->second;
      precursors.irt[i] = value;
      if (std::isnan(value)) { ++unpredicted; }
    }
    // Appending nothing, but the accessor is non-const; the ordering is
    // unchanged, so restore the flag rather than forcing a needless re-sort.
    return unpredicted;
  }


  std::size_t LibraryGenerator::predictFragmentIntensities(
    Library& library, const std::string& ms2_model_path, const DigestParams& params,
    float nce, const std::string& instrument, bool prefer_gpu, unsigned sessions)
  {
    auto& p = library.precursors();
    auto& t = library.transitions();
    const std::size_t n = library.precursorCount();
    if (n == 0) { return 0; }

    PeptDeepPredictor predictor(ms2_model_path, prefer_gpu, 1, inferenceSessions(sessions));
    // State the device: a silent CPU fallback turns a 4-minute stage
    // into a 25-minute one with nothing in the log to explain it.
    writeDeviceLine("MS2", predictor);

    // The new transition arrays are built alongside the old ones and swapped in
    // at the end. Editing in place is not possible: a precursor's fragment
    // count changes, so every later precursor's CSR offset moves.
    Library::TransitionArrays built;
    built.product_mz.reserve(t.product_mz.size());
    built.library_intensity.reserve(t.product_mz.size());
    built.type.reserve(t.product_mz.size());
    built.ordinal.reserve(t.product_mz.size());
    built.charge.reserve(t.product_mz.size());
    built.loss.reserve(t.product_mz.size());

    std::vector<std::uint32_t> begin(n, 0), count(n, 0);

    // Blocked so that peak memory is set by the block, not by the proteome.
    constexpr std::size_t BLOCK = 20000;
    std::size_t unpredicted = 0;
    std::vector<Fragment> fragments;
    std::vector<std::pair<float, Fragment>> ranked;

    for (std::size_t base = 0; base < n; base += BLOCK)
    {
      const std::size_t last = std::min(base + BLOCK, n);

      std::vector<AASequence> peptides;
      std::vector<int> charges;
      peptides.reserve(last - base);
      charges.reserve(last - base);
      for (std::size_t i = base; i < last; ++i)
      {
        peptides.push_back(
          AASequence::fromString(std::string(library.strings().get(p.modified_sequence[i]))));
        charges.push_back(static_cast<int>(p.charge[i]));
      }

      std::vector<PeptDeepPredictor::Failure> failures;
      const auto spectra = predictor.predictMS2(peptides, charges, nce, instrument, &failures);

      for (std::size_t i = base; i < last; ++i)
      {
        const auto& spectrum = spectra[i - base];
        begin[i] = static_cast<std::uint32_t>(built.product_mz.size());

        if (spectrum.positions == 0)
        {
          // Prediction failed for this one. Keep what generate() chose rather
          // than dropping the precursor: an m/z-ranked assay is worse than a
          // predicted one but better than none, and it is counted.
          ++unpredicted;
          for (std::uint32_t k = 0; k < p.transition_count[i]; ++k)
          {
            const std::size_t j = p.transition_begin[i] + k;
            built.product_mz.push_back(t.product_mz[j]);
            built.library_intensity.push_back(t.library_intensity[j]);
            built.type.push_back(t.type[j]);
            built.ordinal.push_back(t.ordinal[j]);
            built.charge.push_back(t.charge[j]);
            built.loss.push_back(t.loss[j]);
          }
          count[i] = p.transition_count[i];
          continue;
        }

        // Re-enumerate without the m/z cap: the cap is what has to be redone.
        enumerateAllFragments(peptides[i - base], params, charges[i - base], fragments);

        const std::size_t residues = peptides[i - base].size();
        // Base-peak relative, so one absolute threshold cannot be right for a
        // strong precursor and wrong for a weak one.
        float peak = 0.0f;
        for (std::size_t q = 0; q < spectrum.positions; ++q)
        {
          for (std::size_t ch = 0; ch < 4; ++ch) { peak = std::max(peak, spectrum.at(q, ch)); }
        }
        const float floor = peak * static_cast<float>(params.min_relative_intensity);
        ranked.clear();
        for (const auto& f : fragments)
        {
          // Position and channel follow alphabase's layout: fragment position
          // q separates prefix q+1 from suffix nAA-q-1, so a b ion of ordinal
          // o sits at q = o-1 and a y ion of ordinal o at q = nAA-1-o.
          std::size_t position = 0;
          std::size_t channel = 0;
          if (f.type == FragmentType::B)
          {
            position = static_cast<std::size_t>(f.ordinal) - 1;
            channel = f.charge == 1 ? 0 : 1;
          }
          else
          {
            position = residues - 1 - static_cast<std::size_t>(f.ordinal);
            channel = f.charge == 1 ? 2 : 3;
          }
          if (position >= spectrum.positions) { continue; }
          const float intensity = spectrum.at(position, channel);
          // A predicted zero is a fragment the model says is not there, and a
          // predicted 1e-7 is the same statement in floating point. Ranking on
          // `> 0` filled spare slots with denormal noise (see
          // min_relative_intensity); the floor is relative to this precursor's
          // own base peak, because the spectrum is base-peak normalised.
          if (!(intensity > floor)) { continue; }
          ranked.emplace_back(intensity, f);
        }

        std::sort(ranked.begin(), ranked.end(),
                  [](const auto& a, const auto& b) {
                    if (a.first != b.first) { return a.first > b.first; }
                    // Ties broken deterministically, so the library does not
                    // depend on sort implementation.
                    if (a.second.mz != b.second.mz) { return a.second.mz > b.second.mz; }
                    if (a.second.type != b.second.type) { return a.second.type < b.second.type; }
                    if (a.second.ordinal != b.second.ordinal)
                    { return a.second.ordinal < b.second.ordinal; }
                    return a.second.charge < b.second.charge;
                  });
        // A quota for doubly-charged fragments, when asked for. Applied by
        // promoting the best unselected z2 ions over the weakest selected z1
        // ones, so the cap is unchanged and only its composition moves.
        if (params.reserved_doubly_charged > 0 && ranked.size() > params.max_fragments)
        {
          const std::size_t cap = params.max_fragments;
          const std::size_t want =
            std::min(params.reserved_doubly_charged, cap);
          std::size_t have = 0;
          for (std::size_t r = 0; r < cap; ++r)
          {
            if (ranked[r].second.charge >= 2) { ++have; }
          }
          for (std::size_t r = cap; r < ranked.size() && have < want; ++r)
          {
            if (ranked[r].second.charge < 2) { continue; }
            // The weakest singly-charged fragment inside the cap makes way.
            std::size_t victim = cap;
            for (std::size_t q = cap; q-- > 0;)
            {
              if (ranked[q].second.charge < 2) { victim = q; break; }
            }
            if (victim == cap) { break; }
            std::swap(ranked[victim], ranked[r]);
            ++have;
          }
        }
        if (ranked.size() > params.max_fragments) { ranked.resize(params.max_fragments); }

        // The same bar `enumerateFragments` applied BEFORE prediction. The
        // base-peak intensity floor above can prune a precursor below it, and
        // until now nothing re-checked, so the count committed here was
        // whatever survived. That is where the S08 library's 366,084 targets
        // with 0-2 fragments came from -- and `appendDecoys` then refused to
        // build decoys for them, because it does apply the bar. The result was
        // 7.3% of targets competing against a null that contained nothing like
        // them. Emptied rather than kept: a precursor below the bar its decoy
        // would have to clear is not searchable, and `-min_library_fragments`
        // drops it symmetrically at load.
        if (ranked.size() < params.min_fragments) { ranked.clear(); }

        for (const auto& [intensity, f] : ranked)
        {
          built.product_mz.push_back(toFixed(f.mz));
          built.library_intensity.push_back(intensity);
          built.type.push_back(f.type);
          built.ordinal.push_back(f.ordinal);
          built.charge.push_back(f.charge);
          built.loss.push_back(LossType::None);
        }
        count[i] = static_cast<std::uint32_t>(ranked.size());
      }
    }

    t = std::move(built);
    p.transition_begin = std::move(begin);
    p.transition_count = std::move(count);
    library.shrinkToFit();
    return unpredicted;
  }



  LibraryGenerator::IrtCalibration
  LibraryGenerator::fitIrtCalibration(const std::string& rt_model_path,
                                      const std::string& standards_file,
                                      bool prefer_gpu)
  {
    std::vector<AASequence> peptides;
    std::vector<double> known;
    std::ifstream in(standards_file);
    if (!in) { throw std::runtime_error("cannot read iRT standards: " + standards_file); }
    std::string line;
    while (std::getline(in, line))
    {
      if (line.empty() || line[0] == '#') { continue; }
      const auto tab = line.find('\t');
      if (tab == std::string::npos) { continue; }
      try
      {
        peptides.push_back(AASequence::fromString(line.substr(0, tab)));
        known.push_back(std::stod(line.substr(tab + 1)));
      }
      catch (const std::exception&) { continue; }
    }
    if (peptides.size() < 3)
    {
      throw std::runtime_error("iRT standards file has fewer than 3 usable peptides: " +
                               standards_file);
    }

    PeptDeepPredictor predictor(rt_model_path, prefer_gpu);
    // State the device: a silent CPU fallback turns a 4-minute stage
    // into a 25-minute one with nothing in the log to explain it.
    writeDeviceLine("IRT-FIT", predictor);
    // Tolerate a standard that cannot be encoded rather than abandoning the
    // calibration: with eleven points, losing one to an exotic modification
    // should cost precision, not the whole line. Their entries come back NaN
    // and are skipped below, which is what makes the "at least three survived"
    // check reachable.
    std::vector<PeptDeepPredictor::Failure> failures;
    const auto raw = predictor.predictRT(peptides, &failures);

    // Ordinary least squares of known iRT on the model's raw output. Eleven
    // points and one line; nothing here warrants a solver.
    double sx = 0, sy = 0, sxx = 0, sxy = 0;
    std::size_t n = 0;
    for (std::size_t i = 0; i < raw.size(); ++i)
    {
      if (std::isnan(raw[i])) { continue; }   // a standard we could not encode
      const double x = raw[i], y = known[i];
      sx += x; sy += y; sxx += x * x; sxy += x * y;
      ++n;
    }
    if (n < 3)
    {
      throw std::runtime_error("only " + std::to_string(n) + " of " +
                               std::to_string(peptides.size()) +
                               " iRT standards could be predicted; at least 3 are "
                               "needed to fit a line with a residual worth checking");
    }
    const double denom = static_cast<double>(n) * sxx - sx * sx;
    if (!(std::abs(denom) > 0.0))
    {
      // Every standard predicted the same value, so there is no line. This is
      // what a saturated or broken model looks like, and it must not silently
      // produce a slope of infinity.
      throw std::runtime_error("iRT standards all predicted the same retention time; "
                               "the model cannot be calibrated");
    }

    IrtCalibration c;
    c.slope = (static_cast<double>(n) * sxy - sx * sy) / denom;
    c.intercept = (sy - c.slope * sx) / static_cast<double>(n);
    c.peptides = n;
    for (std::size_t i = 0; i < raw.size(); ++i)
    {
      if (std::isnan(raw[i])) { continue; }
      c.max_abs_error = std::max(c.max_abs_error, std::abs(c.apply(raw[i]) - known[i]));
    }
    return c;
  }

  void LibraryGenerator::applyIrtCalibration(Library& library,
                                             const IrtCalibration& calibration)
  {
    auto& p = library.precursors();
    for (auto& v : p.irt)
    {
      // NaN means "not predicted" and must stay that way; a calibrated NaN is
      // still NaN, but going through the multiply would be a lie about having
      // a value.
      if (!std::isnan(v)) { v = static_cast<float>(calibration.apply(v)); }
    }
  }

  std::size_t LibraryGenerator::predictCollisionCrossSections(
    Library& library, const std::string& ccs_model_path, bool prefer_gpu,
    unsigned sessions, bool derive_mobility)
  {
    auto& p = library.precursors();
    const std::size_t n = library.precursorCount();
    p.ccs.assign(n, std::numeric_limits<float>::quiet_NaN());
    if (derive_mobility) { p.im.assign(n, std::numeric_limits<float>::quiet_NaN()); }
    if (n == 0) { return 0; }

    PeptDeepPredictor predictor(ccs_model_path, prefer_gpu, 1, inferenceSessions(sessions));
    // State the device: a silent CPU fallback turns a 4-minute stage
    // into a 25-minute one with nothing in the log to explain it.
    writeDeviceLine("CCS", predictor);

    // Blocked, as the MS2 pass is: the whole proteome at once would hold every
    // parsed AASequence live alongside the library.
    constexpr std::size_t BLOCK = 200000;
    std::size_t unpredicted = 0;
    for (std::size_t base = 0; base < n; base += BLOCK)
    {
      const std::size_t last = std::min(base + BLOCK, n);
      std::vector<AASequence> peptides;
      std::vector<int> charges;
      peptides.reserve(last - base);
      charges.reserve(last - base);
      for (std::size_t i = base; i < last; ++i)
      {
        peptides.push_back(
          AASequence::fromString(std::string(library.strings().get(p.modified_sequence[i]))));
        charges.push_back(static_cast<int>(p.charge[i]));
      }

      std::vector<PeptDeepPredictor::Failure> failures;
      const auto ccs = predictor.predictCCS(peptides, charges, &failures);
      for (std::size_t i = base; i < last; ++i)
      {
        p.ccs[i] = ccs[i - base];
        if (std::isnan(p.ccs[i])) { ++unpredicted; continue; }
        // Derive 1/K0 alongside. CCS stays the authoritative value -- this is
        // the nominal mobility a timsTOF would report for it, which is what a
        // diaPASEF consumer needs to use the mobility dimension at all.
        if (derive_mobility)
        {
          p.im[i] = static_cast<float>(
            mobilityFromCCS(p.ccs[i], fromFixed(p.mz[i]), p.charge[i]));
        }
      }
    }
    return unpredicted;
  }

  std::size_t LibraryGenerator::appendDecoys(Library& library, DecoyMethod method,
                                             std::size_t* skipped_out,
                                             std::size_t min_fragments,
                                             bool recompute_decoy_mz)
  {
    if (method == DecoyMethod::None) { return 0; }

    const std::size_t n_targets = library.precursorCount();

    // Targets that already have a decoy, keyed as the decoy stores them.
    // Without this a second call appends a duplicate set: the copies are
    // byte-identical rows sharing a synthesised Precursor.Id, so they merge
    // silently on reload and a precursor plus half its transitions disappear.
    std::set<std::pair<std::uint32_t, std::uint8_t>> already;
    for (std::size_t i = 0; i < n_targets; ++i)
    {
      if (library.precursors().decoy[i])
      {
        already.emplace(library.precursors().modified_sequence[i],
                        library.precursors().charge[i]);
      }
    }
    std::size_t made = 0;
    std::size_t skipped = 0;

    for (std::size_t i = 0; i < n_targets; ++i)
    {
      if (library.precursors().decoy[i]) { continue; }
      if (already.count({library.precursors().modified_sequence[i],
                         library.precursors().charge[i]})) { continue; }

      const std::string sequence(library.strings().get(library.precursors().modified_sequence[i]));

      // Tokenise into residues, each carrying its own modification suffix, so
      // mutation and reversal preserve modifications instead of dropping them.
      struct Token { std::string text; char residue; bool modified; };
      std::vector<Token> tokens;
      for (std::size_t k = 0; k < sequence.size(); ++k)
      {
        if (sequence[k] == '(' || sequence[k] == '[') { continue; }  // leading N-term mod
        Token tok{std::string(1, sequence[k]), sequence[k], false};
        while (k + 1 < sequence.size() && (sequence[k + 1] == '(' || sequence[k + 1] == '['))
        {
          const char close = sequence[k + 1] == '(' ? ')' : ']';
          std::size_t j = k + 1;
          while (j < sequence.size() && sequence[j] != close) { tok.text.push_back(sequence[j++]); }
          if (j < sequence.size()) { tok.text.push_back(sequence[j]); }
          tok.modified = true;
          k = j;
        }
        tokens.push_back(std::move(tok));
      }
      if (tokens.size() < 4) { ++skipped; continue; }

      std::string decoy_sequence;

      if (method == DecoyMethod::Mutate)
      {
        // DIA-NN picks positions near each terminus that are not modified;
        // mutating a modified residue would silently invalidate the modification.
        std::size_t n_pos = tokens.size();
        for (std::size_t k = 1; k + 1 < tokens.size(); ++k)
        {
          if (!tokens[k].modified) { n_pos = k; break; }
        }
        std::size_t c_pos = tokens.size();
        for (std::size_t k = tokens.size() - 2; k >= 1; --k)
        {
          if (!tokens[k].modified && k != n_pos) { c_pos = k; break; }
          if (k == 1) { break; }
        }
        if (n_pos == tokens.size() || c_pos == tokens.size()) { ++skipped; continue; }

        tokens[n_pos].text = std::string(1, mutateResidue(tokens[n_pos].residue));
        tokens[c_pos].text = std::string(1, mutateResidue(tokens[c_pos].residue));
        for (const auto& tok : tokens) { decoy_sequence += tok.text; }
      }
      else if (method == DecoyMethod::Reverse)
      {
        // The whole sequence, termini included. mProphet's original recipe.
        // Cheap and always defined, but a palindromic peptide maps to itself,
        // which the caller catches as a collision below.
        for (auto it = tokens.rbegin(); it != tokens.rend(); ++it)
        { decoy_sequence += it->text; }
      }
      else if (method == DecoyMethod::PseudoReverse)
      {
        // Interior reversed, BOTH termini fixed -- the OpenSWATH and DIA-NN
        // reading of "pseudo-reverse". Holding the C-terminus matters for a
        // tryptic library, where every peptide ends in K or R and reversing it
        // would make the decoy population trivially separable from the targets.
        const std::size_t keep_n = DECOY_KEEP_NTERM, keep_c = DECOY_KEEP_CTERM;
        if (tokens.size() <= keep_n + keep_c + 1) { ++skipped; continue; }
        std::vector<Token> reordered(tokens.begin(), tokens.begin() + keep_n);
        std::vector<Token> mid(tokens.begin() + keep_n, tokens.end() - keep_c);
        std::reverse(mid.begin(), mid.end());
        reordered.insert(reordered.end(), mid.begin(), mid.end());
        reordered.insert(reordered.end(), tokens.end() - keep_c, tokens.end());
        for (const auto& tok : reordered) { decoy_sequence += tok.text; }
      }
      else // Shuffle
      {
        // Permute the interior, termini fixed, and keep permuting until the
        // sequence actually changed. DIA-NN 2.x's Generic mode is this family.
        //
        // The permutation is DETERMINISTIC -- seeded from the sequence itself --
        // because a library is a cache key. A std::random_device here would make
        // two runs of the same config produce different decoys, so the
        // fingerprint would promise reproducibility the file cannot deliver.
        const std::size_t keep_n = DECOY_KEEP_NTERM, keep_c = DECOY_KEEP_CTERM;
        if (tokens.size() <= keep_n + keep_c + 1) { ++skipped; continue; }
        std::vector<Token> mid(tokens.begin() + keep_n, tokens.end() - keep_c);

        std::uint64_t seed = 1469598103934665603ull;      // FNV-1a of the sequence
        for (const char ch : sequence) { seed = (seed ^ static_cast<std::uint8_t>(ch)) * 1099511628211ull; }
        std::mt19937_64 rng(seed);

        std::vector<Token> best;
        for (int attempt = 0; attempt < 20; ++attempt)
        {
          std::shuffle(mid.begin(), mid.end(), rng);
          std::string candidate;
          for (const auto& tok : mid) { candidate += tok.text; }
          std::string original;
          for (auto it = tokens.begin() + keep_n; it != tokens.end() - keep_c; ++it)
          { original += it->text; }
          if (candidate != original) { best = mid; break; }
        }
        // A run of identical residues cannot be shuffled into anything else.
        // DIA-NN falls back to mutation in that case; so do we, rather than
        // emitting a decoy identical to its target.
        if (best.empty())
        {
          std::size_t n_pos = tokens.size();
          for (std::size_t k = keep_n; k + keep_c < tokens.size(); ++k)
          { if (!tokens[k].modified) { n_pos = k; break; } }
          if (n_pos == tokens.size()) { ++skipped; continue; }
          tokens[n_pos].text = std::string(1, mutateResidue(tokens[n_pos].residue));
          for (const auto& tok : tokens) { decoy_sequence += tok.text; }
        }
        else
        {
          std::vector<Token> reordered(tokens.begin(), tokens.begin() + keep_n);
          reordered.insert(reordered.end(), best.begin(), best.end());
          reordered.insert(reordered.end(), tokens.end() - keep_c, tokens.end());
          for (const auto& tok : reordered) { decoy_sequence += tok.text; }
        }
      }

      // Recompute every fragment from the decoy sequence rather than shifting
      // the target's.
      //
      // Shifting b ions by the N-terminal mass difference and y ions by the
      // C-terminal one is wrong whenever an ion spans *both* mutated residues,
      // which the long b and y ions always do -- and those are exactly the ions
      // enumerateFragments keeps, since it caps by descending m/z. Measured on
      // the human proteome before this change: 15.3% of decoy fragment m/z did
      // not correspond to the decoy sequence stored beside them, worst case
      // 76 Th out. The failure is worse than it sounds because it is asymmetric:
      // recomputing from the sequence disagrees on ~15% of decoy transitions and
      // 0% of target ones, and a criterion that behaves differently for the two
      // classes is the anti-conservative FDR mode D7 exists to avoid.
      AASequence decoy;
      try { decoy = AASequence::fromString(decoy_sequence); }
      catch (const std::exception&) { ++skipped; continue; }

      auto& p = library.precursors();
      auto& t = library.transitions();

      const std::uint32_t begin = p.transition_begin[i];
      const std::uint32_t count = p.transition_count[i];
      const std::uint32_t new_begin = static_cast<std::uint32_t>(t.product_mz.size());

      for (std::uint32_t k = 0; k < count; ++k)
      {
        const std::uint32_t s = begin + k;
        const auto charge = t.charge[s] == 0 ? std::int8_t{1} : t.charge[s];
        const Size ord = t.ordinal[s];
        if (ord == 0 || ord >= decoy.size()) { continue; }

        double mz = 0.0;
        try
        {
          switch (t.type[s])
          {
            case FragmentType::A: mz = decoy.getPrefix(ord).getMZ(charge, Residue::AIon); break;
            case FragmentType::B: mz = decoy.getPrefix(ord).getMZ(charge, Residue::BIon); break;
            case FragmentType::C: mz = decoy.getPrefix(ord).getMZ(charge, Residue::CIon); break;
            case FragmentType::X: mz = decoy.getSuffix(ord).getMZ(charge, Residue::XIon); break;
            case FragmentType::Y: mz = decoy.getSuffix(ord).getMZ(charge, Residue::YIon); break;
            case FragmentType::Z: mz = decoy.getSuffix(ord).getMZ(charge, Residue::ZIon); break;
            default:
              // Precursor and unrecognised types cannot be recomputed. Copying
              // the target's value would give the decoy a transition identical
              // to its target's, so drop it instead.
              continue;
          }
        }
        catch (const std::exception&) { continue; }

        const double loss = lossMass(t.loss[s]);
        if (std::isnan(loss)) { continue; }   // cannot be reproduced faithfully
        mz -= loss / charge;

        t.product_mz.push_back(toFixed(mz));
        t.library_intensity.push_back(t.library_intensity[s]);
        t.type.push_back(t.type[s]);
        t.ordinal.push_back(t.ordinal[s]);
        // Store the charge actually used, not the file's 0 placeholder, or the
        // decoy would carry a singly-charged mass under a charge-0 label.
        t.charge.push_back(charge);
        t.loss.push_back(t.loss[s]);
      }

      // The same fragment-count bar the target had to clear. A gate applied to
      // one class and not the other is the anti-conservative FDR mode D7 rule 2
      // names: decoys admitted below their target's bar score lower, so the
      // decoy distribution is not the null the targets were drawn against.
      const auto written = static_cast<std::uint32_t>(t.product_mz.size() - new_begin);
      if (written == 0 || written < min_fragments)
      {
        t.product_mz.resize(new_begin);
        t.library_intensity.resize(new_begin);
        t.type.resize(new_begin);
        t.ordinal.resize(new_begin);
        t.charge.resize(new_begin);
        t.loss.resize(new_begin);
        ++skipped;
        continue;
      }

      // THE PRECURSOR m/z. Inherited from the target by default, and that is a
      // real defect rather than an approximation like the ones below it.
      //
      // The fragment m/z above are recomputed from the decoy peptide; this one
      // is copied. Under `mutate` the decoy's composition genuinely differs, so
      // the stored value is the mass of a different molecule. Two consequences:
      // any MS1-envelope-only quantity is numerically IDENTICAL for a decoy and
      // its target and therefore has zero discriminative power, and the decoy is
      // extracted from the isolation window where its target's real signal lives,
      // so it acquires genuine chromatographic evidence it has no right to.
      // NOT an ODIA peculiarity, and an earlier version of this comment claimed
      // the opposite. Inheriting the precursor m/z is the FIELD CONVENTION:
      // DIA-NN moves fragment m/z only, and OpenSWATH's MRMDecoy ships
      // product_mz_shift = +20 Th with precursor_mz_shift = 0. Both keep the
      // target's precursor m/z deliberately.
      //
      // So the consequence -- MS1-envelope-only quantities carrying no
      // discriminative power -- is shared by every library-based DIA engine, and
      // is presumably why DIA-NN's MS1 evidence is nine CORRELATION channels
      // rather than an envelope match. That is an argument for making ODIA's MS1
      // evidence co-elution-based, not obviously an argument for recomputing.
      //
      // Recomputing also cuts the other way on the null: a decoy extracted where
      // its target's real signal lives is a decoy that must beat real
      // interference, while one relocated to an arbitrary window may face less.
      // Whether that is conservative or anti-conservative is an entrapment
      // question that nothing here measures.
      //
      // The 1/K0 has to move WITH it. On diaPASEF the extraction window is
      // two-dimensional, so a decoy whose m/z changes while its mobility does
      // not is a physically inconsistent ion that can fall outside the window
      // and be DELETED rather than relocated -- and deleting decoys
      // non-uniformly is an FDR change wearing a bug fix's clothes.
      //
      // `mobilityFromCCS` already exists for exactly this and is what fills the
      // targets' own 1/K0 at prediction time; the collision cross-section is a
      // property of the ion, so it is the quantity that legitimately carries
      // over to a rearranged sequence, and 1/K0 is then re-derived from it at
      // the new mass. A library that carries 1/K0 but no CCS is inverted first
      // with `ccsFromMobility`, so both library provenances behave the same.
      double new_mz = fromFixed(p.mz[i]);
      if (recompute_decoy_mz)
      {
        const double dz = decoy.getMZ(static_cast<int>(p.charge[i]));
        if (std::isfinite(dz) && dz > 0.0) { new_mz = dz; }
      }
      p.mz.push_back(toFixed(new_mz));
      p.irt.push_back(p.irt[i]);

      float new_im = p.im[i];
      if (recompute_decoy_mz && std::isfinite(p.im[i]))
      {
        const double ccs_i =
          (i < p.ccs.size() && std::isfinite(p.ccs[i]))
            ? static_cast<double>(p.ccs[i])
            : ccsFromMobility(p.im[i], fromFixed(p.mz[i]), p.charge[i]);
        const double k0 = mobilityFromCCS(ccs_i, new_mz, p.charge[i]);
        // Only on a finite, positive result: a decoy silently losing its
        // mobility would be distinguishable by a missing value, which is the
        // same FDR leak the CCS copy below exists to avoid.
        if (std::isfinite(k0) && k0 > 0.0) { new_im = static_cast<float>(k0); }
      }
      p.im.push_back(new_im);
      // A decoy has the target's composition rearranged, so its cross-section
      // is close to the target's but not identical. Copying is the same
      // approximation already made for iRT and the intensity pattern, and it
      // keeps decoys from being distinguishable by a missing value -- which
      // would be an FDR leak, not a cosmetic gap.
      p.ccs.push_back(i < p.ccs.size() ? p.ccs[i] : std::nanf(""));
      p.charge.push_back(p.charge[i]);
      p.decoy.push_back(1);
      // The decoy keeps the target's sequence, which is DIA-NN's own convention
      // and fixes two things at once. The mutation table is many-to-one, so
      // distinct targets collide on a decoy sequence (18 collisions per 33,075
      // targets in the fixture); with the mutated sequence stored, two decoys
      // synthesised the same Precursor.Id and the tool could not reload its own
      // output. It also makes the row self-consistent: the decoy already
      // inherits the target's precursor m/z by design (D7), so storing the
      // mutated sequence left the two disagreeing by up to 76 Th.
      p.modified_sequence.push_back(p.modified_sequence[i]);
      p.protein_group.push_back(p.protein_group[i]);
      p.transition_begin.push_back(new_begin);
      p.transition_count.push_back(written);
      ++made;
    }
    if (made) { library.markUnsorted(); }
    if (skipped_out) { *skipped_out = skipped; }
    return made;
  }

} // namespace ODIA

namespace ODIA
{
  std::string LibraryGenerator::fingerprintParams(const DigestParams& p,
                                                  const std::string& rt_model,
                                                  const std::string& ms2_model,
                                                  const std::string& ccs_model,
                                                  double nce,
                                                  const std::string& instrument,
                                                  bool irt_rescale)
  {
    std::ostringstream ps;
    ps.imbue(std::locale::classic());          // a cache key must not follow the locale
    ps.setf(std::ios::fixed); ps.precision(3);
    // Round-trip exact, so no two distinct doubles ever share a token.
    auto exact = [](double v) {
      std::ostringstream o;
      o.imbue(std::locale::classic());
      o << std::setprecision(17) << v;
      return o.str();
    };
    auto join = [](const std::vector<std::string>& v) {
      std::string j;
      for (const auto& m : v) { j += m; j += "."; }
      return j.empty() ? std::string("none") : j;
    };
    // v2: v1 omitted the enzyme and the variable-modification IDENTITIES (it
    // carried only their maximum count), so two libraries digested differently
    // shared a key. Bumping misses every v1 cache once, which is the safe
    // direction.
    ps << "v3"
       << ";enz=" << p.enzyme
       << ";len=" << p.min_length << "-" << p.max_length
       << ";mc=" << p.missed_cleavages
       << ";z=";
    for (const int z : p.charges) { ps << z << "."; }
    ps << ";pmz=" << p.precursor_mz_min << "-" << p.precursor_mz_max
       << ";fmz=" << p.fragment_mz_min << "-" << p.fragment_mz_max
       << ";fz=" << p.max_fragment_charge
       << ";frag=" << p.min_fragments << "-" << p.max_fragments
       << ";varmod=" << p.max_variable_modifications
       << ";varmods=" << join(p.variable_modifications)
       // Fixed modifications change every precursor and fragment mass. Without
       // this, flipping the alkylation silently reuses a library built with the
       // other one -- the stale-cache collision that would make a CAM-free run
       // reproduce the CAM-inclusive result (doc/27).
       << ";fixmod=" << join(p.fixed_modifications)
       << ";nme=" << (p.n_terminal_methionine_excision ? 1 : 0)
       // Changes the stored RT of every free-cysteine peptide, so a library
       // built with it is not the same library.
       << ";fcys=" << (p.free_cysteine_rt_correction ? 1 : 0)
       << ";im=" << (p.derive_ion_mobility ? 1 : 0)
       // NOT the stream's fixed(3): min_relative_intensity is routinely 1e-4,
       // which rendered as "0.000" -- the same token as 0.0, so a library that
       // kept every fragment and one that dropped the weak ones shared a key.
       << ";minint=" << exact(p.min_relative_intensity)
       << ";rdc=" << p.reserved_doubly_charged
       << ";rt=" << rt_model
       << ";frgmodel=" << ms2_model
       << ";ccsmodel=" << ccs_model
       << ";nce=" << nce
       << ";inst=" << instrument
       // The RT column's DOMAIN, not just its values: raw 0..1 model output and
       // iRT are different scales of the same numbers, and without this token a
       // raw library and an iRT library of identical content fingerprinted the
       // same. A consumer that reads the wrong one extracts from nowhere.
       << ";rtdomain=" << (irt_rescale ? "irt" : "raw");
    return ps.str();
  }
}
