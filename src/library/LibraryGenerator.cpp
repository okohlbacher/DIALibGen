// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <odia/LibraryGenerator.h>
#include <odia/Hardware.h>
#include <odia/IonLadder.h>
#include <odia/PeptDeepEncoder.h>
#include <odia/PeptDeepPredictor.h>
#include <odia/PeptideParser.h>
#include <odia/Trace.h>

#include "Parallel.h"

#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ModifiedPeptideGenerator.h>
#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/FORMAT/FASTAFile.h>

#include <iostream>
#include <algorithm>
#include <iomanip>
#include <locale>
#include <sstream>
#include <thread>
#include <fstream>
#include <future>
#include <string_view>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <map>
#include <random>
#include <set>
#include <stdexcept>
#include <system_error>
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
      const IonLadder ladder(peptide);
      for (int z = 1; z <= max_z; ++z)
      {
        for (Size i = 1; i < n; ++i)
        {
          const double b = ladder.b(peptide, i, z);
          if (b >= p.fragment_mz_min && b <= p.fragment_mz_max)
          {
            out.push_back({b, FragmentType::B, static_cast<std::uint8_t>(i),
                           static_cast<std::int8_t>(z)});
          }
          const double y = ladder.y(peptide, i, z);
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

  namespace
  {
    /// One protein's digest, before it is merged into the peptide map.
    struct DigestedProtein
    {
      bool invalid = false;               ///< the protein sequence did not parse
      std::size_t ambiguous = 0;          ///< peptides dropped for X/B/Z/J or no mass
      std::vector<std::string> peptides;  ///< accepted, in digest order
      std::vector<std::uint8_t> shard;    ///< the peptide map shard of each
      std::vector<std::uint8_t> sites;    ///< per peptide: ProteinSite bits of this occurrence
    };

    /// Where an occurrence of a peptide sits in its protein, for
    /// protein-terminal modifications. A peptide accumulates the bits of all
    /// its occurrences.
    enum ProteinSite : std::uint8_t
    {
      PROTEIN_N = 1,    ///< starts at the protein N-terminus, or right after an excised initiator Met
      PROTEIN_C = 2,    ///< ends at the protein C-terminus
      INTERNAL_N = 4,   ///< starts anywhere else
      INTERNAL_C = 8    ///< ends anywhere else
    };

    /// The peptide map, one shard of it: peptide -> ';'-separated accessions.
    using PeptideMap = std::unordered_map<std::string, std::string>;
    using PeptideItem = const PeptideMap::value_type*;

    /// Every item of @p shards, in key order: the order one std::map holding
    /// them all would iterate in.
    ///
    /// Bucketed by the first two bytes (unsigned, as std::string compares, a
    /// missing byte as 0, which sorts a prefix first), which is the key order
    /// between buckets; each bucket is then sorted on its own. Keys are
    /// unique across the shards, so the result does not depend on which shard
    /// held what or on the thread count.
    std::vector<PeptideItem> inKeyOrder(const std::vector<PeptideMap>& shards, unsigned threads)
    {
      constexpr std::size_t BUCKETS = std::size_t(1) << 16;
      const auto bucket = [](const std::string& key) {
        const std::size_t c0 = key.size() > 0 ? static_cast<unsigned char>(key[0]) : 0u;
        const std::size_t c1 = key.size() > 1 ? static_cast<unsigned char>(key[1]) : 0u;
        return (c0 << 8) | c1;
      };
      std::vector<std::vector<std::size_t>> at(shards.size(), std::vector<std::size_t>(BUCKETS, 0));
      parallel::forEach(shards.size(), threads, [&](std::size_t sh)
      { for (const auto& item : shards[sh]) { ++at[sh][bucket(item.first)]; } });
      std::vector<std::size_t> bucket_begin(BUCKETS + 1, 0);
      std::size_t total = 0;
      for (std::size_t b = 0; b < BUCKETS; ++b)
      {
        bucket_begin[b] = total;
        for (auto& counts : at)
        {
          const std::size_t c = counts[b];
          counts[b] = total;   // from now on: where this shard's next item of bucket b goes
          total += c;
        }
      }
      bucket_begin[BUCKETS] = total;
      std::vector<PeptideItem> items(total);
      parallel::forEach(shards.size(), threads, [&](std::size_t sh)
      { for (const auto& item : shards[sh]) { items[at[sh][bucket(item.first)]++] = &item; } });
      constexpr std::size_t TASK = 16;   // buckets; one first byte spans 256 of them
      parallel::forEach(BUCKETS / TASK, threads, [&](std::size_t task)
      {
        for (std::size_t b = task * TASK; b < (task + 1) * TASK; ++b)
        {
          std::sort(items.begin() + static_cast<std::ptrdiff_t>(bucket_begin[b]),
                    items.begin() + static_cast<std::ptrdiff_t>(bucket_begin[b + 1]),
                    [](PeptideItem a, PeptideItem b) { return a->first < b->first; });
        }
      });
      return items;
    }

    /// One precursor of one modified form, ready to append: its fragments
    /// are [first, first + count) of the peptide's `fragments`.
    struct EnumeratedPrecursor
    {
      double mz;
      int charge;
      std::uint32_t first;
      std::uint32_t count;
    };

    struct EnumeratedForm
    {
      std::string sequence;
      std::uint64_t variable_mods = 0;   ///< bit k: carries requested variable modification k
      std::uint64_t fixed_mods = 0;      ///< bit k: carries protein-terminal fixed modification k
      std::vector<EnumeratedPrecursor> precursors;
    };

    /// One peptide's forms and precursors. `parsed` false: skipped entirely,
    /// its protein group not even interned -- exactly what the serial loop did.
    struct EnumeratedPeptide
    {
      bool parsed = false;
      std::size_t dropped_precursor_mz = 0;
      std::size_t dropped_too_few_fragments = 0;
      std::size_t precursors = 0;          ///< over all forms
      std::vector<EnumeratedForm> forms;
      std::vector<Fragment> fragments;     ///< every precursor's, one after the other
    };

    /// Items per parallel task and per wave. A wave is computed in parallel
    /// and then merged in input order, so memory is bounded by the wave and
    /// the merge sees exactly the serial order.
    constexpr std::size_t TASK_ITEMS = 256;
    constexpr std::size_t TASKS_PER_THREAD = 8;

    /// Sequences per forEach grain when parsing (parsePeptide): a parse is
    /// about a microsecond, too little to hand out one index at a time.
    constexpr std::size_t PARSE_GRAIN = 256;

    /// Runs @p task on a pool worker, or on this thread, without the
    /// overlap, when no worker can be started (a process or cgroup thread
    /// limit). The result is the same either way.
    ///
    /// The handle waits for the task when it goes, exceptions or not: the
    /// task refers to the caller's locals.
    class Overlapped
    {
    public:
      Overlapped() = default;
      explicit Overlapped(std::future<void> f) : f_(std::move(f)) {}
      Overlapped(Overlapped&&) = default;
      Overlapped& operator=(Overlapped&& other)
      {
        if (f_.valid()) { f_.wait(); }
        f_ = std::move(other.f_);
        return *this;
      }
      ~Overlapped() { if (f_.valid()) { f_.wait(); } }
      bool valid() const { return f_.valid(); }
      void get() { f_.get(); }

    private:
      std::future<void> f_;
    };

    template <class Task>
    Overlapped overlap(const Task& task)
    {
      return Overlapped(parallel::async(task));
    }
  } // namespace

  LibraryGenerator::Stats LibraryGenerator::generate(const std::string& fasta_file,
                                                     const DigestParams& params,
                                                     Library& library,
                                                     unsigned threads)
  {
    threads = parallel::workers(threads);
    if (params.max_fragment_charge < 1 || params.max_fragment_charge > 2)
    { throw std::invalid_argument("max_fragment_charge must be between 1 and 2"); }
    std::set<int> charges(params.charges.begin(), params.charges.end());
    if (charges.empty() || *charges.begin() < 1 || *charges.rbegin() > 8 ||
        charges.size() != params.charges.size())
    { throw std::invalid_argument("precursor_charges must be distinct integers between 1 and 8"); }
    Stats stats;

    // Started on the way out (below) and declared before every local of the
    // digest, so it ends after their destructors have run: freeing the
    // peptide map and the parsed sequences is a phase of its own.
    trace::Scope release_phase;
    trace::Scope fasta_phase("fasta_read", 0);
    std::vector<FASTAFile::FASTAEntry> entries;
    FASTAFile().load(fasta_file, entries);
    stats.proteins = entries.size();
    fasta_phase.fields().rows = static_cast<std::int64_t>(entries.size());
    fasta_phase.end();

    trace::Scope digest_phase("digest", 0);

    // Peptides are deduplicated across proteins; a peptide seen in several
    // proteins accumulates them, which is what protein inference later needs.
    //
    // Digested in parallel and merged in parallel too: the peptide map is
    // split into shards by a hash of the peptide, every shard visits the
    // proteins in FASTA order, so each accession list is built exactly as by
    // one serial pass, and the keys are put in order at the end (inKeyOrder).
    //
    // At most 16 digest threads. The digest is cheap (about 9 CPU-s for a
    // whole human proteome) and beyond 16 threads it got no faster but burnt
    // CPU waiting: at 96 threads 98 CPU-s for the same 2.3 s (2.2 s and 22
    // CPU-s at 16).
    const unsigned digest_threads = std::min(threads, 16u);
    const std::size_t shards = digest_threads;
    std::vector<PeptideMap> peptide_to_proteins(shards);
    const auto digestOne = [&](const ProteaseDigestion& digestion, const FASTAFile::FASTAEntry& entry,
                               DigestedProtein& out)
    {
      AASequence protein;
      std::string sequence = entry.sequence;
      if (!sequence.empty() && sequence.back() == '*') { sequence.pop_back(); }
      try { protein = AASequence::fromString(sequence); }
      catch (const std::exception&) { out.invalid = true; return; }

      // The products as start/end positions: the same peptides, in the same
      // order, as the AASequence overload (which is built on this one), plus
      // where each sits in the protein.
      std::vector<AASequence> peptides;
      std::vector<std::uint8_t> sites;
      const auto digestSpans = [&](const AASequence& chain)
      {
        std::vector<std::pair<std::size_t, std::size_t>> spans;
        digestion.digest(chain, spans, params.min_length, params.max_length);
        for (const auto& [begin, end] : spans)
        {
          peptides.push_back(chain.getSubsequence(begin, static_cast<UInt>(end - begin)));
          sites.push_back(static_cast<std::uint8_t>((begin == 0 ? PROTEIN_N : INTERNAL_N) |
                                                    (end == chain.size() ? PROTEIN_C : INTERNAL_C)));
        }
      };
      digestSpans(protein);

      // Include the Met-excised N-terminus at every allowed missed-cleavage
      // count. Internal peptides repeat the original digest and are removed
      // by the peptide-to-protein map below. The excised chain's first
      // residue is the mature protein's N-terminus.
      if (params.n_terminal_methionine_excision && !entry.sequence.empty() &&
          entry.sequence.front() == 'M')
      {
        try { digestSpans(AASequence::fromString(sequence.substr(1))); }
        catch (const std::exception&) {}   // non-standard residues, as above
      }

      for (std::size_t j = 0; j < peptides.size(); ++j)
      {
        const auto& pep = peptides[j];
        // Ambiguous residues have no unique elemental composition. Keep the
        // other peptides from that protein; never guess a mass for X/B/Z/J.
        std::string residues = pep.toUnmodifiedString();
        if (residues.find_first_of("XBZJ") != std::string::npos)
        { ++out.ambiguous; continue; }
        try
        {
          if (!std::isfinite(pep.getMonoWeight()))
          { ++out.ambiguous; continue; }
        }
        catch (const std::exception&)
        { ++out.ambiguous; continue; }
        out.shard.push_back(static_cast<std::uint8_t>(std::hash<std::string>{}(residues) % shards));
        out.peptides.push_back(std::move(residues));
        out.sites.push_back(sites[j]);
      }
    };

    // Protein-terminal modifications need each peptide's position in its
    // proteins; collected only when one is requested, per shard of the
    // peptide map like the accessions (sitesOf() finds a peptide's).
    const auto proteinTerminal = [](const std::string& name) {
      const auto term = ModificationsDB::getInstance()->getModification(name)->getTermSpecificity();
      return term == ResidueModification::PROTEIN_N_TERM || term == ResidueModification::PROTEIN_C_TERM;
    };
    const bool track_sites =
      std::any_of(params.fixed_modifications.begin(), params.fixed_modifications.end(), proteinTerminal) ||
      std::any_of(params.variable_modifications.begin(), params.variable_modifications.end(), proteinTerminal);
    std::vector<std::unordered_map<std::string, std::uint8_t>> peptide_sites(track_sites ? shards : 0);
    const auto sitesOf = [&](const std::string& peptide) -> std::uint8_t {
      if (!track_sites) { return 0; }
      const auto& map = peptide_sites[std::hash<std::string>{}(peptide) % shards];
      const auto found = map.find(peptide);
      return found == map.end() ? 0 : found->second;
    };

    const std::size_t wave = TASK_ITEMS * TASKS_PER_THREAD * threads;
    std::vector<DigestedProtein> digested;
    for (std::size_t first = 0; first < entries.size(); first += wave)
    {
      const std::size_t count = std::min(wave, entries.size() - first);
      digested.assign(count, DigestedProtein{});
      parallel::forEach((count + TASK_ITEMS - 1) / TASK_ITEMS, digest_threads, [&](std::size_t task)
      {
        // One digestion object per task: it is not shared between threads.
        ProteaseDigestion digestion;
        digestion.setEnzyme(params.enzyme);
        digestion.setMissedCleavages(params.missed_cleavages);
        const std::size_t end = std::min(count, (task + 1) * TASK_ITEMS);
        for (std::size_t k = task * TASK_ITEMS; k < end; ++k)
        { digestOne(digestion, entries[first + k], digested[k]); }
      });

      for (std::size_t k = 0; k < count; ++k)
      {
        if (digested[k].invalid) { ++stats.dropped_proteins; continue; }
        stats.dropped_ambiguous_peptides += digested[k].ambiguous;
      }
      parallel::forEach(shards, digest_threads, [&](std::size_t sh)
      {
        auto& map = peptide_to_proteins[sh];
        for (std::size_t k = 0; k < count; ++k)
        {
          auto& protein = digested[k];
          if (protein.invalid) { continue; }
          const auto& accession = entries[first + k].identifier;
          for (std::size_t q = 0; q < protein.peptides.size(); ++q)
          {
            if (protein.shard[q] != sh) { continue; }
            if (track_sites) { peptide_sites[sh][protein.peptides[q]] |= protein.sites[q]; }
            // Each peptide is visited by its own shard only, so it can be moved.
            auto& proteins = map[std::move(protein.peptides[q])];
            // Compare on ';' boundaries, not by substring: a plain find() drops P1
            // when P12 is already present, which is routine with isoform accessions
            // (P12345 / P12345-2) and with the sp|X|Y form we emit.
            if (proteins.empty()) { proteins = accession; }
            else if (!containsAccession(proteins, accession))
            {
              proteins += ';';
              proteins += accession;
            }
          }
        }
      });
    }
    const std::vector<PeptideItem> ordered = inKeyOrder(peptide_to_proteins, threads);
    stats.peptides = ordered.size();
    digest_phase.fields().rows = static_cast<std::int64_t>(stats.peptides);
    digest_phase.end();

    // Precursors and their m/z-capped placeholder fragments.
    trace::Scope enumerate_phase("enumerate", 0);

    // StringList is std::vector<OpenMS::String>, NOT std::vector<std::string>:
    // the two are distinct types and no implicit conversion exists between the
    // vectors even though String converts to and from std::string element-wise.
    // Built explicitly so this compiles against a stock OpenMS release and not
    // only against a tree where the alias happens to be spelled differently.
    auto toStringList = [](const std::vector<std::string>& v) {
      return StringList(v.begin(), v.end());
    };
    // Protein-terminal modifications are applied here, from the positions the
    // digest recorded: OpenMS's ModifiedPeptideGenerator does not know where
    // a peptide sits in its protein and applies them nowhere, fixed or
    // variable. The residue and peptide-terminal ones go to OpenMS.
    std::vector<std::string> residue_fixed, residue_variable;
    std::vector<const ResidueModification*> protein_fixed, protein_variable;
    std::vector<const ResidueModification*> requested_variable;
    std::vector<std::string> protein_fixed_names;
    {
      const auto* db = ModificationsDB::getInstance();
      for (const auto& name : params.fixed_modifications)
      {
        if (!proteinTerminal(name)) { residue_fixed.push_back(name); continue; }
        protein_fixed.push_back(db->getModification(name));
        protein_fixed_names.push_back(name);
      }
      for (const auto& name : params.variable_modifications)
      {
        const auto* mod = db->getModification(name);
        requested_variable.push_back(mod);
        if (proteinTerminal(name)) { protein_variable.push_back(mod); }
        else { residue_variable.push_back(name); }
      }
    }
    ModifiedPeptideGenerator::MapToResidueType fixed_map =
      ModifiedPeptideGenerator::getModifications(toStringList(residue_fixed));
    ModifiedPeptideGenerator::MapToResidueType variable_map =
      ModifiedPeptideGenerator::getModifications(toStringList(residue_variable));
    if (requested_variable.size() > 64)
    { throw std::invalid_argument("at most 64 variable modifications can be requested"); }
    // Residue-less peptide-terminal variable modifications (origin 'X', such
    // as 'Acetyl (N-term)'). OpenMS's one-modification path places a variable
    // modification only on a residue it names, so at max_variable_modifications
    // 1 these were dropped; they are added here, one per form, after the
    // forms OpenMS made. Its path for 2 or more handles them.
    std::vector<const ResidueModification*> terminal_variable;
    for (const auto* mod : requested_variable)
    {
      const auto term = mod->getTermSpecificity();
      if ((term == ResidueModification::N_TERM || term == ResidueModification::C_TERM) && mod->getOrigin() == 'X' &&
          std::find(terminal_variable.begin(), terminal_variable.end(), mod) == terminal_variable.end())
      { terminal_variable.push_back(mod); }
    }
    const auto variableBits = [&](const AASequence& form) {
      std::uint64_t bits = 0;
      const auto mark = [&](const ResidueModification* mod) {
        if (mod == nullptr) { return; }
        for (std::size_t k = 0; k < requested_variable.size(); ++k)
        {
          if (requested_variable[k] == mod || requested_variable[k]->getFullId() == mod->getFullId())
          { bits |= std::uint64_t{1} << k; }
        }
      };
      mark(form.getNTerminalModification());
      mark(form.getCTerminalModification());
      for (Size k = 0; k < form.size(); ++k) { mark(form[k].getModification()); }
      return bits;
    };

    // A protein-terminal modification fits a peptide that has the site, is
    // free at that terminus and, for a residue-specific one, starts (ends)
    // with its residue.
    const auto fits = [](const ResidueModification* mod, const AASequence& form, std::uint8_t sites) {
      const bool n_term = mod->getTermSpecificity() == ResidueModification::PROTEIN_N_TERM;
      if (!(sites & (n_term ? PROTEIN_N : PROTEIN_C)) || form.empty()) { return false; }
      if (n_term ? form.hasNTerminalModification() : form.hasCTerminalModification()) { return false; }
      const char origin = mod->getOrigin();
      const auto& residue = form[n_term ? 0 : form.size() - 1];
      return origin == 'X' || (!residue.isModified() && residue.getOneLetterCode()[0] == origin);
    };
    const auto withTerminal = [](AASequence form, const ResidueModification* mod) {
      if (mod->getTermSpecificity() == ResidueModification::PROTEIN_N_TERM) { form.setNTerminalModification(mod); }
      else { form.setCTerminalModification(mod); }
      return form;
    };
    // Variable modifications a form carries beyond its base.
    const auto variableCount = [](const AASequence& form, const AASequence& base) {
      std::size_t n = (form.getNTerminalModification() != base.getNTerminalModification()) +
                      (form.getCTerminalModification() != base.getCTerminalModification());
      for (Size k = 0; k < form.size(); ++k) { n += form[k].getModification() != base[k].getModification(); }
      return n;
    };

    // Reject the estimate before multiplying: an extreme fragment cap must
    // not wrap size_t and bypass the 32-bit CSR limit in reserve().
    constexpr auto transition_limit = std::numeric_limits<std::uint32_t>::max();
    Library::checkTransitionCapacity(ordered.size());
    const auto estimated_precursors = ordered.size() * params.charges.size();
    if (params.max_fragments && estimated_precursors > transition_limit / params.max_fragments)
    { throw std::length_error("library exceeds the 32-bit transition limit (4294967295); split the input FASTA"); }
    library.reserve(estimated_precursors, estimated_precursors * params.max_fragments);

    // Enumerate in parallel -- parsing, modifications, precursor m/z and the
    // fragment ladders are the expensive part -- and append in map order:
    // the output offsets of a wave are a prefix sum over its peptides, the
    // arrays are filled in parallel at those offsets, and only the interning
    // of strings stays serial, in the order of the serial loop, on a helper
    // thread that overlaps the next wave's enumeration.
    const auto enumerateOne = [&](const std::string& sequence, EnumeratedPeptide& out)
    {
      AASequence base;
      try { base = AASequence::fromString(sequence); }
      catch (const std::exception&) { return; }
      out.parsed = true;

      ModifiedPeptideGenerator::applyFixedModifications(fixed_map, base);

      // Fixed protein-terminal modifications: the modified base where the
      // peptide is protein-terminal, and the unmodified one as well where it
      // also occurs elsewhere.
      const std::uint8_t sites = sitesOf(sequence);
      std::vector<AASequence> bases{base};
      std::vector<std::uint64_t> base_fixed{0};
      for (std::size_t m = 0; m < protein_fixed.size(); ++m)
      {
        const auto* mod = protein_fixed[m];
        const bool n_term = mod->getTermSpecificity() == ResidueModification::PROTEIN_N_TERM;
        const bool elsewhere = sites & (n_term ? INTERNAL_N : INTERNAL_C);
        for (std::size_t b = 0, end = bases.size(); b < end; ++b)
        {
          if (!fits(mod, bases[b], sites)) { continue; }
          if (elsewhere)
          {
            bases.push_back(withTerminal(bases[b], mod));
            base_fixed.push_back(base_fixed[b] | std::uint64_t{1} << m);
          }
          else
          {
            bases[b] = withTerminal(bases[b], mod);
            base_fixed[b] |= std::uint64_t{1} << m;
          }
        }
      }

      std::vector<AASequence> forms;
      std::vector<std::uint64_t> form_fixed;
      std::vector<std::uint64_t> merged;   ///< variable bits of a form that duplicated another
      for (std::size_t b = 0; b < bases.size(); ++b)
      {
        const AASequence& from = bases[b];
        const std::size_t first = forms.size();
        if (params.variable_modifications.empty()) { forms.push_back(from); }
        else
        {
          ModifiedPeptideGenerator::applyVariableModifications(
            variable_map, from, params.max_variable_modifications, forms, true);
          if (params.max_variable_modifications == 1)
          {
            for (const auto* mod : terminal_variable)
            {
              const bool n_term = mod->getTermSpecificity() == ResidueModification::N_TERM;
              if (n_term ? from.hasNTerminalModification() : from.hasCTerminalModification()) { continue; }
              AASequence form = from;
              if (n_term) { form.setNTerminalModification(mod); } else { form.setCTerminalModification(mod); }
              forms.push_back(std::move(form));
            }
          }
        }
        // Variable protein-terminal modifications, on every form of this base
        // that has room for one more variable modification.
        merged.resize(forms.size(), 0);
        for (const auto* mod : protein_variable)
        {
          for (std::size_t k = first, end = forms.size(); k < end; ++k)
          {
            if (!fits(mod, forms[k], sites) ||
                variableCount(forms[k], from) >= params.max_variable_modifications) { continue; }
            AASequence form = withTerminal(forms[k], mod);
            // 'Acetyl (N-term)' and 'Acetyl (Protein N-term)' write the same
            // form: keep one, credited to both.
            const std::string text = form.toString();
            std::size_t same = first;
            while (same < forms.size() && forms[same].toString() != text) { ++same; }
            if (same < forms.size()) { merged[same] |= variableBits(form); continue; }
            forms.push_back(std::move(form));
            merged.push_back(0);
          }
        }
        form_fixed.resize(forms.size(), base_fixed[b]);
      }

      std::vector<Fragment> fragments;
      for (std::size_t index = 0; index < forms.size(); ++index)
      {
        const auto& form = forms[index];
        EnumeratedForm& f = out.forms.emplace_back();
        f.sequence = form.toString();
        if (!requested_variable.empty()) { f.variable_mods = variableBits(form) | merged[index]; }
        f.fixed_mods = form_fixed[index];
        for (int z : params.charges)
        {
          const double precursor_mz = form.getMZ(z);
          if (precursor_mz < params.precursor_mz_min || precursor_mz > params.precursor_mz_max)
          {
            ++out.dropped_precursor_mz;
            continue;
          }

          enumerateFragments(form, params, z, fragments);
          if (fragments.size() < params.min_fragments)
          {
            ++out.dropped_too_few_fragments;
            continue;
          }
          f.precursors.push_back({precursor_mz, z, static_cast<std::uint32_t>(out.fragments.size()),
                                  static_cast<std::uint32_t>(fragments.size())});
          out.fragments.insert(out.fragments.end(), fragments.begin(), fragments.end());
          ++out.precursors;
        }
      }
    };


    std::vector<std::size_t> variable_counts(requested_variable.size(), 0);
    std::vector<std::size_t> fixed_counts(protein_fixed.size(), 0);
    auto& p = library.precursors();
    auto& t = library.transitions();
    auto& strings = library.strings();
    std::vector<EnumeratedPeptide> enumerated[2];   // the wave being interned, and the next
    std::vector<std::size_t> precursor_at, transition_at;
    std::size_t precursor_total = p.mz.size(), transition_total = t.product_mz.size();
    Overlapped interning;   // after the waves: waited for before they go
    for (std::size_t first = 0; first < ordered.size(); first += wave)
    {
      const std::size_t count = std::min(wave, ordered.size() - first);
      auto& items = enumerated[first / wave % 2];
      items.resize(count);
      parallel::forEach((count + TASK_ITEMS - 1) / TASK_ITEMS, threads, [&](std::size_t task)
      {
        const std::size_t end = std::min(count, (task + 1) * TASK_ITEMS);
        for (std::size_t k = task * TASK_ITEMS; k < end; ++k)
        {
          items[k] = EnumeratedPeptide{};   // the wave before last, freed in parallel
          enumerateOne(ordered[first + k]->first, items[k]);
        }
      });

      // Where every peptide's precursors and transitions go.
      precursor_at.resize(count);
      transition_at.resize(count);
      std::size_t precursors = precursor_total, transitions = transition_total;
      for (std::size_t k = 0; k < count; ++k)
      {
        const auto& peptide = items[k];
        precursor_at[k] = precursors;
        transition_at[k] = transitions;
        if (!peptide.parsed) { continue; }
        stats.dropped_precursor_mz += peptide.dropped_precursor_mz;
        stats.dropped_too_few_fragments += peptide.dropped_too_few_fragments;
        for (const auto& form : peptide.forms)
        {
          for (std::size_t v = 0; v < variable_counts.size(); ++v)
          { if (form.variable_mods >> v & 1) { variable_counts[v] += form.precursors.size(); } }
          for (std::size_t v = 0; v < fixed_counts.size(); ++v)
          { if (form.fixed_mods >> v & 1) { fixed_counts[v] += form.precursors.size(); } }
        }
        precursors += peptide.precursors;
        transitions += peptide.fragments.size();
      }
      Library::checkTransitionCapacity(transition_total, transitions - transition_total);

      // The previous wave's strings are interned: the arrays can grow.
      if (interning.valid()) { interning.get(); }
      p.mz.resize(precursors); p.irt.resize(precursors); p.im.resize(precursors); p.ccs.resize(precursors);
      p.charge.resize(precursors); p.decoy.resize(precursors);
      p.modified_sequence.resize(precursors); p.protein_group.resize(precursors);
      p.transition_begin.resize(precursors); p.transition_count.resize(precursors);
      t.product_mz.resize(transitions); t.library_intensity.resize(transitions); t.type.resize(transitions);
      t.ordinal.resize(transitions); t.charge.resize(transitions); t.loss.resize(transitions);

      // The string handles, serially and in order, beside everything else;
      // they are the only columns it writes.
      const auto* wave_items = &items;
      const std::size_t wave_first = first, wave_precursor = precursor_total;
      interning = overlap([&, wave_items, wave_first, count, wave_precursor]() {
        std::size_t i = wave_precursor;
        for (std::size_t k = 0; k < count; ++k)
        {
          const auto& peptide = (*wave_items)[k];
          if (!peptide.parsed) { continue; }
          const std::uint32_t protein_handle = strings.intern(ordered[wave_first + k]->second);
          for (const auto& form : peptide.forms)
          {
            const std::uint32_t sequence_handle = strings.intern(form.sequence);
            for (std::size_t q = 0; q < form.precursors.size(); ++q, ++i)
            {
              p.modified_sequence[i] = sequence_handle;
              p.protein_group[i] = protein_handle;
            }
          }
        }
      });

      parallel::forEach((count + TASK_ITEMS - 1) / TASK_ITEMS, threads, [&](std::size_t task)
      {
        const std::size_t end = std::min(count, (task + 1) * TASK_ITEMS);
        for (std::size_t k = task * TASK_ITEMS; k < end; ++k)
        {
          const auto& peptide = items[k];
          if (!peptide.parsed) { continue; }
          std::size_t i = precursor_at[k], j = transition_at[k];
          for (const auto& form : peptide.forms)
          {
            for (const auto& precursor : form.precursors)
            {
              p.mz[i] = toFixed(precursor.mz);
              p.irt[i] = std::nanf("");   // filled in by prediction
              p.im[i] = std::nanf("");
              p.ccs[i] = std::nanf("");   // likewise
              p.charge[i] = static_cast<std::uint8_t>(precursor.charge);
              p.decoy[i] = 0;
              p.transition_begin[i] = static_cast<std::uint32_t>(j);
              p.transition_count[i] = precursor.count;
              for (std::uint32_t q = 0; q < precursor.count; ++q, ++j)
              {
                const auto& f = peptide.fragments[precursor.first + q];
                t.product_mz[j] = toFixed(f.mz);
                t.library_intensity[j] = 1.0f;  // placeholder until prediction
                t.type[j] = f.type;
                t.ordinal[j] = f.ordinal;
                t.charge[j] = f.charge;
                t.loss[j] = LossType::None;
              }
              ++i;
            }
          }
        }
      });
      precursor_total = precursors;
      transition_total = transitions;
    }
    if (interning.valid()) { interning.get(); }

    stats.precursors = library.precursorCount();
    stats.transitions = library.transitionCount();
    for (std::size_t k = 0; k < variable_counts.size(); ++k)
    {
      stats.variable_modification_precursors.emplace_back(params.variable_modifications[k], variable_counts[k]);
      if (variable_counts[k] > 0) { continue; }
      // A modification that modifies nothing is a request the library does
      // not honour; until 0.11.1 that went by with exit code 0.
      std::string hint;
      if (params.max_variable_modifications == 0) { hint = " (max_variable_modifications is 0)"; }
      else { hint = " (does any peptide in the m/z and length windows have its site?)"; }
      throw std::invalid_argument("variable modification '" + params.variable_modifications[k] +
                                  "' produced no modified precursor in the whole library" + hint);
    }
    for (std::size_t k = 0; k < fixed_counts.size(); ++k)
    {
      if (fixed_counts[k] > 0) { continue; }
      throw std::invalid_argument("fixed modification '" + protein_fixed_names[k] +
                                  "' produced no modified precursor in the whole library (does any "
                                  "protein-terminal peptide in the m/z and length windows have its site?)");
    }
    enumerate_phase.fields().rows = static_cast<std::int64_t>(stats.precursors);
    enumerate_phase.end();
    release_phase.start("release", 0);
    // The peptide map and the last waves, freed in parallel; the rest goes on return.
    parallel::forEach(peptide_to_proteins.size(), digest_threads, [&](std::size_t sh)
    { PeptideMap().swap(peptide_to_proteins[sh]); });
    for (auto& items : enumerated)
    {
      parallel::forEach((items.size() + TASK_ITEMS - 1) / TASK_ITEMS, threads, [&](std::size_t task)
      {
        for (std::size_t k = task * TASK_ITEMS; k < std::min(items.size(), (task + 1) * TASK_ITEMS); ++k)
        { items[k] = EnumeratedPeptide{}; }
      });
    }
    return stats;
  }


  namespace
  {
    void writeDeviceLine(const char* stage, const PeptDeepPredictor& p, const GpuOptions& gpu)
    {
      const bool coreml = p.provider() == PeptDeepPredictor::Provider::CoreML;
      std::cerr << "[" << stage << "] provider "
                << (coreml ? "CoreML" : p.provider() == PeptDeepPredictor::Provider::CUDA ? "CUDA" : "CPU");
      if (coreml) { std::cerr << " units " << gpu.coreml_units << " concurrent " << p.sessionCount(); }
      else { std::cerr << " device " << p.device() << " sessions " << p.sessionCount(); }
      if (!p.providerNote().empty()) { std::cerr << " (" << p.providerNote() << ")"; }
      // Name every device when the sessions span several.
      const auto& devices = p.sessionDevices();
      if (std::any_of(devices.begin(), devices.end(), [&](int d) { return d != p.device(); }))
      {
        std::cerr << " on devices";
        for (std::size_t k = 0; k < devices.size(); ++k) { std::cerr << (k ? "," : " ") << devices[k]; }
      }
      std::cerr << "\n";
      // Loud, and recorded: the provider line alone was the only sign that RT
      // and CCS had silently left the GPU.
      if (!p.fallbackReason().empty())
      {
        const std::string warning = std::string(stage) + ": " + p.fallbackReason();
        std::cerr << "WARNING: " << warning << "\n";
        if (gpu.warnings != nullptr) { gpu.warnings->push_back(warning); }
      }
    }

    /// CPU inference sessions, from -threads: one per thread, 0 meaning one
    /// per physical core, and never more than the physical cores.
    ///
    /// Until 0.11.1 this was capped at 16. Each session then held an ONNX
    /// arena sized for a 2048-row batch (about 1.9 GB for MS2), and a call
    /// waited for its largest batch, so sessions beyond 16 bought little speed
    /// for a lot of memory. With 64-row batches (PeptDeepPredictor::
    /// CPU_CHUNK_ROWS) an MS2 arena is about 53 MB and every session stays
    /// busy: on 2 x 96 cores, 96 sessions kept 94 busy through MS2. More
    /// sessions than physical cores share cores and are slower (Hardware.h).
    int inferenceSessions(unsigned requested)
    {
      const unsigned cores = physicalCores();
      const unsigned want = requested == 0 ? cores : std::min(requested, cores);
      return static_cast<int>(std::max(1u, want));
    }

    /// Host threads for the work that runs beside a prediction call: parsing
    /// the next block, and selecting and appending the previous one. On the
    /// CPU every session already holds a core, and the helper has a whole
    /// block's prediction to finish in, so it gets an eighth of the threads
    /// instead of a second full team competing with the sessions. On CUDA the
    /// cores are mostly free while the device computes.
    unsigned overlapThreads(const PeptDeepPredictor& predictor, unsigned threads)
    {
      if (predictor.provider() != PeptDeepPredictor::Provider::CPU) { return std::max(1u, threads / 2); }
      return std::max(1u, threads / 8);
    }
  } // namespace

  float LibraryGenerator::freeCysteineRtOffset(std::size_t free_cysteines)
  {
    // Pooled medians of (observed - predicted)/slope over the two benchmarks.
    // Nothing beyond three: 3+ carried 58 peptides on Astral and none on IH1,
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
                                                      bool free_cysteine_rt_correction,
                                                      const GpuOptions& gpu)
  {
    // One prediction per distinct sequence. The RT model has no charge input,
    // so predicting per precursor would repeat identical work for every charge
    // state of the same peptide.
    trace::Scope rt_phase("rt", 0);
    trace::Scope unique_step("unique_parse", 1);
    const auto& p = library.precursors();
    const unsigned threads = parallel::workers(sessions);
    // Distinct sequences in order of first appearance, then parsed in
    // parallel and compacted in that order: the same peptides in the same
    // order as one serial pass, which the batching depends on. Handles are
    // dense indices into the string arena, so flat arrays replace the maps.
    std::vector<std::uint8_t> seen(library.strings().size(), 0);
    std::vector<std::uint32_t> candidates;
    for (std::size_t i = 0; i < library.precursorCount(); ++i)
    {
      const auto handle = p.modified_sequence[i];
      if (seen[handle]) { continue; }
      seen[handle] = 1;
      candidates.push_back(handle);
    }
    std::vector<AASequence> parsed(candidates.size());
    std::vector<std::uint8_t> parsed_ok(candidates.size(), 0);
    parallel::forEach(candidates.size(), threads, [&](std::size_t k)
    {
      try
      {
        parsed[k] = parsePeptide(library.strings().get(candidates[k]));
        parsed_ok[k] = 1;
      }
      catch (const std::exception&) {}   // left NaN, and counted below
    }, PARSE_GRAIN);
    std::vector<AASequence> unique_peptides;
    std::vector<std::uint32_t> slot_handle;
    unique_peptides.reserve(candidates.size());
    slot_handle.reserve(candidates.size());
    for (std::size_t k = 0; k < candidates.size(); ++k)
    {
      if (!parsed_ok[k]) { continue; }
      unique_peptides.push_back(std::move(parsed[k]));
      slot_handle.push_back(candidates[k]);
    }
    parsed.clear(); parsed.shrink_to_fit();
    unique_step.fields().rows = static_cast<std::int64_t>(unique_peptides.size());
    unique_step.end();
    if (unique_peptides.empty()) { return library.precursorCount(); }

    trace::Scope sessions_step("sessions", 1);
    PeptDeepPredictor predictor(rt_model_path, prefer_gpu, 1, inferenceSessions(sessions), -1, gpu);
    sessions_step.fields().rows = static_cast<std::int64_t>(predictor.sessionCount());
    sessions_step.end();
    // State the device: a silent CPU fallback turns a 4-minute stage
    // into a 25-minute one with nothing in the log to explain it.
    writeDeviceLine("RT", predictor, gpu);
    std::vector<PeptDeepPredictor::Failure> failures;
    trace::Scope predict_step("predict", 1);
    predict_step.fields().rows = static_cast<std::int64_t>(unique_peptides.size());
    const auto predicted = predictor.predictRT(unique_peptides, &failures);
    predict_step.end();

    trace::Scope assign_step("assign", 1);
    std::vector<float> by_handle(library.strings().size(), std::nanf(""));
    parallel::forEach(predicted.size(), threads, [&](std::size_t slot)
    {
      float value = predicted[slot];
      // Applied here, in the model's own units, so that any downstream affine
      // rescale to an iRT gauge carries it through unchanged.
      if (free_cysteine_rt_correction && !std::isnan(value))
      {
        value += freeCysteineRtOffset(freeCysteines(unique_peptides[slot]));
      }
      by_handle[slot_handle[slot]] = value;
    }, 4096);

    std::size_t unpredicted = 0;
    auto& precursors = library.precursors();
    for (std::size_t i = 0; i < library.precursorCount(); ++i)
    {
      const float value = by_handle[precursors.modified_sequence[i]];
      precursors.irt[i] = value;
      if (!std::isfinite(value)) { ++unpredicted; }
    }
    // Appending nothing, but the accessor is non-const; the ordering is
    // unchanged, so restore the flag rather than forcing a needless re-sort.
    return unpredicted;
  }


  namespace
  {
    using RankedFragments = std::vector<std::pair<float, Fragment>>;

    /// The fragments a predicted spectrum keeps for one precursor, best
    /// first: every fragment in range, ranked by predicted intensity, capped
    /// at max_fragments (with the doubly-charged quota), emptied below
    /// min_fragments. Pure per precursor, so precursors can run in parallel.
    void selectFragments(const AASequence& peptide, int charge,
                         const PeptDeepPredictor::Spectrum& spectrum,
                         const DigestParams& params,
                         std::vector<Fragment>& fragments, RankedFragments& ranked)
    {
      // Re-enumerate without the m/z cap: the cap is what has to be redone.
      enumerateAllFragments(peptide, params, charge, fragments);

      const std::size_t residues = peptide.size();
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
      // whatever survived. That is where the IH1 library's 366,084 targets
      // with 0-2 fragments came from -- and `appendDecoys` then refused to
      // build decoys for them, because it does apply the bar. The result was
      // 7.3% of targets competing against a null that contained nothing like
      // them. Emptied rather than kept: a precursor below the bar its decoy
      // would have to clear is not searchable, and `-min_library_fragments`
      // drops it symmetrically at load.
      if (ranked.size() < params.min_fragments) { ranked.clear(); }
    }
  } // namespace

  std::unique_ptr<PeptDeepPredictor> LibraryGenerator::openMs2Predictor(
    const std::string& ms2_model_path, bool prefer_gpu, unsigned sessions, const GpuOptions& gpu)
  {
    trace::Scope sessions_step("ms2_sessions", 0);
    auto predictor = std::make_unique<PeptDeepPredictor>(ms2_model_path, prefer_gpu, 1, inferenceSessions(sessions), -1, gpu);
    sessions_step.fields().rows = static_cast<std::int64_t>(predictor->sessionCount());
    // State the device: a silent CPU fallback turns a 4-minute stage
    // into a 25-minute one with nothing in the log to explain it.
    writeDeviceLine("MS2", *predictor, gpu);
    return predictor;
  }

  bool LibraryGenerator::overlapsCpuPhases(PeptDeepPredictor::Provider ms2, bool lstm_on_cpu)
  {
    return ms2 == PeptDeepPredictor::Provider::CoreML || (ms2 == PeptDeepPredictor::Provider::CUDA && lstm_on_cpu);
  }

  bool LibraryGenerator::predictsIonMobility(const std::string& request, const std::string& canonical_instrument)
  {
    if (request == "true") { return true; }
    if (request == "false") { return false; }
    if (request != "auto") { throw std::invalid_argument("predict_ion_mobility must be auto, true or false"); }
    return PeptDeepEncoder::hasIonMobility(canonical_instrument);
  }

  std::size_t LibraryGenerator::predictFragmentIntensities(
    Library& library, const std::string& ms2_model_path, const DigestParams& params,
    float nce, const std::string& instrument, bool prefer_gpu, unsigned sessions,
    const GpuOptions& gpu, Ms2Sink* sink)
  {
    if (params.max_fragment_charge < 1 || params.max_fragment_charge > 2)
    { throw std::invalid_argument("max_fragment_charge must be between 1 and 2"); }
    if (library.precursorCount() == 0) { return 0; }
    const auto predictor = openMs2Predictor(ms2_model_path, prefer_gpu, sessions, gpu);
    return predictFragmentIntensities(library, *predictor, params, nce, instrument, sessions, sink);
  }

  std::size_t LibraryGenerator::predictFragmentIntensities(
    Library& library, PeptDeepPredictor& predictor, const DigestParams& params,
    float nce, const std::string& instrument, unsigned sessions, Ms2Sink* sink)
  {
    if (params.max_fragment_charge < 1 || params.max_fragment_charge > 2)
    { throw std::invalid_argument("max_fragment_charge must be between 1 and 2"); }
    auto& p = library.precursors();
    auto& t = library.transitions();
    const std::size_t n = library.precursorCount();
    if (n == 0) { return 0; }

    // The predictor's construction is the ms2_sessions phase before this one.
    trace::Scope ms2_phase("ms2", 0);
    trace::Scope reserve_step("reserve", 1);

    // The new transition arrays are built alongside the old ones and swapped in
    // at the end. Editing in place is not possible: a precursor's fragment
    // count changes, so every later precursor's CSR offset moves.
    Library::TransitionArrays built;
    Library::checkTransitionCapacity(t.product_mz.size());
    built.product_mz.reserve(t.product_mz.size());
    built.library_intensity.reserve(t.product_mz.size());
    built.type.reserve(t.product_mz.size());
    built.ordinal.reserve(t.product_mz.size());
    built.charge.reserve(t.product_mz.size());
    built.loss.reserve(t.product_mz.size());

    std::vector<std::uint32_t> begin(n, 0), count(n, 0);
    reserve_step.end();

    // Blocked so that peak memory is set by the block, not by the proteome.
    //
    // Pipelined: while the model predicts block k (on the GPU, or on the CPU
    // sessions), a helper thread parses block k+1 and then selects and
    // appends block k-1, so the sessions do not wait for host work between
    // two calls. Selection is parallel per precursor and the append is serial
    // in precursor order, one block after the other, so the arrays come out
    // exactly as from the serial loop, whatever the thread count.
    //
    // The block is also the unit every session waits at: 100,000 precursors
    // on the CPU, so that a call is long against its last batch; 20,000 on
    // CUDA, as measured there. The library does not depend on it.
    // ODIA_MS2_BLOCK overrides it, for tests (several blocks from a small
    // input) and measurements.
    const std::size_t BLOCK = [&]() -> std::size_t {
      if (const char* v = std::getenv("ODIA_MS2_BLOCK"); v != nullptr && *v != '\0')
      {
        char* end = nullptr;
        const long long x = std::strtoll(v, &end, 10);
        if (end != v && x > 0) { return static_cast<std::size_t>(x); }
      }
      // CoreML pads every batch up to a compiled shape (PeptDeepPredictor::
      // padded_), so a block twice the CPU's halves the padding.
      if (predictor.provider() == PeptDeepPredictor::Provider::CoreML) { return std::size_t{200000}; }
      return predictor.provider() == PeptDeepPredictor::Provider::CUDA ? 20000 : 100000;
    }();
    const unsigned threads = parallel::workers(sessions);
    std::size_t unpredicted = 0;

    struct Block
    {
      std::size_t base = 0, last = 0;
      std::vector<AASequence> peptides;
      std::vector<int> charges;
      std::vector<PeptDeepPredictor::Spectrum> spectra;
      std::vector<RankedFragments> ranked;
    };

    const auto parse = [&](std::size_t base, unsigned team)
    {
      trace::Scope parse_step("parse", 2);
      parse_step.fields().idx = static_cast<std::int64_t>(base / BLOCK);
      parse_step.fields().rows = static_cast<std::int64_t>(std::min(base + BLOCK, n) - base);
      Block block;
      block.base = base;
      block.last = std::min(base + BLOCK, n);
      const std::size_t size = block.last - base;
      block.peptides.resize(size);
      block.charges.resize(size);
      parallel::forEach(size, team, [&](std::size_t k)
      {
        block.peptides[k] = parsePeptide(library.strings().get(p.modified_sequence[base + k]));
        block.charges[k] = static_cast<int>(p.charge[base + k]);
      }, PARSE_GRAIN);
      return block;
    };

    const auto select = [&](Block& block, unsigned team)
    {
      trace::Scope post_step("post", 2);
      post_step.fields().idx = static_cast<std::int64_t>(block.base / BLOCK);
      post_step.fields().rows = static_cast<std::int64_t>(block.last - block.base);
      const std::size_t size = block.last - block.base;
      block.ranked.assign(size, RankedFragments{});
      constexpr std::size_t TASK = 256;
      parallel::forEach((size + TASK - 1) / TASK, team, [&](std::size_t task)
      {
        std::vector<Fragment> fragments;
        for (std::size_t k = task * TASK; k < std::min(size, (task + 1) * TASK); ++k)
        {
          if (block.spectra[k].positions == 0) { continue; }   // kept as generated, below
          selectFragments(block.peptides[k], block.charges[k], block.spectra[k], params,
                          fragments, block.ranked[k]);
        }
      });
    };

    const auto append = [&](const Block& block)
    {
      trace::Scope emit_step("emit", 2);
      emit_step.fields().idx = static_cast<std::int64_t>(block.base / BLOCK);
      emit_step.fields().rows = static_cast<std::int64_t>(block.last - block.base);
      for (std::size_t i = block.base; i < block.last; ++i)
      {
        const auto& spectrum = block.spectra[i - block.base];
        begin[i] = static_cast<std::uint32_t>(built.product_mz.size());

        if (spectrum.positions == 0)
        {
          Library::checkTransitionCapacity(built.product_mz.size(), p.transition_count[i]);
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

        const auto& ranked = block.ranked[i - block.base];
        Library::checkTransitionCapacity(built.product_mz.size(), ranked.size());
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
      if (sink == nullptr) { return; }
      // The block's kept precursors, handed on with copies of their final
      // transitions: `built` itself grows on this thread meanwhile.
      Library::PrecursorRows rows;
      std::size_t total = 0;
      for (std::size_t i = block.base; i < block.last; ++i) { total += count[i]; }
      rows.precursor.reserve(block.last - block.base);
      rows.begin.reserve(block.last - block.base);
      rows.count.reserve(block.last - block.base);
      auto& out = rows.transitions;
      out.product_mz.reserve(total); out.library_intensity.reserve(total); out.type.reserve(total);
      out.ordinal.reserve(total); out.charge.reserve(total); out.loss.reserve(total);
      for (std::size_t i = block.base; i < block.last; ++i)
      {
        if (count[i] == 0) { continue; }
        rows.precursor.push_back(static_cast<std::uint32_t>(i));
        rows.begin.push_back(static_cast<std::uint32_t>(out.product_mz.size()));
        rows.count.push_back(count[i]);
        const auto first = built.product_mz.begin() + begin[i];
        out.product_mz.insert(out.product_mz.end(), first, first + count[i]);
        const auto b = static_cast<std::ptrdiff_t>(begin[i]), e = b + static_cast<std::ptrdiff_t>(count[i]);
        out.library_intensity.insert(out.library_intensity.end(), built.library_intensity.begin() + b, built.library_intensity.begin() + e);
        out.type.insert(out.type.end(), built.type.begin() + b, built.type.begin() + e);
        out.ordinal.insert(out.ordinal.end(), built.ordinal.begin() + b, built.ordinal.begin() + e);
        out.charge.insert(out.charge.end(), built.charge.begin() + b, built.charge.begin() + e);
        out.loss.insert(out.loss.end(), built.loss.begin() + b, built.loss.begin() + e);
      }
      trace::Scope hand_step("hand_over", 2);
      hand_step.fields().rows = static_cast<std::int64_t>(rows.precursor.size());
      sink->rows(std::move(rows));
    };

    const unsigned helper = overlapThreads(predictor, threads);
    Block next = parse(0, threads);   // nothing to overlap with yet: every thread
    Block previous;                   // predicted, not yet selected and appended
    for (std::size_t base = 0; base < n; base += BLOCK)
    {
      Block current = std::move(next);
      next = Block{};
      // Declared after `current`, so that it is waited for before `current`
      // goes, should the prediction throw.
      Overlapped host = overlap([&]() {
        if (base + BLOCK < n) { next = parse(base + BLOCK, helper); }
        if (previous.last > previous.base)
        {
          select(previous, helper);
          append(previous);
          previous = Block{};   // freed here, not between two calls
        }
      });
      trace::Scope block_step("block", 1);
      block_step.fields().idx = static_cast<std::int64_t>(base / BLOCK);
      block_step.fields().rows = static_cast<std::int64_t>(current.last - current.base);
      std::vector<PeptDeepPredictor::Failure> failures;
      trace::Scope predict_step("predict", 2);
      predict_step.fields().rows = static_cast<std::int64_t>(current.last - current.base);
      current.spectra = predictor.predictMS2(current.peptides, current.charges, nce, instrument, &failures);
      predict_step.end();
      host.get();
      previous = std::move(current);
    }
    if (previous.last > previous.base)
    {
      // The sessions are done: every thread.
      trace::Scope drain_step("drain", 1);
      select(previous, threads);
      append(previous);
    }

    if (sink != nullptr)
    {
      // Everything handed over is consumed before the library is rewritten.
      trace::Scope drain_sink("drain_sink", 1);
      sink->drain();
    }
    trace::Scope finalize_step("finalize", 1);
    t = std::move(built);
    p.transition_begin = std::move(begin);
    p.transition_count = std::move(count);
    std::vector<std::size_t> retained;
    for (std::size_t i = 0; i < n; ++i)
    { if (p.transition_count[i] > 0) { retained.push_back(i); } }
    if (retained.size() != n) { library = library.subsetByIndex(retained); }
    library.shrinkToFit();
    return unpredicted;
  }



  LibraryGenerator::IrtCalibration
  LibraryGenerator::fitIrtCalibration(const std::string& rt_model_path,
                                      const std::string& standards_file,
                                      bool prefer_gpu, const GpuOptions& gpu)
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
      try
      {
        if (tab == std::string::npos) { throw std::runtime_error("missing tab"); }
        auto peptide = AASequence::fromString(line.substr(0, tab));
        const auto value = line.substr(tab + 1);
        std::size_t consumed = 0;
        const double irt = std::stod(value, &consumed);
        if (!std::isfinite(irt) || value.find_first_not_of(" \t\r", consumed) != std::string::npos)
        { throw std::runtime_error("expected a finite iRT value"); }
        peptides.push_back(std::move(peptide));
        known.push_back(irt);
      }
      catch (const std::exception& e)
      { throw std::runtime_error("invalid iRT standard in " + standards_file + ": " + line + " (" + e.what() + ")"); }
    }
    if (peptides.size() < 3)
    {
      throw std::runtime_error("iRT standards file has fewer than 3 usable peptides: " +
                               standards_file);
    }

    // Eleven peptides: one session is plenty, and each CUDA session costs a
    // cuDNN set-up.
    trace::Scope irt_phase("irt_fit", 0);
    GpuOptions one_session = gpu;
    one_session.sessions_per_device = 1;
    PeptDeepPredictor predictor(rt_model_path, prefer_gpu, 0, 1, -1, one_session);
    // State the device: a silent CPU fallback turns a 4-minute stage
    // into a 25-minute one with nothing in the log to explain it.
    writeDeviceLine("IRT-FIT", predictor, gpu);
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
      if (!std::isfinite(raw[i])) { continue; }   // a standard we could not encode
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
    if (!std::isfinite(c.slope) || !std::isfinite(c.intercept) || c.slope <= 0)
    { throw std::runtime_error("iRT calibration must have a finite positive slope"); }
    c.peptides = n;
    for (std::size_t i = 0; i < raw.size(); ++i)
    {
      if (!std::isfinite(raw[i])) { continue; }
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
    unsigned sessions, bool derive_mobility, const GpuOptions& gpu,
    std::vector<std::uint16_t>* lengths, bool* on_cpu)
  {
    trace::Scope ccs_phase("ccs", 0);
    auto& p = library.precursors();
    const std::size_t n = library.precursorCount();
    p.ccs.assign(n, std::numeric_limits<float>::quiet_NaN());
    if (derive_mobility) { p.im.assign(n, std::numeric_limits<float>::quiet_NaN()); }
    if (n == 0) { return 0; }

    trace::Scope sessions_step("sessions", 1);
    PeptDeepPredictor predictor(ccs_model_path, prefer_gpu, 1, inferenceSessions(sessions), -1, gpu);
    sessions_step.fields().rows = static_cast<std::int64_t>(predictor.sessionCount());
    sessions_step.end();
    // State the device: a silent CPU fallback turns a 4-minute stage
    // into a 25-minute one with nothing in the log to explain it.
    writeDeviceLine("CCS", predictor, gpu);
    if (on_cpu != nullptr) { *on_cpu = predictor.provider() == PeptDeepPredictor::Provider::CPU; }
    if (lengths != nullptr) { lengths->assign(n, 0); }

    // Blocked, as the MS2 pass is: the whole proteome at once would hold every
    // parsed AASequence live alongside the library. Pipelined like it too:
    // while block k predicts, a helper thread parses block k+1 and stores
    // block k-1's values. Every row is written by exactly one block, so the
    // result does not depend on the overlap.
    const std::size_t BLOCK = ccsBlock();
    const unsigned threads = parallel::workers(sessions);
    const unsigned helper = overlapThreads(predictor, threads);
    std::size_t unpredicted = 0;

    struct Block
    {
      std::size_t base = 0, last = 0;
      std::vector<AASequence> peptides;
      std::vector<int> charges;
      std::vector<float> ccs;
    };

    const auto parse = [&](std::size_t base, unsigned team)
    {
      Block block;
      block.base = base;
      block.last = std::min(base + BLOCK, n);
      trace::Scope parse_step("parse", 2);
      parse_step.fields().idx = static_cast<std::int64_t>(base / BLOCK);
      parse_step.fields().rows = static_cast<std::int64_t>(block.last - base);
      block.peptides.resize(block.last - base);
      block.charges.resize(block.last - base);
      parallel::forEach(block.last - base, team, [&](std::size_t k)
      {
        block.peptides[k] = parsePeptide(library.strings().get(p.modified_sequence[base + k]));
        block.charges[k] = static_cast<int>(p.charge[base + k]);
      }, PARSE_GRAIN);
      if (lengths != nullptr)
      {
        for (std::size_t k = 0; k < block.peptides.size(); ++k)
        { (*lengths)[base + k] = static_cast<std::uint16_t>(block.peptides[k].size()); }
      }
      return block;
    };

    const auto post = [&](const Block& block)
    {
      trace::Scope post_step("post", 2);
      post_step.fields().idx = static_cast<std::int64_t>(block.base / BLOCK);
      post_step.fields().rows = static_cast<std::int64_t>(block.last - block.base);
      for (std::size_t i = block.base; i < block.last; ++i)
      {
        p.ccs[i] = block.ccs[i - block.base];
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
    };

    Block next = parse(0, threads);   // nothing to overlap with yet: every thread
    Block previous;                   // predicted, not yet stored
    for (std::size_t base = 0; base < n; base += BLOCK)
    {
      Block current = std::move(next);
      next = Block{};
      // Declared after `current`: waited for before `current` goes.
      Overlapped host = overlap([&]() {
        if (base + BLOCK < n) { next = parse(base + BLOCK, helper); }
        if (previous.last > previous.base)
        {
          post(previous);
          previous = Block{};   // freed here, not between two calls
        }
      });
      trace::Scope block_step("block", 1);
      block_step.fields().idx = static_cast<std::int64_t>(base / BLOCK);
      block_step.fields().rows = static_cast<std::int64_t>(current.last - current.base);
      std::vector<PeptDeepPredictor::Failure> failures;
      trace::Scope predict_step("predict", 2);
      predict_step.fields().rows = static_cast<std::int64_t>(current.last - current.base);
      current.ccs = predictor.predictCCS(current.peptides, current.charges, &failures);
      predict_step.end();
      host.get();
      previous = std::move(current);
    }
    if (previous.last > previous.base) { post(previous); }
    return unpredicted;
  }

  namespace
  {
    /// Which of @p count rows -- the i-th of residue count length(i) -- a CPU
    /// CCS call of rows [0, count) runs as a batch of one: per residue count,
    /// the predictor cuts the rows into MAX_BATCH_ROWS chunks in order, and
    /// only a last chunk of exactly one row is a one-row batch (the CPU cut
    /// into CPU_CHUNK_ROWS pieces never leaves one; PeptDeepPredictor).
    template <class Length>
    void markAlone(std::size_t count, const Length& length, std::vector<std::uint8_t>& alone)
    {
      std::map<std::size_t, std::pair<std::size_t, std::size_t>> groups;   // length -> (rows, last row)
      for (std::size_t k = 0; k < count; ++k)
      {
        auto& g = groups[length(k)];
        ++g.first;
        g.second = k;
      }
      alone.assign(count, 0);
      for (const auto& [len, g] : groups)
      {
        if (g.first % PeptDeepPredictor::MAX_BATCH_ROWS == 1) { alone[g.second] = 1; }
      }
    }
  } // namespace

  std::size_t LibraryGenerator::ccsBlock()
  {
    static const std::size_t block = [] {
      if (const char* v = std::getenv("ODIA_CCS_BLOCK"); v != nullptr && *v != '\0')
      {
        char* end = nullptr;
        const long long x = std::strtoll(v, &end, 10);
        if (end != v && x > 0) { return static_cast<std::size_t>(x); }
      }
      return std::size_t{200000};
    }();
    return block;
  }

  LibraryGenerator::CcsBatchGuard::CcsBatchGuard(std::string ccs_model_path, bool derive_mobility,
                                                 std::vector<std::uint16_t> lengths, bool on_cpu)
    : model_(std::move(ccs_model_path)), derive_mobility_(derive_mobility), lengths_(std::move(lengths)),
      on_cpu_(on_cpu)
  {
    if (!on_cpu_) { return; }
    alone_.assign(lengths_.size(), 0);
    std::vector<std::uint8_t> part;
    for (std::size_t base = 0; base < lengths_.size(); base += ccsBlock())
    {
      const std::size_t count = std::min(ccsBlock(), lengths_.size() - base);
      markAlone(count, [&](std::size_t k) { return lengths_[base + k]; }, part);
      std::copy(part.begin(), part.end(), alone_.begin() + static_cast<std::ptrdiff_t>(base));
    }
  }

  LibraryGenerator::CcsBatchGuard::~CcsBatchGuard() = default;

  std::size_t LibraryGenerator::CcsBatchGuard::settle(Library& library, const std::vector<std::uint32_t>& kept,
                                                      std::size_t from, std::size_t to)
  {
    if (!on_cpu_ || from >= to) { return 0; }
    trace::Scope settle_step("ccs_settle", 2);
    settle_step.fields().rows = static_cast<std::int64_t>(to - from);
    auto& p = library.precursors();
    std::vector<std::uint8_t> old_alone;
    std::size_t redone = 0;
    for (std::size_t base = from; base < to; base += ccsBlock())
    {
      const std::size_t count = std::min(ccsBlock(), to - base);
      markAlone(count, [&](std::size_t k) { return lengths_[kept[base + k]]; }, old_alone);
      for (std::size_t k = 0; k < count; ++k)
      {
        const std::uint32_t i = kept[base + k];
        if (old_alone[k] == alone_[i]) { continue; }
        // Predicted again the old way: alone, or as a batch of two rows (any
        // batch of two or more gives a row the same value).
        if (!predictor_) { predictor_ = std::make_unique<PeptDeepPredictor>(model_, false, 1, 1); }
        const auto peptide = parsePeptide(library.strings().get(p.modified_sequence[i]));
        const int charge = static_cast<int>(p.charge[i]);
        const std::size_t rows = old_alone[k] ? 1 : 2;
        const auto ccs = predictor_->predictCCS(std::vector<AASequence>(rows, peptide), std::vector<int>(rows, charge));
        p.ccs[i] = ccs.front();
        if (derive_mobility_ && !std::isnan(p.ccs[i]))
        { p.im[i] = static_cast<float>(mobilityFromCCS(p.ccs[i], fromFixed(p.mz[i]), p.charge[i])); }
        ++redone;
      }
    }
    settle_step.fields().idx = static_cast<std::int64_t>(redone);
    return redone;
  }

  std::size_t LibraryGenerator::appendDecoys(Library& library, DecoyMethod method,
                                             std::size_t* skipped_out,
                                             std::size_t min_fragments,
                                             bool recompute_decoy_mz)
  {
    if (method == DecoyMethod::None) { return 0; }
    trace::Scope decoy_phase("decoys", 0);

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

      // Parse with OpenMS so terminal and nested modification names cannot
      // become bogus residue tokens. Terminal modifications stay terminal.
      AASequence target;
      try { target = parsePeptide(sequence); }
      catch (const std::exception&) { ++skipped; continue; }
      struct Token { std::string text; char residue; bool modified; };
      std::vector<Token> tokens;
      for (Size k = 0; k < target.size(); ++k)
      {
        const auto& residue = target[k];
        auto single = AASequence::fromString(residue.getOneLetterCode());
        if (residue.isModified()) { single.setModification(0, residue.getModification()->getId()); }
        tokens.push_back({single.toString(), residue.getOneLetterCode()[0], residue.isModified()});
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

        std::uint64_t seed = 14695981039346656037ull;      // FNV-1a of the sequence
        for (const char ch : sequence) { seed = (seed ^ static_cast<std::uint8_t>(ch)) * 1099511628211ull; }
        std::mt19937_64 rng(seed);

        std::vector<Token> best;
        for (int attempt = 0; attempt < 20; ++attempt)
        {
          // Fisher-Yates with rejection sampling: std::shuffle and its
          // distribution are not specified across standard libraries.
          for (std::size_t k = mid.size(); k > 1; --k)
          {
            const std::uint64_t bound = k;
            const std::uint64_t threshold = (std::uint64_t{0} - bound) % bound;
            std::uint64_t draw;
            do { draw = rng(); } while (draw < threshold);
            std::swap(mid[k - 1], mid[draw % bound]);
          }
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
      try
      {
        decoy = AASequence::fromString(decoy_sequence);
        if (target.hasNTerminalModification())
        { decoy.setNTerminalModification(target.getNTerminalModification()->getId()); }
        if (target.hasCTerminalModification())
        { decoy.setCTerminalModification(target.getCTerminalModification()->getId()); }
      }
      catch (const std::exception&) { ++skipped; continue; }
      if (decoy == target) { ++skipped; continue; }

      auto& p = library.precursors();
      auto& t = library.transitions();

      const std::uint32_t begin = p.transition_begin[i];
      const std::uint32_t count = p.transition_count[i];
      Library::checkTransitionCapacity(t.product_mz.size());
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

        Library::checkTransitionCapacity(t.product_mz.size(), 1);
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
    // v5 removes duplicate variable-modification forms and empty assays,
    // rejects ambiguous masses, and pins shuffle across standard libraries.
    // v6: residue-less peptide-terminal variable modifications are applied at
    // max_variable_modifications 1 (they were silently dropped before).
    ps << "v6"
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
       // reproduce the CAM-inclusive result.
       << ";fixmod=" << join(p.fixed_modifications)
       << ";nme=" << (p.n_terminal_methionine_excision ? 1 : 0)
       // Changes the stored RT of every free-cysteine peptide, so a library
       // built with it is not the same library.
       << ";fcys=" << (p.free_cysteine_rt_correction ? 1 : 0)
       // "none": no CCS and no 1/K0 at all (an instrument without ion
       // mobility). A library with them keeps 0.11.5's token.
       << ";im=" << (!p.predict_ion_mobility ? "none" : p.derive_ion_mobility ? "1" : "0")
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
