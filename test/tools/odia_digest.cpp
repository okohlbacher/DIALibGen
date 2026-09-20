// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <odia/LibraryGenerator.h>
#include <odia/DIANNLibraryFile.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <cmath>

#include <fstream>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>

int main(int argc, char** argv)
{
  if (argc != 2) { return 2; }
  std::ofstream(argv[1]) << ">full\nMPEPTIDEKPEPTIDERACDEFGK\n"
                         ">short\nMAKPEPTIDERACDEFGK\n"
                         ">internal\nPEPTIDERACDEFGK\n";
  try
  {
    ODIA::DigestParams params;
    params.charges = {2};
    params.fixed_modifications.clear();
    params.precursor_mz_min = params.fragment_mz_min = 0.0;
    params.precursor_mz_max = params.fragment_mz_max = 10000.0;
    params.min_fragments = 0;

    for (const bool excise : {false, true})
    {
      params.n_terminal_methionine_excision = excise;
      for (std::size_t missed = 0; missed <= 2; ++missed)
      {
        params.missed_cleavages = missed;
        std::set<std::string> expected{"MPEPTIDEK", "PEPTIDER", "ACDEFGK"};
        if (excise) { expected.insert("PEPTIDEK"); }
        if (missed >= 1)
        {
          expected.insert({"MPEPTIDEKPEPTIDER", "PEPTIDERACDEFGK", "MAKPEPTIDER"});
          if (excise) { expected.insert({"PEPTIDEKPEPTIDER", "AKPEPTIDER"}); }
        }
        if (missed >= 2)
        {
          expected.insert({"MPEPTIDEKPEPTIDERACDEFGK", "MAKPEPTIDERACDEFGK"});
          if (excise)
          {
            expected.insert({"PEPTIDEKPEPTIDERACDEFGK", "AKPEPTIDERACDEFGK"});
          }
        }

        ODIA::Library library;
        const auto stats = ODIA::LibraryGenerator::generate(argv[1], params, library);
        std::set<std::string> actual;
        const auto& precursors = library.precursors();
        for (std::size_t i = 0; i < library.precursorCount(); ++i)
        {
          const std::string sequence(library.strings().get(precursors.modified_sequence[i]));
          actual.insert(sequence);
          if (sequence == "PEPTIDER" &&
              library.strings().get(precursors.protein_group[i]) != "full;short;internal")
          {
            throw std::runtime_error("internal peptide lost or duplicated a protein accession");
          }
        }
        if (actual != expected || stats.peptides != expected.size() ||
            library.precursorCount() != expected.size())
        {
          std::cerr << "Met excision=" << excise << ", missed cleavages=" << missed << '\n';
          for (const auto& sequence : expected)
          {
            if (!actual.count(sequence)) { std::cerr << "missing: " << sequence << '\n'; }
          }
          throw std::runtime_error("digest does not match the expected unique peptides");
        }
      }
    }
    // OpenMS appends the base form itself when keep_unmodified=true.
    std::ofstream(argv[1]) << ">oxidation\nMPEPTIDEK\n";
    params.n_terminal_methionine_excision = false;
    params.missed_cleavages = 0;
    params.variable_modifications = {"Oxidation (M)"};
    for (const std::size_t maximum : {0u, 1u, 2u})
    {
      params.max_variable_modifications = maximum;
      ODIA::Library library;
      ODIA::LibraryGenerator::generate(argv[1], params, library);
      if (library.precursorCount() != (maximum == 0 ? 1u : 2u))
      { throw std::runtime_error("variable modifications duplicated or lost the base form"); }
      const std::string tsv = std::string(argv[1]) + ".tsv";
      ODIA::DIANNLibraryFile::storeTSV(tsv, library);
      ODIA::Library roundtrip;
      ODIA::DIANNLibraryFile::loadTSV(tsv, roundtrip);
      if (roundtrip.precursorCount() != library.precursorCount() ||
          roundtrip.transitionCount() != library.transitionCount())
      { throw std::runtime_error("variable modification TSV roundtrip changed assays"); }
    }
    params.variable_modifications.clear();
    std::ofstream(argv[1]) << ">ambiguous\nPEPTIDEXKPEPTIDER\n"
                             ">terminal_stop\nACDEFGK*\n"
                             ">selenium\nACDEUGK\n";
    ODIA::Library unusual;
    const auto unusual_stats = ODIA::LibraryGenerator::generate(argv[1], params, unusual);
    std::set<std::string> expected_unusual{"PEPTIDER", "ACDEFGK", "ACDEUGK"}, actual_unusual;
    for (const auto handle : unusual.precursors().modified_sequence)
    { actual_unusual.emplace(unusual.strings().get(handle)); }
    if (actual_unusual != expected_unusual || unusual_stats.dropped_ambiguous_peptides != 1)
    { throw std::runtime_error("ambiguous peptide affected unrelated valid peptides or terminal stop handling"); }
    for (const auto method : {ODIA::DecoyMethod::Reverse, ODIA::DecoyMethod::PseudoReverse})
    {
      std::ofstream(argv[1]) << ">palindrome\nAPEPEPA\n";
      ODIA::Library palindrome;
      ODIA::LibraryGenerator::generate(argv[1], params, palindrome);
      std::size_t skipped = 0;
      if (ODIA::LibraryGenerator::appendDecoys(palindrome, method, &skipped) != 0 || skipped != 1)
      { throw std::runtime_error("palindromic target was emitted as its own decoy"); }
    }
    std::ofstream(argv[1]) << ">shuffle\nACDEFGHK\n";
    ODIA::Library shuffled;
    ODIA::LibraryGenerator::generate(argv[1], params, shuffled);
    if (ODIA::LibraryGenerator::appendDecoys(shuffled, ODIA::DecoyMethod::Shuffle) != 1)
    { throw std::runtime_error("shuffle produced no decoy"); }
    const auto expected_shuffle = OpenMS::AASequence::fromString("AGHDCFEK");
    const auto& sp = shuffled.precursors();
    const auto& st = shuffled.transitions();
    for (std::size_t j = sp.transition_begin[1]; j < sp.transition_begin[1] + sp.transition_count[1]; ++j)
    {
      const auto ion = st.type[j] == ODIA::FragmentType::B ? OpenMS::Residue::BIon : OpenMS::Residue::YIon;
      const auto part = st.type[j] == ODIA::FragmentType::B ? expected_shuffle.getPrefix(st.ordinal[j]) : expected_shuffle.getSuffix(st.ordinal[j]);
      if (st.product_mz[j] != ODIA::toFixed(part.getMZ(st.charge[j], ion)))
      { throw std::runtime_error("shuffle differs from the pinned cross-platform decoy AGHDCFEK"); }
    }
    std::vector<ODIA::MzFixed> terminal_masses;
    for (const auto* sequence : {".(Acetyl)ACDEFGHK", "(UniMod:1)ACDEFGHK"})
    {
      ODIA::Library terminal;
      ODIA::LibraryGenerator::generate(argv[1], params, terminal);
      terminal.precursors().modified_sequence[0] = terminal.strings().intern(sequence);
      if (ODIA::LibraryGenerator::appendDecoys(terminal, ODIA::DecoyMethod::PseudoReverse) != 1)
      { throw std::runtime_error("terminal modification prevented decoy generation"); }
      const auto begin = terminal.precursors().transition_begin[1];
      std::vector<ODIA::MzFixed> masses(terminal.transitions().product_mz.begin() + begin, terminal.transitions().product_mz.end());
      if (!terminal_masses.empty() && masses != terminal_masses)
      { throw std::runtime_error("terminal modification spelling changed decoy masses"); }
      terminal_masses = std::move(masses);
    }
    std::cout << "Met excision includes all missed-cleavage peptides and preserves protein mapping\n";
  }
  catch (const std::exception& e)
  {
    std::cerr << e.what() << '\n';
    return 1;
  }
  return 0;
}
