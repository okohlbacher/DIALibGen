// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <odia/LibraryGenerator.h>

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
    std::cout << "Met excision includes all missed-cleavage peptides and preserves protein mapping\n";
  }
  catch (const std::exception& e)
  {
    std::cerr << e.what() << '\n';
    return 1;
  }
  return 0;
}
