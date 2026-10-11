// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

// Terminal modifications on a FASTA whose protein termini are known: every
// modified form is counted, by sequence, against the forms written out here.
//
//   one    MSAAAGK DDDDDR FFFFFH   initiator Met: MSAAAGK and, excised, SAAAGK
//   two    GGGGGK DDDDDR HHHHHK    no Met; DDDDDR also in one (internal both times)
//   three  WWWWWK SAAAGK YYYYYR    SAAAGK internal here, protein N-terminal in one
//
// Trypsin/P, no missed cleavages, length 5..30, charge 2, windows open: the
// eight peptides are MSAAAGK SAAAGK DDDDDR FFFFFH GGGGGK HHHHHK WWWWWK YYYYYR.

#include <odia/LibraryGenerator.h>
#include <odia/ModifiedSequence.h>

#include <fstream>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
  std::multiset<std::string> forms(const std::string& fasta, const ODIA::DigestParams& params,
                                   ODIA::LibraryGenerator::Stats* stats = nullptr)
  {
    ODIA::Library library;
    const auto st = ODIA::LibraryGenerator::generate(fasta, params, library);
    if (stats != nullptr) { *stats = st; }
    std::multiset<std::string> out;
    for (const auto handle : library.precursors().modified_sequence)
    { out.emplace(library.strings().get(handle)); }
    return out;
  }

  std::size_t count(const std::multiset<std::string>& all, const std::string& token)
  {
    std::size_t n = 0;
    for (const auto& s : all) { n += s.find(token) != std::string::npos; }
    return n;
  }

  void expect(const std::multiset<std::string>& actual, const std::multiset<std::string>& expected,
              const std::string& what)
  {
    if (actual == expected) { std::cout << what << ": " << actual.size() << " forms\n"; return; }
    std::cerr << what << ":\n  expected:";
    for (const auto& s : expected) { std::cerr << ' ' << s; }
    std::cerr << "\n  actual:  ";
    for (const auto& s : actual) { std::cerr << ' ' << s; }
    std::cerr << '\n';
    throw std::runtime_error(what + ": wrong modified forms");
  }

  bool refused(const std::string& fasta, const ODIA::DigestParams& params)
  {
    try { forms(fasta, params); }
    catch (const std::invalid_argument& e)
    { return std::string(e.what()).find("produced no modified precursor") != std::string::npos; }
    return false;
  }
}

int main(int argc, char** argv)
{
  if (argc != 2) { return 2; }
  const std::string fasta = argv[1];
  std::ofstream(fasta) << ">one\nMSAAAGKDDDDDRFFFFFH\n"
                          ">two\nGGGGGKDDDDDRHHHHHK\n"
                          ">three\nWWWWWKSAAAGKYYYYYR\n";
  try
  {
    ODIA::DigestParams params;
    params.charges = {2};
    params.fixed_modifications.clear();
    params.missed_cleavages = 0;
    params.min_length = 5;
    params.precursor_mz_min = params.fragment_mz_min = 0.0;
    params.precursor_mz_max = params.fragment_mz_max = 10000.0;
    params.min_fragments = 0;

    const std::multiset<std::string> plain{"MSAAAGK", "SAAAGK", "DDDDDR", "FFFFFH",
                                           "GGGGGK", "HHHHHK", "WWWWWK", "YYYYYR"};
    expect(forms(fasta, params), plain, "unmodified");

    // Variable protein N-terminal acetylation: the four protein N-termini,
    // SAAAGK through the excised Met of 'one'; at the default maximum of 1.
    auto expected = plain;
    for (const char* s : {".(Acetyl)MSAAAGK", ".(Acetyl)SAAAGK", ".(Acetyl)GGGGGK", ".(Acetyl)WWWWWK"})
    { expected.insert(s); }
    ODIA::DigestParams p = params;
    p.variable_modifications = {"Acetyl (Protein N-term)"};
    ODIA::LibraryGenerator::Stats stats;
    const auto protein_n = forms(fasta, p, &stats);
    expect(protein_n, expected, "variable Acetyl (Protein N-term)");
    if (stats.variable_modification_precursors.front().second != 4)
    { throw std::runtime_error("variable Acetyl (Protein N-term) not counted"); }

    // Without excision, SAAAGK is internal everywhere.
    p.n_terminal_methionine_excision = false;
    expected = plain;   // SAAAGK still digested from 'three'
    for (const char* s : {".(Acetyl)MSAAAGK", ".(Acetyl)GGGGGK", ".(Acetyl)WWWWWK"}) { expected.insert(s); }
    expect(forms(fasta, p), expected, "variable Acetyl (Protein N-term), no Met excision");

    // Fixed: a peptide only ever at a protein N-terminus is only modified;
    // SAAAGK, N-terminal in 'one' and internal in 'three', is both.
    p = params;
    p.fixed_modifications = {"Acetyl (Protein N-term)"};
    expect(forms(fasta, p),
           {".(Acetyl)MSAAAGK", ".(Acetyl)SAAAGK", "SAAAGK", "DDDDDR", "FFFFFH",
            ".(Acetyl)GGGGGK", "HHHHHK", ".(Acetyl)WWWWWK", "YYYYYR"},
           "fixed Acetyl (Protein N-term)");

    // Residue-less peptide-terminal modifications at the default maximum 1.
    p = params;
    p.variable_modifications = {"Acetyl (N-term)"};
    const auto peptide_n = forms(fasta, p);
    if (peptide_n.size() != 16 || count(peptide_n, ".(Acetyl)") != 8)
    { throw std::runtime_error("Acetyl (N-term) is not on every peptide N-terminus"); }
    std::cout << "variable Acetyl (N-term): " << peptide_n.size() << " forms\n";
    p.variable_modifications = {"Amidated (C-term)"};
    const auto peptide_c = forms(fasta, p);
    if (peptide_c.size() != 16 || count(peptide_c, ".(Amidated)") != 8)
    { throw std::runtime_error("Amidated (C-term) is not on every peptide C-terminus"); }
    std::cout << "variable Amidated (C-term): " << peptide_c.size() << " forms\n";

    // Protein C-terminal: FFFFFH, HHHHHK, YYYYYR.
    p.variable_modifications = {"Amidated (Protein C-term)"};
    expected = plain;
    for (const char* s : {"FFFFFH.(Amidated)", "HHHHHK.(Amidated)", "YYYYYR.(Amidated)"}) { expected.insert(s); }
    expect(forms(fasta, p), expected, "variable Amidated (Protein C-term)");

    // Both acetylations write the same form once, and each is credited.
    p.variable_modifications = {"Acetyl (N-term)", "Acetyl (Protein N-term)"};
    const auto both = forms(fasta, p, &stats);
    if (both.size() != 16 || stats.variable_modification_precursors[0].second != 8 ||
        stats.variable_modification_precursors[1].second != 4)
    { throw std::runtime_error("peptide and protein N-terminal acetylation together: duplicated or uncredited"); }
    std::cout << "Acetyl (N-term) + Acetyl (Protein N-term): " << both.size() << " forms\n";

    // With a residue modification: one variable modification per form at
    // maximum 1, both together at 2.
    p.variable_modifications = {"Oxidation (M)", "Acetyl (Protein N-term)"};
    for (const std::size_t maximum : {1u, 2u})
    {
      p.max_variable_modifications = maximum;
      const auto mixed = forms(fasta, p);
      std::multiset<std::string> methionine;
      for (const auto& s : mixed) { if (s.find('M') != std::string::npos) { methionine.insert(s); } }
      std::multiset<std::string> want{"MSAAAGK", "M(Oxidation)SAAAGK", ".(Acetyl)MSAAAGK"};
      if (maximum == 2) { want.insert(".(Acetyl)M(Oxidation)SAAAGK"); }
      expect(methionine, want, "Oxidation (M) + Acetyl (Protein N-term), max " + std::to_string(maximum));
    }

    // DIA-NN's notation for the TSV and Parquet outputs.
    for (const auto& [stored, written] : std::vector<std::pair<std::string, std::string>>{
           {".(Acetyl)MSAAAGK", "(UniMod:1)MSAAAGK"},
           {".(Acetyl)M(Oxidation)SAAAGK", "(UniMod:1)M(UniMod:35)SAAAGK"},
           {"FFFFFH.(Amidated)", "FFFFFH(UniMod:2)"}})
    {
      const auto canonical = ODIA::diannModifiedSequence(stored);
      if (canonical != written)
      { throw std::runtime_error(stored + " is written " + canonical + ", not " + written); }
    }

    // A requested modification that still matches nothing is refused: no
    // protein N-terminal peptide passes the length window here.
    std::ofstream(fasta) << ">short_n_terminus\nMKAAAAAKGGGGGR\n";
    p = params;
    p.variable_modifications = {"Acetyl (Protein N-term)"};
    if (!refused(fasta, p)) { throw std::runtime_error("variable protein N-terminal modification matching nothing was not refused"); }
    p.variable_modifications.clear();
    p.fixed_modifications = {"Acetyl (Protein N-term)"};
    if (!refused(fasta, p)) { throw std::runtime_error("fixed protein N-terminal modification matching nothing was not refused"); }
    std::cout << "protein-terminal modifications matching nothing are refused\n";
  }
  catch (const std::exception& e)
  {
    std::cerr << e.what() << '\n';
    return 1;
  }
  return 0;
}
