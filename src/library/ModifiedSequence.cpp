// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <odia/ModifiedSequence.h>

#include <OpenMS/CHEMISTRY/ModificationsDB.h>

#include <cctype>
#include <set>
#include <string>

namespace ODIA
{
  namespace
  {
    /// UniMod accession <-> interchange name. Canonical form is the accession:
    /// it is unambiguous, whereas the names differ between tools.
    struct ModAlias { const char* name; const char* unimod; };
    constexpr ModAlias MOD_ALIASES[] = {
      {"Acetyl", "UniMod:1"},
      {"Carbamidomethyl", "UniMod:4"},
      {"Carbamyl", "UniMod:5"},
      {"Deamidated", "UniMod:7"},
      {"Phospho", "UniMod:21"},
      {"Pyro-carbamidomethyl", "UniMod:26"},
      {"Glu->pyro-Glu", "UniMod:27"},
      {"Gln->pyro-Glu", "UniMod:28"},
      {"Methyl", "UniMod:34"},
      {"Oxidation", "UniMod:35"},
      {"Dimethyl", "UniMod:36"},
      {"Trimethyl", "UniMod:37"},
      {"GG", "UniMod:121"},
      {"Label:13C(6)15N(2)", "UniMod:259"},
      {"Label:13C(6)15N(4)", "UniMod:267"},
    };

    std::size_t tokenEnd(std::string_view seq, std::size_t start)
    {
      const char open = seq[start], close = open == '(' ? ')' : ']';
      std::size_t depth = 0;
      for (std::size_t i = start; i < seq.size(); ++i)
      {
        if (seq[i] == open) { ++depth; }
        else if (seq[i] == close && --depth == 0) { return i; }
      }
      return std::string_view::npos;
    }

    char firstResidue(std::string_view seq)
    {
      for (std::size_t i = 0; i < seq.size(); ++i)
      {
        if (seq[i] == '(' || seq[i] == '[')
        {
          i = tokenEnd(seq, i);
          if (i == std::string_view::npos) { return 0; }
        }
        else if (std::isupper(static_cast<unsigned char>(seq[i]))) { return seq[i]; }
      }
      return 0;
    }

    std::string namedAccession(std::string_view token, char residue, bool n_term, bool c_term)
    {
      using Mod = OpenMS::ResidueModification;
      const OpenMS::String name{std::string(token)};
      const auto* database = OpenMS::ModificationsDB::getInstance();
      // Avoid the database's per-token warning for unknown names; the caller counts them.
      if (!database->has(name)) { return {}; }
      std::set<const Mod*> matches;
      database->searchModifications(matches, name, residue ? std::string(1, residue) : std::string{});
      std::string accession;
      for (const auto* mod : matches)
      {
        const auto term = mod->getTermSpecificity();
        const bool same_site = n_term ? (term == Mod::N_TERM || term == Mod::PROTEIN_N_TERM)
                             : c_term ? (term == Mod::C_TERM || term == Mod::PROTEIN_C_TERM)
                                      : term == Mod::ANYWHERE;
        if (!same_site) { continue; }
        if (mod->getUniModRecordId() <= 0) { return {}; }
        const std::string candidate = mod->getUniModAccession();
        if (!accession.empty() && accession != candidate) { return {}; } // ambiguous synonym
        accession = candidate;
      }
      return accession;
    }

  }

  std::string canonicalModifiedSequence(std::string_view seq, std::size_t* unknown)
  {
    // OpenMS marks N-terminal modifications with a leading dot; DIA-NN does not.
    if (seq.size() > 1 && seq[0] == '.' && (seq[1] == '(' || seq[1] == '[')) { seq.remove_prefix(1); }
    std::string out;
    out.reserve(seq.size());
    char residue = 0;
    bool c_term = false;
    for (std::size_t i = 0; i < seq.size();)
    {
      const char c = seq[i];
      if (c != '(' && c != '[')
      {
        out.push_back(c);
        if (std::isupper(static_cast<unsigned char>(c))) { residue = c; c_term = false; }
        else if (c == '.' && residue) { c_term = true; }
        ++i;
        continue;
      }

      // Isotope labels contain nested brackets, which belong to the same token.
      const std::size_t end = tokenEnd(seq, i);
      if (end == std::string_view::npos)
      {
        if (unknown) { ++*unknown; }
        out.append(seq.substr(i));
        break;
      }

      const std::string_view body = seq.substr(i + 1, end - i - 1);
      std::string canonical(body);
      bool known = false;
      for (const auto& a : MOD_ALIASES)
      { if (body == a.name) { canonical = a.unimod; known = true; break; } }
      if (!known && body.size() > 7 && body.substr(0, 7) == "UniMod:") { known = true; }
      if (!known)
      {
        const auto accession = namedAccession(body, residue ? residue : firstResidue(seq.substr(end + 1)),
                                             residue == 0, c_term);
        if (!accession.empty()) { canonical = accession; known = true; }
      }
      // Anything else -- a bare mass shift like (+57.0215), a bare number, an
      // unknown name -- passes through verbatim so it can never silently
      // collide with a known modification, and is counted.
      if (!known && unknown) { ++*unknown; }
      out.push_back('(');
      out.append(canonical);
      out.push_back(')');
      i = end + 1;
    }
    return out;
  }

} // namespace ODIA
