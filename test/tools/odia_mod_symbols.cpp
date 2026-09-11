// Copyright (c) 2026, Oliver Kohlbacher and the DIALibraryGenerator authors.
// SPDX-License-Identifier: BSD-3-Clause
//
// Walks every modification ModificationsDB knows and reports which element
// symbols fall outside AlphaPeptDeep's 109-element list after translation.
//
// This is the only check that can cover the isotope path. The Python reference
// cannot: its parser splits on the first ')', so a name like
// "Label:13C(6)15N(2)" is mangled before it is ever encoded -- meaning the
// rename table could be deleted entirely and the encoder comparison would not
// notice. ModificationsDB also loads PSI-MOD and XLMOD, not only UniMod, so the
// symbol set is wider than unimod.xml alone suggests.

#include <odia/PeptDeepEncoder.h>

#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/Element.h>

#include <iostream>
#include <map>
#include <set>
#include <string>

int main()
{
  auto* db = OpenMS::ModificationsDB::getInstance();
  std::set<std::string> unmapped;
  std::size_t total = 0, empty_formula = 0;

  for (OpenMS::Size i = 0; i < db->getNumberOfModifications(); ++i)
  {
    const auto* mod = db->getModification(i);
    if (mod == nullptr) { continue; }
    ++total;
    const auto& formula = mod->getDiffFormula();
    if (formula.isEmpty()) { ++empty_formula; continue; }
    for (const auto& [element, count] : formula)
    {
      (void)count;
      // Mirrors the encoder's translation, so a symbol reported here is one the
      // encoder would silently drop into the '?' bucket.
      const auto symbol = element->getSymbol();
      static const std::map<std::string, std::string> renames{
        {"(2)H", "2H"}, {"(13)C", "13C"}, {"(15)N", "15N"}, {"(18)O", "18O"},
        {"(1)H", "H"}, {"(12)C", "C"}, {"(14)N", "N"}, {"(16)O", "O"},
      };
      const auto it = renames.find(symbol);
      const auto translated = it == renames.end() ? symbol : it->second;
      if (!ODIA::peptDeepElementIndex(translated)) { unmapped.insert(symbol); }
    }
  }

  std::cout << "modifications: " << total
            << ", without a composition: " << empty_formula << "\n";
  std::cout << "unmapped symbols:";
  for (const auto& s : unmapped) { std::cout << " " << s; }
  std::cout << "\n";
  return 0;
}
