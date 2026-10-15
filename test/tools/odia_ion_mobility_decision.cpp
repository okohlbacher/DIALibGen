// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

// Whether generate predicts CCS and 1/K0, per instrument name and
// -generation:predict_ion_mobility: auto only for an instrument with ion
// mobility (timsTOF), true and false whatever the instrument.

#include <odia/LibraryGenerator.h>
#include <odia/PeptDeepEncoder.h>

#include <iostream>
#include <string>

int main()
{
  for (const char* request : {"auto", "true", "false"})
  {
    for (const char* name : {"QE", "Lumos", "Astral", "Exploris", "timsTOF", "timsTOF Pro", "SciexTOF", "ThermoTOF"})
    {
      const std::string canonical = ODIA::PeptDeepEncoder::canonicalInstrument(name);
      std::cout << request << " " << name << " -> "
                << (ODIA::LibraryGenerator::predictsIonMobility(request, canonical) ? "CCS" : "none") << "\n";
    }
  }
  try { (void)ODIA::LibraryGenerator::predictsIonMobility("yes", "timsTOF"); std::cout << "yes accepted\n"; }
  catch (const std::invalid_argument&) { std::cout << "yes refused\n"; }
  return 0;
}
