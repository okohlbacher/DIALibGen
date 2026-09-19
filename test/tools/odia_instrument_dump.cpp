// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

// Instrument name -> canonical name -> model slot, and the NCE that name
// defaults to. Printed rather than asserted here so the test that reads it can
// pin the whole table: the failure this guards against is an alias quietly
// resolving to the UNTRAINED slot, which nothing downstream can detect.

#include <odia/PeptDeepEncoder.h>

#include <iostream>
#include <string>
#include <vector>

int main(int argc, const char** argv)
{
  std::vector<std::string> names;
  for (int i = 1; i < argc; ++i) { names.emplace_back(argv[i]); }
  if (names.empty())
  { names = {"QE", "Lumos", "timsTOF", "SciexTOF", "ThermoTOF", "Astral", "timsTOF Pro", "Exploris", "nonsense"}; }

  for (const std::string& n : names)
  {
    const std::string c = ODIA::PeptDeepEncoder::canonicalInstrument(n);
    std::cout << n << " -> ";
    if (c.empty()) { std::cout << "REFUSED\n"; continue; }
    std::cout << c << " slot " << ODIA::PeptDeepEncoder::instrumentIndex(c)
              << " nce " << ODIA::PeptDeepEncoder::defaultNce(c) << "\n";
  }
  return 0;
}
