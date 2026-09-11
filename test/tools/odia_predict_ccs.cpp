// Copyright (c) 2026, Oliver Kohlbacher and the DIALibraryGenerator authors.
// SPDX-License-Identifier: BSD-3-Clause
//
// Dumps predicted collision cross-sections so they can be compared against
// test/peptdeep_reference.py. Takes several peptides in one call, because the
// length grouping and the row write-back are only exercised by a batch.

#include <odia/PeptDeepPredictor.h>

#include <OpenMS/CHEMISTRY/AASequence.h>

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

int main(int argc, char** argv)
{
  if (argc < 3)
  {
    std::cerr << "usage: odia_predict_ccs <model.onnx> <sequence:charge> ...\n";
    return 2;
  }
  try
  {
    std::vector<OpenMS::AASequence> peptides;
    std::vector<int> charges;
    for (int i = 2; i < argc; ++i)
    {
      const std::string arg = argv[i];
      const auto colon = arg.rfind(':');
      if (colon == std::string::npos)
      {
        std::cerr << "error: expected sequence:charge, got " << arg << "\n";
        return 2;
      }
      peptides.push_back(OpenMS::AASequence::fromString(arg.substr(0, colon)));
      charges.push_back(std::atoi(arg.c_str() + colon + 1));
    }

    int threads = 0;
    if (const char* env = std::getenv("ODIA_ORT_THREADS")) { threads = std::atoi(env); }

    ODIA::PeptDeepPredictor predictor(argv[1], true, threads);
    const auto ccs = predictor.predictCCS(peptides, charges);

    std::cout << "[";
    for (std::size_t i = 0; i < ccs.size(); ++i)
    {
      std::cout << (i ? ", " : "");
      // JSON has no bare nan; the reference emits NaN and Python accepts it.
      if (std::isnan(ccs[i])) { std::cout << "NaN"; }
      else { std::cout << std::fixed << std::setprecision(6) << ccs[i]; }
    }
    std::cout << "]\n";
  }
  catch (const std::exception& e)
  {
    std::cerr << "error: " << e.what() << "\n";
    return 1;
  }
  return 0;
}
