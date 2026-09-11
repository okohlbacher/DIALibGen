// Copyright (c) 2026, Oliver Kohlbacher and the DIALibraryGenerator authors.
// SPDX-License-Identifier: BSD-3-Clause
//
// Dumps predicted fragment spectra as JSON so they can be compared against
// test/peptdeep_reference.py.
//
// It takes several peptides in one call on purpose. Predicting one at a time
// leaves the whole batch path unexercised, and that path is where the row
// striding, the length grouping and the per-peptide charge live -- a version
// that gave every row the first peptide's charge, or read every spectrum from
// output row 0, passed a suite built on single-peptide calls.

#include <odia/PeptDeepPredictor.h>

#include <OpenMS/CHEMISTRY/AASequence.h>

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace
{
  /// JSON has no non-finite literals; Python's json accepts these spellings,
  /// and C++'s "nan"/"inf" would abort the parse before the comparison could
  /// report which value was bad.
  void printValue(float v)
  {
    if (std::isnan(v)) { std::cout << "NaN"; }
    else if (std::isinf(v)) { std::cout << (v > 0 ? "Infinity" : "-Infinity"); }
    else { std::cout << std::fixed << std::setprecision(6) << v; }
  }
}

int main(int argc, char** argv)
{
  if (argc < 5)
  {
    std::cerr << "usage: odia_predict_ms2 <model.onnx> <nce> <instrument> "
                 "<sequence:charge> [sequence:charge ...]\n";
    return 2;
  }
  try
  {
    const float nce = static_cast<float>(std::atof(argv[2]));
    const std::string instrument = argv[3];

    std::vector<OpenMS::AASequence> peptides;
    std::vector<int> charges;
    std::vector<std::string> labels;
    for (int i = 4; i < argc; ++i)
    {
      const std::string arg = argv[i];
      const auto colon = arg.rfind(':');
      if (colon == std::string::npos)
      {
        std::cerr << "error: expected sequence:charge, got " << arg << "\n";
        return 2;
      }
      labels.push_back(arg.substr(0, colon));
      peptides.push_back(OpenMS::AASequence::fromString(arg.substr(0, colon)));
      charges.push_back(std::atoi(arg.c_str() + colon + 1));
    }

    // ODIA_ORT_THREADS pins the intra-op thread count. Reduction order depends
    // on it, so leaving it to the core-count heuristic makes a cross-process
    // comparison noisy.
    int threads = 0;
    if (const char* env = std::getenv("ODIA_ORT_THREADS")) { threads = std::atoi(env); }

    ODIA::PeptDeepPredictor predictor(argv[1], true, threads);

    // With ODIA_MS2_TOLERATE_FAILURES the unencodable peptides are reported
    // rather than thrown, which is the only way to test that one bad peptide
    // does not take its whole chunk with it -- the behaviour predictRT
    // documents and predictMS2 did not have.
    std::vector<ODIA::PeptDeepPredictor::Failure> failures;
    const bool tolerate = std::getenv("ODIA_MS2_TOLERATE_FAILURES") != nullptr;
    const auto spectra = predictor.predictMS2(peptides, charges, nce, instrument,
                                              tolerate ? &failures : nullptr);

    std::cout << "{\"failed\": [";
    for (std::size_t k = 0; k < failures.size(); ++k)
    {
      std::cout << (k ? ", " : "") << failures[k].index;
    }
    std::cout << "], \"spectra\": [";
    for (std::size_t k = 0; k < spectra.size(); ++k)
    {
      const auto& s = spectra[k];
      std::cout << (k ? ", " : "") << "{\"peptide\": \"" << labels[k] << "\", \"charge\": "
                << charges[k] << ", \"shape\": [" << s.positions << ", "
                << ODIA::PeptDeepPredictor::Spectrum::CHANNELS << "], \"values\": [";
      for (std::size_t p = 0; p < s.positions; ++p)
      {
        std::cout << (p ? ", " : "") << "[";
        for (std::size_t c = 0; c < ODIA::PeptDeepPredictor::Spectrum::CHANNELS; ++c)
        {
          if (c) { std::cout << ", "; }
          printValue(s.at(p, c));
        }
        std::cout << "]";
      }
      std::cout << "]}";
    }
    std::cout << "]}\n";
  }
  catch (const std::exception& e)
  {
    std::cerr << "error: " << e.what() << "\n";
    return 1;
  }
  return 0;
}
