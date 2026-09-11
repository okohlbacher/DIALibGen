// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause
//
// Dumps the raw-to-iRT line ODIA fits, and each standard's calibrated value,
// so they can be checked against OpenMS's own fixture rather than against a
// number this project wrote down.

#include <odia/LibraryGenerator.h>
#include <odia/PeptDeepPredictor.h>

#include <OpenMS/CHEMISTRY/AASequence.h>

#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

int main(int argc, char** argv)
{
  if (argc < 3)
  {
    std::fprintf(stderr, "usage: odia_irt_calibration <rt_model.onnx> <standards.tsv>\n");
    return 2;
  }
  try
  {
    const auto cal = ODIA::LibraryGenerator::fitIrtCalibration(argv[1], argv[2], false);
    std::printf("slope\t%.9f\nintercept\t%.9f\npeptides\t%zu\nmax_abs_error\t%.9f\n",
                cal.slope, cal.intercept, cal.peptides, cal.max_abs_error);

    // Each standard's raw prediction and its calibrated value, so a mismatch
    // can be attributed to the model or to the line rather than to "somewhere".
    std::vector<OpenMS::AASequence> peptides;
    std::vector<std::string> names;
    std::ifstream in(argv[2]);
    std::string line;
    while (std::getline(in, line))
    {
      if (line.empty() || line[0] == '#') { continue; }
      const auto tab = line.find('\t');
      if (tab == std::string::npos) { continue; }
      names.push_back(line.substr(0, tab));
      peptides.push_back(OpenMS::AASequence::fromString(names.back()));
    }
    ODIA::PeptDeepPredictor predictor(argv[1], false);
    const auto raw = predictor.predictRT(peptides);
    for (std::size_t i = 0; i < raw.size(); ++i)
    {
      std::printf("P\t%s\t%.9f\t%.9f\n", names[i].c_str(), double(raw[i]),
                  cal.apply(double(raw[i])));
    }
  }
  catch (const std::exception& e)
  {
    std::fprintf(stderr, "error: %s\n", e.what());
    return 1;
  }
  return 0;
}
