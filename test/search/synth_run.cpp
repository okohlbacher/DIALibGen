// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

// Writes the synthetic -run fixture (synthetic_run.h) for the command-line
// tests: <dir>/library.tsv (a DIA-NN TSV library of targets), <dir>/run.mzML
// and <dir>/truth.tsv (the planted precursors).
//
//   identify_synth_run <dir> [peptides] [planted_fraction] [seed] [im|im-scatter]
//
// "im": a diaPASEF run -- two 1/K0 bands per isolation window, the library's
// 1/K0 reading 0.08 high (identify_ion_mobility's fixture) -- and the true
// 1/K0 in truth.tsv. "im-scatter": the same, but the library's 1/K0 is also
// off by a per-precursor error (SD 0.025 around the offset, as a predictor's
// is), so the calibrated library 1/K0 is NOT where each precursor is: the
// fixture on which a 1/K0 read inside the extraction window would come out
// pulled towards the prediction.

#include "synthetic_run.h"

#include <odia/DIANNLibraryFile.h>

#include <exception>
#include <filesystem>
#include <iostream>
#include <string>

int main(int argc, char** argv)
{
  if (argc < 2)
  {
    std::cerr << "usage: identify_synth_run <dir> [peptides] [planted_fraction] [seed] [im|im-scatter]\n";
    return 2;
  }
  try
  {
    const std::filesystem::path dir = argv[1];
    std::filesystem::create_directories(dir);
    synthrun::Spec spec;
    if (argc > 2) { spec.peptides = std::stoul(argv[2]); }
    if (argc > 3) { spec.planted_fraction = std::stod(argv[3]); }
    if (argc > 4) { spec.seed = std::stoull(argv[4]); }
    const std::string im = argc > 5 ? argv[5] : "";
    if (im == "im" || im == "im-scatter")
    {
      spec.im_bands = 2;
      spec.im_library_offset = 0.08;
      if (im == "im-scatter") { spec.im_library_noise = 0.025; }
    }
    else if (!im.empty())
    {
      std::cerr << "identify_synth_run: unknown run kind '" << im << "' (im or im-scatter)\n";
      return 2;
    }
    const synthrun::Fixture fx = synthrun::make(spec);
    ODIA::DIANNLibraryFile::storeTSV((dir / "library.tsv").string(), fx.library);
    synthrun::writeMzML(spec, fx, (dir / "run.mzML").string());
    synthrun::writeTruth(fx, (dir / "truth.tsv").string(), &spec);
    std::cout << "synthetic run: " << fx.library.precursorCount() << " library precursors, " << fx.planted.size()
              << " planted, " << spec.windows << " windows" << (spec.im_bands ? " x " + std::to_string(spec.im_bands) + " 1/K0 bands" : std::string())
              << ", " << spec.run_s << " s at " << spec.cycle_s << " s per cycle\n";
    return 0;
  }
  catch (const std::exception& e)
  {
    std::cerr << "identify_synth_run: " << e.what() << "\n";
    return 1;
  }
}
