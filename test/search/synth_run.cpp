// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

// Writes the synthetic -run fixture (synthetic_run.h) for the command-line
// tests: <dir>/library.tsv (a DIA-NN TSV library of targets), <dir>/run.mzML
// and <dir>/truth.tsv (the planted precursors).
//
//   identify_synth_run <dir> [peptides] [planted_fraction] [seed]

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
    std::cerr << "usage: identify_synth_run <dir> [peptides] [planted_fraction] [seed]\n";
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
    const synthrun::Fixture fx = synthrun::make(spec);
    ODIA::DIANNLibraryFile::storeTSV((dir / "library.tsv").string(), fx.library);
    synthrun::writeMzML(spec, fx, (dir / "run.mzML").string());
    synthrun::writeTruth(fx, (dir / "truth.tsv").string());
    std::cout << "synthetic run: " << fx.library.precursorCount() << " library precursors, " << fx.planted.size()
              << " planted, " << spec.windows << " windows, " << spec.run_s << " s at " << spec.cycle_s << " s per cycle\n";
    return 0;
  }
  catch (const std::exception& e)
  {
    std::cerr << "identify_synth_run: " << e.what() << "\n";
    return 1;
  }
}
