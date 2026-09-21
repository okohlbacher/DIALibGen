// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

// The run-facing stages of Identifier: loading, calibration and extraction on
// stock OpenMS 3.5.0 (SwathFile::loadMzML, OpenSwathCalibrationWorkflow::
// performRTNormalization, OpenSwathWorkflow::performExtraction). Their
// contracts are documented in include/odia/search/Identifier.h.

#include <odia/search/Identifier.h>

namespace ODIA::search
{
  RunData Identifier::loadRun(const std::string& path)
  {
    (void)path;
    throw NotImplemented("search: loading the run is not implemented yet (built-in identification is in development)");
  }

  Calibration Identifier::calibrate(const Library& library, const SearchSet& set, RunData& run)
  {
    (void)library; (void)set; (void)run;
    throw NotImplemented("search: RT calibration is not implemented yet (built-in identification is in development)");
  }

  void Identifier::extract(const SearchSet& set, RunData& run, const Calibration& calibration, PeakGroups& out)
  {
    (void)set; (void)run; (void)calibration; (void)out;
    throw NotImplemented("search: extraction is not implemented yet (built-in identification is in development)");
  }
}
