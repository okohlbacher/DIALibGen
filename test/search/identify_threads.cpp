// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

// -threads for -run: a positive count is taken as given; 0 is the effective
// CPU budget (ODIA::physicalCores(), as generate counts it), so on a machine
// with several cores the search predicts with more than one ONNX Runtime
// session (it used to become 1). With the PeptDeep MS2 model as argument the
// model itself is built at -threads 0 and its sessions are counted.
//
//   identify_threads [peptdeep_ms2_dynamic.onnx]

#include <odia/Hardware.h>
#include <odia/search/PredictedAssays.h>
#include <odia/search/SearchParams.h>

#include <filesystem>
#include <iostream>
#include <string>
#include <thread>

namespace
{
  int failures = 0;

  void check(bool ok, const std::string& what)
  {
    if (!ok) { ++failures; std::cerr << "FAIL: " << what << '\n'; }
    else { std::cout << "ok   " << what << '\n'; }
  }
}

int main(int argc, char** argv)
{
  using ODIA::search::SearchParams;
  const int cores = static_cast<int>(ODIA::physicalCores());
  std::cout << "     physical cores: " << cores << " (" << ODIA::physicalCoresReason() << "); hardware threads: "
            << std::thread::hardware_concurrency() << '\n';

  check(SearchParams::resolveThreads(4) == 4 && SearchParams::resolveThreads(1) == 1, "a positive -threads is taken as given");
  check(SearchParams::resolveThreads(0) == cores && SearchParams::resolveThreads(-1) == cores,
        "-threads 0 is the physical cores: " + std::to_string(SearchParams::resolveThreads(0)));
  check(SearchParams::predictionSessions(1) == 1 && SearchParams::predictionSessions(8) == 8 &&
          SearchParams::predictionSessions(0) == 1 &&
          SearchParams::predictionSessions(100000) == cores,
        "one prediction session per thread, at least 1, at most the " + std::to_string(cores) + " physical cores");
  const int sessions = SearchParams::predictionSessions(SearchParams::resolveThreads(0));
  if (std::thread::hardware_concurrency() >= 4)
  { check(sessions > 1, "-threads 0 predicts with more than one session on this multi-core machine: " + std::to_string(sessions)); }
  else { std::cout << "     " << std::thread::hardware_concurrency() << " hardware threads: the multi-core check is skipped\n"; }

  if (argc > 1 && std::filesystem::exists(argv[1]))
  {
    const ODIA::search::PeptDeepFragmentModel model(argv[1], "timsTOF", 0.0, SearchParams::resolveThreads(0));
    check(model.sessions() == sessions, "the MS2 model built at -threads 0 has " + std::to_string(model.sessions()) + " sessions");
  }

  if (failures) { std::cerr << failures << " check(s) failed\n"; return 1; }
  return 0;
}
