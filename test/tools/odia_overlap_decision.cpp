// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

// Which MS2 providers generate overlaps RT and CCS with. 0.11.4 overlapped
// them with any GPU: on CUDA, where RT and CCS run on the GPU too, the three
// models then contended for it and a whole-proteome run was 9-12 % slower.
// On CUDA without cuDNN, RT and CCS run on the CPU and overlap MS2 again.

#include <odia/LibraryGenerator.h>

#include <iostream>

int main()
{
  using P = ODIA::PeptDeepPredictor::Provider;
  const auto say = [](bool overlap) { return overlap ? "overlap" : "in order"; };
  std::cout << "CPU -> " << say(ODIA::LibraryGenerator::overlapsCpuPhases(P::CPU)) << "\n"
            << "CUDA -> " << say(ODIA::LibraryGenerator::overlapsCpuPhases(P::CUDA)) << "\n"
            << "CoreML -> " << say(ODIA::LibraryGenerator::overlapsCpuPhases(P::CoreML)) << "\n"
            << "CUDA, RT and CCS on the CPU -> " << say(ODIA::LibraryGenerator::overlapsCpuPhases(P::CUDA, true)) << "\n"
            << "CPU, RT and CCS on the CPU -> " << say(ODIA::LibraryGenerator::overlapsCpuPhases(P::CPU, true)) << "\n";
  return 0;
}
