// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

// Which MS2 providers generate overlaps RT and CCS with. 0.11.4 overlapped
// them with any GPU: on CUDA, where RT and CCS run on the GPU too, the three
// models then contended for it and a whole-proteome run was 9-12 % slower.
// Where RT and CCS run on the CPU (CoreML; CUDA without cuDNN) they overlap
// MS2 only when asked for (ODIA_OVERLAP=1): on CoreML the overlap slowed the
// GPU's host side more than it gained (in order by default from 0.13.1 on
// CoreML, and after 0.13.1 on CUDA without cuDNN too).

#include <odia/LibraryGenerator.h>

#include <iostream>

int main()
{
  using P = ODIA::PeptDeepPredictor::Provider;
  const auto say = [](bool overlap) { return overlap ? "overlap" : "in order"; };
  std::cout << "CPU -> " << say(ODIA::LibraryGenerator::overlapsCpuPhases(P::CPU)) << "\n"
            << "CUDA -> " << say(ODIA::LibraryGenerator::overlapsCpuPhases(P::CUDA)) << "\n"
            << "CoreML -> " << say(ODIA::LibraryGenerator::overlapsCpuPhases(P::CoreML)) << "\n"
            << "CoreML, ODIA_OVERLAP=1 -> " << say(ODIA::LibraryGenerator::overlapsCpuPhases(P::CoreML, true, true)) << "\n"
            << "CUDA, ODIA_OVERLAP=1 -> " << say(ODIA::LibraryGenerator::overlapsCpuPhases(P::CUDA, false, true)) << "\n"
            << "CUDA, RT and CCS on the CPU -> " << say(ODIA::LibraryGenerator::overlapsCpuPhases(P::CUDA, true)) << "\n"
            << "CUDA, RT and CCS on the CPU, ODIA_OVERLAP=1 -> " << say(ODIA::LibraryGenerator::overlapsCpuPhases(P::CUDA, true, true)) << "\n"
            << "CPU, RT and CCS on the CPU -> " << say(ODIA::LibraryGenerator::overlapsCpuPhases(P::CPU, true)) << "\n"
            << "CPU, ODIA_OVERLAP=1 -> " << say(ODIA::LibraryGenerator::overlapsCpuPhases(P::CPU, true, true)) << "\n";
  return 0;
}
