// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause
//
// Prints the length groups PeptDeepEncoder would batch, and whether a
// mixed-length batch is refused. Both were untested, and they are each other's
// only backstop: padding is not inert, so a peptide encoded in the wrong group
// gets a materially different prediction.

#include <odia/PeptDeepEncoder.h>

#include <OpenMS/CHEMISTRY/AASequence.h>

#include <iostream>
#include <vector>

int main(int argc, char** argv)
{
  std::vector<OpenMS::AASequence> peptides;
  for (int i = 1; i < argc; ++i)
  {
    peptides.push_back(OpenMS::AASequence::fromString(argv[i]));
  }

  const auto groups = ODIA::PeptDeepEncoder::groupByLength(peptides);
  std::cout << "groups:";
  for (const auto& g : groups)
  {
    std::cout << " [";
    for (std::size_t i = 0; i < g.size(); ++i) { std::cout << (i ? "," : "") << g[i]; }
    std::cout << "]";
  }
  std::cout << "\n";

  // A mixed-length batch must be refused rather than silently padded.
  bool refused = false;
  try { (void)ODIA::PeptDeepEncoder::encode(peptides); }
  catch (const std::exception&) { refused = true; }
  std::cout << "mixed_batch_refused: " << (groups.size() > 1 ? (refused ? "yes" : "NO")
                                                             : "n/a") << "\n";
  return 0;
}
