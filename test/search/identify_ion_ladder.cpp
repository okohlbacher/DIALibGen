// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

// IonLadder (odia/IonLadder.h) gives every a, b, c, x, y and z ion bit for bit
// what OpenMS's getPrefix(i).getMZ(z, type) / getSuffix(i).getMZ(z, type)
// give: generation's fragments, the search's predicted assays and its decoys'
// fragment m/z all read it, and the search report must not move by a bit.
// Peptides with residue, N- and C-terminal modifications, every ordinal,
// charges 1 to 4. A peptide with an unknown residue X keeps OpenMS's path
// (which throws), as do ordinals outside 1 .. n - 1.

#include <odia/IonLadder.h>

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/Residue.h>

#include <cstring>
#include <iostream>
#include <string>
#include <vector>

namespace
{
  int failures = 0;

  void check(bool ok, const std::string& what)
  {
    std::cout << (ok ? "ok   " : "FAIL ") << what << '\n';
    if (!ok) { ++failures; }
  }

  bool sameBits(double a, double b) { return std::memcmp(&a, &b, sizeof a) == 0; }
}

int main()
{
  using OpenMS::AASequence;
  using OpenMS::Residue;
  const std::vector<std::string> peptides = {
    "PEPTIDEK", "ACDEFGHIKLMNPQRSTVWY", "C(Carbamidomethyl)PEPM(Oxidation)TIDEK",
    ".(Acetyl)SEQUENC(Carbamidomethyl)EK", "PEPS(Phospho)TIDE.(Amidated)", ".(Acetyl)M(Oxidation)PEPTIDER.(Amidated)",
    "GK", "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAR"};
  const std::vector<std::pair<Residue::ResidueType, bool>> types = {
    {Residue::AIon, true}, {Residue::BIon, true}, {Residue::CIon, true},
    {Residue::XIon, false}, {Residue::YIon, false}, {Residue::ZIon, false}};
  std::size_t compared = 0, differing = 0;
  for (const auto& text : peptides)
  {
    const AASequence peptide = AASequence::fromString(text);
    const ODIA::IonLadder ladder(peptide);
    for (OpenMS::Size i = 1; i < peptide.size(); ++i)
    {
      for (int z = 1; z <= 4; ++z)
      {
        for (const auto& [type, prefix] : types)
        {
          const double expected = prefix ? peptide.getPrefix(i).getMZ(z, type) : peptide.getSuffix(i).getMZ(z, type);
          ++compared;
          if (!sameBits(ladder.ion(peptide, i, z, type), expected))
          {
            ++differing;
            std::cerr << text << " " << type << " " << i << " " << z << ": " << ladder.ion(peptide, i, z, type) << " vs " << expected << '\n';
          }
        }
        if (!sameBits(ladder.b(peptide, i, z), peptide.getPrefix(i).getMZ(z, Residue::BIon)) ||
            !sameBits(ladder.y(peptide, i, z), peptide.getSuffix(i).getMZ(z, Residue::YIon)))
        { ++differing; }
      }
    }
  }
  check(differing == 0, "a, b, c, x, y, z ions bit-identical to OpenMS: " + std::to_string(compared) + " compared, " +
                          std::to_string(differing) + " differ");

  // An unknown residue X has no mass: OpenMS throws, and so does the ladder.
  const AASequence unknown = AASequence::fromString("PEPXIDEK");
  const ODIA::IonLadder x(unknown);
  bool threw = false;
  try { (void)x.ion(unknown, 5, 1, Residue::YIon); } catch (...) { threw = true; }
  check(threw, "a peptide with X keeps OpenMS's path (it throws)");

  if (failures) { std::cerr << failures << " check(s) failed\n"; return 1; }
  return 0;
}
