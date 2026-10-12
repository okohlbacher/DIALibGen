// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

// A synthetic in-memory library for the search tests: tryptic-looking peptides
// drawn from a seeded generator, charges 2 and 3, six y ions each with m/z from
// OpenMS, several peptides per protein group. A fixture, not biology.
#pragma once

#include <odia/Library.h>

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/Residue.h>

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <string>
#include <vector>

namespace synth
{
  inline int failures = 0;

#define CHECK(cond)                                                                   \
  do {                                                                                \
    if (!(cond)) { std::fprintf(stderr, "%s:%d: CHECK failed: %s\n", __FILE__, __LINE__, #cond); ++synth::failures; } \
  } while (0)

  struct Precursor
  {
    std::string sequence;
    std::string protein;
    int charge = 2;
    float rt = 0.0f;
    float im = std::nanf("");
    bool decoy = false;
    std::size_t fragments = 6;
  };

  /// Peptides of length 8..16 ending in K or R, no identical-residue runs that
  /// would defeat a shuffle; @p per_protein peptides per protein group.
  inline std::vector<Precursor> peptides(std::size_t n, std::uint64_t seed, std::size_t per_protein = 3)
  {
    static const std::string aa = "ADEFGHILMNPQSTVWY";
    std::mt19937_64 rng(seed);
    std::vector<Precursor> out;
    for (std::size_t k = 0; k < n; ++k)
    {
      const std::size_t len = 8 + rng() % 9;
      std::string s;
      while (s.size() + 1 < len)
      {
        const char c = aa[rng() % aa.size()];
        if (!s.empty() && s.back() == c) { continue; }
        s.push_back(c);
      }
      s.push_back(rng() % 2 ? 'K' : 'R');
      const float rt = static_cast<float>(rng() % 100000) / 1000.0f;   // 0 .. 100
      for (int z : {2, 3})
      {
        Precursor p;
        p.sequence = s;
        p.protein = "SYNPROT_" + std::to_string(k / per_protein);
        p.charge = z;
        p.rt = rt;
        out.push_back(p);
      }
    }
    return out;
  }

  inline void add(ODIA::Library& library, const Precursor& p)
  {
    auto& pre = library.precursors();
    auto& tr = library.transitions();
    const OpenMS::AASequence seq = OpenMS::AASequence::fromString(p.sequence);
    pre.mz.push_back(ODIA::toFixed(seq.getMZ(p.charge)));
    pre.irt.push_back(p.rt);
    pre.im.push_back(p.im);
    pre.ccs.push_back(std::nanf(""));
    pre.charge.push_back(static_cast<std::uint8_t>(p.charge));
    pre.decoy.push_back(p.decoy ? 1 : 0);
    pre.modified_sequence.push_back(library.strings().intern(p.sequence));
    pre.protein_group.push_back(library.strings().intern(p.protein));
    pre.transition_begin.push_back(static_cast<std::uint32_t>(tr.product_mz.size()));
    std::uint32_t written = 0;
    for (std::size_t ord = 2; written < p.fragments && ord + 1 < seq.size(); ++ord)
    {
      double mz = seq.getSuffix(ord).getMZ(1, OpenMS::Residue::YIon);
      if (p.decoy) { mz += 7.0; }   // a file decoy: shifted fragments, target's everything else
      tr.product_mz.push_back(ODIA::toFixed(mz));
      tr.library_intensity.push_back(1.0f / static_cast<float>(1 + written));
      tr.type.push_back(ODIA::FragmentType::Y);
      tr.ordinal.push_back(static_cast<std::uint8_t>(ord));
      tr.charge.push_back(1);
      tr.loss.push_back(ODIA::LossType::None);
      ++written;
    }
    pre.transition_count.push_back(written);
  }

  inline ODIA::Library library(const std::vector<Precursor>& precursors)
  {
    ODIA::Library lib;
    for (const auto& p : precursors) { add(lib, p); }
    return lib;
  }
}
