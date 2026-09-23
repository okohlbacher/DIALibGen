// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

// A deterministic stand-in for the PeptDeep MS2 model in the search tests
// (search:intensities predicted without model files). One real-looking rule
// carries the point: a y ion that starts with P is the most intense of its
// peptide (the proline effect), so a predictor picks the same short,
// Pro-directed fragments for every peptide that has them -- fragments whose
// m/z recur across a proteome. Everything else gets a pseudo-random
// intensity that depends only on the two residues around the cleavage.
#pragma once

#include <odia/search/PredictedAssays.h>

#include <OpenMS/CHEMISTRY/AASequence.h>

#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

namespace toy
{
  inline double unit(std::uint64_t z)
  {
    z += 0x9E3779B97F4A7C15ull;
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
    z ^= z >> 31;
    return static_cast<double>(z >> 11) / 9007199254740992.0;
  }

  class ProlineModel : public ODIA::search::FragmentModel
  {
  public:
    std::size_t calls = 0;
    std::size_t peptides_seen = 0;

    std::vector<ODIA::PeptDeepPredictor::Spectrum> predict(const std::vector<OpenMS::AASequence>& peptides,
                                                           const std::vector<int>& charges) override
    {
      (void)charges;
      ++calls;
      peptides_seen += peptides.size();
      std::vector<ODIA::PeptDeepPredictor::Spectrum> out(peptides.size());
      for (std::size_t k = 0; k < peptides.size(); ++k)
      {
        const OpenMS::AASequence& p = peptides[k];
        const std::size_t n = p.size();
        if (n < 2) { continue; }
        auto& s = out[k];
        s.positions = n - 1;
        s.intensities.assign(s.positions * ODIA::PeptDeepPredictor::Spectrum::CHANNELS, 0.0f);
        for (std::size_t q = 0; q + 1 < n; ++q)
        {
          const char left = p[q].getOneLetterCode()[0], right = p[q + 1].getOneLetterCode()[0];
          const std::uint64_t key = (static_cast<std::uint64_t>(static_cast<unsigned char>(left)) << 8) |
                                    static_cast<unsigned char>(right);
          const float y = right == 'P' ? 1.0f : static_cast<float>(0.1 + 0.3 * unit(key));
          const float b = static_cast<float>(0.05 + 0.2 * unit(key + 7919));
          float* row = s.intensities.data() + q * ODIA::PeptDeepPredictor::Spectrum::CHANNELS;
          row[0] = b;            // b z1
          row[1] = 0.02f * b;    // b z2
          row[2] = y;            // y z1
          row[3] = 0.02f * y;    // y z2
        }
      }
      return out;
    }

    std::string describe() const override { return "toy proline model"; }
  };

  /// The same Pro rule on a STEEP intensity scale: a fragment's intensity
  /// spans twelve orders of magnitude, so most of a peptide's fragments fall
  /// below the floor (PredictedAssays::predicted_floor of its base peak), as
  /// a real MS2 model's do -- on a timsTOF library about 17 of a target's 40
  /// fragments are above it. A library written with its own fragment cap then
  /// holds MORE transitions for a target than the model puts above the floor
  /// for it, which is the case the pair's count rule must not read out of the
  /// library (identify_decoy_exchangeability, section 3b).
  class SteepProlineModel : public ODIA::search::FragmentModel
  {
  public:
    std::size_t peptides_seen = 0;

    std::vector<ODIA::PeptDeepPredictor::Spectrum> predict(const std::vector<OpenMS::AASequence>& peptides,
                                                           const std::vector<int>& charges) override
    {
      (void)charges;
      peptides_seen += peptides.size();
      std::vector<ODIA::PeptDeepPredictor::Spectrum> out(peptides.size());
      for (std::size_t k = 0; k < peptides.size(); ++k)
      {
        const OpenMS::AASequence& p = peptides[k];
        const std::size_t n = p.size();
        if (n < 2) { continue; }
        auto& s = out[k];
        s.positions = n - 1;
        s.intensities.assign(s.positions * ODIA::PeptDeepPredictor::Spectrum::CHANNELS, 0.0f);
        for (std::size_t q = 0; q + 1 < n; ++q)
        {
          const char left = p[q].getOneLetterCode()[0], right = p[q + 1].getOneLetterCode()[0];
          const std::uint64_t key = (static_cast<std::uint64_t>(static_cast<unsigned char>(left)) << 8) |
                                    static_cast<unsigned char>(right);
          const float y = right == 'P' ? 1.0f : static_cast<float>(std::pow(10.0, -12.0 * unit(key)));
          const float b = static_cast<float>(0.3 * std::pow(10.0, -12.0 * unit(key + 7919)));
          float* row = s.intensities.data() + q * ODIA::PeptDeepPredictor::Spectrum::CHANNELS;
          row[0] = b;            // b z1
          row[1] = 0.02f * b;    // b z2
          row[2] = y;            // y z1
          row[3] = 0.02f * y;    // y z2
        }
      }
      return out;
    }

    std::string describe() const override { return "toy steep proline model"; }
  };
}
