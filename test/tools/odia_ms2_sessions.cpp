// Spreading batches across sessions must not change a single bit. The chunking
// is built in serial order precisely so it cannot, but "cannot" is a claim
// about code that has to be checked against the model's actual output.
#include <odia/PeptDeepPredictor.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <vector>

int main(int argc, char** argv)
{
  if (argc != 2)
  {
    std::fprintf(stderr, "usage: odia_ms2_sessions <model.onnx>\n");
    return 2;
  }
  std::vector<OpenMS::AASequence> peptides;
  std::vector<int> charges;
  // One length must cross the batch boundary. Distinct sequences, interleaved
  // lengths, modifications and charges expose lost ordering or repeated rows.
  const std::string amino_acids = "ACDEFGHIKLMNPQRSTVWY";
  for (std::size_t i = 0; i < ODIA::PeptDeepPredictor::MAX_BATCH_ROWS + 17; ++i)
  {
    std::string sequence = "AAAPEPTIK";
    auto encoded = i;
    for (std::size_t j = 0; j < 3; ++j)
    {
      sequence[j] = amino_acids[encoded % amino_acids.size()];
      encoded /= amino_acids.size();
    }
    peptides.push_back(OpenMS::AASequence::fromString(sequence));
    charges.push_back(2 + static_cast<int>(i % 3));
    if (i % 64 == 0)
    {
      const auto length = 7 + i / 64;
      peptides.push_back(OpenMS::AASequence::fromString(
        "AC(Carbamidomethyl)M(Oxidation)" + std::string(length - 4, 'A') + "K"));
      charges.push_back(2 + static_cast<int>(i % 2));
    }
  }

  // ODIA_MS2_GPU asks for CUDA. Off by default so the determinism test stays
  // a CPU test wherever it runs; a GPU box would otherwise silently measure
  // something else than the machine that reported the baseline.
  const bool want_gpu = std::getenv("ODIA_MS2_GPU") != nullptr;

  auto run = [&](int sessions)
  {
    ODIA::PeptDeepPredictor p(argv[1], want_gpu, 1, sessions);
    const bool on_cuda = p.provider() == ODIA::PeptDeepPredictor::Provider::CUDA;
    if (want_gpu && !on_cuda)
    {
      // Falling back is the failure this whole exercise exists to catch: it
      // looks like a slow GPU rather than an absent one.
      std::fprintf(stderr, "ODIA_MS2_GPU was set but the provider is CPU\n");
      std::exit(3);
    }
    std::vector<ODIA::PeptDeepPredictor::Failure> f;
    const auto t0 = std::chrono::steady_clock::now();
    auto out = p.predictMS2(peptides, charges, 30.0f, "timsTOF", &f);
    if (!f.empty() || out.size() != peptides.size())
    {
      throw std::runtime_error("prediction failed or lost peptides");
    }
    for (std::size_t i = 0; i < out.size(); ++i)
    {
      const auto& spectrum = out[i];
      if (spectrum.positions != peptides[i].size() - 1 ||
          spectrum.intensities.size() != spectrum.positions * spectrum.CHANNELS)
      {
        throw std::runtime_error("prediction has an incorrect spectrum shape");
      }
      bool nonzero = false;
      for (const auto intensity : spectrum.intensities)
      {
        if (!std::isfinite(intensity)) { throw std::runtime_error("non-finite intensity"); }
        nonzero = nonzero || intensity > 0;
      }
      if (!nonzero) { throw std::runtime_error("empty predicted spectrum"); }
    }
    const double dt = std::chrono::duration<double>(
                        std::chrono::steady_clock::now() - t0).count();
    std::printf("%2d session(s) on %-4s: %6.2f s  %8.1f peptides/s  (%zu live)\n",
                sessions, on_cuda ? "CUDA" : "CPU", dt, peptides.size() / dt,
                p.sessionCount());
    return out;
  };

  const auto reference = run(1);
  for (const int s : {2, 8, 32})
  {
    const auto other = run(s);
    if (other.size() != reference.size()) { std::printf("  SIZE MISMATCH\n"); return 1; }
    std::size_t differing = 0, compared = 0;
    for (std::size_t i = 0; i < reference.size(); ++i)
    {
      if (other[i].intensities.size() != reference[i].intensities.size())
      {
        std::printf("  length mismatch at %zu\n", i);
        return 1;
      }
      for (std::size_t k = 0; k < reference[i].intensities.size(); ++k)
      {
        ++compared;
        // Bit-for-bit, not a tolerance: any difference at all means the split
        // is leaking into the arithmetic.
        if (other[i].intensities[k] != reference[i].intensities[k]) { ++differing; }
      }
    }
    std::printf("  vs 1 session: %zu of %zu values differ%s\n", differing, compared,
                differing ? "  <-- NOT IDENTICAL" : "  (identical)");
    if (differing) { return 1; }
  }
  return 0;
}
