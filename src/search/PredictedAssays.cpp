// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <odia/search/PredictedAssays.h>

#include <odia/PeptDeepEncoder.h>

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/Residue.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

namespace ODIA::search
{
  namespace
  {
    using OpenMS::AASequence;
    using OpenMS::Residue;

    /// LibraryGenerator::predictFragmentIntensities's order: intensity, then
    /// m/z (both descending), then type, ordinal, charge.
    bool before(const PredictedFragment& a, const PredictedFragment& b)
    {
      if (a.intensity != b.intensity) { return a.intensity > b.intensity; }
      if (a.mz != b.mz) { return a.mz > b.mz; }
      if (a.type != b.type) { return a.type < b.type; }
      if (a.ordinal != b.ordinal) { return a.ordinal < b.ordinal; }
      return a.charge < b.charge;
    }

    /// Whether every fragment of @p decoy lies within @p ppm of one of @p target's.
    bool copies(const std::vector<PredictedFragment>& target, const std::vector<PredictedFragment>& decoy, double ppm)
    {
      std::vector<double> sorted;
      sorted.reserve(target.size());
      for (const auto& f : target) { sorted.push_back(f.mz); }
      std::sort(sorted.begin(), sorted.end());
      for (const auto& f : decoy)
      {
        const double tol = f.mz * ppm * 1e-6;
        const auto it = std::lower_bound(sorted.begin(), sorted.end(), f.mz - tol);
        if (it == sorted.end() || *it > f.mz + tol) { return false; }
      }
      return true;
    }
  }

  // ---- the PeptDeep model ----------------------------------------------------------

  PeptDeepFragmentModel::PeptDeepFragmentModel(const std::string& model_path, const std::string& instrument, double nce,
                                               int sessions)
  {
    instrument_ = PeptDeepEncoder::canonicalInstrument(instrument);
    if (instrument_.empty()) { throw std::invalid_argument("search:instrument: unknown instrument '" + instrument + "'"); }
    nce_ = nce > 0 ? static_cast<float>(nce) : PeptDeepEncoder::defaultNce(instrument_);
    // CPU only: the search runs where the extraction runs. One intra-op thread
    // per session (the recurrent model does not spread; PeptDeepPredictor.h).
    predictor_ = std::make_unique<PeptDeepPredictor>(model_path, false, 1, std::max(1, std::min(sessions, 16)));
  }

  PeptDeepFragmentModel::~PeptDeepFragmentModel() = default;

  std::vector<PeptDeepPredictor::Spectrum> PeptDeepFragmentModel::predict(const std::vector<AASequence>& peptides,
                                                                          const std::vector<int>& charges)
  {
    std::vector<PeptDeepPredictor::Failure> failures;   // unencodable peptides come back empty
    return predictor_->predictMS2(peptides, charges, nce_, instrument_, &failures);
  }

  std::string PeptDeepFragmentModel::describe() const
  {
    std::ostringstream s;
    // Nothing that depends on -threads: the report embeds this line.
    s << "PeptDeep MS2, instrument " << instrument_ << ", NCE " << nce_;
    return s.str();
  }

  // ---- one member ------------------------------------------------------------------

  void PredictedAssays::rank(const AASequence& peptide, int precursor_charge, const PeptDeepPredictor::Spectrum& spectrum,
                             double mz_min, double mz_max, std::vector<PredictedFragment>& out, std::size_t& above_floor)
  {
    out.clear();
    above_floor = 0;
    const std::size_t n = peptide.size();
    if (n < 2 || spectrum.positions == 0) { return; }
    float peak = 0.0f;
    for (std::size_t q = 0; q < spectrum.positions; ++q)
    {
      for (std::size_t ch = 0; ch < 4; ++ch) { peak = std::max(peak, spectrum.at(q, ch)); }
    }
    const float floor = peak * predicted_floor;
    const int max_z = std::min(2, std::max(1, precursor_charge));
    // The range test in fixed point, as the library stores m/z and as
    // searchDecoy tests it: a double bound from fromFixed() can sit a
    // rounding above the library's own extreme fragment and drop it.
    const MzFixed lo = toFixed(mz_min), hi = toFixed(mz_max);
    auto inRange = [&](double mz) { const MzFixed f = toFixed(mz); return f != MZ_INVALID && f >= lo && f <= hi; };
    for (std::size_t o = 1; o < n; ++o)
    {
      // Position q separates prefix q + 1 from suffix n - q - 1 (alphabase's
      // layout, as in LibraryGenerator): b_o at o - 1, y_o at n - 1 - o.
      const std::size_t b_position = o - 1, y_position = n - 1 - o;
      const AASequence prefix = peptide.getPrefix(o), suffix = peptide.getSuffix(o);
      for (int z = 1; z <= max_z; ++z)
      {
        const double b = prefix.getMZ(z, Residue::BIon);
        if (inRange(b) && b_position < spectrum.positions)
        {
          out.push_back({b, std::max(0.0f, spectrum.at(b_position, z == 1 ? 0 : 1)), FragmentType::B,
                         static_cast<std::uint8_t>(o), static_cast<std::int8_t>(z)});
        }
        const double y = suffix.getMZ(z, Residue::YIon);
        if (inRange(y) && y_position < spectrum.positions)
        {
          out.push_back({y, std::max(0.0f, spectrum.at(y_position, z == 1 ? 2 : 3)), FragmentType::Y,
                         static_cast<std::uint8_t>(o), static_cast<std::int8_t>(z)});
        }
      }
    }
    std::sort(out.begin(), out.end(), before);
    for (const auto& f : out) { above_floor += f.intensity > floor ? 1 : 0; }
  }

  std::size_t PredictedAssays::pairCount(std::size_t library_count, std::size_t target_above, std::size_t decoy_above,
                                         std::size_t target_all, std::size_t decoy_all, std::size_t min_fragments)
  {
    const std::size_t wanted = std::min(library_count, std::max(min_fragments, std::min(target_above, decoy_above)));
    const std::size_t count = std::min({wanted, target_all, decoy_all});
    return count >= min_fragments ? count : 0;
  }

  // ---- pairs -----------------------------------------------------------------------

  void PredictedAssays::predict(const Library& library, const std::vector<std::size_t>& targets, const std::vector<std::string>& decoys,
                                const DecoyRules& rules, FragmentModel& model, std::vector<PairPrediction>& out,
                                const std::vector<std::uint8_t>* counts)
  {
    if (decoys.size() != targets.size() || (counts && counts->size() != targets.size()))
    { throw std::invalid_argument("PredictedAssays::predict: one decoy sequence (and count) per target"); }
    out.resize(targets.size());
    const auto& pre = library.precursors();
    const double mz_min = fromFixed(rules.fragment_min), mz_max = fromFixed(rules.fragment_max);

    // Members of the pairs with a decoy, target then decoy: 2 * m peptides.
    // Both are predicted, always: the pair's rule reads each member through
    // its own prediction, and a target read out of the library instead brings
    // its library's own fragment count into a rule that must not see it.
    std::vector<std::size_t> pairs;
    for (std::size_t k = 0; k < targets.size(); ++k) { if (!decoys[k].empty()) { pairs.push_back(k); } }
    const auto m = static_cast<std::ptrdiff_t>(pairs.size());
    std::vector<AASequence> peptides(2 * pairs.size());
    std::vector<int> charges(2 * pairs.size(), 0);
    std::vector<char> parsed(pairs.size(), 1);
#pragma omp parallel for schedule(dynamic, 256)
    for (std::ptrdiff_t s = 0; s < m; ++s)
    {
      const auto j = static_cast<std::size_t>(s);
      const std::size_t k = pairs[j], i = targets[k];
      try
      {
        peptides[2 * j] = AASequence::fromString(std::string(library.strings().get(pre.modified_sequence[i])));
        peptides[2 * j + 1] = AASequence::fromString(decoys[k]);
      }
      // An unparsable pair is predicted as a stub and discarded; the model
      // sees the same input whatever failed, so the others do not move.
      catch (...) { parsed[j] = 0; peptides[2 * j] = peptides[2 * j + 1] = AASequence::fromString("PEPTIDEK"); }
      charges[2 * j] = charges[2 * j + 1] = std::max(1, static_cast<int>(pre.charge[i]));
    }
    const std::vector<PeptDeepPredictor::Spectrum> spectra = model.predict(peptides, charges);
    if (spectra.size() != peptides.size())
    { throw std::runtime_error("fragment model returned " + std::to_string(spectra.size()) + " spectra for " + std::to_string(peptides.size()) + " peptides"); }

#pragma omp parallel for schedule(dynamic, 256)
    for (std::ptrdiff_t s = 0; s < m; ++s)
    {
      const auto j = static_cast<std::size_t>(s);
      const std::size_t k = pairs[j], i = targets[k];
      PairPrediction& p = out[k];
      p = PairPrediction();
      if (!parsed[j] || spectra[2 * j].positions == 0 || spectra[2 * j + 1].positions == 0)
      { p.outcome = DecoyOutcome::Unpredictable; continue; }
      std::size_t above_t = 0, above_d = 0;
      rank(peptides[2 * j], charges[2 * j], spectra[2 * j], mz_min, mz_max, p.target, above_t);
      rank(peptides[2 * j + 1], charges[2 * j + 1], spectra[2 * j + 1], mz_min, mz_max, p.decoy, above_d);
      std::size_t count = 0;
      if (counts) { count = std::min<std::size_t>({(*counts)[k], p.target.size(), p.decoy.size()}); }
      else
      {
        count = pairCount(pre.transition_count[i], above_t, above_d, p.target.size(), p.decoy.size(), rules.min_fragments);
      }
      if (count < std::max<std::size_t>(1, rules.min_fragments))
      {
        p.outcome = DecoyOutcome::TooFewFragments;
        p.target.clear();
        p.decoy.clear();
        continue;
      }
      p.target.resize(count);
      p.decoy.resize(count);
      if (!counts && copies(p.target, p.decoy, rules.copy_ppm))
      {
        p.outcome = DecoyOutcome::Copy;
        p.target.clear();
        p.decoy.clear();
      }
    }
  }

  void PredictedAssays::apply(SearchSet& set, const Library& input, const DecoyRules& rules, FragmentModel& model,
                              const std::vector<std::uint8_t>& counts)
  {
    const std::size_t P = set.pairs();
    if (counts.size() != P) { throw std::invalid_argument("PredictedAssays::apply: one count per pair"); }
    Library& lib = set.library;
    auto& pre = lib.precursors();
    const auto& old = lib.transitions();
    Library::TransitionArrays built;
    std::vector<std::uint32_t> begin(2 * P, 0), count(2 * P, 0);
    std::vector<std::vector<PredictedFragment>> assays(2 * P);

    for (std::size_t base = 0; base < P; base += block_pairs)
    {
      const std::size_t last = std::min(P, base + block_pairs);
      std::vector<std::size_t> targets;
      std::vector<std::string> decoys(last - base);
      std::vector<std::uint8_t> fixed(counts.begin() + static_cast<std::ptrdiff_t>(base), counts.begin() + static_cast<std::ptrdiff_t>(last));
      for (std::size_t k = base; k < last; ++k) { targets.push_back(set.source[k]); }
      std::string failure;
      const auto n = static_cast<std::ptrdiff_t>(last - base);
#pragma omp parallel for schedule(dynamic, 64)
      for (std::ptrdiff_t s = 0; s < n; ++s)
      {
        try
        {
          const DecoyAssay a = searchDecoy(input, targets[static_cast<std::size_t>(s)], rules);
          if (a.outcome == DecoyOutcome::Made) { decoys[static_cast<std::size_t>(s)] = a.sequence; }
        }
        catch (const std::exception& e)
        {
#pragma omp critical(odia_predicted_failure)
          { if (failure.empty()) { failure = e.what(); } }
        }
      }
      if (!failure.empty()) { throw std::runtime_error("search decoys: " + failure); }
      std::vector<PairPrediction> predicted;
      predict(input, targets, decoys, rules, model, predicted, &fixed);
      for (std::size_t k = base; k < last; ++k)
      {
        PairPrediction& p = predicted[k - base];
        if (decoys[k - base].empty() || p.outcome != DecoyOutcome::Made || p.count() != counts[k])
        {
          throw std::logic_error("search: pair " + std::to_string(k) + " (" + std::string(set.modifiedSequence(k)) +
                                 ") could not be predicted with its " + std::to_string(counts[k]) + " fragments again");
        }
        assays[k] = std::move(p.target);
        assays[P + k] = std::move(p.decoy);
      }
    }

    std::size_t total = 0;
    for (const auto& a : assays) { total += a.size(); }
    Library::checkTransitionCapacity(0, total);
    built.product_mz.reserve(total);
    built.library_intensity.reserve(total);
    built.type.reserve(total);
    built.ordinal.reserve(total);
    built.charge.reserve(total);
    built.loss.reserve(total);
    for (std::size_t i = 0; i < 2 * P; ++i)
    {
      begin[i] = static_cast<std::uint32_t>(built.product_mz.size());
      for (const auto& f : assays[i])
      {
        built.product_mz.push_back(toFixed(f.mz));
        built.library_intensity.push_back(f.intensity);
        built.type.push_back(f.type);
        built.ordinal.push_back(f.ordinal);
        built.charge.push_back(f.charge);
        built.loss.push_back(LossType::None);
      }
      count[i] = static_cast<std::uint32_t>(assays[i].size());
      std::vector<PredictedFragment>().swap(assays[i]);
    }
    (void)old;
    lib.transitions() = std::move(built);
    pre.transition_begin = std::move(begin);
    pre.transition_count = std::move(count);
    lib.shrinkToFit();
  }
}
