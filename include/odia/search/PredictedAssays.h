// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

/// search:intensities predicted -- both members of a search pair get their
/// fragments and intensities from ONE fragment model, each from its own
/// sequence, by ONE rule.
///
/// Why. Paired target-decoy competition assumes a decoy wins against a null
/// (absent) target as often as it loses. A decoy that keeps its target's
/// fragment slots and intensities (search:intensities library) breaks that:
/// the slots were chosen by a predictor for the TARGET's sequence -- the
/// compositions real peptides fragment into (Pro-directed y ions, say) -- and
/// the decoy puts random compositions into them. Measured on entrapment
/// libraries (known-absent targets and their own decoys): targets reached the
/// prefilter depth 1.18x (timsTOF) and 1.11x (Astral) as often as their
/// decoys, and won 240 : 168 of the pairs identified at q <= 0.01, so q-values
/// came out about 1.4x too optimistic. Predicting the decoy from its own
/// sequence with the target's model makes it what a null target is: a
/// peptide nobody measured, seen through the same predictor.
///
/// The rule, the same for both members of every pair:
///   * the universe: b and y ions of ordinal 1 .. n-1, fragment charge 1 ..
///     min(2, precursor charge), no neutral loss, m/z within the library's
///     TARGET fragment range (DecoyRules::fragment_min/max);
///   * the intensity: the model's, base-peak normalised; "above the floor"
///     when above predicted_floor of the member's base peak (the generator's
///     min_relative_intensity);
///   * the rank: intensity, then m/z (descending), then ion type, ordinal and
///     charge -- LibraryGenerator::predictFragmentIntensities's order;
///   * the count, per pair: min(the target's library transition count,
///     max(min_assay_fragments, fewer of the two members' above-floor counts),
///     each member's universe). A symmetric function of the pair: both
///     members get the same number, their own top fragments;
///   * a pair whose decoy's chosen fragments all lie within copy_ppm of the
///     target's is a copy and leaves the search whole, as does one the model
///     cannot predict (Unpredictable) or with fewer than min_assay_fragments.
/// For a library DIALibGen generated with the same model, instrument and NCE,
/// the target's assay is its library assay (up to the count). That is checked,
/// not assumed (libraryIsModel: a sample of targets re-predicted), and then
/// only the decoys are predicted -- half the model's work; a library from
/// another predictor, or one whose intensities were rewritten, has both
/// members predicted.
#pragma once

#include <odia/Library.h>
#include <odia/PeptDeepPredictor.h>
#include <odia/search/CandidateSelector.h>
#include <odia/search/SearchDecoys.h>

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace OpenMS { class AASequence; }

namespace ODIA::search
{
  /// A fragment intensity model: one spectrum per peptide, PeptDeep's layout
  /// ((n - 1) positions x 8 channels; channels 0-3 are b z1, b z2, y z1, y z2).
  /// A peptide it cannot predict gets an empty spectrum (positions 0). Must be
  /// deterministic in its input and safe to call from one thread at a time.
  class FragmentModel
  {
  public:
    virtual ~FragmentModel() = default;
    virtual std::vector<PeptDeepPredictor::Spectrum> predict(const std::vector<OpenMS::AASequence>& peptides,
                                                             const std::vector<int>& charges) = 0;
    /// One line for the log and the report: nothing that depends on the
    /// thread count (the report must not).
    virtual std::string describe() const = 0;
  };

  /// PeptDeep's MS2 model (PeptDeepPredictor::predictMS2) on the CPU.
  class PeptDeepFragmentModel : public FragmentModel
  {
  public:
    /// @param nce <= 0: the instrument's default. Throws on an unknown instrument.
    PeptDeepFragmentModel(const std::string& model_path, const std::string& instrument, double nce, int sessions);
    ~PeptDeepFragmentModel() override;
    std::vector<PeptDeepPredictor::Spectrum> predict(const std::vector<OpenMS::AASequence>& peptides,
                                                     const std::vector<int>& charges) override;
    std::string describe() const override;
    const std::string& instrument() const { return instrument_; }
    float nce() const { return nce_; }

  private:
    std::unique_ptr<PeptDeepPredictor> predictor_;
    std::string instrument_;
    float nce_ = 30.0f;
  };

  /// One predicted fragment of one pair member.
  struct PredictedFragment
  {
    double mz = 0.0;
    float intensity = 0.0f;
    FragmentType type = FragmentType::Y;
    std::uint8_t ordinal = 0;
    std::int8_t charge = 1;
  };

  /// One pair's predicted assays: both members' fragments in rank order,
  /// already cut to the pair's count (empty unless outcome is Made).
  struct PairPrediction
  {
    DecoyOutcome outcome = DecoyOutcome::Made;
    std::vector<PredictedFragment> target, decoy;
    std::size_t count() const { return target.size(); }
  };

  class PredictedAssays
  {
  public:
    /// Base-peak fraction a predicted fragment must exceed to count as
    /// "above the floor" (LibraryGenerator's min_relative_intensity).
    static constexpr float predicted_floor = 1e-4f;
    /// Pairs predicted per model call.
    static constexpr std::size_t block_pairs = 8192;

    /// @p peptide's fragment universe with @p spectrum's intensities, in rank
    /// order (see the file comment); @p above_floor receives how many lead
    /// the list above the floor. Deterministic.
    static void rank(const OpenMS::AASequence& peptide, int precursor_charge, const PeptDeepPredictor::Spectrum& spectrum,
                     double mz_min, double mz_max, std::vector<PredictedFragment>& out, std::size_t& above_floor);

    /// The pair's count (see the file comment); 0 when it cannot reach
    /// @p min_fragments.
    static std::size_t pairCount(std::size_t library_count, std::size_t target_above, std::size_t decoy_above,
                                 std::size_t target_all, std::size_t decoy_all, std::size_t min_fragments);

    /// Predict and rank the pairs of @p targets (input library indices) whose
    /// decoy sequences are @p decoys (DecoyAssay::sequence; an empty string
    /// marks a pair without a decoy, left with its outcome untouched in
    /// @p out). With @p counts, each pair takes that many fragments (a count
    /// fixed by an earlier call, so a later one on a different batch cannot
    /// change it); without, pairCount decides. @p out is resized to
    /// targets.size(); entries with an empty decoy keep their outcome.
    ///
    /// With @p library_targets, the targets are not predicted: their library
    /// assays are the model's own (libraryIsModel), ranked by the same rule.
    static void predict(const Library& library, const std::vector<std::size_t>& targets, const std::vector<std::string>& decoys,
                        const DecoyRules& rules, FragmentModel& model, std::vector<PairPrediction>& out,
                        const std::vector<std::uint8_t>* counts = nullptr, bool library_targets = false);

    /// A target's library assay as predicted fragments, in the rule's rank
    /// order; empty when a transition is not a b/y ion without loss.
    static std::vector<PredictedFragment> libraryAssay(const Library& library, std::size_t i);

    /// Whether @p library's target assays are @p model's predictions under
    /// the search's rule: of @p sample (target indices), how many (@p matched)
    /// have exactly their library fragments as the model's top fragments,
    /// intensities within library_match_tolerance; true when at least
    /// library_match_share of them do.
    static bool libraryIsModel(const Library& library, const std::vector<std::size_t>& sample, const DecoyRules& rules,
                               FragmentModel& model, std::size_t* matched = nullptr);
    /// Relative intensity tolerance and share of matching targets of libraryIsModel.
    /// Batches of another composition move a prediction in its last bits, which
    /// can swap two nearly tied fragments (1 % of the targets of a test).
    static constexpr float library_match_tolerance = 1e-3f;
    static constexpr double library_match_share = 0.95;
    /// Targets libraryIsModel re-predicts.
    static constexpr std::size_t library_check_sample = 2000;

    /// Replace the transitions of every precursor of @p set (a search set:
    /// targets at [0, P), decoys at [P, 2P)) by its member's predicted assay,
    /// @p counts[k] fragments for pair k. The decoys are rebuilt from @p input
    /// (searchDecoy), exactly as the set was built. Throws std::logic_error
    /// when a pair cannot be predicted now although it could before.
    static void apply(SearchSet& set, const Library& input, const DecoyRules& rules, FragmentModel& model,
                      const std::vector<std::uint8_t>& counts, bool library_targets = false);
  };
}
