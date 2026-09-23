// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

/// search:candidates evidence -- which target-decoy pairs a search extracts,
/// chosen by fragment evidence in the run instead of at random.
///
///   pairs     every eligible target (CandidateSelector::eligible) that gets a
///             search decoy (searchDecoy, the decoy that is later searched),
///             with the m/z of each member's top prefilter_fragments
///             fragments: with search:intensities predicted, each member's
///             own most intense by the fragment model (PredictedAssays.h);
///             with library, the target's by library intensity and its
///             decoy's in the SAME slots;
///   sweep     one pass over the run's MS2 spectra: per isolation window, one
///             m/z-sorted index of both members' fragments; per spectrum, its
///             search:prefilter_top_peaks most intense peaks are matched at
///             search:prefilter_ppm. Per member: the most distinct indexed
///             fragments matched in one spectrum (depth), the spectra reaching
///             search:prefilter_depth, and the run time of the best spectrum
///             (deepest, then most matched intensity, then earliest);
///   choose    a pair is kept when its target OR its decoy reaches the depth
///             (pair-union); the kept set is checked to hold as many decoys
///             as targets (ratio guard) and capped to search:max_pairs by a
///             stratified, label-blind rank;
///   seeds     the calibration's seeds: targets that reach the depth
///             THEMSELVES and clearly beat their own decoy, one per peptide,
///             the best of each bin of the central library-RT range.
///
/// Label symmetry is the design rule: the index, the sweep and the matching
/// treat a decoy fragment exactly as a target fragment, and every decision
/// about a pair reads its two members through a symmetric function (the
/// larger depth, the larger evidence), never through the label. A selection
/// that depended on the label would make the target-decoy FDR estimate
/// anti-conservative. Seeds are the one target-only choice: calibration
/// decides no identification.
#pragma once

#include <odia/Library.h>
#include <odia/search/CandidateSelector.h>
#include <odia/search/PredictedAssays.h>
#include <odia/search/SearchParams.h>

#include <OpenMS/OPENSWATHALGO/DATAACCESS/SwathMap.h>

#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace ODIA::search
{
  /// What the sweep found for one pair member.
  struct MemberEvidence
  {
    std::uint8_t depth = 0;      ///< most distinct indexed fragments matched in one spectrum
    std::uint32_t spectra = 0;   ///< spectra in which it matched search:prefilter_depth or more
    float rt = std::numeric_limits<float>::quiet_NaN();   ///< run time of its best spectrum, seconds
    float intensity = 0.0f;      ///< matched peak intensity in that spectrum
  };

  /// The pairs the prefilter can choose from: every eligible target that gets
  /// a search decoy, in draw order. Member 2k is pair k's target, 2k + 1 its decoy.
  struct PrefilterPairs
  {
    std::vector<std::pair<std::uint64_t, std::size_t>> targets;   ///< (draw key, input library index)
    std::vector<std::uint8_t> fragments;   ///< per pair: fragments indexed per member (<= prefilter_fragments)
    /// Per pair: fragments per member of the searched assays (predicted: the
    /// pair's count, which PredictedAssays::apply must reproduce).
    std::vector<std::uint8_t> assay_fragments;
    /// Per member, prefilter_fragments slots (unused slots NaN): the target's
    /// top predicted fragments, and the decoy's fragments in the same slots.
    std::vector<float> fragment_mz;
    std::vector<double> precursor_mz;      ///< per pair
    std::vector<float> library_rt;         ///< per pair, the input library's RT
    std::vector<std::uint8_t> charge;      ///< per pair
    SelectionStats stats;                  ///< eligibility and decoy counts
    std::size_t processed = 0;             ///< targets of the ordered input examined (all unless stopped early)
    double decoy_seconds = 0.0;

    std::size_t size() const { return targets.size(); }
    const float* mz(std::size_t member) const { return fragment_mz.data() + member * SearchParams::prefilter_fragments; }
  };

  /// What the sweep read, for the provenance.
  struct SweepStats
  {
    std::size_t windows = 0;       ///< distinct isolation windows indexed
    std::size_t maps = 0;          ///< MS2 maps swept
    std::size_t spectra = 0;
    std::size_t peaks = 0;         ///< peaks matched against the index (at most top_peaks per spectrum)
    std::size_t entries = 0;       ///< index entries, both classes, every window
    std::size_t index_bytes = 0;
    double index_seconds = 0.0;
    double sweep_seconds = 0.0;
  };

  /// What choose() kept, and why.
  struct PrefilterSelection
  {
    std::vector<std::size_t> kept;                 ///< pair indices, ascending (draw order)
    std::vector<std::size_t> depth_targets, depth_decoys;   ///< members by depth 0 .. prefilter_fragments
    std::size_t targets_passing = 0;               ///< targets at search:prefilter_depth or more
    std::size_t decoys_passing = 0;
    std::size_t both_passing = 0;
    std::size_t union_pairs = 0;                   ///< pairs with either member passing
    std::size_t capped = 0;                        ///< union pairs removed by search:max_pairs
    std::size_t strata = 0;                        ///< non-empty strata of the cap (0: no cap applied)
    std::size_t kept_targets_passing = 0;          ///< among the kept pairs
    std::size_t kept_decoys_passing = 0;
  };

  class EvidencePrefilter
  {
  public:
    using Log = std::function<void(const std::string&)>;

    /// The pair universe: eligible targets (CandidateSelector::eligible with
    /// @p windows), each with its search decoy (searchDecoy with
    /// CandidateSelector::decoyRules); targets without a decoy are counted in
    /// stats.no_decoy and its reasons. Deterministic whatever the thread count.
    /// With search:intensities predicted, @p model predicts both members
    /// (required; std::invalid_argument without it).
    /// @p progress, when given, hears about the prediction every 30 s or so
    /// (it takes minutes on a whole-proteome library).
    static PrefilterPairs pairs(const Library& library, const SearchParams& params, const std::vector<IsolationWindow>& windows,
                                FragmentModel* model = nullptr, const Log& progress = Log());

    /// pairs() over targets already chosen and ordered (draw order), with the
    /// counts of the choice in @p stats. With @p stop_after, stops after the
    /// block in which that many pairs got a decoy (processed says how far).
    static PrefilterPairs pairsOf(const Library& library, const SearchParams& params,
                                  const std::vector<std::pair<std::uint64_t, std::size_t>>& ordered,
                                  const SelectionStats& stats, FragmentModel* model, std::size_t stop_after,
                                  const Log& progress = Log());

    /// One sweep over the MS2 maps of @p maps (MS1 maps are skipped): the
    /// evidence of every member, 2 * pairs.size() entries. A member is
    /// indexed in every map whose isolation window holds its precursor m/z;
    /// its evidence is the best over those maps, its spectra their sum.
    /// Deterministic whatever the thread count. @p swap_classes indexes every
    /// target fragment under its decoy's member and vice versa -- a test hook
    /// for label symmetry: the result must then be the exact swap.
    static std::vector<MemberEvidence> sweep(const PrefilterPairs& pairs, const std::vector<OpenSwath::SwathMap>& maps,
                                             const SearchParams& params, SweepStats* stats = nullptr, bool swap_classes = false);

    /// Pair-union at search:prefilter_depth, the ratio guard, and the
    /// stratified label-blind cap to search:max_pairs: strata are isolation
    /// window (the first of @p windows holding the precursor; 0 without
    /// windows) x library-RT decile of the union x charge; each gets its
    /// largest-remainder share of the cap, and within a stratum pairs rank by
    /// the better member's (depth, spectra), then by draw key.
    static PrefilterSelection choose(const PrefilterPairs& pairs, const std::vector<MemberEvidence>& evidence,
                                     const SearchParams& params, const std::vector<IsolationWindow>& windows);

    /// Calibration seeds: targets whose OWN depth reaches search:prefilter_depth
    /// and that beat their own decoy clearly (at least seed_min_spectra
    /// spectra at that depth and seed_decoy_factor times the decoy's), the
    /// best (depth, spectra, draw key) per stripped peptide sequence, at most
    /// calibration_seeds / calibration_seed_bins from each of
    /// calibration_seed_bins bins of @p scale (the central library-RT range,
    /// SearchSet::rt_robust); a target outside it is not a seed. Best first.
    static std::vector<SeedHint> seeds(const Library& library, const PrefilterPairs& pairs,
                                       const std::vector<MemberEvidence>& evidence, const SearchParams& params,
                                       const RtScale& scale);

    /// The ratio guard: the decoy:target count ratio of a candidate set must
    /// lie in [candidate_ratio_low, candidate_ratio_high]. Throws SearchAbort
    /// with the counts. Pair-union keeps it at exactly 1; this catches a
    /// selection that keeps members instead of pairs.
    static void checkRatio(std::size_t targets, std::size_t decoys);

    /// Everything: the search set of the kept pairs (CandidateSelector::
    /// fromTargets), its seeds, and the record in prefilter_json. Timings go
    /// to @p seconds (a JSON object), never into the set.
    /// With search:intensities predicted, @p model predicts both members of
    /// every pair and the searched set's assays (PredictedAssays::apply).
    static SearchSet select(const Library& library, const SearchParams& params, const std::vector<IsolationWindow>& windows,
                            const std::vector<OpenSwath::SwathMap>& maps, const Log& info, std::string* seconds = nullptr,
                            FragmentModel* model = nullptr);

    /// search:candidates random with search:intensities predicted: the
    /// random draw of CandidateSelector::select (search:subset lowest draw
    /// keys, then the first search:max_pairs with a decoy), whose pairs are
    /// predicted like the evidence prefilter's.
    static SearchSet selectRandom(const Library& library, const SearchParams& params, const std::vector<IsolationWindow>& windows,
                                  FragmentModel& model);
  };
}
