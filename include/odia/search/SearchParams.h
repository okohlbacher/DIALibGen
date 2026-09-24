// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

/// Settings of the built-in identification step (`-run`, TOPP subsection
/// `search:`). See docs/design/built-in-identification.md.
///
/// One struct, filled from the TOPP options in src/Refine.cpp and read by every
/// stage of ODIA::search. Nothing here is a user option unless it has a
/// `search:` twin; the guard band below is a constant of the method.
#pragma once

#include <odia/LibraryGenerator.h>

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>

namespace ODIA::search
{
  /// How the run is held while it is searched (`search:readoptions`).
  enum class ReadMode
  {
    Auto,     ///< cache for runs above 3 GB, normal otherwise (the loader decides)
    Normal,   ///< the whole run in memory
    Cache     ///< per-window cache files in search:cache_dir, read back on demand
  };

  const char* toString(ReadMode m);
  ReadMode parseReadMode(const std::string& s);   ///< throws std::invalid_argument

  /// Where the fragments and intensities of both members of a search pair come
  /// from (`search:intensities`).
  enum class Intensities
  {
    /// Both the target and its decoy are predicted by the same fragment model
    /// from their OWN sequences, and each takes its own most intense fragments
    /// (PredictedAssays.h). A decoy is then built the way a null target is,
    /// which paired competition needs.
    Predicted,
    /// The library's target assay; the decoy re-uses its target's fragment
    /// slots and intensities with recomputed m/z. Known to make decoys weaker
    /// than null targets (entrapment targets beat their own decoys about
    /// 1.4 : 1 at q <= 0.01): a comparison and test hook, not a way to search.
    /// Kept because the exchangeability test measures Predicted against it,
    /// and only so: every run that selects it is warned about, in the log and
    /// in the report's search.warnings.
    Library
  };

  const char* toString(Intensities i);
  Intensities parseIntensities(const std::string& s);   ///< throws std::invalid_argument

  /// The decoy methods the search accepts: shuffle and pseudo_reverse, both
  /// of which keep the target's termini. `mutate` is refused (its substitution
  /// table is DIA-NN's), `reverse` (it moves the C-terminal K/R, which makes
  /// decoys separable from tryptic targets) and `none`.
  DecoyMethod parseSearchDecoyMethod(const std::string& s);   ///< throws std::invalid_argument
  const char* searchDecoyName(DecoyMethod m);

  struct SearchParams
  {
    // ---- candidates ----------------------------------------------------------
    /// Candidate selection:
    ///   * "evidence" (default): target-decoy pairs with fragment evidence in
    ///     the run (EvidencePrefilter.h): a pair is kept when its target OR its
    ///     decoy has prefilter_depth of its top predicted fragments among the
    ///     prefilter_top_peaks most intense peaks of one MS2 spectrum, within
    ///     prefilter_ppm -- the same rule for both classes -- and the kept
    ///     pairs are capped to max_pairs label-blind;
    ///   * "random": a deterministic, paired, label-blind random subset of the
    ///     library's targets (subset), capped to max_pairs.
    std::string candidates = "evidence";
    /// Targets drawn by the random selection before the cap; 0 = every
    /// eligible target. Random selection only.
    std::size_t subset = 0;
    /// Cap on target-decoy pairs after selection; 0 = no cap. Search time is
    /// linear in pairs, identifications in the share of the library present
    /// in the run. Measured on a 6.4 GB Orbitrap Astral run with a library of
    /// which DIA-NN identified 2.9 % (8 threads, random selection): 40,000
    /// pairs gave 1,045 identifications (tuning's validation cohort 93 units,
    /// below its floor of 100), 100,000 gave 2,804 (250) in 9 min, 200,000
    /// gave 5,821 (454) in 20 min at 4.5 GB peak.
    std::size_t max_pairs = 200000;
    /// Evidence prefilter: the co-occurrence depth a pair member must reach
    /// (distinct fragments of its top prefilter_fragments seen in one spectrum).
    int prefilter_depth = 5;
    /// Evidence prefilter: the most intense peaks of each MS2 spectrum it reads.
    std::size_t prefilter_top_peaks = 1000;
    /// Evidence prefilter: fragment match tolerance, ppm either side.
    double prefilter_ppm = 10.0;
    /// How the in-memory search decoys are built (decoys already in the library
    /// file are ignored by the search and left untouched).
    DecoyMethod decoys = DecoyMethod::Shuffle;
    /// Where both members' fragments and intensities come from (Intensities).
    Intensities intensities = Intensities::Predicted;
    /// search:intensities predicted: the PeptDeep MS2 model file, its content
    /// hash (DIANNLibraryFile::hashFile; what the report records instead of
    /// the path), the instrument and the NCE it is run with. An empty
    /// instrument means "timsTOF on an ion-mobility run, QE otherwise"; an
    /// NCE <= 0 the instrument's default (PeptDeepEncoder::defaultNce).
    std::string ms2_model;
    std::string ms2_model_hash;
    std::string instrument;
    double nce = -1.0;
    /// Where instrument and NCE came from, for the provenance ("search:instrument",
    /// "library", "run").
    std::string instrument_source;
    /// Salt of the candidate draw. Changes WHICH pairs are searched, never how.
    std::uint64_t seed = 42;
    /// Extraction passes. 1 in this version.
    int passes = 1;

    // ---- extraction ----------------------------------------------------------
    /// Full width of the RT extraction window in seconds; 0 = from the calibration.
    double rt_window = 0.0;
    /// Full width of the fragment m/z extraction window in ppm; 0 = automatic.
    double mz_ppm = 0.0;
    /// Ion mobility on a diaPASEF run (window 1/K0 limits AND a per-peak 1/K0
    /// array): full width of the 1/K0 extraction window. 0 = automatic, from
    /// the ion-mobility calibration (im_window_padding x its residuals,
    /// clamped to [im_window_min, im_window_max]); > 0 = this width; -1 = off:
    /// the run is searched by m/z and RT only, as before M2, and the report
    /// carries no 1/K0. Ignored on a run without ion mobility.
    double im_window = 0.0;
    /// MS1 traces and MS1 sub-scores.
    bool ms1 = true;
    /// Score the RT and 1/K0 deviation sub-scores. false removes them from the
    /// fit (an ablation switch for tuning, whose purpose is to correct them).
    bool rt_im_scores = true;
    /// Precursors per extraction batch within one isolation window, inside one
    /// OpenSWATH call (OpenSwathWorkflow's batchSize counts compounds).
    std::size_t batch_size = 500;
    /// Precursors per OpenSWATH extraction call; target-decoy pairs stay together.
    std::size_t chunk = 20000;

    // ---- calibration ---------------------------------------------------------
    double calibration_min_rsq = 0.70;
    double calibration_min_coverage = 0.30;
    /// Fall back to a linear map over the library RT range when the calibration
    /// fails. A test hook, recorded in the provenance; off = fail closed.
    bool allow_bootstrap = false;

    // ---- run input -----------------------------------------------------------
    ReadMode readoptions = ReadMode::Auto;
    /// Where cache files go; empty = the directory of -out_ids.
    std::string cache_dir;

    // ---- scoring and run-level guards -----------------------------------------
    /// Abort when fewer target precursors than this are identified at q <= 0.01.
    std::size_t min_ids = 200;
    /// Abort when more than this fraction of the scored target precursors is
    /// identified at q <= 0.01: no honest null looks like that.
    double max_target_fraction = 0.5;
    /// Precursors (targets AND decoys) up to this precursor q-value go into the report.
    double report_max_q = 0.10;
    /// Protein-group prefix marking entrapment proteins; empty = no entrapment estimate.
    std::string entrapment_tag;
    /// Also run the label-swap and random-label self-checks, aborting on
    /// failure (validation gate (e): on every run). They catch a classifier
    /// that leaks labels. Decoys built weaker than null targets they cannot
    /// see, and neither can the null-pair balance below in general: that
    /// asymmetry grows with the score and sits in the tail where the 1 % cut
    /// falls. Only known nulls see it: search:entrapment_tag's winner test.
    bool selftest = true;

    /// Threads for extraction and the classifier's folds (TOPP -threads).
    int threads = 1;

    /// Why the caller needs observed 1/K0 values (e.g. "-write_im"); empty =
    /// it does not. When set, a run without ion mobility is refused right
    /// after it is loaded, before anything is searched. How the tool was
    /// asked, not how the search runs: not part of toJson.
    std::string require_ion_mobility;

    // ---- constants of the method, not options --------------------------------
    /// The q-value an identification is counted at in every guard and log line.
    static constexpr double identification_q = 0.01;
    /// Decoy precursors with peak groups per target precursor with peak groups.
    /// A target and its decoy share precursor m/z, RT, 1/K0 and isolation window,
    /// so extraction should find candidates for both about equally often; outside
    /// this band the pairing or the extraction is broken and the FDR meaningless.
    static constexpr double decoy_ratio_low = 0.8;
    static constexpr double decoy_ratio_high = 1.25;
    /// Pooled vs paired identification count ratio ABOVE which the provenance
    /// carries a warning. Below 1 is expected: the pooled estimator also counts
    /// decoys that light up with their present target and lose their pair.
    static constexpr double pooled_vs_paired_warn = 2.0;
    /// The null-pair balance (a bulk diagnostic, NOT an exchangeability test:
    /// it cannot see an asymmetry confined to the high-scoring tail) looks at
    /// the lowest this share of pair winners ...
    static constexpr double null_balance_fraction = 0.25;
    /// ... and warns when targets and decoys win there more unevenly than this
    /// many binomial standard deviations.
    static constexpr double null_balance_warn_z = 3.0;
    /// Fewest transitions a search assay may have, for targets and decoys alike
    /// (the generator's own floor).
    static constexpr std::size_t min_assay_fragments = 3;
    /// The most OpenMP threads a stock OpenSWATH call (calibration, extraction)
    /// gets, whatever -threads says: stock 3.5.0 serialises feature scoring on
    /// a process-wide lock, and more threads only spin (see RunStages.cpp).
    static constexpr int openswath_max_threads = 8;
    /// Evidence prefilter: the fragments indexed per pair member, its own most
    /// intense (search:intensities predicted: each member by its own
    /// prediction; library: by the target's library intensity, in the same
    /// slots for the decoy).
    static constexpr std::size_t prefilter_fragments = 6;
    /// With search:entrapment_tag: pairs whose target is an entrapment peptide
    /// (known absent) must be won by the target and by its decoy equally
    /// often. A one-sided binomial z above this at q <= 0.01 or q <= 0.10, at
    /// any level, is warned about: the decoys are weaker than null targets.
    static constexpr double entrapment_winner_warn_z = 3.0;
    /// Evidence-seeded calibration: at most this many seeds, spread over the
    /// library RT range in calibration_seed_bins bins. The range is the
    /// central calibration_rt_quantile .. 1 - calibration_rt_quantile of the
    /// eligible targets' library RT, so a few out-of-range library rows
    /// cannot squeeze the bins; a seed outside it (a kit precursor included)
    /// is not used.
    static constexpr std::size_t calibration_seeds = 2000;
    static constexpr std::size_t calibration_seed_bins = 100;
    static constexpr double calibration_rt_quantile = 0.001;
    /// An evidence seed must beat its own decoy clearly: at least
    /// seed_min_spectra spectra at search:prefilter_depth, and at least
    /// seed_decoy_factor times its decoy's. A library-RT bin without such a
    /// target stays empty.
    static constexpr std::uint32_t seed_min_spectra = 3;
    static constexpr std::uint32_t seed_decoy_factor = 3;
    /// A candidate set whose decoy:target count ratio leaves this band is a
    /// broken selection (pairs are kept whole, so it is exactly 1).
    static constexpr double candidate_ratio_low = 0.8;
    static constexpr double candidate_ratio_high = 1.25;

    // ---- ion mobility (diaPASEF), recorded with the calibration --------------
    /// search:im_window 0: the full 1/K0 extraction width is 2 x
    /// im_window_padding x the larger of the 0.99 quantile of the calibration
    /// residuals and im_window_normal_quantile robust SDs of them (the RT
    /// window's rule), clamped to this range. The floor keeps a whole
    /// ion-mobility peak (FWHM about 0.02 on a timsTOF) inside the window
    /// whatever the calibration says; the cap keeps the window below the
    /// width of one diaPASEF isolation window's 1/K0 range.
    static constexpr double im_window_padding = 1.3;
    static constexpr double im_window_normal_quantile = 2.5758;
    static constexpr double im_window_min = 0.04;
    static constexpr double im_window_max = 0.16;
    /// The ion-mobility calibration (library 1/K0 -> run 1/K0, a robust line
    /// through the RT calibration's seeds, each measured at its apex): at
    /// least this many seeds must be measured, with at least
    /// im_seed_min_fragments fragments each co-locating in 1/K0 ...
    static constexpr std::size_t im_calibration_min_points = 20;
    static constexpr std::size_t im_seed_min_fragments = 3;
    /// ... and the line's slope must lie within this band: library and run
    /// 1/K0 differ by a calibration, not by a factor.
    static constexpr double im_slope_low = 0.8;
    static constexpr double im_slope_high = 1.25;
    /// OpenSWATH's 1/K0 of a peak group (reportedMobility, a diagnostic since
    /// the M2 review) is its MS2 value (im_drift), NaN when the MS1 value
    /// (im_ms1_drift) exists and differs by more than this. The REPORT's 1/K0
    /// is re-measured over the windows' whole 1/K0 range instead (mobilityAt),
    /// NaN with fewer than im_seed_min_fragments fragments at its apex.
    static constexpr double im_ms1_agreement = 0.02;
    /// An ion-mobility search needs a library 1/K0 (an IM value or a CCS) for
    /// most targets: each precursor is extracted at its calibrated library
    /// 1/K0, and the 1/K0 calibration measures seeds that have one. Below this
    /// share of the library's targets the search is refused before anything
    /// is searched (search:im_window -1 searches without ion mobility).
    static constexpr double im_library_min_share = 0.5;
    /// The loader reads every spectrum of a diaPASEF run once: a peak whose
    /// 1/K0 lies more than im_limits_tolerance outside its own isolation
    /// window's 1/K0 limits counts against those limits, and above
    /// im_limits_outside_warn of the MS2 peaks the limits are not on the
    /// per-peak values' calibration (8.9 % in files of the converter that
    /// wrote a linear approximation), which window assignment relies on: a
    /// warning, and the share is recorded.
    static constexpr double im_limits_tolerance = 0.001;
    static constexpr double im_limits_outside_warn = 0.01;

    /// Throws std::invalid_argument naming the first bad setting.
    void validate() const;

    /// The settings as a JSON object. With @p execution false it leaves out the
    /// settings that decide how the search runs but must not change what it
    /// finds -- chunk, batch_size, readoptions (the mode actually used is
    /// recorded with the run) and cache_dir (a path of this machine) -- which
    /// is what the report embeds; the provenance sidecar gets everything.
    std::string toJson(bool execution = true) const;
  };

  /// Thrown by a run-level guard. The message always carries the counts.
  struct SearchAbort : std::runtime_error
  {
    using std::runtime_error::runtime_error;
  };

  /// Thrown by a stage that does not exist yet.
  struct NotImplemented : std::logic_error
  {
    using std::logic_error::logic_error;
  };
}
