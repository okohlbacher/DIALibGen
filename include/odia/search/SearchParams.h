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

  /// The decoy methods the search accepts: shuffle and pseudo_reverse, both
  /// of which keep the target's termini. `mutate` is refused (its substitution
  /// table is DIA-NN's), `reverse` (it moves the C-terminal K/R, which makes
  /// decoys separable from tryptic targets) and `none`.
  DecoyMethod parseSearchDecoyMethod(const std::string& s);   ///< throws std::invalid_argument
  const char* searchDecoyName(DecoyMethod m);

  struct SearchParams
  {
    // ---- candidates ----------------------------------------------------------
    /// Candidate selection. Only "random" exists in this version: a
    /// deterministic, paired, label-blind random subset of the library's targets.
    std::string candidates = "random";
    /// Targets drawn by the random selection; 0 = every eligible target.
    std::size_t subset = 100000;
    /// Cap on target-decoy pairs after selection; 0 = no cap. Search time is
    /// linear in pairs.
    std::size_t max_pairs = 40000;
    /// How the in-memory search decoys are built (decoys already in the library
    /// file are ignored by the search and left untouched).
    DecoyMethod decoys = DecoyMethod::Shuffle;
    /// Salt of the candidate draw. Changes WHICH pairs are searched, never how.
    std::uint64_t seed = 42;
    /// Extraction passes. 1 in this version.
    int passes = 1;

    // ---- extraction ----------------------------------------------------------
    /// Full width of the RT extraction window in seconds; 0 = from the calibration.
    double rt_window = 0.0;
    /// Full width of the fragment m/z extraction window in ppm; 0 = automatic.
    double mz_ppm = 0.0;
    /// Full width of the 1/K0 extraction window; 0 = automatic, -1 = off.
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
    /// that leaks labels; decoys built weaker than null targets they cannot
    /// see -- that is what the null-pair balance diagnostic is for.
    bool selftest = true;

    /// Threads for extraction and the classifier's folds (TOPP -threads).
    int threads = 1;

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
    /// The null-pair balance looks at the lowest this share of pair winners ...
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
