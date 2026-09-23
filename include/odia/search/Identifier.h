// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

/// The built-in identification step behind `-run`: from a library and one DIA
/// run to an identification report, which then serves as `-ids`.
///
///   loadRun                          SwathFile::loadMzML, IM-limit fix, diaPASEF detection   [seam]
///   candidates                       evidence prefilter (default) or paired random subset inside
///                                    the run's isolation windows + in-memory decoys; both members'
///                                    assays predicted by one model (search:intensities predicted)
///   calibrate                        performRTNormalization on seed assays; windows          [seam]
///   extract                          chunked performExtraction -> PeakGroups                 [seam]
///   checkExtraction                  decoy:target band
///   scorePeakGroups, checkGuards     odia-core LDA + competition + roll-ups; run-level guards
///   ReportWriter                     -out_ids
///
/// The three seams are virtual so tests can replace the run with a synthetic
/// one; their product implementations live in src/search/RunStages.cpp. Every
/// guard aborts with counts (SearchAbort); nothing falls back silently.
#pragma once

#include <odia/Library.h>
#include <odia/search/CandidateSelector.h>
#include <odia/search/PredictedAssays.h>
#include <odia/search/ReportWriter.h>
#include <odia/search/Scoring.h>
#include <odia/search/SearchParams.h>

#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>
#include <OpenMS/METADATA/ExperimentalSettings.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/SwathMap.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <functional>
#include <limits>
#include <memory>
#include <string>
#include <vector>

namespace ODIA::search
{
  /// A loaded run: what loadRun() returns and extraction reads.
  struct RunData
  {
    std::string path;
    /// The report's Run value. Identifier fills in the file stem when empty.
    std::string name;
    /// One map per isolation window, MS1 maps flagged ms1 (SwathFile::loadMzML).
    /// Their lower/upper bounds decide which candidates can be searched.
    std::vector<OpenSwath::SwathMap> maps;
    /// SwathFile's metadata; identify() releases it as soon as loadRun returns.
    std::shared_ptr<OpenMS::ExperimentalSettings> meta;
    /// diaPASEF: MS2 windows carry 1/K0 limits AND spectra carry a per-peak 1/K0 array.
    bool ion_mobility = false;
    /// ... and the MS1 spectra carry a per-peak 1/K0 array too (MS1 traces
    /// and MS1 scores are then read within the precursor's 1/K0 range).
    bool ms1_ion_mobility = false;
    /// What the loader actually used ("normal" or "cache").
    std::string read_mode;
    /// Bytes the cache files took after loading (0 in memory).
    std::uintmax_t cache_bytes = 0;
    /// A directory the loader created for cache files. identify() removes it
    /// (recursively) once extraction is over, also when a later stage throws,
    /// and always after the maps (which hold its files open) are released.
    std::filesystem::path scratch;
    /// Anything else the loader wants in the provenance, as a JSON object; may be empty.
    std::string provenance_json;
  };

  /// A calibrated run: what calibrate() returns and extraction reads.
  struct Calibration
  {
    /// Run RT (seconds) -> assay RT (SearchSet::rt_scale, [0, 100]): the
    /// orientation OpenSwathWorkflow::performExtraction takes.
    OpenMS::TransformationDescription rt;
    /// Library 1/K0 -> run 1/K0, a line fitted to the seeds' measured 1/K0
    /// (IonMobility.h); the assays map every library 1/K0 through it
    /// (AssayOptions::im_map). Empty (no data points) without ion mobility.
    OpenMS::TransformationDescription im;
    double rt_window = 0.0;      ///< full RT extraction width, seconds
    double mz_ppm = 0.0;         ///< full MS2 m/z extraction width, ppm
    double ms1_mz_ppm = 0.0;     ///< full MS1 m/z extraction width, ppm
    double im_window = -1.0;     ///< full 1/K0 extraction width; -1 = no ion mobility (none, or search:im_window -1)
    std::size_t seeds = 0;       ///< seed assays searched
    std::size_t points = 0;      ///< points the RT model was fitted to
    double rsq = std::numeric_limits<double>::quiet_NaN();
    /// points / seeds, recorded as points_per_seed. Not what
    /// search:calibration_min_coverage bounds: that is the share of the FOUND
    /// seeds the stock outlier removal must keep.
    double coverage = std::numeric_limits<double>::quiet_NaN();
    bool bootstrap = false;      ///< search:allow_bootstrap was USED (calibration had failed)
    /// Anything else for the provenance, as a JSON object; may be empty.
    std::string provenance_json;
  };

  struct IdentificationResult
  {
    std::string report;                 ///< the -out_ids path written
    std::size_t report_rows = 0;
    std::size_t report_targets = 0;
    std::size_t report_decoys = 0;
    std::size_t identified = 0;         ///< target precursors at q <= 0.01
    std::string run_name;               ///< the report's Run value
    /// The provenance sidecar's "search" block, as a JSON object.
    std::string provenance_json;
  };

  class Identifier
  {
  public:
    using Log = std::function<void(const std::string&)>;

    explicit Identifier(SearchParams params, Log info = Log(), Log warn = Log());
    virtual ~Identifier();
    Identifier(const Identifier&) = delete;
    Identifier& operator=(const Identifier&) = delete;

    /// Search @p run_path for @p library's precursors and write the report to
    /// @p out_ids (refused if it exists). @p library is not modified. Throws
    /// SearchAbort on a guard, NotImplemented from a missing stage, and
    /// std::exception otherwise; the report is written only on success.
    IdentificationResult identify(const Library& library, const std::string& run_path,
                                  const std::string& out_ids, const std::string& tool_version);

    const SearchParams& params() const { return params_; }

    /// The fragment model search:intensities predicted uses instead of the
    /// PeptDeep model named by search's ms2_model (a test seam: a synthetic
    /// run needs a model that predicts what it planted).
    void setFragmentModel(std::shared_ptr<FragmentModel> model) { fragment_model_ = std::move(model); }

    /// The report rows of a scored search: every competition winner, target or
    /// decoy, at precursor q <= search:report_max_q, in SearchSet order.
    static std::vector<ReportRow> reportRows(const SearchSet& set, const PeakGroups& groups,
                                             const ScoringOutcome& outcome, const SearchParams& params);

    /// The MS2 isolation windows of a loaded run (its non-MS1 maps), by lower bound.
    static std::vector<IsolationWindow> isolationWindows(const RunData& run);

    /// The RT calibration's second model: LOWESS (stock, span by
    /// cross-validation, its own iterative outlier cut on @p stock) replaces
    /// @p line, fitted to @p points (the line's inliers), when at least 200
    /// points remain, it increases over the run and the assay range, and its
    /// 5-fold cross-validated error on ONE point set (the union of both
    /// inlier sets, residuals capped at 3 line scales) is below 0.97x the
    /// line's. Points are (run seconds, assay RT). record_json says why.
    struct RtModelChoice
    {
      bool lowess = false;
      OpenMS::TransformationDescription fit;           ///< the LOWESS fit, when used
      std::vector<std::pair<double, double>> points;   ///< its inliers, when used
      std::string record_json;
    };
    static RtModelChoice chooseRtModel(const std::vector<std::pair<double, double>>& stock,
                                       const std::vector<std::pair<double, double>>& points,
                                       const OpenMS::TransformationDescription& line, double run_first, double run_span);

  protected:
    // ---- seams: src/search/RunStages.cpp ------------------------------------
    /// Load @p path (search:readoptions, search:cache_dir): swap reversed 1/K0
    /// window limits, detect diaPASEF, record the scratch directory. Auto
    /// caches runs above 3 GB; the cache directory is created under
    /// search:cache_dir, else under scratchParent(). On a diaPASEF run it
    /// checks that sampled MS2 spectra carry a 1/K0 array (unless
    /// search:im_window -1) and whether the MS1 spectra do.
    virtual RunData loadRun(const std::string& path);
    /// Calibrate RT (and 1/K0) on seed assays -- targets of @p set, or kit
    /// peptides found in @p library (which is the unmodified input) -- with the
    /// thresholds search:calibration_min_rsq / _min_coverage; size the
    /// extraction windows (search:rt_window, mz_ppm, im_window; 0 = from the
    /// calibration). A failed calibration throws SearchAbort unless
    /// search:allow_bootstrap. On a diaPASEF run (unless search:im_window -1)
    /// each seed of the RT calibration is measured in 1/K0 at its apex
    /// (mobilityApex) and a robust line maps library to run 1/K0; its failure
    /// throws SearchAbort, whatever search:allow_bootstrap says.
    virtual Calibration calibrate(const Library& library, const SearchSet& set, RunData& run);
    /// Extract every precursor of @p set in chunks (AssayBuilder::chunks with
    /// search:chunk, search:batch_size) and append one row per reported peak
    /// group to @p out (PeakGroups::add), freeing each chunk's features.
    virtual void extract(const SearchSet& set, RunData& run, const Calibration& calibration, PeakGroups& out);

    void info(const std::string& message) const;
    void warn(const std::string& message) const;

    /// calibrate()'s ion-mobility step: measure the seeds behind the RT
    /// model's @p points (run seconds, assay RT) in 1/K0 at their apex, fit
    /// library -> run 1/K0 robustly, validate, and size the window; sets
    /// @p cal's im and im_window and returns the provenance record (JSON).
    /// @p seed holds the seed assays, @p picked their input library indices
    /// in the same order. Throws SearchAbort when the calibration fails.
    std::string calibrateMobility(const Library& library, const SearchSet& set, const RunData& run,
                                  const OpenSwath::LightTargetedExperiment& seed, const std::vector<std::size_t>& picked,
                                  const std::vector<std::pair<double, double>>& points, Calibration& cal);

    /// Where a cache-mode load puts its scratch directory when search:cache_dir
    /// is empty: the directory of -out_ids (set by identify()), which is where
    /// the user already expects large files; the system temporary directory
    /// when the loader runs outside identify().
    std::filesystem::path scratchParent() const;

  private:
    SearchParams params_;
    Log info_, warn_;
    std::filesystem::path output_dir_;
    std::shared_ptr<FragmentModel> fragment_model_;
  };
}
