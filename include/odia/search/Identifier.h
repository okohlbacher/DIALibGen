// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

/// The built-in identification step behind `-run`: from a library and one DIA
/// run to an identification report, which then serves as `-ids`.
///
///   candidates (CandidateSelector)   paired random subset + in-memory decoys
///   loadRun                          SwathFile::loadMzML, IM-limit fix, diaPASEF detection   [seam]
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
#include <odia/search/ReportWriter.h>
#include <odia/search/Scoring.h>
#include <odia/search/SearchParams.h>

#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>
#include <OpenMS/METADATA/ExperimentalSettings.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/SwathMap.h>

#include <cstddef>
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
    std::vector<OpenSwath::SwathMap> maps;
    std::shared_ptr<OpenMS::ExperimentalSettings> meta;
    /// diaPASEF: MS2 windows carry 1/K0 limits AND spectra carry a per-peak 1/K0 array.
    bool ion_mobility = false;
    /// What the loader actually used ("normal" or "cache").
    std::string read_mode;
    /// A directory the loader created for cache files. identify() removes it
    /// (recursively) once extraction is over, also when a later stage throws.
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
    /// 1/K0 correction as performRTNormalization returns it; the assays apply
    /// it inverted to library 1/K0 (AssayOptions::im_map).
    OpenMS::TransformationDescription im;
    double rt_window = 0.0;      ///< full RT extraction width, seconds
    double mz_ppm = 0.0;         ///< full MS2 m/z extraction width, ppm
    double ms1_mz_ppm = 0.0;     ///< full MS1 m/z extraction width, ppm
    double im_window = -1.0;     ///< full 1/K0 extraction width; -1 = no ion mobility
    std::size_t seeds = 0;       ///< seed assays searched
    std::size_t points = 0;      ///< points the RT model was fitted to
    double rsq = std::numeric_limits<double>::quiet_NaN();
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

    /// The report rows of a scored search: every competition winner, target or
    /// decoy, at precursor q <= search:report_max_q, in SearchSet order.
    static std::vector<ReportRow> reportRows(const SearchSet& set, const PeakGroups& groups,
                                             const ScoringOutcome& outcome, const SearchParams& params);

  protected:
    // ---- seams: src/search/RunStages.cpp ------------------------------------
    /// Load @p path (search:readoptions, search:cache_dir): swap reversed 1/K0
    /// window limits, detect diaPASEF, record the scratch directory.
    virtual RunData loadRun(const std::string& path);
    /// Calibrate RT (and 1/K0) on seed assays -- targets of @p set, or kit
    /// peptides found in @p library (which is the unmodified input) -- with the
    /// thresholds search:calibration_min_rsq / _min_coverage; size the
    /// extraction windows (search:rt_window, mz_ppm, im_window; 0 = from the
    /// calibration). A failed calibration throws SearchAbort unless
    /// search:allow_bootstrap.
    virtual Calibration calibrate(const Library& library, const SearchSet& set, RunData& run);
    /// Extract every precursor of @p set in chunks (AssayBuilder::chunks with
    /// search:chunk, search:batch_size) and append one row per reported peak
    /// group to @p out (PeakGroups::add), freeing each chunk's features.
    virtual void extract(const SearchSet& set, RunData& run, const Calibration& calibration, PeakGroups& out);

    void info(const std::string& message) const;
    void warn(const std::string& message) const;

  private:
    SearchParams params_;
    Log info_, warn_;
  };
}
