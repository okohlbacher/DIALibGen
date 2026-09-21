// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

// The run-facing stages of Identifier: loading, calibration and extraction on
// stock OpenMS 3.5.0 (SwathFile::loadMzML, OpenSwathCalibrationWorkflow::
// performRTNormalization, OpenSwathWorkflow::performExtraction). Their
// contracts are documented in include/odia/search/Identifier.h.
//
// Everything here is the stock API as the bioconda OpenMS 3.5.0 ships it. The
// traps it works around were measured on that build (a throwaway probe on an
// Astral slice, a full Astral run and a diaPASEF slice); each is named where
// it is handled:
//   * OpenSwathHelper::sampleExperiment de-duplicates by compound.sequence, so
//     it is fed the stripped sequence while sampling and the seeds are rebuilt
//     sequence-blind afterwards; its kit-priority path is not used -- every
//     kit precursor in the library is a seed anyway;
//   * performRTNormalization "succeeds" on a handful of points with a
//     nonsense slope, and its outlier removal stops as soon as r^2 reaches
//     min_rsq, so its points are refit robustly (RobustLine.h) and the result
//     is validated here;
//   * MRMRTNormalizer is not a DefaultParamHandler: every key it reads must be
//     present (outlierMethod included);
//   * the m/z window estimate exists only when an m/z correction function is
//     set, and that correction wraps every map in a transforming accessor
//     (+70 % extraction CPU measured) -- so the calibration runs on a COPY of
//     the maps and extraction reads the originals;
//   * Scores:use_elution_model_score defaults to true and Param::setValue on a
//     misspelled key only warns, so every key is set explicitly and checked;
//   * MS1 isotope traces above the monoisotopic one are placed at
//     +k * 1.00336 Th whatever the charge, so only the monoisotopic trace is
//     extracted (the MS1 isotope scores come from full scans regardless);
//   * features are appended to the FeatureMap in thread-completion order, so
//     each precursor's peak groups are ranked by apex RT to get an id that
//     does not depend on scheduling or chunking.

#include <odia/search/AssayBuilder.h>
#include <odia/search/Identifier.h>
#include <odia/search/RobustLine.h>

#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFinderScoring.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathHelper.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathOSWWriter.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathWorkflow.h>
#include <OpenMS/ANALYSIS/OPENSWATH/SwathMapMassCorrection.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataWritingConsumer.h>
#include <OpenMS/FORMAT/SwathFile.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>
#include <OpenMS/METADATA/MetaInfoRegistry.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/ISpectrumAccess.h>
#include <OpenMS/SYSTEM/File.h>

#include <nlohmann/json.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <streambuf>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace ODIA::search
{
  namespace
  {
    namespace fs = std::filesystem;
    using json = nlohmann::json;
    using OpenMS::Param;

    // ---- constants of the method (recorded in the provenance) -----------------

    /// search:readoptions auto holds a run in memory up to this size and caches
    /// it per window on disk above it (bytes of the mzML file).
    constexpr std::uintmax_t cache_above_bytes = 3000000000ull;
    /// Fragment and precursor m/z extraction width (full, ppm) when nothing
    /// narrower is known; the calibration's estimate may only narrow it.
    constexpr double default_mz_ppm = 20.0;
    constexpr double mz_floor_ppm = 4.0;
    constexpr double mz_cap_ppm = 50.0;
    /// m/z width (full, ppm) the seed assays are searched with, over the whole run.
    constexpr double calibration_mz_ppm = 20.0;
    /// Seed sampling over the search's targets, after every kit precursor:
    /// 100 RT bins from the most intense 40 %, as OpenSwathWorkflow's auto_irt
    /// linear calibration, but 20 per bin instead of 5 -- a predicted library
    /// has far fewer detectable precursors than the empirical libraries that
    /// recipe was made for (an Astral slice: 0 accepted points from 496 seeds).
    constexpr std::size_t seed_bins = 100;
    constexpr std::size_t seeds_per_bin = 20;
    constexpr double seed_top_fraction = 0.4;
    constexpr unsigned seed_sampling_seed = 5489;
    /// A calibration is accepted only with at least this many points ...
    constexpr std::size_t calibration_min_points = 20;
    /// ... covering at least this fraction of the seeds' library-RT span ...
    constexpr double calibration_min_span = 0.5;
    /// ... and a slope that maps the library RT range onto a run span within
    /// this factor of the run's actual length.
    constexpr double calibration_span_factor = 10.0;
    /// The robust refit of the stock points (RobustLine.h): inliers within
    /// this many robust SDs, the SD never taken below this (library units on
    /// the [0, 100] assay scale).
    constexpr double robust_cut = 3.0;
    constexpr double robust_min_scale = 0.1;
    /// RT window: quantile and padding of the seed residuals (TOPP's recipe),
    /// never narrower than this (seconds).
    constexpr double rt_window_quantile = 0.99;
    constexpr double rt_window_padding = 1.3;
    constexpr double rt_window_floor_s = 30.0;
    /// Monoisotopic MS1 trace only: stock 3.5.0 places the isotope traces at
    /// +k * 1.00336 Th regardless of charge (ChromatogramExtractor::prepare_coordinates).
    constexpr int ms1_isotope_traces = 0;

    /// A standard stream silenced while a stock call runs: OpenSWATH prints
    /// one line per window and batch to std::cout ("Thread 3_0 will analyze
    /// ..."), and extraction one "Detected N empty chromatograms" warning per
    /// batch to std::cerr -- thousands of lines on an Astral run. Scoped and
    /// restored on exceptions; exceptions themselves are not affected.
    class Quiet
    {
    public:
      explicit Quiet(std::ostream& stream) : stream_(stream), previous_(stream.rdbuf(&null_)) {}
      ~Quiet() { stream_.rdbuf(previous_); }
      Quiet(const Quiet&) = delete;
      Quiet& operator=(const Quiet&) = delete;

    private:
      struct Null : std::streambuf
      {
        int overflow(int c) override { return traits_type::not_eof(c); }
        std::streamsize xsputn(const char*, std::streamsize n) override { return n; }
      } null_;
      std::ostream& stream_;
      std::streambuf* previous_;
    };

    /// Both standard streams: stock extraction and calibration print per-window
    /// progress to std::cout. Only std::cout and std::cerr are silenced: lines
    /// that stock OpenMS writes through its own log streams ("Detected N empty
    /// chromatograms", "N spectra and 0 chromatograms stored", "Read
    /// chromatogram while reading SWATH files") or with printf (the m/z
    /// regression parameters) still reach the console, and the printf cannot
    /// be silenced portably.
    struct QuietStreams
    {
      Quiet out{std::cout};
      Quiet err{std::cerr};
    };

    /// The OpenMP team for one stock OpenSWATH call, restored afterwards. Stock
    /// 3.5.0 scores every feature through MetaInfo::setValue/getValue, which
    /// take a process-wide critical section (MetaInfoRegistry); above about
    /// eight threads the team mostly spins on it. Measured on an Astral slice
    /// (40,000 pairs): extraction 108 s at 4 threads, 62-64 s at 8, 57 s at 16
    /// and 94-106 s at 32, with 3x the CPU of 16. The result does not depend on it.
    class ThreadCap
    {
    public:
      explicit ThreadCap(int cap)
      {
#ifdef _OPENMP
        previous_ = omp_get_max_threads();
        used_ = std::max(1, std::min(previous_, cap));
        omp_set_num_threads(used_);
#else
        (void)cap;
#endif
      }
      ~ThreadCap()
      {
#ifdef _OPENMP
        omp_set_num_threads(previous_);
#endif
      }
      ThreadCap(const ThreadCap&) = delete;
      ThreadCap& operator=(const ThreadCap&) = delete;
      int used() const { return used_; }
      int requested() const { return previous_; }

    private:
      int previous_ = 1;
      int used_ = 1;
    };

    std::string fixed(double v, int digits)
    {
      std::ostringstream s;
      s.setf(std::ios::fixed);
      s.precision(digits);
      s << v;
      return s.str();
    }

    json num(double v) { return std::isfinite(v) ? json(v) : json(nullptr); }

    /// Set @p key, which must already exist: Param::setValue on a misspelled
    /// key adds it with a warning and changes nothing the algorithm reads.
    void set(Param& p, const std::string& key, const OpenMS::ParamValue& value)
    {
      if (!p.exists(key))
      { throw std::logic_error("search: OpenMS parameter '" + key + "' does not exist in this OpenMS build (stock 3.5.0 is required)"); }
      p.setValue(key, value);
    }

    /// OpenSWATH's peak picking and scoring: the TOPP tool's "Scoring" block,
    /// every key set explicitly. Elution-model and ion-series scores are off
    /// (the former defaults to ON; the latter needs sequences, which the
    /// assays do not carry).
    Param featureFinderParam(bool ms1, bool ion_mobility)
    {
      Param ff = OpenMS::MRMFeatureFinderScoring().getDefaults();
      set(ff, "stop_report_after_feature", 5);
      set(ff, "rt_normalization_factor", 100.0);   // assay RT lives on [0, 100]
      set(ff, "Scores:use_ms1_correlation", ms1 ? "true" : "false");
      set(ff, "Scores:use_ms1_fullscan", ms1 ? "true" : "false");
      set(ff, "Scores:use_ms1_mi", ms1 ? "true" : "false");
      set(ff, "Scores:use_mi_score", "true");
      set(ff, "Scores:use_elution_model_score", "false");
      set(ff, "Scores:use_ionseries_scores", "false");
      set(ff, "Scores:use_ion_mobility_scores", ion_mobility ? "true" : "false");
      set(ff, "use_ms1_ion_mobility", ion_mobility ? "true" : "false");
      set(ff, "TransitionGroupPicker:min_peak_width", -1.0);
      set(ff, "TransitionGroupPicker:recalculate_peaks", "true");
      set(ff, "TransitionGroupPicker:compute_peak_quality", "false");
      set(ff, "TransitionGroupPicker:minimal_quality", -1.5);
      set(ff, "TransitionGroupPicker:background_subtraction", "none");
      set(ff, "TransitionGroupPicker:compute_peak_shape_metrics", "false");
      set(ff, "TransitionGroupPicker:recalculate_peaks_max_z", 0.75);
      set(ff, "TransitionGroupPicker:PeakPickerChromatogram:use_gauss", "false");
      set(ff, "TransitionGroupPicker:PeakPickerChromatogram:sgolay_polynomial_order", 3);
      set(ff, "TransitionGroupPicker:PeakPickerChromatogram:sgolay_frame_length", 11);
      set(ff, "TransitionGroupPicker:PeakPickerChromatogram:peak_width", -1.0);
      set(ff, "TransitionGroupPicker:PeakPickerChromatogram:remove_overlapping_peaks", "true");
      set(ff, "TransitionGroupPicker:PeakPickerChromatogram:write_sn_log_messages", "false");
      set(ff, "TransitionGroupPicker:PeakPickerChromatogram:method", "corrected");
      set(ff, "TransitionGroupPicker:PeakPickerChromatogram:signal_to_noise", 0.1);
      set(ff, "TransitionGroupPicker:PeakPickerChromatogram:gauss_width", 30.0);
      set(ff, "uis_threshold_sn", -1);
      set(ff, "uis_threshold_peak_area", 0);
      set(ff, "EMGScoring:max_iteration", 10);
      return ff;
    }

    /// MRMRTNormalizer's settings ("Calibration:RTNormalization"). Not a
    /// DefaultParamHandler: every key it reads must be present. Linear model,
    /// iterative residual outlier removal, and estimateBestPeptides ON -- the
    /// TOPP default (off) stops outlier removal at r^2 ~ min_rsq by
    /// construction; on, the probe measured r^2 0.977 and residuals near an
    /// oracle linear fit.
    Param irtDetectionParam()
    {
      Param p;
      p.setValue("alignmentMethod", "linear");
      p.setValue("lowess:auto_span", "true");
      p.setValue("lowess:span", 0.05);
      p.setValue("lowess:auto_span_min", 0.15);
      p.setValue("lowess:auto_span_max", 0.80);
      p.setValue("lowess:auto_span_grid", "0.005,0.01,0.05,0.15,0.25,0.30,0.50,0.70,0.90");
      p.setValue("b_spline:num_nodes", 5);
      p.setValue("outlierMethod", "iter_residual");
      p.setValue("useIterativeChauvenet", "false");
      p.setValue("RANSACMaxIterations", 1000);
      p.setValue("RANSACMaxPercentRTThreshold", 3);
      p.setValue("RANSACSamplingSize", 10);
      p.setValue("estimateBestPeptides", "true");
      p.setValue("InitialQualityCutoff", 0.5);
      p.setValue("OverallQualityCutoff", 5.5);
      p.setValue("NrRTBins", 10);
      p.setValue("MinPeptidesPerBin", 1);
      p.setValue("MinBinsFilled", 8);
      return p;
    }

    /// SwathMapMassCorrection for the calibration. The correction function is
    /// set only so the m/z windows get ESTIMATED (with "none" the estimator
    /// returns before measuring anything); the corrected maps it produces are
    /// a copy and are discarded.
    Param massCorrectionParam()
    {
      Param p = OpenMS::SwathMapMassCorrection().getDefaults();
      set(p, "mz_extraction_window", calibration_mz_ppm);
      set(p, "mz_extraction_window_ppm", "true");
      set(p, "im_extraction_window", -1.0);
      set(p, "mz_correction_function", "regression_delta_ppm");
      set(p, "im_correction_function", "none");
      return p;
    }

    /// ChromExtractParams is a POD without a constructor: every field is set.
    OpenMS::ChromExtractParams chromParams(double mz_ppm, double im_window, double rt_window)
    {
      OpenMS::ChromExtractParams cp;
      cp.min_upper_edge_dist = 0.0;
      cp.mz_extraction_window = mz_ppm;
      cp.ppm = true;
      cp.im_extraction_window = im_window;
      cp.extraction_function = "tophat";
      cp.rt_extraction_window = rt_window;
      cp.extra_rt_extract = 0.0;
      return cp;
    }

    /// The sub-scores handed to the classifier, in column order. A fixed list:
    /// columns a run cannot fill stay NaN and are dropped label-blind by the
    /// classifier. Deliberately a non-collinear subset of OpenSWATH's var_*
    /// scores (the library_* and xcorr_* families each carry near-duplicates),
    /// plus the MS1 scores with search:ms1. The RT deviation score is
    /// var_norm_rt_score, the one search:rt_im_scores false removes.
    ///
    /// MS1-MS2 co-elution: stock 3.5.0 computes var_ms1_xcorr_coelution,
    /// var_ms1_xcorr_shape and var_ms1_mi_score only from TWO or more precursor
    /// isotope traces (OpenSwathScoring, `precursor_ids.size() > 1`), and this
    /// search extracts the monoisotopic trace only (ms1_isotope_traces), so
    /// those three were NaN on every run. The *_contrast variants correlate
    /// that one MS1 trace with the fragment traces and are always computed.
    std::vector<std::string> scoreColumns(bool ms1)
    {
      std::vector<std::string> c = {
        "var_xcorr_coelution_weighted", "var_xcorr_shape_weighted", "var_xcorr_coelution", "var_xcorr_shape",
        "var_library_corr", "var_library_rmsd", "var_library_sangle", "var_library_manhattan",
        "var_intensity_score", "var_log_sn_score", "var_massdev_score_weighted",
        "var_isotope_correlation_score", "var_isotope_overlap_score", "var_mi_weighted_score", "var_norm_rt_score"};
      if (ms1)
      {
        for (const char* s : {"var_ms1_ppm_diff", "var_ms1_isotope_correlation", "var_ms1_isotope_overlap",
                              "var_ms1_xcorr_coelution_contrast", "var_ms1_xcorr_shape_contrast", "var_ms1_mi_contrast_score"})
        { c.emplace_back(s); }
      }
      return c;
    }

    /// Uppercase residues outside (...) and [...]: "AC(UniMod:4)K" -> "ACK".
    std::string strippedSequence(std::string_view modified)
    {
      std::string out;
      int depth = 0;
      for (const char c : modified)
      {
        if (c == '(' || c == '[') { ++depth; continue; }
        if (c == ')' || c == ']') { depth = std::max(0, depth - 1); continue; }
        if (depth == 0 && c >= 'A' && c <= 'Z') { out.push_back(c); }
      }
      return out;
    }

    /// The iRT and CiRT kit sequences OpenMS ships (share/OpenMS/CHEMISTRY),
    /// read from their PeptideSequence column. Missing files are reported, not fatal.
    std::unordered_set<std::string> kitSequences(std::vector<std::string>& found, std::vector<std::string>& missing)
    {
      std::unordered_set<std::string> out;
      for (const char* name : {"CHEMISTRY/cirtkit.tsv", "CHEMISTRY/irtkit.tsv"})
      {
        std::string path;
        try { path = OpenMS::File::find(name); }
        catch (const std::exception&) { missing.emplace_back(name); continue; }
        std::ifstream in(path);
        std::string line;
        if (!in || !std::getline(in, line)) { missing.emplace_back(name); continue; }
        auto split = [](const std::string& s) {
          std::vector<std::string> f;
          std::size_t start = 0;
          for (;;)
          {
            const std::size_t tab = s.find('\t', start);
            f.push_back(s.substr(start, tab == std::string::npos ? std::string::npos : tab - start));
            if (tab == std::string::npos) { break; }
            start = tab + 1;
          }
          return f;
        };
        const auto header = split(line);
        const auto col = std::find(header.begin(), header.end(), "PeptideSequence");
        if (col == header.end()) { missing.emplace_back(std::string(name) + " (no PeptideSequence column)"); continue; }
        const auto k = static_cast<std::size_t>(col - header.begin());
        while (std::getline(in, line))
        {
          if (!line.empty() && line.back() == '\r') { line.pop_back(); }
          const auto f = split(line);
          if (k < f.size() && !f[k].empty()) { out.insert(strippedSequence(f[k])); }
        }
        found.emplace_back(name);
      }
      return out;
    }

    /// First and last MS2 spectrum time over all isolation windows, seconds.
    std::pair<double, double> runSpan(const std::vector<OpenSwath::SwathMap>& maps)
    {
      double lo = std::numeric_limits<double>::infinity(), hi = -lo;
      for (const auto& m : maps)
      {
        if (m.ms1 || !m.sptr) { continue; }
        const std::size_t n = m.sptr->getNrSpectra();
        if (n == 0) { continue; }
        lo = std::min(lo, m.sptr->getSpectrumMetaById(0).RT);
        hi = std::max(hi, m.sptr->getSpectrumMetaById(static_cast<int>(n - 1)).RT);
      }
      return {lo, hi};
    }

    double quantile(std::vector<double> v, double q)
    {
      if (v.empty()) { return std::numeric_limits<double>::quiet_NaN(); }
      std::sort(v.begin(), v.end());
      const double pos = q * static_cast<double>(v.size() - 1);
      const auto a = static_cast<std::size_t>(std::floor(pos));
      const std::size_t b = std::min(v.size() - 1, a + 1);
      return v[a] + (pos - static_cast<double>(a)) * (v[b] - v[a]);
    }

    /// Everything the parameter setters above would refuse, checked before
    /// the run is read (which can take minutes): a missing key is a build
    /// problem, not a run problem.
    void checkParameters()
    {
      (void)featureFinderParam(true, true);
      (void)featureFinderParam(false, false);
      (void)massCorrectionParam();
    }

    /// A strict total order on score rows (NaN after every number), so equal
    /// apex times still rank the same way on every run.
    bool rowLess(double apex_a, const float* a, double apex_b, const float* b, std::size_t width)
    {
      if (apex_a != apex_b) { return apex_a < apex_b; }
      for (std::size_t c = 0; c < width; ++c)
      {
        const bool na = std::isnan(a[c]), nb = std::isnan(b[c]);
        if (na != nb) { return nb; }
        if (!na && a[c] != b[c]) { return a[c] < b[c]; }
      }
      return false;
    }
  }

  // ---- loading ----------------------------------------------------------------

  RunData Identifier::loadRun(const std::string& path)
  {
    checkParameters();
    RunData run;
    run.path = path;
    std::error_code ec;
    if (!fs::is_regular_file(path, ec))
    { throw std::runtime_error("search: the run " + path + " is not a readable file"); }
    const std::uintmax_t bytes = fs::file_size(path);
    const bool cache = params().readoptions == ReadMode::Cache ||
                       (params().readoptions == ReadMode::Auto && bytes > cache_above_bytes);
    run.read_mode = cache ? "cache" : "normal";

    // Cache files go to a private directory, which identify() removes once
    // extraction is over. OpenMS never deletes its cache files itself.
    std::string tmp;
    if (cache)
    {
      const fs::path parent = params().cache_dir.empty() ? scratchParent() : fs::path(params().cache_dir);
      fs::create_directories(parent);
      std::random_device random;
      for (int attempt = 0; attempt < 16 && run.scratch.empty(); ++attempt)
      {
        const fs::path candidate = parent / (".dialibgen-run-cache-" + std::to_string(random()) + "-" + std::to_string(random()));
        if (fs::create_directory(candidate)) { run.scratch = candidate; }
      }
      if (run.scratch.empty()) { throw std::runtime_error("search: cannot create a cache directory in " + parent.string()); }
      // SwathFile treats tmp as a directory only when it ends in '/', on every platform.
      tmp = run.scratch.generic_string() + "/";
    }
    // Said before the read, which can take minutes: where the cache goes, so
    // a killed process leaves a directory the user can find.
    info("search run: reading " + path + " (" + std::to_string(bytes / 1000000) + " MB) " +
         (cache ? "into per-window cache files in " + run.scratch.string() + " (about the size of the run; removed after "
                  "extraction)"
                : std::string("into memory")));
    try
    {
      OpenMS::SwathFile file;
      file.setLogType(OpenMS::ProgressLogger::NONE);
      QuietStreams quiet;
      run.maps = file.loadMzML(path, tmp, run.meta, run.read_mode);
    }
    catch (...)
    {
      if (!run.scratch.empty()) { fs::remove_all(run.scratch, ec); }
      throw;
    }

    try
    {
      std::size_t ms1 = 0, ms2 = 0, spectra_ms2 = 0, swapped = 0, with_limits = 0;
      for (auto& m : run.maps)
      {
        if (m.ms1) { ++ms1; continue; }
        ++ms2;
        spectra_ms2 += m.sptr ? m.sptr->getNrSpectra() : 0;
        // Reversed 1/K0 limits (lower > upper), as some converters write them.
        if (m.imLower >= 0 && m.imUpper >= 0)
        {
          ++with_limits;
          if (m.imLower > m.imUpper) { std::swap(m.imLower, m.imUpper); ++swapped; }
        }
      }
      if (ms2 == 0 || spectra_ms2 == 0)
      {
        throw SearchAbort("search: the run " + path + " has no DIA isolation windows with spectra (" + std::to_string(ms2) +
                          " MS2 windows, " + std::to_string(ms1) + " MS1 maps); DIA mzML is required");
      }
      // Ion mobility needs both: window limits AND a per-peak 1/K0 array. A
      // file with limits but no array is usable for RT only.
      bool array = false;
      for (const auto& m : run.maps)
      {
        if (m.ms1 || !m.sptr || m.sptr->getNrSpectra() == 0) { continue; }
        if (m.sptr->getSpectrumById(0)->getDriftTimeArray()) { array = true; break; }
      }
      run.ion_mobility = with_limits > 0 && array;
      const auto [first, last] = runSpan(run.maps);
      if (swapped > 0)
      { info("search run: " + std::to_string(swapped) + " of " + std::to_string(ms2) + " isolation windows had reversed 1/K0 limits (swapped)"); }
      if (run.ion_mobility)
      {
        warn("search: the run has ion mobility (diaPASEF); this version extracts it WITHOUT an ion-mobility window, "
             "so every window is searched by m/z and RT only, and IM in the report is empty (NaN)");
      }
      else if (with_limits > 0)
      { warn("search: the isolation windows carry 1/K0 limits but the spectra no per-peak 1/K0 array; the run is searched by m/z and RT only"); }
      const json loader = {
        {"bytes", bytes}, {"read_mode", run.read_mode}, {"auto_cache_above_bytes", cache_above_bytes},
        {"ms2_windows", ms2}, {"ms1_maps", ms1}, {"ms2_spectra", spectra_ms2},
        {"rt_range_s", {num(first), num(last)}},
        {"im_limits", with_limits}, {"im_limits_swapped", swapped}, {"im_per_peak_array", array},
        {"im_extraction", false}};
      run.provenance_json = loader.dump();
    }
    catch (...)
    {
      run.maps.clear();
      if (!run.scratch.empty()) { fs::remove_all(run.scratch, ec); }
      throw;
    }
    return run;
  }

  // ---- calibration ------------------------------------------------------------

  Calibration Identifier::calibrate(const Library& library, const SearchSet& set, RunData& run)
  {
    Calibration cal;
    json detail = json::object();
    const auto [run_first, run_last] = runSpan(run.maps);
    const double run_span = run_last - run_first;
    if (!(run_span > 0))
    { throw SearchAbort("search: the run's MS2 spectra span no time (" + fixed(run_first, 1) + " .. " + fixed(run_last, 1) + " s); cannot calibrate RT"); }

    // 1. Seeds: EVERY kit precursor the library holds (the iRT and CiRT kits
    //    exist to be found in any run), and a stock-sampled share of the
    //    search's other targets. Input-library indices, targets only, one per
    //    precursor id.
    std::vector<std::string> kit_files, kit_missing;
    const std::unordered_set<std::string> kit = kitSequences(kit_files, kit_missing);
    const auto& pre = library.precursors();
    auto eligible = [&](std::size_t i) {
      return !pre.decoy[i] && pre.mz[i] != MZ_INVALID && pre.charge[i] != 0 && std::isfinite(pre.irt[i]) &&
             pre.transition_count[i] >= SearchParams::min_assay_fragments;
    };
    auto idOf = [&](std::size_t i) {
      std::string id(library.strings().get(pre.modified_sequence[i]));
      id += std::to_string(pre.charge[i]);
      return id;
    };
    std::vector<std::size_t> picked;
    std::unordered_set<std::string> picked_ids;
    std::size_t kit_in_library = 0, kit_outside_windows = 0;
    const std::vector<IsolationWindow> windows = isolationWindows(run);
    if (!kit.empty())
    {
      for (std::size_t i = 0; i < library.precursorCount(); ++i)
      {
        if (!eligible(i) || !kit.count(strippedSequence(library.strings().get(pre.modified_sequence[i])))) { continue; }
        ++kit_in_library;
        // A seed outside every isolation window cannot be found.
        if (!windows.empty() && !CandidateSelector::inWindow(fromFixed(pre.mz[i]), windows)) { ++kit_outside_windows; continue; }
        if (picked_ids.insert(idOf(i)).second) { picked.push_back(i); }   // a duplicated (sequence, charge): first index
      }
    }
    const std::size_t picked_kit = picked.size();

    // 2. The stock sampler over the other targets: RT bins, the most intense
    //    share first. It de-duplicates by compound.sequence, so it sees the
    //    STRIPPED sequence (one seed per peptide, whatever its charge states);
    //    one transition per compound carries the summed library intensity it
    //    ranks by. The seed assays are rebuilt sequence-blind below.
    {
      std::vector<std::size_t> candidates;
      candidates.reserve(set.pairs());
      for (std::size_t k = 0; k < set.pairs(); ++k) { candidates.push_back(set.source[k]); }
      std::sort(candidates.begin(), candidates.end());
      candidates.erase(std::unique(candidates.begin(), candidates.end()), candidates.end());
      OpenSwath::LightTargetedExperiment pool;
      pool.compounds.reserve(candidates.size());
      pool.transitions.reserve(candidates.size());
      std::unordered_map<std::string, std::size_t> index_of;
      for (const std::size_t i : candidates)
      {
        const std::string id = idOf(i);
        if (picked_ids.count(id) || !index_of.emplace(id, i).second) { continue; }
        OpenSwath::LightCompound c;
        c.id = id;
        c.sequence = strippedSequence(library.strings().get(pre.modified_sequence[i]));
        c.charge = pre.charge[i];
        c.rt = set.rt_scale.toAssay(pre.irt[i]);
        c.drift_time = -1.0;
        OpenSwath::LightTransition t;
        t.transition_name = id + "_sum";
        t.peptide_ref = id;
        t.library_intensity = 0.0;
        for (std::uint32_t k = 0; k < pre.transition_count[i]; ++k)
        { t.library_intensity += library.transitions().library_intensity[pre.transition_begin[i] + k]; }
        t.precursor_mz = fromFixed(pre.mz[i]);
        t.decoy = false;
        pool.compounds.push_back(std::move(c));
        pool.transitions.push_back(std::move(t));
      }
      // The sampler keeps the top share first and only then needs three
      // candidates (it throws below that): pools under ceil(3 / top share)
      // contribute no sampled seeds, and the calibration guard below says so.
      if (static_cast<double>(pool.compounds.size()) * seed_top_fraction >= 3.0)
      {
        try
        {
          QuietStreams quiet;
          const OpenSwath::LightTargetedExperiment sampled = OpenMS::OpenSwathHelper::sampleExperiment(
            pool, seed_bins, seeds_per_bin, seed_sampling_seed, true, seed_top_fraction);
          for (const auto& c : sampled.compounds) { picked.push_back(index_of.at(c.id)); }
        }
        catch (const std::exception& e) { detail["seed_sampling_failure"] = e.what(); }
      }
      detail["seed_candidates"] = pool.compounds.size();
    }
    const OpenSwath::LightTargetedExperiment seed = AssayBuilder::buildTargets(library, picked, set.rt_scale);
    cal.seeds = seed.compounds.size();
    detail["kit_files"] = kit_files;
    detail["kit_files_missing"] = kit_missing;
    detail["kit_sequences"] = kit.size();
    detail["kit_precursors_in_library"] = kit_in_library;
    detail["kit_precursors_outside_windows"] = kit_outside_windows;
    detail["seeds_kit"] = picked_kit;
    detail["sampling"] = {{"bins", seed_bins}, {"per_bin", seeds_per_bin}, {"top_fraction", seed_top_fraction},
                          {"seed", seed_sampling_seed}};
    detail["seed_mz_ppm"] = calibration_mz_ppm;
    detail["model"] = "linear";
    detail["outliers"] = "iter_residual, estimateBestPeptides";

    // 3. performRTNormalization on a COPY of the maps (the m/z estimate wraps
    //    the maps it is given; extraction must read the originals).
    std::string failure;
    double mz_estimate = -1.0, ms1_mz_estimate = -1.0;
    if (cal.seeds < 3)
    { failure = std::to_string(cal.seeds) + " seed assays (" + std::to_string(detail["seed_candidates"].get<std::size_t>()) + " candidates)"; }
    else
    {
      try
      {
        std::vector<OpenSwath::SwathMap> maps = run.maps;
        OpenMS::OpenSwathCalibrationWorkflow workflow;
        workflow.setLogType(OpenMS::ProgressLogger::NONE);
        OpenMS::TransformationDescription im_trafo;
        const Param ff = featureFinderParam(false, false);
        const OpenMS::ChromExtractParams cp_irt = chromParams(calibration_mz_ppm, -1.0, -1.0);
        const ThreadCap team(SearchParams::openswath_max_threads);
        QuietStreams quiet;
        // pasef false: this version searches ion-mobility runs by m/z and RT only.
        cal.rt = workflow.performRTNormalization(seed, maps, im_trafo, params().calibration_min_rsq,
                                                 params().calibration_min_coverage, ff, cp_irt, irtDetectionParam(),
                                                 massCorrectionParam(), "", 0, false, false);
        cal.im = im_trafo;
        mz_estimate = workflow.getEstimatedMzWindow();
        ms1_mz_estimate = workflow.getEstimatedMs1MzWindow();
      }
      catch (const std::exception& e)
      {
        failure = std::string("performRTNormalization failed: ") + e.what();
        // The stock binned-coverage check (8 of 10 library-RT bins) fails when
        // too few seeds are found across the library's RT range: a run that
        // covers part of the gradient only (a slice), or a library that does
        // not match the run. Nothing here can tell the two apart.
        if (failure.find("not enough bins") != std::string::npos)
        {
          failure += " -- the seeds found in the run do not cover most of the library RT range: either the run covers "
                     "only part of the gradient (a slice), or the library does not match the run (organism, "
                     "modifications, a different gradient)";
        }
      }
    }

    // 4. Refit robustly, then validate: stock stops removing outliers as soon as
    //    r^2 reaches min_rsq, and accepts a handful of points and any slope.
    std::vector<std::pair<double, double>> points;
    if (failure.empty())
    {
      std::vector<std::pair<double, double>> stock;
      for (const auto& dp : cal.rt.getDataPoints()) { stock.emplace_back(dp.first, dp.second); }
      detail["points_stock"] = stock.size();
      const RobustLine line = robustLine(stock, robust_cut, robust_min_scale);
      if (!line.valid())
      { failure = "no line through the " + std::to_string(stock.size()) + " calibration points"; }
      else
      {
        for (std::size_t k = 0; k < stock.size(); ++k) { if (line.inlier[k]) { points.push_back(stock[k]); } }
        detail["robust"] = {{"cut_sd", robust_cut}, {"inliers", line.inliers}, {"iterations", line.iterations},
                            {"scale_library", num(line.scale)},
                            {"scale_s", num(line.slope > 0 ? line.scale / line.slope : std::numeric_limits<double>::quiet_NaN())}};
        if (points.size() >= 2)
        {
          cal.rt = OpenMS::TransformationDescription();
          cal.rt.setDataPoints(points);
          cal.rt.fitModel("linear", Param());
        }
      }
    }
    cal.points = failure.empty() ? points.size() : 0;
    double slope = std::numeric_limits<double>::quiet_NaN(), intercept = slope;
    if (failure.empty())
    {
      const double n = static_cast<double>(points.size());
      double sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0, ylo = std::numeric_limits<double>::infinity(), yhi = -ylo;
      for (const auto& dp : points)
      {
        sx += dp.first; sy += dp.second; sxx += dp.first * dp.first; syy += dp.second * dp.second; sxy += dp.first * dp.second;
        ylo = std::min(ylo, dp.second); yhi = std::max(yhi, dp.second);
      }
      const double vx = n * sxx - sx * sx, vy = n * syy - sy * sy, cxy = n * sxy - sx * sy;
      cal.rsq = (vx > 0 && vy > 0) ? (cxy * cxy) / (vx * vy) : std::numeric_limits<double>::quiet_NaN();
      slope = vx > 0 ? cxy / vx : std::numeric_limits<double>::quiet_NaN();   // assay units per second
      intercept = n > 0 ? (sy - slope * sx) / n : slope;
      double seed_lo = std::numeric_limits<double>::infinity(), seed_hi = -seed_lo;
      for (const auto& c : seed.compounds) { seed_lo = std::min(seed_lo, c.rt); seed_hi = std::max(seed_hi, c.rt); }
      const double span = (seed_hi > seed_lo && points.size() > 1) ? (yhi - ylo) / (seed_hi - seed_lo) : 0.0;
      cal.coverage = cal.seeds > 0 ? static_cast<double>(points.size()) / static_cast<double>(cal.seeds) : 0.0;
      const double implied_s = (slope > 0) ? 100.0 / slope : std::numeric_limits<double>::infinity();
      detail["seed_rt_span_covered"] = num(span);
      detail["library_range_implied_s"] = num(implied_s);
      detail["run_span_s"] = run_span;
      if (points.size() < calibration_min_points)
      { failure = std::to_string(points.size()) + " calibration points, fewer than " + std::to_string(calibration_min_points); }
      else if (!(cal.rsq >= params().calibration_min_rsq))
      { failure = "r^2 " + fixed(cal.rsq, 3) + " below search:calibration_min_rsq " + fixed(params().calibration_min_rsq, 2); }
      else if (!(slope > 0))
      { failure = "the RT map does not increase (slope " + fixed(slope, 6) + " library units per second)"; }
      else if (implied_s > calibration_span_factor * run_span || implied_s < run_span / calibration_span_factor)
      {
        failure = "the RT map spreads the library RT range over " + fixed(implied_s, 0) + " s of a " + fixed(run_span, 0) +
                  " s run (slope " + fixed(slope, 6) + "); more than " + fixed(calibration_span_factor, 0) + "x off";
      }
      else if (span < calibration_min_span)
      { failure = "the points cover " + fixed(100.0 * span, 1) + " % of the seeds' library RT span, less than " + fixed(100.0 * calibration_min_span, 0) + " %"; }
    }
    detail["slope_library_per_s"] = num(slope);
    detail["intercept_library"] = num(intercept);

    if (!failure.empty())
    {
      detail["failure"] = failure;
      if (!params().allow_bootstrap)
      {
        throw SearchAbort("search: RT calibration failed: " + failure + " (" + std::to_string(cal.seeds) + " seed assays, " +
                          std::to_string(picked_kit) + " of them iRT/CiRT kit precursors, " + std::to_string(cal.points) +
                          " points kept); nothing is searched on a guessed RT scale unless search:allow_bootstrap is set");
      }
      // The test hook: the library RT range laid linearly over the run, and a
      // window as wide as the run unless search:rt_window says otherwise.
      cal.bootstrap = true;
      cal.rt = OpenMS::TransformationDescription();
      cal.rt.setDataPoints(std::vector<std::pair<double, double>>{{run_first, 0.0}, {run_last, 100.0}});
      cal.rt.fitModel("linear", Param());
      cal.im = OpenMS::TransformationDescription();
      cal.rt_window = params().rt_window > 0 ? params().rt_window : 2.0 * run_span;
      cal.mz_ppm = params().mz_ppm > 0 ? params().mz_ppm : default_mz_ppm;
      cal.ms1_mz_ppm = params().mz_ppm > 0 ? params().mz_ppm : default_mz_ppm;
      cal.im_window = -1.0;
      detail["rt_window_rule"] = params().rt_window > 0 ? "search:rt_window" : "bootstrap: twice the run span";
      cal.provenance_json = detail.dump();
      return cal;
    }

    // 5. Windows. RT: the seed residuals (TOPP's quantile and padding), with a
    //    floor; m/z: the calibration's estimate may only NARROW the default.
    {
      OpenMS::TransformationDescription inverse = cal.rt;
      inverse.invert();
      std::vector<double> residual;
      residual.reserve(points.size());
      for (const auto& dp : points) { residual.push_back(std::fabs(inverse.apply(dp.second) - dp.first)); }
      detail["residual_s"] = {{"median", num(quantile(residual, 0.5))}, {"p95", num(quantile(residual, 0.95))},
                              {"p99", num(quantile(residual, 0.99))}, {"max", num(quantile(residual, 1.0))}};
    }
    if (params().rt_window > 0)
    {
      cal.rt_window = params().rt_window;
      detail["rt_window_rule"] = "search:rt_window";
    }
    else
    {
      const double estimate = cal.rt.estimateWindow(rt_window_quantile, true, true, rt_window_padding);
      cal.rt_window = std::min(2.0 * run_span, std::max(rt_window_floor_s, estimate));
      detail["rt_window_estimate_s"] = num(estimate);
      detail["rt_window_rule"] = "2 x " + fixed(rt_window_quantile, 2) + " quantile of the seed residuals x " +
                                 fixed(rt_window_padding, 1) + ", at least " + fixed(rt_window_floor_s, 0) + " s";
    }
    auto mzWindow = [&](double estimate) {
      if (params().mz_ppm > 0) { return params().mz_ppm; }
      if (!(estimate > 0)) { return default_mz_ppm; }
      return std::min(default_mz_ppm, std::min(mz_cap_ppm, std::max(mz_floor_ppm, estimate)));
    };
    cal.mz_ppm = mzWindow(mz_estimate);
    cal.ms1_mz_ppm = mzWindow(ms1_mz_estimate);
    detail["mz_estimate_ppm"] = num(mz_estimate > 0 ? mz_estimate : std::numeric_limits<double>::quiet_NaN());
    detail["ms1_mz_estimate_ppm"] = num(ms1_mz_estimate > 0 ? ms1_mz_estimate : std::numeric_limits<double>::quiet_NaN());
    detail["mz_rule"] = params().mz_ppm > 0 ? std::string("search:mz_ppm")
                                            : "estimate clamped to [" + fixed(mz_floor_ppm, 0) + ", " + fixed(mz_cap_ppm, 0) +
                                                "] ppm, never wider than the default " + fixed(default_mz_ppm, 0) + " ppm";
    cal.im_window = -1.0;
    if (run.ion_mobility && params().im_window > 0)
    { warn("search: search:im_window " + fixed(params().im_window, 3) + " is ignored: this version searches ion-mobility runs without an ion-mobility window"); }
    cal.provenance_json = detail.dump();
    return cal;
  }

  // ---- extraction -------------------------------------------------------------

  void Identifier::extract(const SearchSet& set, RunData& run, const Calibration& calibration, PeakGroups& out)
  {
    bool has_ms1 = false;
    for (const auto& m : run.maps) { has_ms1 = has_ms1 || (m.ms1 && m.sptr && m.sptr->getNrSpectra() > 0); }
    const bool ms1 = params().ms1 && has_ms1;
    if (params().ms1 && !has_ms1) { warn("search: search:ms1 is on but the run has no MS1 spectra; MS1 traces and sub-scores are off"); }

    const std::vector<std::string> columns = scoreColumns(ms1);
    out.setColumns(columns);
    const std::size_t width = columns.size();
    // Meta value names -> registry indices. registerName, not getIndex: a name
    // OpenSWATH has not used yet in this process has no index, and getIndex
    // would answer "unknown" -- the first search would read that sub-score as
    // missing everywhere and the next one would not.
    auto& registry = OpenMS::MetaInfoInterface::metaRegistry();
    std::vector<OpenMS::UInt> column_index;
    column_index.reserve(width);
    for (const auto& name : columns) { column_index.push_back(registry.registerName(name)); }
    const OpenMS::UInt peptide_ref = registry.registerName("PeptideRef");
    const OpenMS::UInt left_width = registry.registerName("leftWidth");
    const OpenMS::UInt right_width = registry.registerName("rightWidth");

    const Param ff = featureFinderParam(ms1, false);
    const OpenMS::ChromExtractParams cp = chromParams(calibration.mz_ppm, -1.0, calibration.rt_window);
    const OpenMS::ChromExtractParams cp_ms1 = chromParams(calibration.ms1_mz_ppm, -1.0, calibration.rt_window);
    // pasef false: this version searches ion-mobility runs by m/z and RT only.
    OpenMS::OpenSwathWorkflow workflow(ms1, false, false, false, -1);
    workflow.setLogType(OpenMS::ProgressLogger::NONE);
    OpenMS::OpenSwathOSWWriter osw("", 0);   // an empty name: inactive
    OpenMS::NoopMSDataWritingConsumer chromatograms("");

    const auto chunks = AssayBuilder::chunks(set, params().chunk);
    const float nan = std::numeric_limits<float>::quiet_NaN();
    std::size_t done = 0;
    for (std::size_t c = 0; c < chunks.size(); ++c)
    {
      const auto started = std::chrono::steady_clock::now();
      const AssayChunk assays = AssayBuilder::build(set, chunks[c], AssayOptions{});
      std::unordered_map<std::string, std::size_t> precursor_of;
      precursor_of.reserve(assays.precursor.size());
      for (std::size_t k = 0; k < assays.precursor.size(); ++k)
      { precursor_of.emplace(assays.experiment.compounds[k].id, assays.precursor[k]); }

      // Extract, flatten, free the features. Rows are ordered by (precursor,
      // apex RT, sub-scores) before they get their ids: the FeatureMap's own
      // order depends on which thread finished first.
      struct Row { std::size_t precursor; double apex, left, right; std::size_t offset; };
      std::vector<Row> rows;
      std::vector<float> values;
      std::size_t n_features = 0;
      {
        OpenMS::FeatureMap features;
        {
          const ThreadCap team(SearchParams::openswath_max_threads);
          if (c == 0 && team.used() < team.requested())
          {
            info("search extraction: OpenSWATH runs on " + std::to_string(team.used()) + " of the " +
                 std::to_string(team.requested()) + " threads; stock OpenSWATH does not scale beyond " +
                 std::to_string(SearchParams::openswath_max_threads));
          }
          QuietStreams quiet;
          workflow.performExtraction(run.maps, calibration.rt, cp, cp_ms1, ff, assays.experiment, features, true, osw,
                                     &chromatograms, static_cast<int>(params().batch_size), ms1_isotope_traces, false);
        }
        n_features = features.size();
        rows.reserve(features.size());
        values.reserve(features.size() * width);
        for (const OpenMS::Feature& f : features)
        {
          const std::string id = f.getMetaValue(peptide_ref).toString();
          const auto it = precursor_of.find(id);
          if (it == precursor_of.end())
          { throw std::logic_error("search: OpenSWATH returned a feature for '" + id + "', which is not in the chunk it was given"); }
          Row r;
          r.precursor = it->second;
          r.apex = f.getRT();
          const OpenMS::DataValue& l = f.getMetaValue(left_width);
          const OpenMS::DataValue& h = f.getMetaValue(right_width);
          r.left = l.isEmpty() ? std::numeric_limits<double>::quiet_NaN() : static_cast<double>(l);
          r.right = h.isEmpty() ? std::numeric_limits<double>::quiet_NaN() : static_cast<double>(h);
          r.offset = values.size();
          for (const OpenMS::UInt k : column_index)
          {
            const OpenMS::DataValue& v = f.getMetaValue(k);
            const bool numeric = v.valueType() == OpenMS::DataValue::DOUBLE_VALUE || v.valueType() == OpenMS::DataValue::INT_VALUE;
            values.push_back(numeric ? static_cast<float>(static_cast<double>(v)) : nan);
          }
          rows.push_back(r);
        }
      }

      std::sort(rows.begin(), rows.end(), [&](const Row& a, const Row& b) {
        if (a.precursor != b.precursor) { return a.precursor < b.precursor; }
        return rowLess(a.apex, values.data() + a.offset, b.apex, values.data() + b.offset, width);
      });
      std::size_t targets = 0, decoys = 0;
      std::int64_t rank = 0;
      for (std::size_t r = 0; r < rows.size(); ++r)
      {
        const bool first = r == 0 || rows[r].precursor != rows[r - 1].precursor;
        rank = first ? 0 : rank + 1;
        if (first) { (set.isDecoy(rows[r].precursor) ? decoys : targets)++; }
        out.add(set, rows[r].precursor, rank, values.data() + rows[r].offset, static_cast<float>(rows[r].apex),
                static_cast<float>(rows[r].left), static_cast<float>(rows[r].right), nan);
      }
      done += assays.precursor.size();
      const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
      info("search extraction: chunk " + std::to_string(c + 1) + "/" + std::to_string(chunks.size()) + ", " +
           std::to_string(assays.precursor.size()) + " precursors: " + std::to_string(n_features) + " peak groups for " +
           std::to_string(targets) + " targets and " + std::to_string(decoys) + " decoys (" + std::to_string(done) + "/" +
           std::to_string(set.size()) + " done, " + fixed(seconds, 1) + " s)");
    }
  }
}
