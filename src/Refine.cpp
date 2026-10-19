// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <odia/DIANNLibraryFile.h>
#include <odia/AtomicFile.h>
#include <odia/DeviceRoute.h>
#include <odia/Hardware.h>
#include <odia/LibraryExport.h>
#include <memory>
#include <odia/Library.h>
#include <odia/LibraryGenerator.h>
#include <odia/LibraryRefiner.h>
#include <odia/Trace.h>
#include <odia/search/Identifier.h>
#include <odia/search/MzPeakRun.h>
#include <odia/search/SearchParams.h>

#ifdef DIALIBGEN_WITH_FINETUNE
#include <odia/tune/Trainer.h>
#include <odia/PeptDeepEncoder.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <unordered_set>
#endif

#include "DIALibGen.h"

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/APPLICATIONS/ParameterInformation.h>

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <fstream>
#include <future>
#include <iomanip>
#include <map>
#include <random>
#include <sstream>
#include <system_error>

using json = nlohmann::json;

namespace
{
  /// The run's fingerprint for the provenance: the file's FNV-1a64, or for an
  /// unpacked mzPeak archive (a directory) the FNV-1a64 of its members'
  /// "relative path<TAB>FNV-1a64" lines in path order.
  std::string runHash(const std::string& run)
  {
    std::error_code ec;
    if (!std::filesystem::is_directory(run, ec)) { return ODIA::DIANNLibraryFile::hashFile(run); }
    std::vector<std::filesystem::path> members;
    for (auto it = std::filesystem::recursive_directory_iterator(run, ec); !ec && it != std::filesystem::recursive_directory_iterator();
         it.increment(ec))
    {
      if (it->is_regular_file(ec)) { members.push_back(it->path()); }
    }
    std::sort(members.begin(), members.end());
    std::uint64_t h = 1469598103934665603ull;
    for (const auto& m : members)
    {
      const std::string line = std::filesystem::relative(m, run).generic_string() + "\t" + ODIA::DIANNLibraryFile::hashFile(m.string()) + "\n";
      for (const unsigned char c : line) { h ^= c; h *= 1099511628211ull; }
    }
    std::ostringstream out;
    out << std::hex << std::setw(16) << std::setfill('0') << h;
    return out.str();
  }

  const char* rtUnitName(ODIA::RefineParams::RtUnit u)
  { return u == ODIA::RefineParams::RtUnit::MinMax ? "minmax" : "observed"; }
  const char* intensityNormName(ODIA::RefineParams::IntensityNorm n)
  {
    switch (n)
    {
      case ODIA::RefineParams::IntensityNorm::BasePeak: return "base_peak";
      case ODIA::RefineParams::IntensityNorm::Sum: return "sum";
      case ODIA::RefineParams::IntensityNorm::Raw: return "raw";
      case ODIA::RefineParams::IntensityNorm::LibraryMax: break;
    }
    return "library_max";
  }

  ODIA::RefineParams::IntensityNorm intensityNormFrom(const std::string& s)
  {
    if (s == "library_max") { return ODIA::RefineParams::IntensityNorm::LibraryMax; }
    if (s == "base_peak") { return ODIA::RefineParams::IntensityNorm::BasePeak; }
    if (s == "sum") { return ODIA::RefineParams::IntensityNorm::Sum; }
    if (s == "raw") { return ODIA::RefineParams::IntensityNorm::Raw; }
    throw std::runtime_error("intensity_norm must be library_max, base_peak, sum or raw; got '" + s + "'");
  }

  const char* dedupName(ODIA::RefineParams::Dedup d)
  { return d == ODIA::RefineParams::Dedup::HighestEvidence ? "highest_evidence" : "lowest_q"; }

  json effectiveConfig(const ODIA::RefineParams& p)
  {
    return json{
      {"schema_version", 1},
      {"filter", p.filter},
      {"q_precursor", p.q_precursor}, {"q_global", p.q_global}, {"q_protein", p.q_protein},
      {"require_gates", p.require_gates},
      {"min_fragments", p.min_fragments},
      {"write_rt", p.write_rt}, {"write_im", p.write_im}, {"write_intensity", p.write_intensity},
      {"rt_unit", rtUnitName(p.rt_unit)},
      {"dedup", dedupName(p.dedup)},
      {"im_min_charge", p.im_min_charge},
      {"im_ramp_top", p.im_ramp_top}, {"im_ramp_margin", p.im_ramp_margin},
      {"min_match_fraction", p.min_match_fraction},
      {"intensity_min_correlation", p.intensity_min_correlation},
      {"intensity_restrict", p.intensity_restrict}, {"intensity_rerank", p.intensity_rerank},
      {"intensity_min_fragments", p.intensity_min_fragments},
      {"intensity_norm", intensityNormName(p.intensity_norm)},
      {"intensity_min_relative", p.intensity_min_relative},
      {"intensity_mz_tol_ppm", p.intensity_mz_tol_ppm},
      {"intensity_max_mz_mismatch", p.intensity_max_mz_mismatch},
      {"allow_mixed_intensity", p.allow_mixed_intensity}};
  }

  json num(double v) { return std::isfinite(v) ? json(v) : json(nullptr); }
}

  static void applyJson_(const json& j, ODIA::RefineParams& p)
  {
    if (!j.is_object()) { throw std::runtime_error("refinement config must be a JSON object"); }
    const json ref = effectiveConfig(p);
    for (const auto& [k, v] : j.items())
    {
      if (k == "tuning") { continue; }   // the training settings: refine_ reads them
      if (!ref.contains(k)) { throw std::runtime_error("unknown config key: " + k); }
    }
    if (j.contains("schema_version") && j["schema_version"] != 1)
    { throw std::runtime_error("unsupported schema_version " + j["schema_version"].dump() + "; this tool writes 1"); }

    auto get_bool = [&](const char* k, bool& dst) { if (j.contains(k)) { dst = j.at(k).get<bool>(); } };
    auto get_num = [&](const char* k, double& dst, double lo, double hi)
    {
      if (!j.contains(k)) { return; }
      const double v = j.at(k).get<double>();
      if (!(v >= lo && v <= hi)) { throw std::runtime_error(std::string(k) + " out of range"); }
      dst = v;
    };
    get_bool("filter", p.filter); get_bool("require_gates", p.require_gates);
    get_bool("write_rt", p.write_rt); get_bool("write_im", p.write_im); get_bool("write_intensity", p.write_intensity);
    get_num("q_precursor", p.q_precursor, 0.0, 1.0); get_num("q_global", p.q_global, 0.0, 1.0); get_num("q_protein", p.q_protein, 0.0, 1.0);
    get_num("im_ramp_top", p.im_ramp_top, 0.0, 10.0); get_num("im_ramp_margin", p.im_ramp_margin, 0.0, 1.0);
    get_num("min_match_fraction", p.min_match_fraction, 0.0, 1.0);
    get_bool("intensity_restrict", p.intensity_restrict); get_bool("intensity_rerank", p.intensity_rerank);
    get_bool("allow_mixed_intensity", p.allow_mixed_intensity);
    get_num("intensity_min_correlation", p.intensity_min_correlation, -1.0, 1.0);
    get_num("intensity_min_relative", p.intensity_min_relative, 0.0, 1.0);
    get_num("intensity_mz_tol_ppm", p.intensity_mz_tol_ppm, 0.0, 1000.0);
    get_num("intensity_max_mz_mismatch", p.intensity_max_mz_mismatch, 0.0, 1.0);
    auto get_size = [&](const char* name, std::size_t& value)
    {
      if (!j.contains(name)) { return; }
      const auto& v = j.at(name);
      if (!v.is_number_integer() || (v.is_number_integer() && !v.is_number_unsigned() && v.get<std::int64_t>() < 0))
      { throw std::runtime_error(std::string(name) + " must be a nonnegative integer"); }
      value = v.get<std::size_t>();
    };
    get_size("intensity_min_fragments", p.intensity_min_fragments);
    if (j.contains("intensity_norm")) { p.intensity_norm = intensityNormFrom(j.at("intensity_norm").get<std::string>()); }
    get_size("min_fragments", p.min_fragments);
    if (j.contains("im_min_charge"))
    {
      if (!j.at("im_min_charge").is_number_integer() || j.at("im_min_charge") < 1 || j.at("im_min_charge") > 8)
      { throw std::runtime_error("im_min_charge must be an integer in [1, 8]"); }
      p.im_min_charge = j.at("im_min_charge").get<int>();
    }
    if (j.contains("rt_unit"))
    {
      const auto s = j.at("rt_unit").get<std::string>();
      if (s == "observed") { p.rt_unit = ODIA::RefineParams::RtUnit::Observed; }
      else if (s == "minmax") { p.rt_unit = ODIA::RefineParams::RtUnit::MinMax; }
      else { throw std::runtime_error("rt_unit must be observed or minmax, not '" + s + "'"); }
    }
    if (j.contains("dedup"))
    {
      const auto s = j.at("dedup").get<std::string>();
      if (s == "lowest_q") { p.dedup = ODIA::RefineParams::Dedup::LowestQ; }
      else if (s == "highest_evidence") { p.dedup = ODIA::RefineParams::Dedup::HighestEvidence; }
      else { throw std::runtime_error("dedup must be lowest_q or highest_evidence, not '" + s + "'"); }
    }
  }


void DIALibGen::registerRefinementOptions_()
  {
    // Shorthands: every option's advanced flag comes from advanced_ (the
    // visible list in DIALibGen.cpp), so registration order is help order.
    const auto flag = [this](const char* name, const std::string& text) { registerFlag_(name, text, advanced_(name)); };
    const auto deprecated = [](const char* use) { return std::string("Deprecated (0.12, removed in 0.13): use ") + use; };
    registerInputFile_("ids", "<file>", "", "Refine, tune: the DIA-NN report.parquet of one run (or, with "
                       "-empirical_library, a pre-filtered library)", false, advanced_("ids"));
    setValidFormats_("ids", {"parquet"}, false);
    flag("no_filter", "Refine: keep the precursors the run did not identify (needs -tune or -no_write_rt, so that "
                      "RT stays in one unit; JSON key filter)");
    flag("no_write_rt", "Refine: keep the predicted RT (JSON key write_rt)");
    flag("write_im", "Refine: also write the observed 1/K0 (ion-mobility runs), for charges >= -im_min_charge");
    flag("write_intensity", "Refine: also write the observed fragment intensities (a report exported with "
                            "--report-lib-info (DIA-NN 1.9) or --export-quant (2.x)); fragments are cross-checked on "
                            "m/z, and a run that replaces nothing is an error");
    registerDoubleOption_("q_precursor", "<q>", 0.01, "Refine: precursor q-value cut-off (>= 1 disables); -q_global "
                          "and -q_protein add gates (--helphelp)", false, advanced_("q_precursor"));
    // Registered in EVERY build, including one without libtorch, so --help,
    // -write_ini and -write_ctd describe the same tool on every platform. -tune
    // itself is what fails when the stage was not compiled in.
    flag("tune", "Refine: fine-tune RT (and CCS) on -ids and re-predict the whole library first, so precursors the run "
                 "did not identify are corrected too (a training run plus whole-library prediction)");
    registerStringOption_("tune_heads", "<which>", "auto", "Tune: models to fine-tune; auto = RT, plus CCS when the "
                          "report has enough ion-mobility data", false, advanced_("tune_heads"));
    setValidStrings_("tune_heads", {"auto", "rt", "ccs", "both"});
    registerStringOption_("tune_out_models", "<dir>", "", "Tune: also keep the tuned models (ONNX, with their .tune.json "
                          "and .trajectory.tsv) in this directory", false, advanced_("tune_out_models"));

    registerInputFile_("run", "<file>", "",
                       "Experimental: a centroided DIA run (mzML, or an mzPeak archive or its unpacked directory) that "
                       "the built-in identification searches instead of reading -ids; not validated yet (-ids is). "
                       "Settings: search:", false, advanced_("run"));
    // Not checked by TOPP (OpenMS 3.5 has no mzPeak type): refine_ checks the extension.
    setValidFormats_("run", {"mzML", "mzpeak"}, false);
    registerOutputFile_("out_ids", "<file>", "",
                        "With -run: the identification report it writes (Parquet, DIA-NN column names), then used "
                        "as -ids. Default: <out>.ids.parquet. Never overwritten", false, advanced_("out_ids"));
    setValidFormats_("out_ids", {"parquet"}, false);
    registerOutputFile_("out_report", "<file>", "",
                        "Refine, tune: per-axis residuals before the overwrite (TSV)", false, advanced_("out_report"));
    setValidFormats_("out_report", {"tsv"}, false);
    flag("empirical_library", "Refine: -ids is a pre-filtered empirical library, not a report: gates whose columns are "
                              "missing are bypassed and each bypass is recorded (without it a missing gate column "
                              "is an error; JSON key require_gates)");
    registerDoubleOption_("q_global", "<q>", 1.0, "Refine: global (peptide-level) q-value cut-off; 1 = off",
                          false, advanced_("q_global"));
    registerDoubleOption_("q_protein", "<q>", 1.0, "Refine: protein q-value cut-off; 1 = off", false, advanced_("q_protein"));
    registerIntOption_("min_fragments", "<n>", 0, "Refine: least distinct fragments the report identified per "
                       "precursor (needs fragment identities in -ids); 0 = off", false, advanced_("min_fragments"));
    registerStringOption_("rt_unit", "<unit>", "observed", "Refine: observed = the run's own RT units; minmax = "
                          "rescaled to 0..100 over the matched set (not with -no_filter or tuning)", false,
                          advanced_("rt_unit"));
    setValidStrings_("rt_unit", {"observed", "minmax"});
    registerStringOption_("dedup", "<rule>", "lowest_q", "Refine: which observation wins for a precursor seen more "
                          "than once: lowest_q (ties: lower PEP, then first seen) or highest_evidence (needs an "
                          "Evidence column)", false, advanced_("dedup"));
    setValidStrings_("dedup", {"lowest_q", "highest_evidence"});
    registerDoubleOption_("min_match_fraction", "<f>", 0.0, "Refine: refuse when fewer than this fraction of the passing "
                          "report precursors are in the library; 0 = only when none is", false, advanced_("min_match_fraction"));
    registerIntOption_("im_min_charge", "<z>", 2, "Refine, with -write_im: charges below this keep their predicted 1/K0 "
                       "(z1 is censored at the ramp top on timsTOF diaPASEF)", false, advanced_("im_min_charge"));
    registerDoubleOption_("im_ramp_top", "<1/K0>", 0.0, "Refine, with -write_im: the mobility ramp's top, if known; "
                          "observations within -im_ramp_margin of it are censored. 0 = unknown", false,
                          advanced_("im_ramp_top"));
    registerDoubleOption_("im_ramp_margin", "<1/K0>", 0.02, "Refine, with -write_im: see -im_ramp_top", false,
                          advanced_("im_ramp_margin"));
    registerDoubleOption_("intensity_min_correlation", "<r>", 0.0, "Refine, with -write_intensity: least fragment "
                          "quality (correlation in DIA-NN 1.9, Score in 2.x); -1 disables the gate", false,
                          advanced_("intensity_min_correlation"));
    flag("intensity_no_restrict", "Refine, with -write_intensity: replace a precursor only when every transition is "
                                  "trusted, else keep its predictions (neutral-loss transitions never match; JSON key "
                                  "intensity_restrict)");
    flag("intensity_no_rerank", "Refine, with -write_intensity: keep a replaced precursor's transitions in their order "
                                "(JSON key intensity_rerank)");
    registerIntOption_("intensity_min_fragments", "<n>", 3, "Refine, with -write_intensity: a replaced precursor keeps "
                       "at least this many transitions, or keeps its predictions", false, advanced_("intensity_min_fragments"));
    registerStringOption_("intensity_norm", "<rule>", "library_max", "Refine, with -write_intensity: how an observed area "
                          "becomes an intensity; library_max scales to the precursor's previous maximum (the file's own "
                          "convention)", false, advanced_("intensity_norm"));
    setValidStrings_("intensity_norm", {"library_max", "base_peak", "sum", "raw"});
    registerDoubleOption_("intensity_min_relative", "<fraction>", 1e-4, "Refine, with -write_intensity: drop observed "
                          "intensities below this fraction of the kept maximum", false, advanced_("intensity_min_relative"));
    registerDoubleOption_("intensity_mz_tol_ppm", "<ppm>", 20.0, "Refine, with -write_intensity: a fragment matches only "
                          "when the report's m/z agrees with the library's within this", false,
                          advanced_("intensity_mz_tol_ppm"));
    registerDoubleOption_("intensity_max_mz_mismatch", "<frac>", 0.01, "Refine, with -write_intensity: refuse when more "
                          "than this fraction of matches fail the m/z check; 1 surveys a suspect pairing", false,
                          advanced_("intensity_max_mz_mismatch"));
    flag("allow_mixed_intensity", "Refine: permit -write_intensity with -no_filter, which leaves observed and predicted "
                                  "intensities in one library (recorded)");

    registerStringOption_("tune_models", "<dir>", "", deprecated("-models"), false, advanced_("tune_models"));
    flag("tune_predict_gpu", deprecated("-machine:device (re-prediction follows it)"));
    registerIntOption_("tune_predict_sessions", "<n>", 0, "Tune: ONNX Runtime sessions for the re-prediction; 0 = "
                       "-threads", false, advanced_("tune_predict_sessions"));
    flag("tune_keep_free_cysteine_offset", "Tune: keep the free-cysteine RT offset when re-predicting with the tuned RT "
                                           "model (off: the offset was fitted to the stock model, and a model tuned on "
                                           "the run has learned the effect)");

    registerTOPPSubsection_("filter", "Fine-tuning: the training set (not the library: that is -q_precursor and the other gates)");
    registerDoubleOption_("filter:q_value", "<q>", 0.01, "Precursor q-value cut-off of the training set", false,
                          advanced_("filter:q_value"));
    registerIntOption_("filter:min_charge", "<z>", 2, "CCS: lowest precursor charge (RT pools all charges; z1 is "
                       "censored at the ramp top on timsTOF)", false, advanced_("filter:min_charge"));
    setMinInt_("filter:min_charge", 1); setMaxInt_("filter:min_charge", 8);
    flag("filter:allow_z1", "CCS: permit filter:min_charge 1 (censored observations then train)");
    registerDoubleOption_("filter:rt_spread_max", "<min>", 0.2, "Drop a sequence whose charge states' RTs span more "
                          "than this (minutes)", false, advanced_("filter:rt_spread_max"));
    registerDoubleOption_("filter:rt_max_minutes", "<min>", 0.0, "RT normalisation; 0 = the run's latest observed RT",
                          false, advanced_("filter:rt_max_minutes"));

    registerTOPPSubsection_("cohort", "Fine-tuning: protein-level cohorts (frozen before subsampling)");
    registerIntOption_("cohort:train_size", "<n>", 0, "Training units (sequences for RT, sequence x charge for CCS); "
                       "0 = the whole pool", false, advanced_("cohort:train_size"));
    registerDoubleOption_("cohort:train_frac", "<f>", 0.0, "Instead of train_size: the fraction of the pool", false,
                          advanced_("cohort:train_frac"));
    flag("cohort:full_fit", "Train on every unit of the run, test and validation included: validation and test are then "
                            "in-sample, so judge the result on a different run");
    flag("cohort:no_inner_val", "No inner validation cohort: its units rejoin the pool and checkpoints are selected on "
                                "the test cohort, which is then no longer held out");

    registerTOPPSubsection_("train", "Fine-tuning: recipe");
    registerIntOption_("train:epochs", "<n>", 100, "Length of the cosine schedule (and the most epochs)", false,
                       advanced_("train:epochs"));
    registerIntOption_("train:warmup", "<n>", 10, "Linear warm-up epochs", false, advanced_("train:warmup"));
    registerDoubleOption_("train:lr", "<lr>", 1e-4, "Peak learning rate (Adam)", false, advanced_("train:lr"));
    registerIntOption_("train:batch_size", "<n>", 1024, "Batch size within a length group", false,
                       advanced_("train:batch_size"));

    registerTOPPSubsection_("stop", "Fine-tuning: convergence");
    registerIntOption_("stop:eval_every", "<n>", 1, "Validate every n epochs", false, advanced_("stop:eval_every"));
    registerIntOption_("stop:min_epochs", "<n>", 20, "Never stop before this epoch", false, advanced_("stop:min_epochs"));
    registerIntOption_("stop:patience", "<n>", 10, "Stop after this many epochs without progress (not before "
                       "max(min_epochs, warmup))", false, advanced_("stop:patience"));
    registerDoubleOption_("stop:rel_tol", "<f>", 0.005, "Progress: the selection metric beats the best so far by this "
                          "fraction", false, advanced_("stop:rel_tol"));
    registerDoubleOption_("stop:abs_tol", "<f>", 0.0, "Progress: ... by this absolute amount (0 = use rel_tol)", false,
                          advanced_("stop:abs_tol"));
    registerDoubleOption_("stop:max_seconds", "<s>", 0.0, "Time budget per head, checked between epochs (the final "
                          "evaluation may exceed it); 0 = none", false, advanced_("stop:max_seconds"));
    registerStringOption_("stop:select", "<metric>", "calibrated_sd", "Selection metric on the validation cohort", false,
                          advanced_("stop:select"));
    setValidStrings_("stop:select", {"calibrated_sd", "rmse"});

    registerTOPPSubsection_("machine", "Devices");
    registerStringOption_("machine:device", "<dev>", "auto",
                          "Auto: the GPU when one is usable and pays off (NVIDIA: CUDA, RT and CCS with cuDNN 9; Apple: "
                          "CoreML for MS2 of large libraries, MPS for fine-tuning), else the CPU; cpu: CPU only (the "
                          "reference); also cuda[:N], coreml, mps", false, advanced_("machine:device"));
    registerStringOption_("machine:coreml_units", "<units>", "gpu",
                          "CoreML compute units: gpu (CPU and GPU, fastest), all (Core ML's choice, Neural Engine "
                          "included), ane (CPU and Neural Engine), cpu (checks)", false, advanced_("machine:coreml_units"));
    setValidStrings_("machine:coreml_units", {"gpu", "all", "ane", "cpu"});
    registerIntOption_("machine:train_threads", "<n>", 0,
                       "CPU threads for fine-tuning (0: 16 on Linux and Windows, 8 on macOS). On Linux the tuned weights "
                       "depend on the count (recorded); on macOS they do not", false,
                       advanced_("machine:train_threads"));
    registerIntOption_("machine:threads", "<n>", 0, deprecated("-machine:train_threads"), false, advanced_("machine:threads"));
    registerStringOption_("machine:cores", "<which>", "all",
                          "MacOS: which cores -threads 0 counts: all, or performance (saves energy, not time). Does not "
                          "change the output", false, advanced_("machine:cores"));
    setValidStrings_("machine:cores", {"all", "performance"});
    flag("machine:no_cudnn", "CUDA fine-tuning without cuDNN (when only its loader shim is installed, as in "
                             "pytorch.org's libtorch zips); slower");
    registerIntOption_("machine:seed", "<n>", 20260803, "Seed of the training subsample and batch order", false,
                       advanced_("machine:seed"));
    registerStringOption_("machine:tf32", "<true/false>", "false",
                          "TensorFloat-32 on NVIDIA GPUs, for prediction and fine-tuning: faster, but about 1 % of "
                          "precursors select a different fragment set (off: GPU and CPU agree to float rounding)",
                          false, advanced_("machine:tf32"));
    setValidStrings_("machine:tf32", {"true", "false"});
    registerIntOption_("machine:gpu_sessions", "<n>", 4,
                       "CUDA: ONNX Runtime sessions per GPU and model, each with its own stream; CoreML: concurrent "
                       "batches. Does not change the output", false, advanced_("machine:gpu_sessions"));
    setMinInt_("machine:gpu_sessions", 1);
    setMaxInt_("machine:gpu_sessions", 16);
    registerIntList_("machine:gpu_devices", "<ids>", OpenMS::IntList(),
                     "CUDA devices to spread the sessions over (e.g. 0 1). Default: the first that works. Does not "
                     "change the output", false, advanced_("machine:gpu_devices"));
    setMinInt_("machine:gpu_devices", 0);
    flag("machine:require_gpu", "Fail instead of predicting on the CPU when -machine:device asks for a GPU (or auto finds "
                                "one) that a model cannot use");

    registerTOPPSubsection_("search", "Built-in identification with -run (experimental; docs/design/built-in-identification.md)");
    registerStringOption_("search:candidates", "<rule>", "evidence", "Candidate pairs: evidence = target-decoy pairs of "
                          "which either member has fragment evidence in the run (search:prefilter_*); random = a "
                          "deterministic, label-blind random subset (search:subset)", false, advanced_("search:candidates"));
    setValidStrings_("search:candidates", {"evidence", "random"});
    registerIntOption_("search:subset", "<n>", 0, "Random candidates: targets drawn before the cap (0 = every eligible "
                       "target)", false, advanced_("search:subset"));
    registerIntOption_("search:max_pairs", "<n>", 200000, "Most target-decoy pairs searched; time is linear in pairs (0 = "
                       "no cap). Evidence candidates are capped per isolation window, RT decile and charge, keeping the "
                       "most evidence", false, advanced_("search:max_pairs"));
    registerIntOption_("search:prefilter_depth", "<n>", 5, "Evidence candidates: a pair is kept when either member has "
                       "this many of its 6 most intense predicted fragments among one spectrum's top peaks", false,
                       advanced_("search:prefilter_depth"));
    setMinInt_("search:prefilter_depth", 1); setMaxInt_("search:prefilter_depth", 6);
    registerIntOption_("search:prefilter_top_peaks", "<n>", 1000, "Evidence candidates: the most intense peaks per MS2 "
                       "spectrum that the prefilter matches", false, advanced_("search:prefilter_top_peaks"));
    setMinInt_("search:prefilter_top_peaks", 1);
    registerDoubleOption_("search:prefilter_ppm", "<ppm>", 10.0, "Evidence candidates: fragment tolerance, ppm either "
                          "side", false, advanced_("search:prefilter_ppm"));
    setMinFloat_("search:prefilter_ppm", 0.1); setMaxFloat_("search:prefilter_ppm", 100.0);
    registerStringOption_("search:decoys", "<method>", "shuffle", "How the search builds its in-memory decoys (both keep "
                          "the termini); decoys in the library are not searched and stay in the output", false,
                          advanced_("search:decoys"));
    setValidStrings_("search:decoys", {"shuffle", "pseudo_reverse"});
    registerStringOption_("search:intensities", "<source>", "predicted", "Validation only: fragments of both pair "
                          "members. predicted = predicted from each member's own sequence by the MS2 model (-models); "
                          "library = the target's assay, the decoy in its target's fragment slots (weaker decoys)",
                          false, advanced_("search:intensities"));
    setValidStrings_("search:intensities", {"predicted", "library"});
    registerStringOption_("search:instrument", "<name>", "auto", "Predicted intensities: the MS2 model's instrument; auto = "
                          "the one the library records, else timsTOF for an ion-mobility run and QE otherwise", false,
                          advanced_("search:instrument"));
    registerDoubleOption_("search:nce", "<nce>", -1.0, "Predicted intensities: the MS2 model's collision energy; -1 = the "
                          "one the library records with its instrument, else the instrument's", false,
                          advanced_("search:nce"));
    setMinFloat_("search:nce", -1.0); setMaxFloat_("search:nce", 100.0);
    registerIntOption_("search:seed", "<n>", 42, "Salt of the candidate draw: changes which pairs are searched, not how",
                       false, advanced_("search:seed"));
    registerIntOption_("search:passes", "<n>", 1, deprecated("nothing (1 is its only value)"), false, advanced_("search:passes"));
    setMinInt_("search:passes", 1); setMaxInt_("search:passes", 1);
    registerDoubleOption_("search:rt_window", "<s>", 0.0, "Full RT extraction window, seconds (0 = from the calibration)",
                          false, advanced_("search:rt_window"));
    setMinFloat_("search:rt_window", 0.0);
    registerDoubleOption_("search:mz_ppm", "<ppm>", 0.0, "Full fragment m/z extraction window, ppm (0 = automatic)", false,
                          advanced_("search:mz_ppm"));
    setMinFloat_("search:mz_ppm", 0.0);
    registerDoubleOption_("search:im_window", "<1/K0>", 0.0, "Full 1/K0 extraction window on ion-mobility runs: 0 = from "
                          "the run's 1/K0 calibration; > 0 = this width; -1 = off (no observed 1/K0: -write_im and the "
                          "CCS head are refused)", false, advanced_("search:im_window"));
    setMinFloat_("search:im_window", -1.0);
    registerStringOption_("search:ms1", "<true/false>", "true", "Extract MS1 traces and use the MS1 sub-scores", false,
                          advanced_("search:ms1"));
    setValidStrings_("search:ms1", {"true", "false"});
    registerStringOption_("search:rt_im_scores", "<true/false>", "true", "Validation only: false keeps the RT and 1/K0 "
                          "deviation sub-scores out of the classifier (an ablation)", false, advanced_("search:rt_im_scores"));
    setValidStrings_("search:rt_im_scores", {"true", "false"});
    registerDoubleOption_("search:calibration_min_rsq", "<r2>", 0.70, "Least r^2 of the RT calibration on the seed assays",
                          false, advanced_("search:calibration_min_rsq"));
    registerDoubleOption_("search:calibration_min_coverage", "<f>", 0.30, "Least fraction of the seed assays found in the "
                          "run that the RT outlier removal must keep", false, advanced_("search:calibration_min_coverage"));
    for (const char* name : {"search:calibration_min_rsq", "search:calibration_min_coverage"})
    { setMinFloat_(name, 0.0); setMaxFloat_(name, 1.0); }
    flag("search:allow_bootstrap", deprecated("ODIA_SEARCH_ALLOW_BOOTSTRAP=1"));
    registerStringOption_("search:readoptions", "<mode>", "auto", "How the run is held: normal = in memory, cache = "
                          "per-window cache files, auto = cache for runs above 3 GB", false, advanced_("search:readoptions"));
    setValidStrings_("search:readoptions", {"auto", "normal", "cache"});
    registerStringOption_("search:cache_dir", "<dir>", "", "Directory of the cache files (default: that of -out_ids); "
                          "removed after the search", false, advanced_("search:cache_dir"));
    registerStringOption_("search:ms1_in_memory", "<mode>", "auto", "Memory: hold the MS1 map in memory during "
                          "extraction (auto = when the free memory is at least twice its size; on; off). Changes speed "
                          "and memory, never the report", false, advanced_("search:ms1_in_memory"));
    setValidStrings_("search:ms1_in_memory", {"auto", "on", "off"});
    registerIntOption_("search:min_ids", "<n>", 200, "Abort when fewer target precursors pass q <= 0.01", false,
                       advanced_("search:min_ids"));
    registerDoubleOption_("search:max_target_fraction", "<f>", 0.5, "Abort when more than this fraction of the scored "
                          "targets passes q <= 0.01 (no honest decoy set looks like that)", false,
                          advanced_("search:max_target_fraction"));
    setMinFloat_("search:max_target_fraction", 0.01); setMaxFloat_("search:max_target_fraction", 1.0);
    registerDoubleOption_("search:report_max_q", "<q>", 0.10, "Precursors (targets and decoys) up to this q-value go "
                          "into -out_ids", false, advanced_("search:report_max_q"));
    setMinFloat_("search:report_max_q", 0.01); setMaxFloat_("search:report_max_q", 1.0);
    registerStringOption_("search:entrapment_tag", "<prefix>", "", "Validation only: prefix of entrapment protein ids "
                          "(at the start or after a '|'); logs and records the entrapment FDP estimate and winner test",
                          false, advanced_("search:entrapment_tag"));
    registerStringOption_("search:selftest", "<true/false>", "true", "Also score with swapped and with random pair labels "
                          "and abort unless both identify (almost) nothing (catches a classifier that leaks labels)",
                          false, advanced_("search:selftest"));
    setValidStrings_("search:selftest", {"true", "false"});
    registerIntOption_("search:batch_size", "<n>", 100, "Memory: precursors per extraction batch within one isolation "
                       "window (the report does not depend on it)", false, advanced_("search:batch_size"));
    registerIntOption_("search:chunk", "<n>", 20000, "Memory: precursors per extraction call (pairs stay together); "
                       "bounds extraction memory, not the calibration's", false, advanced_("search:chunk"));
    for (const char* name : {"search:subset", "search:max_pairs", "search:seed", "search:min_ids"}) { setMinInt_(name, 0); }
    setMinInt_("search:batch_size", 1);
    setMinInt_("search:chunk", 2);

    for (const char* name : {"min_fragments", "intensity_min_fragments", "tune_predict_sessions", "cohort:train_size",
                            "train:warmup", "stop:min_epochs", "stop:patience", "machine:seed"}) { setMinInt_(name, 0); }
    for (const char* name : {"train:epochs", "train:batch_size", "stop:eval_every"}) { setMinInt_(name, 1); }
    setMinInt_("machine:threads", 0);
    setMinInt_("machine:train_threads", 0);
    // Flat refinement options follow the existing effective-config schema;
    // training options have their own prefixes and are recognized in main_.
    const std::map<std::string, std::string> flags = {{"filter", "no_filter"}, {"require_gates", "empirical_library"},
      {"write_rt", "no_write_rt"}, {"intensity_restrict", "intensity_no_restrict"}, {"intensity_rerank", "intensity_no_rerank"}};
    const json defaults = effectiveConfig(ODIA::RefineParams{});
    for (const auto& [key, value] : defaults.items())
    {
      if (key == "schema_version") { continue; }
      const auto flag = flags.find(key);
      refinement_options_.insert(flag == flags.end() ? key : flag->second);
    }
    refinement_options_.insert({"ids", "out_report", "run", "out_ids"});
  }

ODIA::search::SearchParams DIALibGen::searchParams_()
  {
    auto count = [this](const char* name) { return static_cast<std::size_t>(std::max(0, getIntOption_(name))); };
    ODIA::search::SearchParams s;
    s.candidates = getStringOption_("search:candidates");
    s.subset = count("search:subset");
    s.max_pairs = count("search:max_pairs");
    s.prefilter_depth = getIntOption_("search:prefilter_depth");
    s.prefilter_top_peaks = count("search:prefilter_top_peaks");
    s.prefilter_ppm = getDoubleOption_("search:prefilter_ppm");
    s.decoys = ODIA::search::parseSearchDecoyMethod(getStringOption_("search:decoys"));
    s.intensities = ODIA::search::parseIntensities(getStringOption_("search:intensities"));
    if (getStringOption_("search:instrument") != "auto")
    {
      s.instrument = getStringOption_("search:instrument");
      s.instrument_source = "search:instrument";
    }
    s.nce = getDoubleOption_("search:nce");
    s.seed = static_cast<std::uint64_t>(count("search:seed"));
    s.passes = getIntOption_("search:passes");   // deprecated: 1, its only value
    s.rt_window = getDoubleOption_("search:rt_window");
    s.mz_ppm = getDoubleOption_("search:mz_ppm");
    s.im_window = getDoubleOption_("search:im_window");
    s.ms1 = getStringOption_("search:ms1") == "true";
    s.rt_im_scores = getStringOption_("search:rt_im_scores") == "true";
    s.batch_size = count("search:batch_size");
    s.chunk = count("search:chunk");
    s.calibration_min_rsq = getDoubleOption_("search:calibration_min_rsq");
    s.calibration_min_coverage = getDoubleOption_("search:calibration_min_coverage");
    // ODIA_SEARCH_ALLOW_BOOTSTRAP=1, or the deprecated -search:allow_bootstrap.
    if (const char* e = std::getenv("ODIA_SEARCH_ALLOW_BOOTSTRAP"); e != nullptr && std::string(e) == "1") { s.allow_bootstrap = true; }
    s.allow_bootstrap = s.allow_bootstrap || getFlag_("search:allow_bootstrap");
    s.readoptions = ODIA::search::parseReadMode(getStringOption_("search:readoptions"));
    s.cache_dir = getStringOption_("search:cache_dir");
    s.ms1_in_memory = ODIA::search::parseMs1InMemory(getStringOption_("search:ms1_in_memory"));
    s.min_ids = count("search:min_ids");
    s.max_target_fraction = getDoubleOption_("search:max_target_fraction");
    s.report_max_q = getDoubleOption_("search:report_max_q");
    s.entrapment_tag = getStringOption_("search:entrapment_tag");
    s.selftest = getStringOption_("search:selftest") == "true";
    s.threads = ODIA::search::SearchParams::resolveThreads(getIntOption_("threads"));
    s.validate();
    return s;
  }


namespace
{
  /// The training settings a recipe's "tuning" object holds (and -write_config writes).
  const std::vector<std::string>& trainingSettings()
  {
    static const std::vector<std::string> names = {
      "tune_heads", "tune_keep_free_cysteine_offset",
      "filter:q_value", "filter:min_charge", "filter:allow_z1", "filter:rt_spread_max", "filter:rt_max_minutes",
      "cohort:train_size", "cohort:train_frac", "cohort:full_fit", "cohort:no_inner_val",
      "train:epochs", "train:warmup", "train:lr", "train:batch_size",
      "stop:eval_every", "stop:min_epochs", "stop:patience", "stop:rel_tol", "stop:abs_tol", "stop:max_seconds",
      "stop:select", "machine:device", "machine:tf32", "machine:train_threads", "machine:seed", "machine:no_cudnn"};
    return names;
  }
}

OpenMS::ParamValue DIALibGen::trainingSetting_(const std::string& name) const
{
  // The command line, a deprecated name's value, the INI, the recipe, the default.
  if (argv_supplied_.count(name)) { return getParam_().getValue(name); }
  if (const auto alias = aliased_.find(name); alias != aliased_.end()) { return alias->second; }
  if (supplied_.count(name) || !tuning_recipe_.contains(name)) { return getParam_().getValue(name); }
  const json& v = tuning_recipe_.at(name);
  const auto& info = findEntry_(name);
  using PI = OpenMS::ParameterInformation;
  if (info.type == PI::FLAG)
  {
    if (!v.is_boolean()) { throw std::runtime_error("tuning." + name + " must be true or false"); }
    return OpenMS::ParamValue(v.get<bool>() ? "true" : "false");
  }
  if (info.type == PI::INT)
  {
    if (!v.is_number_integer()) { throw std::runtime_error("tuning." + name + " must be an integer"); }
    const long long x = v.get<long long>();
    if (x < info.min_int || x > info.max_int)
    { throw std::runtime_error("tuning." + name + " must be between " + std::to_string(info.min_int) + " and " + std::to_string(info.max_int)); }
    return OpenMS::ParamValue(static_cast<int>(x));
  }
  if (info.type == PI::DOUBLE)
  {
    if (!v.is_number() || !std::isfinite(v.get<double>()) || v.get<double>() < info.min_float || v.get<double>() > info.max_float)
    { throw std::runtime_error("tuning." + name + " must be a number within its range (see --helphelp)"); }
    return OpenMS::ParamValue(v.get<double>());
  }
  if (!v.is_string()) { throw std::runtime_error("tuning." + name + " must be a string"); }
  const std::string text = v.get<std::string>();
  if (!info.valid_strings.empty() && std::find(info.valid_strings.begin(), info.valid_strings.end(), text) == info.valid_strings.end())
  { throw std::runtime_error("tuning." + name + ": '" + text + "' is not one of its values"); }
  return OpenMS::ParamValue(text);
}

json DIALibGen::trainingRecipe_() const
{
  json out = json::object();
  for (const auto& name : trainingSettings())
  {
    const OpenMS::ParamValue v = trainingSetting_(name);
    const auto type = findEntry_(name).type;
    using PI = OpenMS::ParameterInformation;
    if (type == PI::FLAG) { out[name] = v.toString() == "true"; }
    else if (type == PI::INT) { out[name] = static_cast<int>(v); }
    else if (type == PI::DOUBLE) { out[name] = static_cast<double>(v); }
    else { out[name] = v.toString(); }
  }
  return out;
}

DIALibGen::ExitCodes DIALibGen::refine_(bool tune_only)
  {
    ODIA::RefineParams p;
    p.q_precursor = getDoubleOption_("q_precursor");
    p.q_global = getDoubleOption_("q_global");
    p.q_protein = getDoubleOption_("q_protein");
    p.min_fragments = static_cast<std::size_t>(std::max(0, getIntOption_("min_fragments")));
    p.filter = !getFlag_("no_filter");
    p.require_gates = !getFlag_("empirical_library");
    p.write_rt = !getFlag_("no_write_rt");
    p.write_im = getFlag_("write_im");
    p.write_intensity = getFlag_("write_intensity");
    p.intensity_min_correlation = getDoubleOption_("intensity_min_correlation");
    p.intensity_restrict = !getFlag_("intensity_no_restrict");
    p.intensity_rerank = !getFlag_("intensity_no_rerank");
    p.intensity_min_fragments = static_cast<std::size_t>(std::max(0, getIntOption_("intensity_min_fragments")));
    p.intensity_norm = intensityNormFrom(getStringOption_("intensity_norm"));
    p.intensity_min_relative = getDoubleOption_("intensity_min_relative");
    p.intensity_mz_tol_ppm = getDoubleOption_("intensity_mz_tol_ppm");
    p.intensity_max_mz_mismatch = getDoubleOption_("intensity_max_mz_mismatch");
    p.allow_mixed_intensity = getFlag_("allow_mixed_intensity");
    p.im_min_charge = getIntOption_("im_min_charge");
    p.im_ramp_top = getDoubleOption_("im_ramp_top");
    p.im_ramp_margin = getDoubleOption_("im_ramp_margin");
    p.min_match_fraction = getDoubleOption_("min_match_fraction");
    p.rt_unit = getStringOption_("rt_unit") == "minmax" ? ODIA::RefineParams::RtUnit::MinMax : ODIA::RefineParams::RtUnit::Observed;
    p.dedup = getStringOption_("dedup") == "highest_evidence" ? ODIA::RefineParams::Dedup::HighestEvidence : ODIA::RefineParams::Dedup::LowestQ;

    const bool tune = tune_only || getFlag_("tune");
    if (tune_only) { p.filter = false; p.write_rt = false; }
    const json cli = effectiveConfig(p);

    if (const std::string cfg = getStringOption_("config"); !cfg.empty())
    {
      std::ifstream in(cfg);
      if (!in) { writeLogError_("cannot read config: " + cfg); return INPUT_FILE_NOT_FOUND; }
      try
      {
        const json recipe = json::parse(in, nullptr, true, true);
        applyJson_(recipe, p);
        if (recipe.contains("tuning"))
        {
          const json& t = recipe["tuning"];
          if (!t.is_object()) { throw std::runtime_error("tuning must be a JSON object"); }
          for (const auto& [k, v] : t.items())
          {
            if (std::find(trainingSettings().begin(), trainingSettings().end(), k) == trainingSettings().end())
            { throw std::runtime_error("unknown config key: tuning." + k); }
            (void)v;
          }
          tuning_recipe_ = t;
          if (!tune) { warn_("the recipe's tuning settings have no effect without -tune"); }
        }
        // Explicit TOPP options (CLI or INI) override the JSON compatibility config.
        json overrides = json::object();
        const std::map<std::string, std::string> flags = {{"filter", "no_filter"}, {"require_gates", "empirical_library"},
          {"write_rt", "no_write_rt"}, {"intensity_restrict", "intensity_no_restrict"}, {"intensity_rerank", "intensity_no_rerank"}};
        for (const auto& [key, value] : cli.items())
        {
          const auto f = flags.find(key);
          if (supplied_.count(f == flags.end() ? key : f->second)) { overrides[key] = value; }
        }
        applyJson_(overrides, p);
      }
      catch (const std::exception& e) { writeLogError_(std::string("config: ") + e.what()); return ILLEGAL_PARAMETERS; }
    }
    try { applyJson_(effectiveConfig(p), p); }
    catch (const std::exception& e) { writeLogError_(e.what()); return ILLEGAL_PARAMETERS; }
    if (tune_only && (p.filter || p.write_rt || p.write_im || p.write_intensity))
    { writeLogError_("-mode tune preserves the whole library and writes predictions only; filtering or observed-value replacement requires -mode refine"); return ILLEGAL_PARAMETERS; }
    if (tune && p.rt_unit == ODIA::RefineParams::RtUnit::MinMax)
    { writeLogError_("tuning and rt_unit=minmax are incompatible: predictions use the reference run's minutes"); return ILLEGAL_PARAMETERS; }
    std::string heads;
    try { heads = trainingSetting_("tune_heads").toString(); }
    catch (const std::exception& e) { writeLogError_(std::string("config: ") + e.what()); return ILLEGAL_PARAMETERS; }
    if (p.write_rt && !p.filter && (!tune || heads == "ccs"))
    { writeLogError_("observed RT with the filter off would mix reference-run minutes with library predictions; use -no_write_rt or -tune with the RT head (rt or both)"); return ILLEGAL_PARAMETERS; }
    if (tune_only)
    {
      ODIA::RefineParams defaults;
      defaults.filter = defaults.write_rt = false;
      const json expected = effectiveConfig(defaults), actual = effectiveConfig(p);
      for (const auto& [key, value] : actual.items())
      {
        // An INI from before 0.12 carries the gates' old default, 0.01.
        if ((key == "q_global" || key == "q_protein") && old_ini_ && supplied_.count(key) && !argv_supplied_.count(key) &&
            value == 0.01) { continue; }
        if (value != expected.at(key))
        { writeLogError_(key + " has no effect in -mode tune; use the training options or select -mode refine"); return ILLEGAL_PARAMETERS; }
      }
    }

    // Settings without an effect on this refinement: a WARNING in 0.12,
    // refused from 0.13. Only values other than the default count.
    if (!tune_only)
    {
      const ODIA::RefineParams d;
      const auto inert = [this](const char* option, const char* needs)
      { warn_(std::string("-") + option + " has no effect without -" + needs + "; remove it (refused from 0.13)"); };
      if (!p.write_intensity)
      {
        if (p.intensity_min_correlation != d.intensity_min_correlation) { inert("intensity_min_correlation", "write_intensity"); }
        if (p.intensity_restrict != d.intensity_restrict) { inert("intensity_no_restrict", "write_intensity"); }
        if (p.intensity_rerank != d.intensity_rerank) { inert("intensity_no_rerank", "write_intensity"); }
        if (p.intensity_min_fragments != d.intensity_min_fragments) { inert("intensity_min_fragments", "write_intensity"); }
        if (p.intensity_norm != d.intensity_norm) { inert("intensity_norm", "write_intensity"); }
        if (p.intensity_min_relative != d.intensity_min_relative) { inert("intensity_min_relative", "write_intensity"); }
        if (p.intensity_mz_tol_ppm != d.intensity_mz_tol_ppm) { inert("intensity_mz_tol_ppm", "write_intensity"); }
        if (p.intensity_max_mz_mismatch != d.intensity_max_mz_mismatch) { inert("intensity_max_mz_mismatch", "write_intensity"); }
        if (p.allow_mixed_intensity != d.allow_mixed_intensity) { inert("allow_mixed_intensity", "write_intensity"); }
      }
      if (!p.write_im)
      {
        if (p.im_min_charge != d.im_min_charge) { inert("im_min_charge", "write_im"); }
        if (p.im_ramp_top != d.im_ramp_top) { inert("im_ramp_top", "write_im"); }
        if (p.im_ramp_margin != d.im_ramp_margin) { inert("im_ramp_margin", "write_im"); }
      }
    }

    json eff = effectiveConfig(p);
    if (tune)
    {
      try { eff["tuning"] = trainingRecipe_(); }
      catch (const std::exception& e) { writeLogError_(std::string("config: ") + e.what()); return ILLEGAL_PARAMETERS; }
    }
    if (const std::string wc = getStringOption_("write_config"); !wc.empty())
    {
      try
      {
        ODIA::AtomicFile staged(wc);
        std::ofstream o(staged.temporaryPath());
        o << eff.dump(2) << '\n'; o.close();
        if (!o) { throw std::runtime_error("cannot write " + wc); }
        staged.commit();
      }
      catch (const std::exception& e) { writeLogError_(e.what()); return CANNOT_WRITE_OUTPUT_FILE; }
      writeLogInfo_("wrote effective config to " + wc);
      return EXECUTION_OK;
    }

    // A directory given as an input file -- a Bruker .d above all -- crashes
    // TOPP's file-type probe in getStringOption_ (an uncaught
    // std::ios_base::failure). Refuse it before any getter can run that probe.
    for (const std::string name : {"run", "ids"})
    {
      const std::string raw = getParam_().getValue(name).toString();
      std::error_code ec;
      if (raw.empty() || !std::filesystem::is_directory(raw, ec)) { continue; }
      std::string stem = raw;
      while (stem.size() > 1 && (stem.back() == '/' || stem.back() == '\\')) { stem.pop_back(); }
      const bool bruker = stem.size() > 2 && (stem.ends_with(".d") || stem.ends_with(".D"));
      // An unpacked mzPeak archive is a directory holding mzpeak_index.json.
      if (name == "run" && ODIA::search::mzpeak::isMzPeak(raw)) { continue; }
      if (name == "run" && bruker)
      {
        writeLogError_("-run " + raw + " is a Bruker .d directory, which this version cannot read; convert it to a "
                       "per-window mzML with a per-peak 1/K0 array (mzpeak-convert --to mzml) or an mzPeak archive (mzpeak-convert) and pass that");
      }
      else { writeLogError_("-" + name + " " + raw + " is a directory; it takes a file" + (name == "run" ? " (mzML)" : " (report.parquet)")); }
      return ILLEGAL_PARAMETERS;
    }
    const std::string in = getStringOption_("in");
    std::string ids = getStringOption_("ids");
    // Read as given: TOPP's input-file getter probes the type, which OpenMS 3.5
    // does not know for mzPeak and which fails on a directory.
    const std::string run = getParam_().getValue("run").toString();
    if (!run.empty())
    {
      std::error_code ec;
      std::string name = std::filesystem::path(run).filename().string();
      std::transform(name.begin(), name.end(), name.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
      const bool mzpeak = ODIA::search::mzpeak::isMzPeak(run);
      if (!mzpeak && !name.ends_with(".mzml") && !name.ends_with(".mzml.gz"))
      { writeLogError_("-run " + run + ": give an mzML (.mzML) or an mzPeak archive (.mzpeak)"); return ILLEGAL_PARAMETERS; }
      if (!std::filesystem::exists(run, ec))
      { writeLogError_("-run " + run + " does not exist"); return INPUT_FILE_NOT_FOUND; }
      if (mzpeak && !ODIA::search::mzpeak::available())
      {
        writeLogError_("-run " + run + " is an mzPeak archive, and this build reads mzML only (DIALIBGEN_WITH_MZPEAK off); "
                       "convert it with mzpeak-convert --to mzml");
        return ILLEGAL_PARAMETERS;
      }
    }
    // TOPP's writable-file probe follows and then removes dangling symlinks.
    // Inspect destinations before any output-file getter can run that probe.
    const std::string out = getParam_().getValue("out").toString();
    if (in.empty() || out.empty() || (ids.empty() && run.empty()))
    { writeLogError_("-in, -ids and -out are required (-run can take the place of -ids)"); return ILLEGAL_PARAMETERS; }
    if (!ids.empty() && !run.empty())
    { writeLogError_("give exactly one of -ids and -run: -run writes its own identification report and uses it as -ids"); return ILLEGAL_PARAMETERS; }
#ifndef DIALIBGEN_WITH_FINETUNE
    if (tune)
    { writeLogError_("this build has no fine-tuning stage; configure DIALIBGEN_BUILD_FINETUNE with libtorch"); return ILLEGAL_PARAMETERS; }
#endif
    ODIA::LibraryFormat format;
    try { format = ODIA::LibraryExport::resolve(out, getStringOption_("out_format")); }
    catch (const std::invalid_argument& e) { writeLogError_(e.what()); return ILLEGAL_PARAMETERS; }

    // The built-in identification step (-run): refusals, settings, report path.
    std::string out_ids;
    std::unique_ptr<ODIA::search::SearchParams> search;
    if (!run.empty())
    {
      if (p.write_intensity)
      { writeLogError_("-write_intensity is not available with -run yet: the built-in search does not measure fragment intensities"); return ILLEGAL_PARAMETERS; }
      if (!p.require_gates)
      { writeLogError_("-empirical_library does not apply to -run: the built-in report carries every gate column"); return ILLEGAL_PARAMETERS; }
      if (p.min_fragments > 0)
      { writeLogError_("-min_fragments is not available with -run yet: the built-in report carries no fragment identities"); return ILLEGAL_PARAMETERS; }
      // Observed 1/K0 comes only from an ion-mobility (diaPASEF) run searched
      // with its ion mobility: whether the run has it is known once it is
      // loaded, which the search checks before searching anything
      // (SearchParams::require_ion_mobility); search:im_window -1 is refused here.
      std::vector<std::string> needs_im;
      if (p.write_im) { needs_im.push_back("-write_im"); }
      // auto settles on the heads once the report is there (no ion mobility: RT only).
      if (tune && heads != "rt" && heads != "auto") { needs_im.push_back("-tune_heads " + heads); }
      std::string needs_im_text;
      for (const auto& n : needs_im) { needs_im_text += (needs_im_text.empty() ? "" : " and ") + n; }
      if (!needs_im.empty() && getDoubleOption_("search:im_window") == -1.0)
      {
        writeLogError_(needs_im_text + " needs observed 1/K0 values, and search:im_window -1 searches without ion mobility, so "
                       "the built-in search measures none" + std::string(tune ? "; use -tune_heads rt" : ""));
        return ILLEGAL_PARAMETERS;
      }
      // Everything that would fail AFTER the search (minutes on a full run)
      // is checked before it.
      if (tune)
      {
        const int epochs = static_cast<int>(trainingSetting_("train:epochs")), warmup = static_cast<int>(trainingSetting_("train:warmup"));
        if (epochs < 1 || warmup < 0 || warmup > epochs)
        { writeLogError_("train:warmup must be in [0, train:epochs]"); return ILLEGAL_PARAMETERS; }
        if (static_cast<int>(trainingSetting_("cohort:train_size")) > 0 && static_cast<double>(trainingSetting_("cohort:train_frac")) > 0)
        { writeLogError_("give cohort:train_size or cohort:train_frac, not both"); return ILLEGAL_PARAMETERS; }
#ifdef DIALIBGEN_WITH_FINETUNE
        const std::string models = modelDirectory_();
        if (models.empty())
        { writeLogError_("-tune needs -models (or $DIALIBGEN_MODEL_DIR): the stock peptdeep models to start from"); return ILLEGAL_PARAMETERS; }
        for (const char* head : {"rt", "ccs"})
        {
          // auto needs the RT model; it tunes CCS only when that model is there.
          if (std::string(head) == "rt" ? heads == "ccs" : heads != "ccs" && heads != "both") { continue; }
          const std::string file = std::string("peptdeep_") + head + "_dynamic.onnx";
          if (!std::filesystem::exists(std::filesystem::path(models) / file))
          { writeLogError_("no " + file + " in " + models); return ILLEGAL_PARAMETERS; }
        }
#endif
      }
      try { search = std::make_unique<ODIA::search::SearchParams>(searchParams_()); }
      catch (const std::exception& e) { writeLogError_(e.what()); return ILLEGAL_PARAMETERS; }
      if (getIntOption_("threads") <= 0)
      {
        // -threads 0 is the effective CPU budget, as generate counts it: the
        // search's prediction sessions, its OpenMP loops (TOPPBase passes 0 to
        // omp_set_num_threads, which leaves OpenMP at every logical CPU) and
        // the classifier. OpenSWATH's own team stays capped at
        // SearchParams::openswathThreads() (16 with the nested window loop).
        writeLogInfo_("-threads 0: " + std::to_string(search->threads) + " threads (" + ODIA::physicalCoresReason() + ")");
#ifdef _OPENMP
        omp_set_num_threads(search->threads);
#endif
      }
      search->require_ion_mobility = needs_im_text;
      out_ids = getParam_().getValue("out_ids").toString();
      if (out_ids.empty()) { out_ids = out + ".ids.parquet"; }
      if (!out_ids.ends_with(".parquet"))
      { writeLogError_("-out_ids must end in .parquet"); return ILLEGAL_PARAMETERS; }
      if (search->intensities == ODIA::search::Intensities::Predicted)
      {
        // The MS2 model both members of every pair are predicted with: where
        // -tune finds its stock models.
        const std::string models = modelDirectory_();
        const std::filesystem::path ms2 = models.empty() ? std::filesystem::path() : std::filesystem::path(models) / "peptdeep_ms2_dynamic.onnx";
        if (models.empty() || !std::filesystem::exists(ms2))
        {
          writeLogError_("-run predicts the fragments of every target and decoy it searches (search:intensities predicted) "
                         "and needs peptdeep_ms2_dynamic.onnx: " +
                         (models.empty() ? std::string("give -models or $DIALIBGEN_MODEL_DIR") : "none in " + models));
          return ILLEGAL_PARAMETERS;
        }
        search->ms2_model = ms2.string();
        search->ms2_model_hash = ODIA::DIANNLibraryFile::hashFile(search->ms2_model);
        if (!search->instrument.empty() && ODIA::PeptDeepEncoder::canonicalInstrument(search->instrument).empty())
        {
          writeLogError_("search:instrument: unknown instrument '" + search->instrument +
                         "' (QE, Lumos, timsTOF, SciexTOF, ThermoTOF, or auto)");
          return ILLEGAL_PARAMETERS;
        }
        // The instrument and NCE the library was generated with, when
        // DIALibGen generate recorded them: the targets are then predicted
        // exactly as the library predicted them.
        if (search->instrument.empty())
        {
          try
          {
            const std::string config = ODIA::DIANNLibraryFile::readFingerprint(in, "odia.config_json");
            if (!config.empty())
            {
              const json j = json::parse(config);
              // 0.12 records the instrument as named and the MS2 model's slot.
              const char* key = j.contains("model_instrument") ? "model_instrument" : "instrument";
              if (j.contains(key) && j[key].is_string() && !j[key].get<std::string>().empty())
              {
                search->instrument = j[key].get<std::string>();
                search->instrument_source = "library";
                if (search->nce <= 0 && j.contains("nce") && j["nce"].is_number() && j["nce"].get<double>() > 0)
                { search->nce = j["nce"].get<double>(); }
              }
            }
          }
          catch (const std::exception&) {}   // not a DIALibGen library: the run decides
        }
      }
    }

    const std::string report = getParam_().getValue("out_report").toString();
    std::set<std::filesystem::path> destinations;
    for (const auto& file : {out, out + ".refine.json", report, out_ids})
    {
      if (file.empty()) { continue; }
      if (std::filesystem::exists(file) || std::filesystem::is_symlink(file))
      { writeLogError_("refusing to overwrite existing output: " + file); return CANNOT_WRITE_OUTPUT_FILE; }
      if (!destinations.insert(std::filesystem::weakly_canonical(file)).second)
      { writeLogError_("output paths must be distinct: " + file); return ILLEGAL_PARAMETERS; }
    }
    (void)getStringOption_("out");
    (void)getStringOption_("out_report");
    if (!run.empty()) { (void)getStringOption_("out_ids"); }

    ODIA::Library library;
    ODIA::RefineStats st;
    json tune_prov = json::object();
    json search_prov;
    std::string run_hash;
    // The run's and the library's content hashes (FNV-1a, one byte at a time,
    // ~0.9 GB/s: 36 s for a 30 GB diaPASEF mzML) run on threads of their own
    // from the start, beside the library load and the search, instead of on
    // the main thread afterwards (SearchParams::asyncIo). get() waits for each
    // where its value is needed and rethrows what went wrong there.
    std::future<std::string> run_hash_ahead, library_hash_ahead;
    if (ODIA::search::SearchParams::asyncIo())
    {
      auto hashing = [](const char* name, std::string path, bool run_hash) {
        std::string value;
        {
          ODIA::trace::Scope hash_step(name, 1);
          value = run_hash ? runHash(path) : ODIA::DIANNLibraryFile::hashFile(path);
        }
        ODIA::trace::flushThread();
        return value;
      };
      try
      {
        if (!run.empty()) { run_hash_ahead = std::async(std::launch::async, hashing, "run_hash_bg", run, true); }
        library_hash_ahead = std::async(std::launch::async, hashing, "library_hash_bg", in, false);
      }
      catch (const std::system_error&) {}   // no thread to spare: hashed in order below
    }
    // The search's run cache, removed on a thread while refine and tuning go
    // on; waited for before the provenance is written.
    std::shared_ptr<ODIA::search::CacheRemoval> cache_removal;
    // Once the search has written its report, a later failure must say that the
    // report is still there: -out_ids is never overwritten, so the same command
    // would now be refused.
    bool report_written = false;
    auto keptReport = [&]() {
      if (report_written)
      {
        writeLogError_("the identification report " + out_ids + " was kept; rerun with -ids " + out_ids +
                       " instead of -run to reuse it, or remove it");
      }
    };
    try
    {
      ODIA::trace::Scope load_phase("library_load", 0);
      ODIA::DIANNLibraryFile::load(in, library, static_cast<unsigned>(std::max(0, getIntOption_("threads"))));
      load_phase.end();
      writeLogInfo_("library: " + std::to_string(library.precursorCount()) + " precursors, " +
                    std::to_string(library.transitionCount()) + " transitions");
      if (search)
      {
        // Protein groups: tuning builds its held-out cohorts from them, and the
        // report's PG.Q.Value (the -q_protein gate) is 1 without them. Checked
        // before the search, not after it.
        std::size_t targets = 0, grouped = 0;
        for (std::size_t i = 0; i < library.precursorCount(); ++i)
        {
          if (library.precursors().decoy[i]) { continue; }
          ++targets;
          grouped += library.strings().get(library.precursors().protein_group[i]).empty() ? 0 : 1;
        }
        if (grouped == 0 && (tune || p.q_protein < 1.0))
        {
          throw std::runtime_error("none of the library's " + std::to_string(targets) + " target precursors has a Protein.Group; " +
                                   (tune ? std::string("tuning needs protein groups for its held-out cohorts")
                                         : std::string("every PG.Q.Value would be 1 and -q_protein would reject everything "
                                                       "(pass -q_protein 1 to refine without the protein gate)")));
        }
        writeLogWarn_("-run is EXPERIMENTAL: the built-in identification has not passed its entrapment validation yet");
        if (p.filter)
        { writeLogWarn_("-mode refine with -run keeps only precursors identified among the at most search:max_pairs pairs "
                        "searched; the built-in workflow is -mode refine -tune -no_filter"); }
        ODIA::search::Identifier identifier(*search, [this](const std::string& m) { writeLogInfo_(m); },
                                            [this](const std::string& m) { writeLogWarn_(m); });
        const ODIA::search::IdentificationResult found = identifier.identify(library, run, out_ids, DIALIBGEN_VERSION);
        report_written = true;
        search_prov = json::parse(found.provenance_json);
        cache_removal = found.cache_removal;
        {
          ODIA::trace::Scope hash_phase("run_hash", 0);
          run_hash = run_hash_ahead.valid() ? run_hash_ahead.get() : runHash(run);
        }
        // The hand-off: the report just written IS the reference from here on,
        // read, gated and hashed by exactly the code that reads -ids.
        ids = out_ids;
        writeLogInfo_("identification report " + out_ids + " serves as -ids: " + std::to_string(found.identified) +
                      " target precursors at q <= 0.01");
      }
      ODIA::LibraryRefiner::ObsMap obs;
      if (!tune_only)
      {
        ODIA::trace::Scope read_phase("apply", 0);
        obs = ODIA::LibraryRefiner::readObservations(ids, p, st);
      }

#ifdef DIALIBGEN_WITH_FINETUNE
      if (tune)
      {
        namespace fs = std::filesystem;
        const std::string models = modelDirectory_();
        if (models.empty())
        { throw std::runtime_error("-tune needs -models (or $DIALIBGEN_MODEL_DIR): the stock peptdeep models to start from"); }

        // A foreign library can carry mass-only modifications that the model
        // cannot encode. Reject these before training either head.
        ODIA::trace::Scope check_phase("tune_check", 0);
        std::unordered_set<std::uint32_t> checked;
        for (const auto sequence_id : library.precursors().modified_sequence)
        {
          if (!checked.insert(sequence_id).second) { continue; }
          const std::string sequence(library.strings().get(sequence_id));
          try { (void)ODIA::PeptDeepEncoder::encode(OpenMS::AASequence::fromString(sequence)); }
          catch (const std::exception& error)
          { throw std::runtime_error("tuning cannot encode library precursor '" + sequence + "': " + error.what()); }
        }

        check_phase.end();
        const std::string keep = getStringOption_("tune_out_models");
        struct WorkDirectory
        {
          fs::path path;
          bool temporary = false;
          ~WorkDirectory() { if (temporary) { std::error_code ec; fs::remove_all(path, ec); } }
        } directory;
        if (keep.empty())
        {
          std::random_device random;
          for (int attempt = 0; attempt < 16; ++attempt)
          {
            directory.path = fs::temp_directory_path() / ("dialibgen-tune-" + std::to_string(random()) + "-" + std::to_string(random()));
#ifdef _WIN32
            const bool created = fs::create_directory(directory.path);
#else
            const bool created = ::mkdir(directory.path.c_str(), 0700) == 0;
            if (!created && errno != EEXIST)
            { throw std::runtime_error("cannot create private tuning directory: " + std::error_code(errno, std::generic_category()).message()); }
#endif
            if (created) { directory.temporary = true; break; }
          }
          if (!directory.temporary) { throw std::runtime_error("cannot create a unique tuning directory"); }
        }
        else { directory.path = keep; fs::create_directories(directory.path); }
        const fs::path& work = directory.path;
        ODIA::tune::TuneParams tp;
        const auto setting = [this](const char* name) { return trainingSetting_(name); };
        tp.report = ids;
        tp.q_value = setting("filter:q_value");
        tp.min_charge = static_cast<int>(setting("filter:min_charge"));
        tp.allow_z1 = setting("filter:allow_z1").toString() == "true";
        tp.rt_spread_max = setting("filter:rt_spread_max");
        tp.rt_max_minutes = setting("filter:rt_max_minutes");
        tp.train_size = static_cast<std::size_t>(std::max(0, static_cast<int>(setting("cohort:train_size"))));
        tp.train_frac = setting("cohort:train_frac");
        tp.inner_val = setting("cohort:no_inner_val").toString() != "true";
        tp.full_fit = setting("cohort:full_fit").toString() == "true";
        tp.epochs = static_cast<int>(setting("train:epochs"));
        tp.warmup = static_cast<int>(setting("train:warmup"));
        tp.lr = setting("train:lr");
        tp.batch_size = static_cast<int>(setting("train:batch_size"));
        tp.eval_every = static_cast<int>(setting("stop:eval_every"));
        tp.min_epochs = static_cast<int>(setting("stop:min_epochs"));
        tp.patience = static_cast<int>(setting("stop:patience"));
        tp.rel_tol = setting("stop:rel_tol");
        tp.abs_tol = setting("stop:abs_tol");
        tp.max_seconds = setting("stop:max_seconds");
        tp.select = setting("stop:select").toString() == "rmse" ? ODIA::tune::Select::Rmse : ODIA::tune::Select::CalibratedSd;
        tp.device = setting("machine:device").toString();
        // 0: the platform default. Resolved here so that the trainer, and the
        // provenance, see a concrete count (ODIA::tune::defaultTrainingThreads).
        tp.threads = static_cast<int>(setting("machine:train_threads"));
        if (tp.threads <= 0) { tp.threads = ODIA::tune::defaultTrainingThreads(); }
        tp.cudnn = setting("machine:no_cudnn").toString() != "true";
        tp.tf32 = setting("machine:tf32").toString() == "true";
        tp.seed = static_cast<std::uint32_t>(static_cast<int>(setting("machine:seed")));
        if (tp.epochs < 1 || tp.warmup < 0 || tp.warmup > tp.epochs)
        { throw std::runtime_error("train:warmup must be in [0, train:epochs]"); }
        if (tp.train_size && tp.train_frac > 0)
        { throw std::runtime_error("give cohort:train_size or cohort:train_frac, not both"); }

        // The heads. auto: RT, and CCS when the report has the ion mobility
        // for it -- the same data selection and cohort check as training, so
        // that a report with too few usable 1/K0 values drops CCS instead of
        // failing the run.
        bool want_rt = heads != "ccs", want_ccs = heads == "ccs" || heads == "both";
        std::string heads_note;
        if (heads == "auto")
        {
          const fs::path ccs_model = fs::path(models) / "peptdeep_ccs_dynamic.onnx";
          if (!fs::exists(ccs_model)) { heads_note = "no CCS model in " + models; }
          else if (ODIA::tune::reportIonMobilityRows(ids) == 0) { heads_note = "the report has no ion mobility"; }
          else
          {
            ODIA::tune::TuneParams probe = tp;
            probe.head = ODIA::tune::HeadKind::CCS;
            probe.model_in = ccs_model.string();
            std::ostringstream quiet;
            try { ODIA::tune::validateTrainingReport(probe, quiet); want_ccs = true; }
            catch (const std::exception& e)
            {
              heads_note = std::string("too little usable ion mobility: ") + e.what();
              warn_("-tune_heads auto: CCS not tuned (" + heads_note + "); -tune_heads both makes this an error");
            }
          }
          writeLogInfo_(std::string("-tune_heads auto: ") + (want_ccs ? "RT and CCS (the report has ion mobility)"
                                                                       : "RT only (" + heads_note + ")"));
        }
        tune_prov["heads"] = {{"requested", heads}, {"tuned", want_rt && want_ccs ? "both" : want_ccs ? "ccs" : "rt"},
                              {"note", heads_note}};

        for (const char* head : {"rt", "ccs"})
        {
          if ((std::string(head) == "rt" && !want_rt) || (std::string(head) == "ccs" && !want_ccs)) { continue; }
          for (const char* suffix : {".onnx", ".onnx.tune.json", ".onnx.trajectory.tsv"})
          {
            const fs::path file = work / (std::string("peptdeep_") + head + "_dynamic" + suffix);
            if (fs::exists(file) || fs::is_symlink(file)) { throw std::runtime_error("refusing to overwrite tuned model artifact: " + file.string()); }
          }
        }

        // Validate every requested head before the first expensive training run.
        for (const bool ccs : {false, true})
        {
          if (ccs ? !want_ccs : !want_rt) { continue; }
          const char* file = ccs ? "peptdeep_ccs_dynamic.onnx" : "peptdeep_rt_dynamic.onnx";
          tp.head = ccs ? ODIA::tune::HeadKind::CCS : ODIA::tune::HeadKind::RT;
          tp.model_in = (fs::path(models) / file).string();
          tp.model_out = (work / file).string();
          if (!fs::exists(tp.model_in)) { throw std::runtime_error("no " + std::string(file) + " in " + models); }
          ODIA::tune::validateTrainingReport(tp, std::cout);
        }

        auto run_head = [&](ODIA::tune::HeadKind head, const char* file)
        {
          tp.head = head;
          tp.model_in = (fs::path(models) / file).string();
          tp.model_out = (work / file).string();
          // finetune throws when nothing beat the stock model, the same way it
          // throws on a corrupt report -- TuneResult carries no "exported" flag
          // to tell them apart, so a head that does not improve aborts the run
          // rather than silently refining with stock predictions.
          return ODIA::tune::finetune(tp, std::cout);
        };
        auto training_provenance = [](const ODIA::tune::TuneResult& result) {
          const std::string path = result.model_out + ".tune.json";
          std::ifstream input(path);
          if (!input) { throw std::runtime_error("cannot read training provenance: " + path); }
          return json::parse(input);
        };

        const bool free_cys = setting("tune_keep_free_cysteine_offset").toString() == "true";
        // Re-prediction follows -machine:device as generate does: CUDA when
        // auto finds a usable GPU (RT and CCS need cuDNN), the CPU otherwise
        // (CoreML has no LSTM; mps is the training device). The deprecated
        // -tune_predict_gpu still forces the GPU.
        std::vector<std::string> device_warnings;
        const ODIA::GpuOptions gpu_options = gpuOptions_(&device_warnings);
        bool gpu = getFlag_("tune_predict_gpu");
        {
          const std::string device = setting("machine:device").toString();
          if (gpu && device == "cpu")
          { throw std::runtime_error("-tune_predict_gpu (deprecated) contradicts -machine:device cpu; drop -tune_predict_gpu"); }
          // -machine:require_gpu asks the re-prediction for the GPU where there
          // is a CUDA route at all (a Mac has none for RT and CCS: there it is
          // the training device's matter, MPS).
          const ODIA::DeviceCaps caps = ODIA::probeDeviceCaps();
          const bool require = gpu_options.require && caps.cuda_provider;
          const ODIA::DeviceRoute route = ODIA::chooseDeviceRoute(device == "mps" || device == "coreml" ? "auto" : device,
                                                                  require, caps, library.precursorCount(),
                                                                  static_cast<std::size_t>(-1), want_ccs);
          if (!route.error.empty()) { throw std::runtime_error(route.error); }
          gpu = gpu || require || (route.lstm_gpu && !route.lstm_coreml);
          writeLogInfo_(std::string("re-prediction on the ") + (gpu ? "GPU (CUDA)" : "CPU"));
        }
        // 0 follows -threads, the run's CPU budget; 0.11.2 lifted the cap of
        // 16 sessions, so 0 meaning "all cores" would take a 2 x 96-core node
        // whatever -threads said. The predictor never exceeds the physical
        // cores (Hardware.h).
        const int requested_sessions = getIntOption_("tune_predict_sessions");   // a resource, not a recipe
        const unsigned sessions = static_cast<unsigned>(
          std::max(0, requested_sessions > 0 ? requested_sessions : getIntOption_("threads")));

        if (want_rt)
        {
          ODIA::trace::Scope train_phase("tune_rt", 0);
          const ODIA::tune::TuneResult r = run_head(ODIA::tune::HeadKind::RT, "peptdeep_rt_dynamic.onnx");
          train_phase.end();
          writeLogInfo_("tuned RT: " + r.stop_reason + ", best epoch " + std::to_string(r.best_epoch) +
                        " of " + std::to_string(r.epochs_run));
          // predictRetentionTimes returns what it could NOT predict, not what it
          // did -- undocumented, and it reads exactly the other way round.
          const std::size_t unpredicted =
            ODIA::LibraryGenerator::predictRetentionTimes(library, r.model_out, gpu, sessions, free_cys, gpu_options);
          if (unpredicted) { throw std::runtime_error("tuned RT could not encode " + std::to_string(unpredicted) + " precursors; refusing mixed RT units"); }
          const std::size_t n = library.precursorCount() - unpredicted;
          // The model emits rt_norm = RT / rt_max_minutes, and that denominator
          // lives nowhere else. Multiplying it back puts the WHOLE library in the
          // reference run's minutes -- the same unit refine() writes for the
          // matched precursors, so the column stays one coherent object.
          for (auto& v : library.precursors().irt) { v *= r.rt_max_minutes; }
          p.library_rt_in_minutes = true;
          tune_prov["rt"] = {{"stop_reason", r.stop_reason}, {"best_epoch", r.best_epoch},
                             {"epochs_run", r.epochs_run}, {"rt_max_minutes", r.rt_max_minutes},
                             {"model_sha256", r.model_out_sha256}, {"stock_sha256", r.model_in_sha256},
                             {"training", training_provenance(r)},
                             {"repredicted", n}, {"unpredicted", unpredicted}};
          writeLogInfo_("re-predicted RT for " + std::to_string(n) + " of " + std::to_string(library.precursorCount()) +
                        " precursors, in the run's minutes" +
                        (unpredicted ? " (" + std::to_string(unpredicted) + " could not be encoded)" : ""));
        }
        if (want_ccs)
        {
          ODIA::trace::Scope train_phase("tune_ccs", 0);
          const ODIA::tune::TuneResult r = run_head(ODIA::tune::HeadKind::CCS, "peptdeep_ccs_dynamic.onnx");
          train_phase.end();
          writeLogInfo_("tuned CCS: " + r.stop_reason + ", best epoch " + std::to_string(r.best_epoch) +
                        " of " + std::to_string(r.epochs_run));
          // derive_mobility rewrites the WHOLE 1/K0 column, so it runs only when
          // the CCS head actually tuned -- otherwise a library that arrived with
          // measured mobilities would lose them to stock predictions.
          const std::size_t unpredicted =
            ODIA::LibraryGenerator::predictCollisionCrossSections(library, r.model_out, gpu, sessions, true, gpu_options);
          if (unpredicted) { throw std::runtime_error("tuned CCS could not encode " + std::to_string(unpredicted) + " precursors; refusing a partially re-predicted library"); }
          const std::size_t n = library.precursorCount() - unpredicted;
          tune_prov["ccs"] = {{"stop_reason", r.stop_reason}, {"best_epoch", r.best_epoch},
                              {"epochs_run", r.epochs_run},
                              {"model_sha256", r.model_out_sha256}, {"stock_sha256", r.model_in_sha256},
                              {"training", training_provenance(r)},
                              {"repredicted", n}, {"unpredicted", unpredicted}};
          writeLogInfo_("re-predicted CCS and 1/K0 for " + std::to_string(n) + " of " +
                        std::to_string(library.precursorCount()) + " precursors" +
                        (unpredicted ? " (" + std::to_string(unpredicted) + " could not be encoded)" : ""));
        }

        if (!device_warnings.empty()) { tune_prov["device_warnings"] = device_warnings; }
        tune_prov["full_fit"] = tp.full_fit;   // recorded for either head, not only when RT was tuned
        tune_prov["models_kept"] = keep.empty() ? json(nullptr) : json(fs::absolute(work).string());
      }
#endif
      if (tune_only)
      {
        st.library_before = st.library_after = library.precursorCount();
        writeLogInfo_("tune mode preserved all " + std::to_string(library.precursorCount()) + " precursors; no observed values written");
      }
      else
      {
      writeLogInfo_("reference" + (st.run.empty() ? std::string() : " (run " + st.run + ")") + ": " +
                    std::to_string(st.ids_rows) + " rows, " + std::to_string(st.ids_precursors) + " precursors, " +
                    std::to_string(st.ids_passing) + " passing the gates; rejected: " +
                    std::to_string(st.ids_decoy) + " decoy, " + std::to_string(st.ids_q_invalid) + " invalid q, " +
                    std::to_string(st.ids_q_above) + " above threshold, " + std::to_string(st.ids_charge_invalid) + " bad charge, " +
                    std::to_string(st.ids_sequence_invalid) + " empty sequence; " + std::to_string(st.ids_fragment_invalid) +
                    " invalid fragment identities ignored; " + std::to_string(st.ids_too_few_fragments) +
                    " precursors below min_fragments");
      for (const auto& g : st.gates_bypassed)
      { writeLogWarn_("gate " + g + " BYPASSED: the reference has no such column (-empirical_library); recorded in provenance"); }
      if (st.ids_unknown_mod_tokens)
      { writeLogWarn_(std::to_string(st.ids_unknown_mod_tokens) + " modification tokens were neither a known name nor a UniMod accession; passed through verbatim"); }
      if (st.ids_ramp_censored)
      { writeLogInfo_(std::to_string(st.ids_ramp_censored) + " observed 1/K0 values within " + std::to_string(p.im_ramp_margin) +
                      " of the declared ramp top " + std::to_string(p.im_ramp_top) + " treated as censored"); }

      {
        ODIA::trace::Scope apply_phase("apply", 0);
        ODIA::LibraryRefiner::refine(library, obs, p, st);
      }
      if (st.lib_unknown_mod_tokens)
      { writeLogWarn_(std::to_string(st.lib_unknown_mod_tokens) + " library modification tokens could not be resolved to UniMod; passed through verbatim and may fail the reference join"); }
      }
    }
    catch (const std::exception& e)
    {
      // A search failure names its stage already ("search: ...").
      const std::string what = e.what();
      writeLogError_(what.rfind("search", 0) == 0 ? what : "refine: " + what);
      keptReport();
      return UNEXPECTED_RESULT;
    }

    if (!tune_only)
    {
    std::ostringstream m;
    m.setf(std::ios::fixed); m.precision(1);
    m << "matched " << st.ids_passing - st.ids_unmatched << " of " << st.ids_passing << " reference precursors ("
      << 100.0 * st.match_fraction << "%); " << st.ids_unmatched << " observed but absent from the library";
    writeLogInfo_(m.str());
    auto resid = [&](const char* axis, std::size_t n, double mean, double sd, double p95, int prec)
    {
      if (!n) { return; }
      std::ostringstream r; r.setf(std::ios::fixed); r.precision(prec);
      r << axis << " residual BEFORE refinement (library prediction - observed), n=" << n
        << ": mean " << mean << ", sd " << sd << ", p95 |resid| " << p95;
      writeLogInfo_(r.str());
    };
    resid("RT", st.rt_resid_n, st.rt_resid_mean, st.rt_resid_sd, st.rt_resid_p95, 4);
    resid("1/K0", st.im_resid_n, st.im_resid_mean, st.im_resid_sd, st.im_resid_p95, 5);
    writeLogInfo_("wrote " + std::to_string(st.rt_written) + " RT values" +
                  (st.rt_missing ? " (" + std::to_string(st.rt_missing) + " matched targets had no observed RT and keep their prediction)" : "") +
                  " and " + std::to_string(st.im_written) + " 1/K0 values" +
                  (st.im_charge_excluded ? " (" + std::to_string(st.im_charge_excluded) + " below im_min_charge kept their prediction)" : "") +
                  (st.im_missing ? " (" + std::to_string(st.im_missing) + " had no usable observed 1/K0)" : ""));
    writeLogInfo_("library " + std::to_string(st.library_before) + " -> " + std::to_string(st.library_after) +
                  " precursors (" + std::to_string(st.matched) + " targets + " + std::to_string(st.decoys_kept) + " decoys)" +
                  (p.filter ? "" : " (filter off)"));
    if (p.write_intensity)
    {
      std::ostringstream r; r.setf(std::ios::fixed); r.precision(3);
      r << "replaced fragment intensities on " << st.intensity_replaced_precursors << " of " << st.intensity_candidate_precursors
        << " candidate precursors (" << st.intensity_kept_predicted << " kept their predictions); transitions per replaced precursor "
        << st.intensity_transitions_before << " -> " << st.intensity_transitions_after
        << "; observed base peak was already the library's top no-loss transition in " << 100.0 * st.intensity_rank_agreement << "%";
      writeLogInfo_(r.str());
      writeLogInfo_("fragment fates over " + std::to_string(st.intensity_candidate_transitions) + " candidate transitions: " +
                    std::to_string(st.intensity_matched_transitions) + " matched (" + std::to_string(st.intensity_gated_zero_quant) +
                    " zero area, " + std::to_string(st.intensity_gated_correlation) + " below the correlation gate, " +
                    std::to_string(st.intensity_gated_floor) + " below the floor), " + std::to_string(st.intensity_unmatched_in_library) +
                    " not in the report, " + std::to_string(st.intensity_mz_mismatch) + " m/z disagreements, " +
                    std::to_string(st.intensity_loss_bearing) + " loss-bearing; " + std::to_string(st.intensity_observed_not_in_library) +
                    " report fragments the library never carried were NOT added");
      if (st.intensity_loss_bearing)
      { warn_(std::to_string(st.intensity_loss_bearing) + " neutral-loss transitions have no observed intensity (the "
              "report has no loss field) and were dropped from the replaced precursors"); }
      if (!p.filter)
      { warn_("mixed intensities: matched precursors carry observed, the others predicted intensities "
              "(-allow_mixed_intensity)"); }
    }
    if (p.write_rt && st.rt_written)
    {
      writeLogInfo_("RT now holds the reference run's observed retention times (" + std::string(rtUnitName(p.rt_unit)) +
                    "), not iRT: use this library for that run's gradient");
    }
    }

    // Provenance: the recipe, the inputs by content hash, the run, the units and every
    // count above -- embedded in the Parquet and always written as a sidecar (M13).
    ODIA::trace::Scope provenance_phase("provenance", 0);
    // Everything that ran beside the main thread is done before the provenance.
    if (cache_removal)
    {
      const std::string message = cache_removal->wait();
      if (!message.empty()) { writeLogWarn_(message); }
      cache_removal.reset();
    }
    const std::string library_hash = library_hash_ahead.valid() ? library_hash_ahead.get() : ODIA::DIANNLibraryFile::hashFile(in);
    json prov = {
      {"tool", "DIALibGen"}, {"tool_version", DIALIBGEN_VERSION}, {"mode", tune_only ? "tune" : "refine"},
      {"config", eff},
      {"inputs", {{"library", std::filesystem::absolute(in).string()}, {"library_fnv1a64", library_hash},
                  {"reference", std::filesystem::absolute(ids).string()}, {"reference_fnv1a64", ODIA::DIANNLibraryFile::hashFile(ids)},
                  {"reference_run", st.run}}},
      {"gates_bypassed", st.gates_bypassed},
      {"tune", tune_prov},
      {"reference", {{"rows", st.ids_rows}, {"precursors", st.ids_precursors}, {"passing", st.ids_passing},
                     {"decoy", st.ids_decoy}, {"q_invalid", st.ids_q_invalid}, {"q_above", st.ids_q_above},
                     {"charge_invalid", st.ids_charge_invalid}, {"unknown_mod_tokens", st.ids_unknown_mod_tokens},
                     {"sequence_invalid", st.ids_sequence_invalid}, {"fragment_invalid", st.ids_fragment_invalid},
                     {"too_few_fragments", st.ids_too_few_fragments},
                     {"ramp_censored", st.ids_ramp_censored}, {"unmatched", st.ids_unmatched}}},
      {"library", {{"before", st.library_before}, {"after", st.library_after}, {"matched_targets", st.matched},
                   {"unknown_mod_tokens", st.lib_unknown_mod_tokens},
                   {"rt_repredicted_in_reference_minutes", p.library_rt_in_minutes},
                   {"decoys_kept", st.decoys_kept}, {"match_fraction", st.match_fraction},
                   {"rt_written", st.rt_written}, {"rt_missing", st.rt_missing},
                   {"im_written", st.im_written}, {"im_missing", st.im_missing}, {"im_charge_excluded", st.im_charge_excluded}}},
      {"intensity", p.write_intensity ? json{
                   {"candidate_precursors", st.intensity_candidate_precursors}, {"replaced_precursors", st.intensity_replaced_precursors},
                   {"replaced_decoys", st.intensity_replaced_decoys}, {"kept_predicted", st.intensity_kept_predicted},
                   {"decoy_asymmetry", st.intensity_decoy_asymmetry}, {"duplicate_key", st.intensity_duplicate_key},
                   {"candidate_transitions", st.intensity_candidate_transitions}, {"matched_transitions", st.intensity_matched_transitions},
                   {"mz_mismatch", st.intensity_mz_mismatch}, {"mz_mismatch_fraction", st.intensity_mz_mismatch_fraction},
                   {"unmatched_in_library", st.intensity_unmatched_in_library}, {"loss_bearing", st.intensity_loss_bearing},
                   {"observed_not_in_library", st.intensity_observed_not_in_library},
                   {"gated_zero_quant", st.intensity_gated_zero_quant}, {"gated_correlation", st.intensity_gated_correlation},
                   {"gated_floor", st.intensity_gated_floor}, {"bad_tokens", st.intensity_bad_tokens},
                   {"row_length_mismatch", st.intensity_row_length_mismatch},
                   {"rank_agreement", num(st.intensity_rank_agreement)},
                   {"transitions_before", num(st.intensity_transitions_before)}, {"transitions_after", num(st.intensity_transitions_after)},
                   {"mixed_provenance", !p.filter}} : json(nullptr)},
      {"residual_before", {{"rt", {{"n", st.rt_resid_n}, {"mean", num(st.rt_resid_mean)}, {"sd", num(st.rt_resid_sd)}, {"p95_abs", num(st.rt_resid_p95)}}},
                           {"im", {{"n", st.im_resid_n}, {"mean", num(st.im_resid_mean)}, {"sd", num(st.im_resid_sd)}, {"p95_abs", num(st.im_resid_p95)}}},
                           {"sd_convention", "ddof=1; p95 = lower nearest-rank quantile of |residual|; NaN -> null when n<2"}}},
      {"units", {{"rt", st.rt_written ? (p.rt_unit == ODIA::RefineParams::RtUnit::MinMax ? "0..100 minmax over the matched set" :
                         (!p.filter ? "reference-run minutes: observed where available, tuned predictions otherwise" : "the reference run's own RT units")) :
                         (tune_prov.contains("rt") ? "tuned model predictions in reference-run minutes" : "unchanged (library prediction)")},
                 {"im", p.write_im ? "observed 1/K0 for charges >= im_min_charge; CCS cleared where written" :
                         (tune_prov.contains("ccs") ? "1/K0 derived from tuned CCS predictions" : "unchanged (library prediction)")},
                 {"intensity", p.write_intensity ? std::string("observed fragment areas from the reference run, normalised by ") + intensityNormName(p.intensity_norm) +
                                                   (p.filter ? "" : "; UNMATCHED precursors keep MS2-model predictions") : std::string("unchanged (library prediction)")}}},
      {"warning", tune_only ? "Predictions adapted to one reference run; evaluate transfer on a different run." :
                             "Observed values describe the reference run and its gradient; evaluate transfer on a different run."}};
    if (search)
    {
      prov["inputs"]["run"] = std::filesystem::absolute(run).string();
      prov["inputs"]["run_fnv1a64"] = run_hash;
      prov["search"] = search_prov;
    }
    if (tune_only)
    {
      prov["reference"] = nullptr;
      prov["residual_before"] = nullptr;
      prov["inputs"]["reference_run"] = nullptr;
      for (const char* key : {"matched_targets", "decoys_kept", "match_fraction"}) { prov["library"].erase(key); }
    }

    provenance_phase.end();
    try
    {
      ODIA::trace::Scope write_phase("write", 0);
      // The writers' output does not depend on the thread count.
      const unsigned write_threads = static_cast<unsigned>(std::max(0, getIntOption_("threads")));
      ODIA::AtomicFile output(out), sidecar(out + ".refine.json");
      // The library embeds what is a function of the inputs and settings;
      // the search's timings and thread count stay in the sidecar, so the
      // same -run command writes the same library bytes.
      json embedded = prov;
      if (embedded.contains("search") && embedded["search"].is_object()) { embedded["search"].erase("resources"); }
      if (format == ODIA::LibraryFormat::DIALibGenParquet)
      {
        ODIA::DIANNLibraryFile::storeParquetCompact(output.temporaryPath().string(), library, ODIA::DIANNLibraryFile::Fingerprint{},
                                                    embedded.dump(), {}, write_threads);
      }
      else if (format == ODIA::LibraryFormat::DIANNTSV)
      { ODIA::DIANNLibraryFile::storeTSV(output.temporaryPath().string(), library, write_threads); }
      else
      {
        std::vector<std::pair<std::string, std::string>> metadata;
        if (ODIA::LibraryExport::carriesRecipe(format)) { metadata.emplace_back("odia.config_json", embedded.dump()); }
        ODIA::LibraryExport::store(output.temporaryPath().string(), library, format, metadata, write_threads);
      }
      std::ofstream side(sidecar.temporaryPath());
      side << prov.dump(2) << '\n'; side.close();
      if (!side) { throw std::runtime_error("cannot write " + out + ".refine.json"); }

      std::unique_ptr<ODIA::AtomicFile> staged_report;
      if (!report.empty())
      {
        staged_report = std::make_unique<ODIA::AtomicFile>(report);
        std::ofstream o(staged_report->temporaryPath());
        o << "metric\tvalue\n";
        for (const auto& [k, v] : prov["reference"].items()) { o << "ids_" << k << '\t' << v << '\n'; }
        for (const auto& [k, v] : prov["library"].items()) { o << "library_" << k << '\t' << v << '\n'; }
        o.setf(std::ios::fixed); o.precision(6);
        o << "rt_resid_n\t" << st.rt_resid_n << "\nrt_resid_mean\t" << st.rt_resid_mean << "\nrt_resid_sd\t" << st.rt_resid_sd
          << "\nrt_resid_p95\t" << st.rt_resid_p95 << "\nim_resid_n\t" << st.im_resid_n << "\nim_resid_mean\t" << st.im_resid_mean
          << "\nim_resid_sd\t" << st.im_resid_sd << "\nim_resid_p95\t" << st.im_resid_p95 << '\n';
        o.close();
        if (!o) { throw std::runtime_error("cannot write " + report); }
      }
      // Prepare every file before publishing any. Commit the library last so
      // its appearance means the corresponding provenance is already present.
      sidecar.commit();
      if (staged_report) { staged_report->commit(); }
      output.commit();
      writeLogInfo_("wrote " + out + " and " + out + ".refine.json");
      if (staged_report) { writeLogInfo_("wrote report to " + report); }
    }
    catch (const std::exception& e) { writeLogError_(std::string("write: ") + e.what()); keptReport(); return CANNOT_WRITE_OUTPUT_FILE; }

    return EXECUTION_OK;
  }
