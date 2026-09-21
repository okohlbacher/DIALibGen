// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <odia/DIANNLibraryFile.h>
#include <odia/AtomicFile.h>
#include <memory>
#include <odia/Library.h>
#include <odia/LibraryGenerator.h>
#include <odia/LibraryRefiner.h>
#include <odia/search/Identifier.h>
#include <odia/search/SearchParams.h>

#ifdef DIALIBGEN_WITH_FINETUNE
#include <odia/tune/Trainer.h>
#include <odia/PeptDeepEncoder.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <unordered_set>
#endif

#include "DIALibGen.h"

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <nlohmann/json.hpp>

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <random>
#include <sstream>

using json = nlohmann::json;

namespace
{
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
    { if (!ref.contains(k)) { throw std::runtime_error("unknown config key: " + k); } }
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
    registerInputFile_("ids", "<file>", "",
                       "DIA-NN report.parquet, or a pre-filtered library with -empirical_library. Modification names are canonicalised.", false);
    setValidFormats_("ids", {"parquet"}, false);
    registerInputFile_("run", "<file>", "",
                       "EXPERIMENTAL, in development: a centroided DIA run (mzML) that the built-in identification step "
                       "searches instead of reading -ids. Give exactly one of -ids and -run. Settings: search:.", false);
    setValidFormats_("run", {"mzML"}, false);
    registerOutputFile_("out_ids", "<file>", "",
                        "With -run: the identification report (Parquet, DIA-NN column names), which then serves as -ids. "
                        "Default: <out>.ids.parquet. Never overwritten.", false);
    setValidFormats_("out_ids", {"parquet"}, false);
    registerOutputFile_("out_report", "<file>", "",
                        "Per-axis residual report (TSV), measured BEFORE the overwrite.", false);
    setValidFormats_("out_report", {"tsv"}, false);
    registerFlag_("no_filter", "Keep precursors the reference did not identify; requires -no_write_rt or successful RT tuning so RT units stay consistent.");
    registerFlag_("empirical_library", "Declare -ids a pre-filtered empirical library rather than a report: gates "
                                       "whose columns are absent are BYPASSED and each bypass is recorded. Without "
                                       "this, a missing gate column is an error.");
    registerFlag_("no_write_rt", "Keep predicted RT rather than replacing it with observed RT.");
    registerFlag_("write_im", "Also overwrite 1/K0 with the observed value for charges >= -im_min_charge; off by default.");
    registerFlag_("write_intensity", "Replace predicted fragment intensities with the reference run's observed ones. "
                                     "Needs Fragment.Info/Fragment.Quant.Raw (1.9) or Fr.N.Id/Quantity (2.x), exported with "
                                     "--report-lib-info or --export-quant respectively. Every match is cross-checked on fragment m/z, and a run that "
                                     "replaces nothing is an error.");
    registerDoubleOption_("intensity_min_correlation", "<r>", 0.0, "Require fragment quality greater than this: "
                          "correlation in DIA-NN 1.9, Score in DIA-NN 2.x. -1 disables the quality gate.", false);
    registerFlag_("intensity_no_restrict", "Replace a precursor only when EVERY one of its transitions is trusted, else keep "
                                           "its predictions whole. Neutral-loss transitions cannot match report fragments, so "
                                           "precursors carrying them keep their predictions. Transition counts cannot change.");
    registerFlag_("intensity_no_rerank", "Keep a replaced precursor's transitions in their original order.");
    registerIntOption_("intensity_min_fragments", "<n>", 3, "A replaced precursor keeps at least this many transitions or "
                       "keeps its predictions whole.", false);
    registerStringOption_("intensity_norm", "<rule>", "library_max", "How an observed area becomes a library intensity. "
                          "library_max scales to the maximum that precursor already held, which preserves the file's own "
                          "convention -- base peak = 1 is not an invariant of these libraries.", false);
    setValidStrings_("intensity_norm", {"library_max", "base_peak", "sum", "raw"});
    registerDoubleOption_("intensity_min_relative", "<fraction>", 1e-4,
                          "Drop observed intensities below this fraction of the kept maximum.", false);
    registerDoubleOption_("intensity_mz_tol_ppm", "<ppm>", 20.0, "A fragment matches only when the report's m/z agrees with "
                          "the library's to within this.", false);
    registerDoubleOption_("intensity_max_mz_mismatch", "<frac>", 0.01, "Refuse when more than this fraction of identity "
                          "matches fail the m/z check. 1 surveys a suspect pairing instead of failing.", false);
    registerFlag_("allow_mixed_intensity", "Permit -write_intensity with -no_filter, which leaves observed and predicted "
                                           "intensities in one library. Recorded in the provenance.");
    registerDoubleOption_("q_precursor", "<q>", 0.01, "Precursor q-value gate (>= 1 disables).", false);
    registerDoubleOption_("q_global", "<q>", 0.01, "Global/peptide q-value gate (>= 1 disables).", false);
    registerDoubleOption_("q_protein", "<q>", 0.01, "Protein q-value gate (>= 1 disables).", false);
    registerIntOption_("min_fragments", "<n>", 0,
                       "Minimum DISTINCT reference fragments per precursor; needs fragment identities in -ids. 0 = off.", false);
    registerStringOption_("rt_unit", "<unit>", "observed",
                          "observed = the reference's own units; minmax = rescaled to 0..100 over the matched set "
                          "(refused with -no_filter).", false);
    setValidStrings_("rt_unit", {"observed", "minmax"});
    registerStringOption_("dedup", "<rule>", "lowest_q",
                          "Which observation wins when a precursor was seen more than once: lowest_q (ties: lower "
                          "PEP, then first seen) or highest_evidence (needs an Evidence column).", false);
    setValidStrings_("dedup", {"lowest_q", "highest_evidence"});
    registerIntOption_("im_min_charge", "<z>", 2, "Charges below this never receive an observed 1/K0 (z1 is censored "
                                                  "at the ramp top on timsTOF diaPASEF).", false);
    registerDoubleOption_("im_ramp_top", "<1/K0>", 0.0, "The instrument's mobility ramp top, if known; observations "
                                                        "within -im_ramp_margin of it are treated as censored. 0 = unknown.", false);
    registerDoubleOption_("im_ramp_margin", "<1/K0>", 0.02, "See -im_ramp_top.", false);
    registerDoubleOption_("min_match_fraction", "<f>", 0.0,
                          "Refuse when fewer than this fraction of passing reference precursors match the library. "
                          "0 = refuse only when nothing matches.", false);

    // Registered in EVERY build, including one without libtorch, so --help,
    // -write_ini and -write_ctd describe the same tool on every platform. -tune
    // itself is what fails when the stage was not compiled in.
    registerFlag_("tune", "Fine-tune the RT and CCS models on -ids and re-predict the whole library through them "
                          "BEFORE refining, so precursors the reference never identified are corrected too. "
                          "Needs a build with the fine-tuning stage. Turns a seconds-long refinement into a "
                          "training run plus whole-library inference.");
    registerStringOption_("tune_models", "<dir>", "", "Directory holding the stock peptdeep_{rt,ccs}_dynamic.onnx. "
                                                      "Default: $DIALIBGEN_MODEL_DIR or the bundled models.", false);
    registerStringOption_("tune_heads", "<which>", "both", "Which models to tune.", false);
    setValidStrings_("tune_heads", {"rt", "ccs", "both"});
    registerStringOption_("tune_out_models", "<dir>", "", "Keep the tuned ONNX files and their .tune.json and "
                                                          ".trajectory.tsv sidecars here. Default: a scratch "
                                                          "directory, removed on exit -- the deliverable is the library.", false);
    registerFlag_("tune_predict_gpu", "Use the GPU for the re-prediction pass (the ONNX one, not training).");
    registerIntOption_("tune_predict_sessions", "<n>", 0, "ONNX Runtime sessions for the re-prediction pass; 0 = default.", false);
    registerFlag_("tune_keep_free_cysteine_offset",
                  "Keep DIALibGen's free-cysteine RT offset when re-predicting with a TUNED RT model. Off by "
                  "default: that offset was fitted against the STOCK model, and a model tuned on this run's own "
                  "identifications has had the chance to learn the effect itself -- applying both counts it twice.");

    registerTOPPSubsection_("filter", "Fine-tuning: which observations train the model. NOTE: this is the TRAINING-set filter, not the library filter -- that is -q_precursor and friends.");
    registerDoubleOption_("filter:q_value", "<q>", 0.01, "Precursor Q.Value threshold", false);
    registerIntOption_("filter:min_charge", "<z>", 2, "CCS only: lowest precursor charge used (RT collapses all charges). z1 is censored at the mobility ramp top on timsTOF", false);
    setMinInt_("filter:min_charge", 1); setMaxInt_("filter:min_charge", 8);
    registerFlag_("filter:allow_z1", "CCS: permit filter:min_charge 1 (censored observations enter training)");
    registerDoubleOption_("filter:rt_spread_max", "<min>", 0.2, "Drop a sequence whose charge states' RTs span more than this (minutes)", false);
    registerDoubleOption_("filter:rt_max_minutes", "<min>", 0.0, "rt_norm denominator; 0 = the run's maximum observed RT", false);

    registerTOPPSubsection_("cohort", "Fine-tuning: protein-level cohorts (frozen before subsampling)");
    registerIntOption_("cohort:train_size", "<n>", 0, "Training units (sequences for rt, sequence x charge for ccs); 0 = full pool", false);
    registerDoubleOption_("cohort:train_frac", "<f>", 0.0, "Alternative to train_size: fraction of the pool", false);
    registerFlag_("cohort:full_fit", "Train on EVERY unit of the run, the test and validation cohorts included: the closest "
                                     "fit the data allows. val and TEST are then in-sample; judge this mode by searching a "
                                     "DIFFERENT run with the result, never by its own numbers");
    registerFlag_("cohort:no_inner_val", "No inner validation cohort: its units rejoin the pool and checkpoints are selected on TEST, which is then no longer a held-out number");

    registerTOPPSubsection_("train", "Fine-tuning: recipe");
    registerIntOption_("train:epochs", "<n>", 100, "Horizon of the cosine schedule (and the maximum epochs)", false);
    registerIntOption_("train:warmup", "<n>", 10, "Linear warmup epochs", false);
    registerDoubleOption_("train:lr", "<lr>", 1e-4, "Peak learning rate (Adam)", false);
    registerIntOption_("train:batch_size", "<n>", 1024, "Batch size within a length group", false);

    registerTOPPSubsection_("stop", "Fine-tuning: convergence");
    registerIntOption_("stop:eval_every", "<n>", 1, "Validate every n epochs", false);
    registerIntOption_("stop:min_epochs", "<n>", 20, "Never stop before this epoch", false);
    registerIntOption_("stop:patience", "<n>", 10, "Stop after this many epochs without progress (never before max(min_epochs, warmup))", false);
    registerDoubleOption_("stop:rel_tol", "<f>", 0.005, "Progress = the selection metric beats the anchor by this fraction", false);
    registerDoubleOption_("stop:abs_tol", "<f>", 0.0, "Progress = beats the anchor by this absolute amount (0 = use rel_tol)", false);
    registerDoubleOption_("stop:max_seconds", "<s>", 0.0, "Wall-clock budget per head, checked at epoch boundaries; final evaluation/export may exceed it (0 = none)", false);
    registerStringOption_("stop:select", "<metric>", "calibrated_sd", "Selection metric on the validation cohort", false);
    setValidStrings_("stop:select", {"calibrated_sd", "rmse"});

    registerTOPPSubsection_("machine", "Fine-tuning: device");
    registerStringOption_("machine:device", "<dev>", "cpu", "cpu or cuda[:N]", false);
    registerIntOption_("machine:threads", "<n>", 4, "CPU threads used for training", false);
    registerFlag_("machine:no_cudnn", "CUDA: do not use cuDNN (needed when only its loader shim is installed, as in pytorch.org's libtorch zips); slower");
    registerIntOption_("machine:seed", "<n>", 20260803, "Seed for the training subsample and batch order", false);

    registerTOPPSubsection_("search", "Built-in identification with -run (EXPERIMENTAL): candidates, decoys, extraction, "
                                      "calibration and the run-level guards. See docs/design/built-in-identification.md");
    registerStringOption_("search:candidates", "<rule>", "random", "Candidate selection: random = a deterministic, label-blind "
                          "random subset of target-decoy pairs", false);
    setValidStrings_("search:candidates", {"random"});
    registerIntOption_("search:subset", "<n>", 100000, "Targets drawn from the library (0 = every eligible target)", false);
    registerIntOption_("search:max_pairs", "<n>", 40000, "Cap on the target-decoy pairs searched; time is linear in pairs (0 = no cap)", false);
    registerStringOption_("search:decoys", "<method>", "shuffle", "How the search's in-memory decoys are built from the selected "
                          "targets; both methods keep the termini. Decoys in the library file are not searched and stay in the output", false);
    setValidStrings_("search:decoys", {"shuffle", "pseudo_reverse"});
    registerIntOption_("search:seed", "<n>", 42, "Salt of the candidate draw: changes which pairs are searched, not how", false);
    registerIntOption_("search:passes", "<n>", 1, "Extraction passes (1 in this version)", false);
    setMinInt_("search:passes", 1); setMaxInt_("search:passes", 1);
    registerDoubleOption_("search:rt_window", "<s>", 0.0, "Full RT extraction window, seconds (0 = from the calibration)", false);
    setMinFloat_("search:rt_window", 0.0);
    registerDoubleOption_("search:mz_ppm", "<ppm>", 0.0, "Full fragment m/z extraction window, ppm (0 = automatic)", false);
    setMinFloat_("search:mz_ppm", 0.0);
    registerDoubleOption_("search:im_window", "<1/K0>", 0.0, "Full 1/K0 extraction window on ion-mobility runs (0 = automatic, -1 = off)", false);
    setMinFloat_("search:im_window", -1.0);
    registerStringOption_("search:ms1", "<true/false>", "true", "Extract MS1 traces and use the MS1 sub-scores", false);
    setValidStrings_("search:ms1", {"true", "false"});
    registerStringOption_("search:rt_im_scores", "<true/false>", "true", "Let the RT and 1/K0 deviation sub-scores into the "
                          "classifier; false is an ablation for tuning, which exists to correct those deviations", false);
    setValidStrings_("search:rt_im_scores", {"true", "false"});
    registerDoubleOption_("search:calibration_min_rsq", "<r2>", 0.70, "Least r^2 of the RT calibration on the seed assays", false);
    registerDoubleOption_("search:calibration_min_coverage", "<f>", 0.30, "Least fraction of the seed assays found in the run that "
                          "the RT outlier removal must keep", false);
    for (const char* name : {"search:calibration_min_rsq", "search:calibration_min_coverage"})
    { setMinFloat_(name, 0.0); setMaxFloat_(name, 1.0); }
    registerFlag_("search:allow_bootstrap", "When the RT calibration fails, map the library RT range linearly onto the run "
                                            "instead of aborting. A test hook, recorded in the provenance");
    registerStringOption_("search:readoptions", "<mode>", "auto", "How the run is held: normal = in memory, cache = per-window "
                          "cache files, auto = cache for runs above 3 GB", false);
    setValidStrings_("search:readoptions", {"auto", "normal", "cache"});
    registerStringOption_("search:cache_dir", "<dir>", "", "Directory for cache files (default: the directory of -out_ids); "
                          "they are removed after the search", false);
    registerIntOption_("search:min_ids", "<n>", 200, "Abort when fewer target precursors pass q <= 0.01", false);
    registerDoubleOption_("search:max_target_fraction", "<f>", 0.5, "Abort when more than this fraction of the scored target "
                          "precursors passes q <= 0.01: no honest decoy set looks like that", false);
    setMinFloat_("search:max_target_fraction", 0.01); setMaxFloat_("search:max_target_fraction", 1.0);
    registerDoubleOption_("search:report_max_q", "<q>", 0.10, "Precursors, targets and decoys, up to this precursor q-value "
                          "go into -out_ids", false);
    setMinFloat_("search:report_max_q", 0.01); setMaxFloat_("search:report_max_q", 1.0);
    registerStringOption_("search:entrapment_tag", "<prefix>", "", "Protein-group prefix of entrapment proteins: log and record "
                          "the combined entrapment FDP estimate. A validation aid", false);
    registerStringOption_("search:selftest", "<true/false>", "true", "Also score with swapped and with random pair labels, "
                          "and abort unless both identify (almost) nothing. These catch a classifier that leaks labels, not "
                          "decoys that are built weaker than null targets", false);
    setValidStrings_("search:selftest", {"true", "false"});
    registerIntOption_("search:batch_size", "<n>", 500, "Advanced: precursors per extraction batch within one isolation window", false);
    registerIntOption_("search:chunk", "<n>", 20000, "Advanced: precursors per extraction call; target-decoy pairs stay together", false);
    for (const char* name : {"search:subset", "search:max_pairs", "search:seed", "search:min_ids"}) { setMinInt_(name, 0); }
    setMinInt_("search:batch_size", 1);
    setMinInt_("search:chunk", 2);

    for (const char* name : {"min_fragments", "intensity_min_fragments", "tune_predict_sessions", "cohort:train_size",
                            "train:warmup", "stop:min_epochs", "stop:patience", "machine:seed"}) { setMinInt_(name, 0); }
    for (const char* name : {"train:epochs", "train:batch_size", "stop:eval_every", "machine:threads"}) { setMinInt_(name, 1); }
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
    s.decoys = ODIA::search::parseSearchDecoyMethod(getStringOption_("search:decoys"));
    s.seed = static_cast<std::uint64_t>(count("search:seed"));
    s.passes = getIntOption_("search:passes");
    s.rt_window = getDoubleOption_("search:rt_window");
    s.mz_ppm = getDoubleOption_("search:mz_ppm");
    s.im_window = getDoubleOption_("search:im_window");
    s.ms1 = getStringOption_("search:ms1") == "true";
    s.rt_im_scores = getStringOption_("search:rt_im_scores") == "true";
    s.batch_size = count("search:batch_size");
    s.chunk = count("search:chunk");
    s.calibration_min_rsq = getDoubleOption_("search:calibration_min_rsq");
    s.calibration_min_coverage = getDoubleOption_("search:calibration_min_coverage");
    s.allow_bootstrap = getFlag_("search:allow_bootstrap");
    s.readoptions = ODIA::search::parseReadMode(getStringOption_("search:readoptions"));
    s.cache_dir = getStringOption_("search:cache_dir");
    s.min_ids = count("search:min_ids");
    s.max_target_fraction = getDoubleOption_("search:max_target_fraction");
    s.report_max_q = getDoubleOption_("search:report_max_q");
    s.entrapment_tag = getStringOption_("search:entrapment_tag");
    s.selftest = getStringOption_("search:selftest") == "true";
    s.threads = std::max(1, getIntOption_("threads"));
    s.validate();
    return s;
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
        applyJson_(json::parse(in, nullptr, true, true), p);
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
    if (p.write_rt && !p.filter && (!tune || getStringOption_("tune_heads") == "ccs"))
    { writeLogError_("observed RT with the filter off would mix reference-run minutes with library predictions; use -no_write_rt or -tune with the RT head (rt or both)"); return ILLEGAL_PARAMETERS; }
    if (tune_only)
    {
      ODIA::RefineParams defaults;
      defaults.filter = defaults.write_rt = false;
      const json expected = effectiveConfig(defaults), actual = effectiveConfig(p);
      for (const auto& [key, value] : actual.items())
      {
        if (value != expected.at(key))
        { writeLogError_(key + " has no effect in -mode tune; use the training options or select -mode refine"); return ILLEGAL_PARAMETERS; }
      }
    }

    const json eff = effectiveConfig(p);
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
      if (name == "run" && bruker)
      {
        writeLogError_("-run " + raw + " is a Bruker .d directory, which this version cannot read; convert it to a "
                       "frame-merged mzML with a per-peak 1/K0 array (mzpeak-convert --to mzml) and pass that");
      }
      else { writeLogError_("-" + name + " " + raw + " is a directory; it takes a file" + (name == "run" ? " (mzML)" : " (report.parquet)")); }
      return ILLEGAL_PARAMETERS;
    }
    const std::string in = getStringOption_("in");
    std::string ids = getStringOption_("ids");
    const std::string run = getStringOption_("run");
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
    if (!out.ends_with(".parquet") && !out.ends_with(".tsv"))
    { writeLogError_("-out must end in .parquet or .tsv"); return ILLEGAL_PARAMETERS; }

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
      if (p.write_im)
      { writeLogError_("-write_im is not available with -run yet: the built-in search does not measure 1/K0, so nothing would be written"); return ILLEGAL_PARAMETERS; }
      // Everything that would fail AFTER the search (minutes on a full run)
      // is checked before it.
      if (tune)
      {
        const std::string heads = getStringOption_("tune_heads");
        if (heads != "rt")
        {
          writeLogError_("-tune_heads " + heads + " needs observed 1/K0 values, which the built-in search (-run) does not "
                         "measure yet; use -tune_heads rt");
          return ILLEGAL_PARAMETERS;
        }
        const int epochs = getIntOption_("train:epochs"), warmup = getIntOption_("train:warmup");
        if (epochs < 1 || warmup < 0 || warmup > epochs)
        { writeLogError_("train:warmup must be in [0, train:epochs]"); return ILLEGAL_PARAMETERS; }
        if (getIntOption_("cohort:train_size") > 0 && getDoubleOption_("cohort:train_frac") > 0)
        { writeLogError_("give cohort:train_size or cohort:train_frac, not both"); return ILLEGAL_PARAMETERS; }
#ifdef DIALIBGEN_WITH_FINETUNE
        std::string models = getStringOption_("tune_models");
        if (models.empty()) { if (const char* e = std::getenv("DIALIBGEN_MODEL_DIR"); e && *e) { models = e; } }
        if (models.empty()) { models = bundledModelDir(); }
        if (models.empty())
        { writeLogError_("-tune needs -tune_models (or $DIALIBGEN_MODEL_DIR): the stock peptdeep models to start from"); return ILLEGAL_PARAMETERS; }
        if (!std::filesystem::exists(std::filesystem::path(models) / "peptdeep_rt_dynamic.onnx"))
        { writeLogError_("no peptdeep_rt_dynamic.onnx in " + models); return ILLEGAL_PARAMETERS; }
#endif
      }
      try { search = std::make_unique<ODIA::search::SearchParams>(searchParams_()); }
      catch (const std::exception& e) { writeLogError_(e.what()); return ILLEGAL_PARAMETERS; }
      out_ids = getParam_().getValue("out_ids").toString();
      if (out_ids.empty()) { out_ids = out + ".ids.parquet"; }
      if (!out_ids.ends_with(".parquet"))
      { writeLogError_("-out_ids must end in .parquet"); return ILLEGAL_PARAMETERS; }
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
      ODIA::DIANNLibraryFile::load(in, library);
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
        run_hash = ODIA::DIANNLibraryFile::hashFile(run);
        // The hand-off: the report just written IS the reference from here on,
        // read, gated and hashed by exactly the code that reads -ids.
        ids = out_ids;
        writeLogInfo_("identification report " + out_ids + " serves as -ids: " + std::to_string(found.identified) +
                      " target precursors at q <= 0.01");
      }
      ODIA::LibraryRefiner::ObsMap obs;
      if (!tune_only) { obs = ODIA::LibraryRefiner::readObservations(ids, p, st); }

#ifdef DIALIBGEN_WITH_FINETUNE
      if (tune)
      {
        namespace fs = std::filesystem;
        std::string models = getStringOption_("tune_models");
        if (models.empty())
        { if (const char* e = std::getenv("DIALIBGEN_MODEL_DIR"); e && *e) { models = e; } }
        if (models.empty()) { models = bundledModelDir(); }
        if (models.empty())
        { throw std::runtime_error("-tune needs -tune_models (or $DIALIBGEN_MODEL_DIR): the stock peptdeep models to start from"); }

        const std::string heads = getStringOption_("tune_heads");
        const bool want_rt = heads != "ccs", want_ccs = heads != "rt";

        // A foreign library can carry mass-only modifications that the model
        // cannot encode. Reject these before training either head.
        std::unordered_set<std::uint32_t> checked;
        for (const auto sequence_id : library.precursors().modified_sequence)
        {
          if (!checked.insert(sequence_id).second) { continue; }
          const std::string sequence(library.strings().get(sequence_id));
          try { (void)ODIA::PeptDeepEncoder::encode(OpenMS::AASequence::fromString(sequence)); }
          catch (const std::exception& error)
          { throw std::runtime_error("tuning cannot encode library precursor '" + sequence + "': " + error.what()); }
        }

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
        for (const char* head : {"rt", "ccs"})
        {
          if ((std::string(head) == "rt" && !want_rt) || (std::string(head) == "ccs" && !want_ccs)) { continue; }
          for (const char* suffix : {".onnx", ".onnx.tune.json", ".onnx.trajectory.tsv"})
          {
            const fs::path file = work / (std::string("peptdeep_") + head + "_dynamic" + suffix);
            if (fs::exists(file) || fs::is_symlink(file)) { throw std::runtime_error("refusing to overwrite tuned model artifact: " + file.string()); }
          }
        }

        ODIA::tune::TuneParams tp;
        tp.report = ids;
        tp.q_value = getDoubleOption_("filter:q_value");
        tp.min_charge = getIntOption_("filter:min_charge");
        tp.allow_z1 = getFlag_("filter:allow_z1");
        tp.rt_spread_max = getDoubleOption_("filter:rt_spread_max");
        tp.rt_max_minutes = getDoubleOption_("filter:rt_max_minutes");
        tp.train_size = static_cast<std::size_t>(std::max(0, getIntOption_("cohort:train_size")));
        tp.train_frac = getDoubleOption_("cohort:train_frac");
        tp.inner_val = !getFlag_("cohort:no_inner_val");
        tp.full_fit = getFlag_("cohort:full_fit");
        tp.epochs = getIntOption_("train:epochs");
        tp.warmup = getIntOption_("train:warmup");
        tp.lr = getDoubleOption_("train:lr");
        tp.batch_size = getIntOption_("train:batch_size");
        tp.eval_every = getIntOption_("stop:eval_every");
        tp.min_epochs = getIntOption_("stop:min_epochs");
        tp.patience = getIntOption_("stop:patience");
        tp.rel_tol = getDoubleOption_("stop:rel_tol");
        tp.abs_tol = getDoubleOption_("stop:abs_tol");
        tp.max_seconds = getDoubleOption_("stop:max_seconds");
        tp.select = getStringOption_("stop:select") == "rmse" ? ODIA::tune::Select::Rmse : ODIA::tune::Select::CalibratedSd;
        tp.device = getStringOption_("machine:device");
        tp.threads = getIntOption_("machine:threads");
        tp.cudnn = !getFlag_("machine:no_cudnn");
        tp.seed = static_cast<std::uint32_t>(getIntOption_("machine:seed"));
        if (tp.epochs < 1 || tp.warmup < 0 || tp.warmup > tp.epochs)
        { throw std::runtime_error("train:warmup must be in [0, train:epochs]"); }
        if (tp.train_size && tp.train_frac > 0)
        { throw std::runtime_error("give cohort:train_size or cohort:train_frac, not both"); }
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

        const bool free_cys = getFlag_("tune_keep_free_cysteine_offset");
        const bool gpu = getFlag_("tune_predict_gpu");
        const unsigned sessions = static_cast<unsigned>(std::max(0, getIntOption_("tune_predict_sessions")));

        if (want_rt)
        {
          const ODIA::tune::TuneResult r = run_head(ODIA::tune::HeadKind::RT, "peptdeep_rt_dynamic.onnx");
          writeLogInfo_("tuned RT: " + r.stop_reason + ", best epoch " + std::to_string(r.best_epoch) +
                        " of " + std::to_string(r.epochs_run));
          // predictRetentionTimes returns what it could NOT predict, not what it
          // did -- undocumented, and it reads exactly the other way round.
          const std::size_t unpredicted =
            ODIA::LibraryGenerator::predictRetentionTimes(library, r.model_out, gpu, sessions, free_cys);
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
          const ODIA::tune::TuneResult r = run_head(ODIA::tune::HeadKind::CCS, "peptdeep_ccs_dynamic.onnx");
          writeLogInfo_("tuned CCS: " + r.stop_reason + ", best epoch " + std::to_string(r.best_epoch) +
                        " of " + std::to_string(r.epochs_run));
          // derive_mobility rewrites the WHOLE 1/K0 column, so it runs only when
          // the CCS head actually tuned -- otherwise a library that arrived with
          // measured mobilities would lose them to stock predictions.
          const std::size_t unpredicted =
            ODIA::LibraryGenerator::predictCollisionCrossSections(library, r.model_out, gpu, sessions, true);
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

      ODIA::LibraryRefiner::refine(library, obs, p, st);
      if (st.lib_unknown_mod_tokens)
      { writeLogWarn_(std::to_string(st.lib_unknown_mod_tokens) + " library modification tokens could not be resolved to UniMod; passed through verbatim and may fail the reference join"); }
      }
    }
    catch (const std::exception& e) { writeLogError_(std::string("refine: ") + e.what()); keptReport(); return UNEXPECTED_RESULT; }

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
      { writeLogWarn_(std::to_string(st.intensity_loss_bearing) + " neutral-loss transitions could not be matched: the report has "
                      "no loss field. Under the default restriction they were dropped from replaced precursors."); }
      if (!p.filter)
      { writeLogWarn_("MIXED INTENSITY PROVENANCE: matched precursors carry this run's observed intensities, unmatched ones "
                      "carry MS2-model predictions (-allow_mixed_intensity)."); }
    }
    if (p.write_rt && st.rt_written)
    {
      writeLogInfo_("NOTE: the RT column now holds the REFERENCE RUN's observed retention times (unit: " +
                    std::string(rtUnitName(p.rt_unit)) + "), not iRT. This library is a per-run object.");
    }
    }

    // Provenance: the recipe, the inputs by content hash, the run, the units and every
    // count above -- embedded in the Parquet and always written as a sidecar (M13).
    json prov = {
      {"tool", "DIALibGen"}, {"tool_version", DIALIBGEN_VERSION}, {"mode", tune_only ? "tune" : "refine"},
      {"config", eff},
      {"inputs", {{"library", std::filesystem::absolute(in).string()}, {"library_fnv1a64", ODIA::DIANNLibraryFile::hashFile(in)},
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

    try
    {
      ODIA::AtomicFile output(out), sidecar(out + ".refine.json");
      if (out.ends_with(".parquet"))
      {
        // The library embeds what is a function of the inputs and settings;
        // the search's timings and thread count stay in the sidecar, so the
        // same -run command writes the same library bytes.
        json embedded = prov;
        if (embedded.contains("search") && embedded["search"].is_object()) { embedded["search"].erase("resources"); }
        ODIA::DIANNLibraryFile::storeParquetCompact(output.temporaryPath().string(), library, ODIA::DIANNLibraryFile::Fingerprint{}, embedded.dump());
      }
      else if (out.ends_with(".tsv"))
      { ODIA::DIANNLibraryFile::storeTSV(output.temporaryPath().string(), library); }
      else { writeLogError_("-out must end in .parquet or .tsv"); return ILLEGAL_PARAMETERS; }
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
