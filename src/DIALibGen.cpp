// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

#include "DIALibGen.h"
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/FORMAT/ParamXMLFile.h>
#include <OpenMS/APPLICATIONS/ParameterInformation.h>
#include <odia/AtomicFile.h>
#include <odia/DIANNLibraryFile.h>
#include <odia/DeviceRoute.h>
#include <odia/Hardware.h>
#include <odia/LibraryExport.h>
#include <odia/Trace.h>
#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
namespace fs = std::filesystem;

std::string DIALibGen::help_mode_;

namespace
{
  /// Every option the tool registers (advanced_ records them): unknownOptions
  /// checks the command line against it before TOPP parses it.
  std::set<std::string>& knownOptions()
  {
    static std::set<std::string> names;
    return names;
  }

  /// The options --help shows: per mode for "--help -mode X", and their union
  /// (key "") otherwise. Everything else is registered as advanced, which TOPP
  /// shows with --helphelp and passes on to -write_ini and -write_ctd.
  const std::map<std::string, std::set<std::string>>& visibleOptions()
  {
    static const std::map<std::string, std::set<std::string>> visible = [] {
      std::map<std::string, std::set<std::string>> v;
      v["generate"] = {"mode", "in", "out", "config", "generation:instrument", "generation:nce",
                       "generation:precursor_mz", "generation:fragment_mz", "generation:precursor_charges",
                       "generation:fixed_modifications", "generation:variable_modifications", "generation:enzyme",
                       "generation:missed_cleavages", "machine:device"};
      v["refine"] = {"mode", "in", "out", "config", "ids", "no_filter", "no_write_rt", "write_im", "write_intensity",
                     "q_precursor", "tune", "tune_heads", "tune_out_models", "machine:device"};
      v["tune"] = {"mode", "in", "out", "config", "ids", "tune_heads", "tune_out_models", "machine:device"};
      v["convert"] = {"mode", "in", "out", "out_format"};
      v["prepare"] = {"mode", "generation:peptide_length"};
      for (const char* mode : {"generate", "refine", "tune"}) { v[""].insert(v[mode].begin(), v[mode].end()); }
      return v;
    }();
    return visible;
  }

  /// Defaults that changed in 0.12. An INI written by an older version lists
  /// the old value (-write_ini writes every option), which the mode check must
  /// not mistake for a deliberate setting of another mode's option.
  struct ChangedDefault { const char* option; OpenMS::ParamValue old_value; const char* now; };
  const std::vector<ChangedDefault>& changedDefaults()
  {
    static const std::vector<ChangedDefault> table = {
      {"tune_heads", OpenMS::ParamValue(std::string("both")), "auto"},
      {"q_global", OpenMS::ParamValue(0.01), "1 (off)"},
      {"q_protein", OpenMS::ParamValue(0.01), "1 (off)"},
      {"generation:precursor_mz", OpenMS::ParamValue(std::vector<double>{350.0, 1200.0}), "from the instrument"},
      {"generation:fragment_mz", OpenMS::ParamValue(std::vector<double>{200.0, 1800.0}), "from the instrument"}};
    return table;
  }

  /// "0.11.5" < "0.12": the INI predates the 0.12 defaults. No version means
  /// a hand-written file, which is treated as old.
  bool olderThan012(const std::string& version)
  {
    if (version.empty()) { return true; }
    int major = 0, minor = 0;
    if (std::sscanf(version.c_str(), "%d.%d", &major, &minor) < 1) { return true; }
    return major == 0 && minor < 12;
  }

  std::string lower(std::string s)
  {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return s;
  }

  /// "-instance 01" is instance 1 to TOPP (a number).
  std::string normalInstance(const std::string& s)
  {
    try { std::size_t used = 0; const int n = std::stoi(s, &used); if (used == s.size()) { return std::to_string(n); } }
    catch (const std::exception&) {}
    return s;
  }

  /// A spectral library rather than a FASTA: .parquet or .tsv, unless the
  /// file starts with '>' (after a BOM or white space), as a FASTA does.
  bool isLibraryInput(const std::string& path)
  {
    const std::string name = lower(fs::path(path).filename().string());
    if (!(name.ends_with(".parquet") || name.ends_with(".tsv"))) { return false; }
    std::ifstream in(path, std::ios::binary);
    char c = 0;
    int skipped = 0;
    while (in.get(c) && skipped < 4096)
    {
      const unsigned char u = static_cast<unsigned char>(c);
      if (u == 0xEF || u == 0xBB || u == 0xBF || std::isspace(u)) { ++skipped; continue; }
      return c != '>';
    }
    return true;
  }

  std::size_t editDistance(const std::string& a, const std::string& b)
  {
    std::vector<std::size_t> row(b.size() + 1);
    for (std::size_t j = 0; j <= b.size(); ++j) { row[j] = j; }
    for (std::size_t i = 1; i <= a.size(); ++i)
    {
      std::size_t diagonal = row[0];
      row[0] = i;
      for (std::size_t j = 1; j <= b.size(); ++j)
      {
        const std::size_t up = row[j];
        row[j] = std::min({row[j] + 1, row[j - 1] + 1, diagonal + (a[i - 1] == b[j - 1] ? 0 : 1)});
        diagonal = up;
      }
    }
    return row[b.size()];
  }
}

const std::vector<DIALibGen::Removed>& DIALibGen::removedOptions()
{
  static const std::vector<Removed> table = {
    {"tune_models", "use -models"},
    {"machine:threads", "use -machine:train_threads"},
    {"tune_predict_gpu", "re-prediction follows -machine:device"},
    {"search:passes", "1 was its only value, so delete it"},
    {"search:allow_bootstrap", "set ODIA_SEARCH_ALLOW_BOOTSTRAP=1 instead"}};
  return table;
}

bool DIALibGen::advanced_(const std::string& name)
{
  knownOptions().insert(name);
  const auto& visible = visibleOptions();
  const auto mode = visible.find(help_mode_);
  return !(mode == visible.end() ? visible.at("") : mode->second).count(name);
}

DIALibGen::DIALibGen()
  : TOPPBase("DIALibGen",
             "Predict DIA spectral libraries from a FASTA and adapt them to a run. Uses all physical cores unless "
             "-threads is given. --help -mode <mode> shows one mode's options. Docs: "
             "https://github.com/okohlbacher/DIALibGen/blob/main/docs/usage.md",
             false,
             {{"Zeng WF et al.", "AlphaPeptDeep, whose models DIALibGen runs: a modular deep learning framework to "
               "predict peptide properties for proteomics", "Nat Commun 2022; 13: 7238", "10.1038/s41467-022-34904-3"}})
{
  version_ = DIALIBGEN_VERSION;
  verboseVersion_ = OpenMS::String(DIALIBGEN_VERSION) + " (OpenMS " + OpenMS::VersionInfo::getVersion() + ")";
}

void DIALibGen::registerOptionsAndFlags_()
{
  registerStringOption_("mode", "<mode>", "generate",
                        "Generate: FASTA to library; refine: write a run's observed RT (and 1/K0) into a library; tune: "
                        "fine-tune RT (and CCS) on a run and re-predict the whole library; convert: change the format; "
                        "prepare (macOS): precompile the GPU model. Without -mode: generate for a FASTA, convert for a "
                        "library.", false, advanced_("mode"));
  setValidStrings_("mode", {"generate", "refine", "tune", "convert", "prepare"});
  registerInputFile_("in", "<file>", "", "Protein FASTA (generate) or spectral library (refine, tune, convert)", false,
                     advanced_("in"));
  setValidFormats_("in", {"fasta", "parquet", "tsv"}, false);
  registerOutputFile_("out", "<file>", "", "Library to write: .tsv for DIA-NN; .pqp for OpenSWATH (add "
                      "-generation:decoys); .parquet for DIALibGen (DIA-NN cannot read it). Other layouts: -out_format",
                      false, advanced_("out"));
  setValidFormats_("out", {"parquet", "tsv", "pqp"}, false);
  registerStringOption_("out_format", "<format>", "auto", "Layout of -out. auto: by extension (.parquet dialibgen_parquet, "
                        ".tsv diann_tsv, .pqp openswath_pqp). diann2_parquet (.parquet): DIA-NN 2.0's own Parquet "
                        "library layout. openswath_tsv (.tsv): OpenSWATH's assay TSV, as TargetedFileConverter reads it",
                        false, advanced_("out_format"));
  {
    const auto names = ODIA::LibraryExport::formatNames();
    setValidStrings_("out_format", OpenMS::StringList(names.begin(), names.end()));
  }
  registerInputFile_("config", "<file>", "", "JSON recipe for the mode (-write_config writes one); options on the "
                     "command line or in the INI override it", false, advanced_("config"));
  setValidFormats_("config", {"json"}, false);
  registerOutputFile_("write_config", "<file>", "", "Write the mode's effective recipe (JSON) and exit", false,
                      advanced_("write_config"));
  setValidFormats_("write_config", {"json"}, false);
  registerStringOption_("models", "<dir>", "", "Directory of the stock peptdeep_{rt,ms2,ccs}_dynamic.onnx models for "
                        "generate, -tune and -run. Default: $DIALIBGEN_MODEL_DIR, else the bundled models. "
                        "-generation:{rt,ms2,ccs}_model name single files", false, advanced_("models"));
  registerGenerationOptions_();
  registerRefinementOptions_();
}

bool DIALibGen::given_(const std::string& name) const
{
  if (argv_supplied_.count(name)) { return true; }
  if (!supplied_.count(name)) { return false; }
  const auto& info = findEntry_(name);
  const auto default_value = info.type == OpenMS::ParameterInformation::FLAG ? OpenMS::ParamValue("false") : info.default_value;
  return getParam_().getValue(name) != default_value;
}

std::string DIALibGen::modelDirectory_() const
{
  std::string models = getParam_().getValue("models").toString();
  if (models.empty()) { if (const char* e = std::getenv("DIALIBGEN_MODEL_DIR"); e && *e) { models = e; } }
  if (models.empty()) { models = bundledModelDir(); }
  return models;
}

std::vector<std::pair<std::string, std::string>> DIALibGen::unknownOptions(int argc, const char** argv)
{
  DIALibGen probe;
  probe.registerOptionsAndFlags_();
  std::set<std::string> known = knownOptions();
  // TOPP's own (TOPPBase::main registers them after the tool's).
  for (const char* name : {"ini", "log", "instance", "debug", "threads", "write_ini", "write_ctd", "write_nested_cwl",
                           "write_cwl", "write_nested_json", "write_json", "no_progress", "force", "test", "-help",
                           "-helphelp"})
  { known.insert(name); }
  std::vector<std::pair<std::string, std::string>> unknown;
  for (int i = 1; i < argc; ++i)
  {
    const std::string arg = argv[i];
    const bool option = arg.size() >= 2 && arg[0] == '-' &&
                        (std::isalpha(static_cast<unsigned char>(arg[1])) ||
                         (arg[1] == '-' && arg.size() >= 3 && std::isalpha(static_cast<unsigned char>(arg[2]))));
    if (!option || known.count(arg.substr(1))) { continue; }
    std::string best;
    std::size_t best_distance = 4;   // suggest only what is close
    for (const auto& name : known)
    {
      if (name.starts_with("-")) { continue; }
      const std::size_t d = editDistance(arg.substr(1), name);
      if (d < best_distance) { best_distance = d; best = name; }
    }
    // A known name after the section: -instrument -> -generation:instrument.
    if (best.empty())
    {
      for (const auto& name : known)
      { if (name.ends_with(":" + arg.substr(1))) { best = name; break; } }
    }
    unknown.emplace_back(arg, best);
  }
  return unknown;
}

OpenMS::TOPPBase::ExitCodes DIALibGen::main_(int argc, const char** argv)
{
  // ODIA_TRACE=<dir>: phase events and /proc samples for this run (docs/tracing.md).
  const ODIA::trace::Session trace_session(argc, argv);
  std::string ini, instance = "1", ini_version;
  for (int i = 1; i < argc; ++i)
  {
    const std::string token = argv[i];
    if (token.starts_with("-")) { supplied_.insert(token.substr(1)); argv_supplied_.insert(token.substr(1)); }
    if (i + 1 < argc && token == "-ini") { ini = argv[i + 1]; }
    if (i + 1 < argc && token == "-instance") { instance = normalInstance(argv[i + 1]); }
  }
  // A TOPP JSON parameter file is TOPP's to read: this XML reader cannot.
  if (!ini.empty() && !lower(ini).ends_with(".json"))
  {
    OpenMS::Param parameters;
    try { OpenMS::ParamXMLFile().load(ini, parameters); }
    catch (const std::exception& e) { writeLogError_("cannot read the INI " + ini + ": " + e.what()); return ILLEGAL_PARAMETERS; }
    for (auto it = parameters.begin(); it != parameters.end(); ++it)
    {
      const std::string name = it.getName();
      if (name == "DIALibGen:version") { ini_version = it->value.toString(); continue; }
      for (const auto& prefix : {"DIALibGen:" + instance + ":", std::string("common:DIALibGen:"), std::string("common:")})
      { if (name.starts_with(prefix)) { supplied_.insert(name.substr(prefix.size())); break; } }
    }
  }
  ODIA::setDiagnosticsLevel(debug_level_);
  const auto& parameters = getParam_();
  const bool old_ini = !ini.empty() && olderThan012(ini_version);
  old_ini_ = old_ini;

  // Inspect the raw paths before TOPP's file getters can probe them.
  const std::string write_config = parameters.getValue("write_config").toString();
  const bool with_run = !parameters.getValue("run").toString().empty();
  const bool with_ids = !parameters.getValue("ids").toString().empty();
  std::string mode = getStringOption_("mode");
  // No -mode on the command line (a stock INI's "generate" does not count):
  // a library as -in converts, as a FASTA generates.
  if (!argv_supplied_.count("mode") && mode == "generate" && write_config.empty())
  {
    const std::string in = parameters.getValue("in").toString();
    if (!in.empty() && isLibraryInput(in))
    {
      if (with_ids || with_run)
      {
        writeLogError_("-in " + in + " is a library and generate reads a FASTA: give -mode refine or -mode tune");
        return ILLEGAL_PARAMETERS;
      }
      mode = "convert";
      writeLogInfo_("no -mode: -in is a library, so -mode convert");
    }
  }
  mode_ = mode;

  // Before anything counts the cores: -threads 0 on macOS may mean the
  // performance cores only (every mode; it changes no output).
  ODIA::preferPerformanceCores(getStringOption_("machine:cores") == "performance");
  // Arrow's pools follow -threads, not the hardware-thread count.
  try { ODIA::DIANNLibraryFile::useThreads(static_cast<unsigned>(std::max(0, getIntOption_("threads")))); }
  catch (const std::exception& e)
  { warn_(std::string("cannot apply -threads to the Parquet reader and writer (") + e.what() + "); using their default"); }
  if (!write_config.empty() && (fs::exists(write_config) || fs::is_symlink(write_config)))
  { writeLogError_("refusing to overwrite config output: " + write_config); return CANNOT_WRITE_OUTPUT_FILE; }
  const bool generating = mode == "generate" || mode == "prepare";
  for (const auto& option : supplied_)
  {
    // Which cores -threads 0 counts: every mode.
    if (option == "machine:cores") { continue; }
    // The GPU options steer prediction in generate as well as in tune.
    const bool prediction_device = option == "machine:require_gpu" || option == "machine:tf32" ||
                                   option == "machine:gpu_sessions" || option == "machine:gpu_devices" ||
                                   option == "machine:device" || option == "machine:coreml_units";
    // prepare takes generate's options: the length window and the CoreML settings.
    if (prediction_device && generating) { continue; }
    // The stock models: generate, -tune and -run.
    if (option == "models" && mode != "convert" && (mode != "refine" || getFlag_("tune") || with_run)) { continue; }
    const bool training = option.starts_with("tune") || option.starts_with("filter:") || option == "models" ||
      option.starts_with("cohort:") || option.starts_with("train:") || option.starts_with("stop:") || option.starts_with("machine:");
    // The built-in identification's settings: refine and tune, and only with -run.
    const bool search = option.starts_with("search:");
    const bool refinement = refinement_options_.count(option) || training || search;
    const bool generation = option.starts_with("generation:") || option == "irt_standards" || option == "instrument_file";
    bool wrong_mode = generating ? refinement : generation;
    if (mode == "convert") { wrong_mode = refinement || generation || option == "config" || option == "write_config"; }
    if (mode == "refine" && training && !getFlag_("tune")) { wrong_mode = true; }
    if (mode == "tune" && refinement && !training && !search && option != "ids" && option != "out_report" &&
        option != "run" && option != "out_ids" && option != "no_filter" && option != "no_write_rt") { wrong_mode = true; }
    const bool needs_run = !wrong_mode && (search || option == "out_ids") && !with_run;
    // A generated INI contains every mode's defaults. Only reject an explicit
    // value that would change behavior if it belonged to the selected mode.
    if (!wrong_mode && !needs_run) { continue; }
    const auto& info = findEntry_(option);
    const auto default_value = info.type == OpenMS::ParameterInformation::FLAG ? OpenMS::ParamValue("false") : info.default_value;
    const auto& value = parameters.getValue(option);
    if (value == default_value) { continue; }
    // ... or the default of the version that wrote the INI.
    if (old_ini && !argv_supplied_.count(option) &&
        std::any_of(changedDefaults().begin(), changedDefaults().end(),
                    [&](const ChangedDefault& c) { return option == c.option && value == c.old_value; }))
    { continue; }
    if (needs_run)
    { writeLogError_("-" + option + " has no effect without -run; give -run or remove the option"); }
    else
    {
      writeLogError_("-" + option + " has no effect in -mode " + mode +
                     (mode == "refine" && training ? "; enable -tune or remove the training option" :
                                                    "; select the matching mode or remove the option"));
    }
    return ILLEGAL_PARAMETERS;
  }

  // An INI from before 0.12 carries the old defaults, and its values apply as
  // given: say which ones differ from today's, once.
  if (old_ini)
  {
    std::string listed;
    const bool tuning = mode == "tune" || (mode == "refine" && getFlag_("tune"));
    for (const auto& c : changedDefaults())
    {
      const std::string option = c.option;
      const bool relevant = option.starts_with("generation:") ? generating
                          : option == "tune_heads" ? tuning : mode == "refine";
      if (!relevant || argv_supplied_.count(option) || !supplied_.count(option) || parameters.getValue(option) != c.old_value)
      { continue; }
      listed += (listed.empty() ? "" : ", ") + ("-" + option + " " + c.old_value.toString() + " (now " + c.now + ")");
    }
    if (mode != "convert" && supplied_.count("threads") && !argv_supplied_.count("threads") && getIntOption_("threads") == 1)
    { listed += std::string(listed.empty() ? "" : ", ") + "-threads 1 (now all cores)"; }
    if (!listed.empty())
    {
      writeLogInfo_("the INI is from DIALibGen " + (ini_version.empty() ? std::string("before 0.12") : ini_version) +
                    " and keeps its values: " + listed + "; delete those lines for the 0.12 defaults");
    }
  }

  // A setting without an effect on this run is refused (Generate.cpp,
  // Refine.cpp); machine:* only warns: an INI written on a Mac must run on Linux.
  if (const std::string units = getStringOption_("machine:coreml_units"); units != "gpu" && mode != "convert")
  {
#ifdef __APPLE__
    const std::string device = getStringOption_("machine:device");
    if (device == "cpu" || device.starts_with("cuda"))
    { warn_("-machine:coreml_units " + units + " has no effect with -machine:device " + device); }
#else
    warn_("-machine:coreml_units " + units + " has no effect here: CoreML is macOS only");
#endif
  }

  if (mode == "convert") { return convert_(); }
  if (mode == "prepare") { return prepare_(); }
  return mode == "generate" ? generate_() : refine_(mode == "tune");
}

OpenMS::TOPPBase::ExitCodes DIALibGen::convert_()
{
  // Read the raw path, as refine does: TOPP's writable-file probe would follow
  // and then remove a dangling symlink before it could be refused.
  const std::string in = getParam_().getValue("in").toString();
  const std::string out = getParam_().getValue("out").toString();
  if (in.empty() || out.empty()) { writeLogError_("-mode convert needs -in and -out"); return ILLEGAL_PARAMETERS; }
  if (!fs::exists(in)) { writeLogError_("-in " + in + " does not exist"); return INPUT_FILE_NOT_FOUND; }
  if (fs::exists(out) || fs::is_symlink(out))
  { writeLogError_("refusing to overwrite existing output: " + out); return CANNOT_WRITE_OUTPUT_FILE; }
  ODIA::LibraryFormat format;
  try { format = ODIA::LibraryExport::resolve(out, getStringOption_("out_format")); }
  catch (const std::invalid_argument& e) { writeLogError_(e.what()); return ILLEGAL_PARAMETERS; }
  ODIA::Library library;
  try
  {
    ODIA::DIANNLibraryFile::load(in, library, static_cast<unsigned>(std::max(0, getIntOption_("threads"))));
    // A DIALibGen Parquet keeps its identity and recipe through a conversion:
    // conversion changes neither. A TSV has none to keep.
    ODIA::DIANNLibraryFile::Fingerprint fp;
    std::string config_json;
    std::vector<std::pair<std::string, std::string>> extra;
    if (lower(in).ends_with(".parquet"))
    {
      using F = ODIA::DIANNLibraryFile;
      fp.fasta_hash = F::readFingerprint(in, "odia.fasta_fnv1a64");
      const std::string bytes = F::readFingerprint(in, "odia.fasta_bytes");
      try { fp.fasta_bytes = bytes.empty() ? 0 : std::stoull(bytes); } catch (const std::exception&) { fp.fasta_bytes = 0; }
      fp.params = F::readFingerprint(in, "odia.params");
      fp.decoy_method = F::readFingerprint(in, "odia.decoy_method");
      const std::string target = F::readFingerprint(in, "odia.target_fingerprint");
      const std::string prefix = fp.fasta_hash + ":" + std::to_string(fp.fasta_bytes) + ":";
      fp.target_params = target.starts_with(prefix) ? target.substr(prefix.size()) : std::string();
      config_json = F::readFingerprint(in, "odia.config_json");
      if (const std::string w = F::readFingerprint(in, "odia.warnings"); !w.empty()) { extra.emplace_back("odia.warnings", w); }
    }
    if (format == ODIA::LibraryFormat::DIALibGenParquet)
    { ODIA::DIANNLibraryFile::storeParquetCompact(out, library, fp, config_json, extra); }
    else if (format == ODIA::LibraryFormat::DIANNTSV) { ODIA::DIANNLibraryFile::store(out, library); }
    else
    {
      auto metadata = extra;
      if (ODIA::LibraryExport::carriesRecipe(format) && !config_json.empty())
      { metadata.insert(metadata.begin(), {"odia.config_json", config_json}); }
      ODIA::LibraryExport::store(out, library, format, metadata);
    }
  }
  catch (const std::exception& e) { writeLogError_(e.what()); return INPUT_FILE_CORRUPT; }
  writeLogInfo_("converted " + std::to_string(library.precursorCount()) + " precursors (" +
                std::to_string(library.decoyCount()) + " decoys) and " +
                std::to_string(library.transitionCount()) + " transitions: " + in + " -> " + out +
                " (" + ODIA::LibraryExport::name(format) + ")");
  return EXECUTION_OK;
}

ODIA::GpuOptions DIALibGen::gpuOptions_(std::vector<std::string>* warnings)
{
  ODIA::GpuOptions gpu;
  gpu.require = getFlag_("machine:require_gpu");
  gpu.tf32 = getStringOption_("machine:tf32") == "true";
  gpu.sessions_per_device = getIntOption_("machine:gpu_sessions");
  // CoreML when asked for (generate decides per model, ODIA::chooseDeviceRoute).
  // mps names the training backend; its re-prediction goes through CoreML too.
  const std::string device = getStringOption_("machine:device");
  gpu.coreml = device == "coreml" || device == "mps";
  const std::string units = getStringOption_("machine:coreml_units");
  gpu.coreml_units = units == "all" ? "ALL" : units == "ane" ? "CPUAndNeuralEngine" : units == "cpu" ? "CPUOnly" : "CPUAndGPU";
  for (const int device : getIntList_("machine:gpu_devices")) { gpu.devices.push_back(device); }
  // cuda:N names the device for prediction too, unless -machine:gpu_devices does.
  if (const int n = ODIA::cudaDeviceIndex(device); gpu.devices.empty() && n >= 0) { gpu.devices.push_back(n); }
  gpu.warnings = warnings;
  return gpu;
}

namespace
{
  enum class DescriptorRequest { None, Ini, Ctd, CwlOrJson };

  bool given(int argc, const char** argv, const char* opt)
  {
    for (int i = 1; i < argc; ++i)
    { if (std::strcmp(argv[i], opt) == 0) { return true; } }
    return false;
  }

  std::string valueOf(int argc, const char** argv, const char* opt)
  {
    for (int i = 1; i + 1 < argc; ++i)
    { if (std::strcmp(argv[i], opt) == 0) { return argv[i + 1]; } }
    return {};
  }

  DescriptorRequest requestedDescriptor(int argc, const char** argv)
  {
    if (given(argc, argv, "-write_ini")) { return DescriptorRequest::Ini; }
    if (given(argc, argv, "-write_ctd")) { return DescriptorRequest::Ctd; }
    if (given(argc, argv, "-write_cwl") || given(argc, argv, "-write_nested_cwl") ||
        given(argc, argv, "-write_json") || given(argc, argv, "-write_nested_json"))
    { return DescriptorRequest::CwlOrJson; }
    return DescriptorRequest::None;
  }

  /// Whether the INI sets -threads (instance or common section): its value
  /// is then respected, 1 included.
  bool iniSetsThreads(const std::string& ini, const std::string& instance)
  {
    if (ini.empty()) { return false; }
    // A JSON parameter file: TOPP reads it; its -threads (or TOPP's 1) applies.
    if (lower(ini).ends_with(".json")) { return true; }
    try
    {
      OpenMS::Param parameters;
      OpenMS::ParamXMLFile().load(ini, parameters);
      return parameters.exists("DIALibGen:" + instance + ":threads") || parameters.exists("common:DIALibGen:threads") ||
             parameters.exists("common:threads");
    }
    catch (...) { return false; }   // TOPP reports the broken INI itself
  }

  /// -write_ini writes TOPP's registered defaults, -threads 1 among them. The
  /// written file gets 0 (all cores, the real default since 0.12) unless an
  /// -ini it was merged from set the count.
  void finishWrittenIni(const std::string& path, const std::string& merged_from, const std::string& instance)
  {
    OpenMS::Param parameters;
    OpenMS::ParamXMLFile().load(path, parameters);
    const std::string base = "DIALibGen:" + instance + ":";
    if (!iniSetsThreads(merged_from, instance) && parameters.exists(base + "threads"))
    {
      parameters.setValue(base + "threads", 0, "Number of threads; 0 = all physical cores (the default when the "
                          "command line and the INI do not set it)");
    }
    OpenMS::ParamXMLFile().store(path, parameters);
  }

  /// The removed names (DIALibGen::removedOptions) that @p argv or the INI
  /// gives, one message each.
  std::vector<std::string> removedIn(int argc, const char** argv, const std::string& ini, const std::string& instance)
  {
    std::vector<std::string> lines;
    const auto& removed = DIALibGen::removedOptions();
    for (int i = 1; i < argc; ++i)
    {
      for (const auto& r : removed)
      {
        if (argv[i] == "-" + std::string(r.old_name))
        { lines.push_back("-" + std::string(r.old_name) + " was removed in 0.13: " + r.use); }
      }
    }
    // A JSON parameter file is TOPP's to read; a broken INI TOPP reports itself.
    if (ini.empty() || lower(ini).ends_with(".json")) { return lines; }
    OpenMS::Param parameters;
    try { OpenMS::ParamXMLFile().load(ini, parameters); }
    catch (...) { return lines; }
    for (const auto& r : removed)
    {
      for (const std::string& key : {"DIALibGen:" + instance + ":" + r.old_name, "common:DIALibGen:" + std::string(r.old_name)})
      {
        if (!parameters.exists(key)) { continue; }
        lines.push_back(ini + ": " + r.old_name + " was removed in 0.13: " + r.use);
        break;
      }
    }
    return lines;
  }

  class ToolHandlerRegistration
  {
  public:
    ToolHandlerRegistration()
    {
      const char* existing = std::getenv("OPENMS_TTD_INTERNAL_PATH");
      if (existing != nullptr && *existing != '\0') { return; }

      std::error_code ec;
      const fs::path base = fs::temp_directory_path(ec);
      if (ec) { return; }
      try { staged_ = std::make_unique<ODIA::AtomicFile>(base / "DIALibGen.ttd"); }
      catch (const std::exception&) { return; }
      const fs::path target = staged_->temporaryPath();
      const fs::path dir = target.parent_path();

      static const char kTtd[] = "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n"
                                 "<ttd><tool status=\"internal\">"
                                 "<name>DIALibGen</name>"
                                 "<category/><type/></tool></ttd>\n";
      constexpr std::size_t kTtdLen = sizeof(kTtd) - 1;
      {
        std::ofstream os(target, std::ios::binary | std::ios::trunc);
        os.write(kTtd, static_cast<std::streamsize>(kTtdLen));
        os.close();
        if (!os) { staged_.reset(); return; }
      }

      const std::string dir_str = dir.string();
#ifdef _WIN32
      const bool set_ok = (_putenv_s("OPENMS_TTD_INTERNAL_PATH", dir_str.c_str()) == 0);
#else
      const bool set_ok = (::setenv("OPENMS_TTD_INTERNAL_PATH", dir_str.c_str(), 1) == 0);
#endif
      if (!set_ok) { staged_.reset(); }
    }

    ToolHandlerRegistration(const ToolHandlerRegistration&) = delete;
    ToolHandlerRegistration& operator=(const ToolHandlerRegistration&) = delete;

  private:
    std::unique_ptr<ODIA::AtomicFile> staged_;
  };
}

int main(int argc, const char** argv)
{
#ifdef _WIN32
  size_t sz = 0;
  if (getenv_s(&sz, nullptr, 0, "OPENMS_DISABLE_UPDATE_CHECK") != 0 || sz == 0)
  { _putenv_s("OPENMS_DISABLE_UPDATE_CHECK", "ON"); }
#else
  ::setenv("OPENMS_DISABLE_UPDATE_CHECK", "ON", 0);
#endif

  // The Core ML probe child (PeptDeepPredictor::probeCoremlMain): no TOPP
  // parsing, one line on stdout, the exit code is the verdict.
  if (argc == 4 && std::strcmp(argv[1], "--probe-coreml") == 0)
  { return ODIA::PeptDeepPredictor::probeCoremlMain(argv[2], argv[3]); }

  const bool help = given(argc, argv, "--help") || given(argc, argv, "--helphelp");
  // --help -mode X: X's options only (--helphelp: everything).
  if (given(argc, argv, "--help") && !given(argc, argv, "--helphelp")) { DIALibGen::setHelpMode(valueOf(argc, argv, "-mode")); }
  const DescriptorRequest descriptor = requestedDescriptor(argc, argv);
  const bool want_descriptor =
    !help && (descriptor == DescriptorRequest::Ctd || descriptor == DescriptorRequest::CwlOrJson);

  // A name removed in 0.13, on the command line or in the INI: one line
  // naming the replacement (TOPP would only call an INI's name unknown).
  if (!help && argc > 1)
  {
    const std::string ini = valueOf(argc, argv, "-ini");
    const std::string instance = normalInstance(valueOf(argc, argv, "-instance").empty() ? "1" : valueOf(argc, argv, "-instance"));
    std::vector<std::string> removed = removedIn(argc, argv, ini, instance);
    if (!removed.empty())
    {
      for (const auto& line : removed) { std::cerr << "Error: " << line << std::endl; }
      return OpenMS::TOPPBase::ILLEGAL_PARAMETERS;
    }
  }

  // An unknown option: one line, not the whole help that TOPP prints after it.
  if (!help && argc > 1)
  {
    const auto unknown = DIALibGen::unknownOptions(argc, argv);
    if (!unknown.empty())
    {
      for (const auto& [option, suggestion] : unknown)
      {
        std::cerr << "Error: Unknown option '" << option << "'"
                  << (suggestion.empty() ? std::string() : " (did you mean -" + suggestion + "?)") << std::endl;
      }
      std::cerr << "DIALibGen --help lists the main options, --help -mode <mode> one mode's, --helphelp all." << std::endl;
      return OpenMS::TOPPBase::ILLEGAL_PARAMETERS;
    }
  }

#ifndef ENABLE_TDL
  if (want_descriptor && descriptor == DescriptorRequest::CwlOrJson)
  {
    std::cerr << "DIALibGen: -write_cwl and -write_json need an OpenMS built with ENABLE_TDL=ON (this one is not); "
                 "use -write_ctd" << std::endl;
    return OpenMS::TOPPBase::UNKNOWN_ERROR;
  }
#endif

  // -threads defaults to every physical core (TOPP's default is 1): when
  // neither the command line nor the INI sets it, -threads 0 is appended
  // (a command-line value outranks the INI's, so the INI is looked at first).
  std::vector<const char*> args(argv, argv + argc);
  const std::string ini = valueOf(argc, argv, "-ini");
  const std::string instance = normalInstance(valueOf(argc, argv, "-instance").empty() ? "1" : valueOf(argc, argv, "-instance"));
  if (!help && descriptor == DescriptorRequest::None && argc > 1 && !given(argc, argv, "-threads") &&
      !iniSetsThreads(ini, instance))
  {
    args.push_back("-threads");
    args.push_back("0");
  }

  std::unique_ptr<ToolHandlerRegistration> ttd;
  if (want_descriptor) { ttd = std::make_unique<ToolHandlerRegistration>(); }

  DIALibGen tool;
  OpenMS::TOPPBase::ExitCodes rc;
  if (ttd)
  {
    try { rc = tool.main(static_cast<int>(args.size()), args.data()); }
    catch (const std::exception& e)
    {
      std::cerr << "DIALibGen: tool description export failed: " << e.what() << std::endl;
      rc = OpenMS::TOPPBase::UNKNOWN_ERROR;
    }
  }
  else
  {
    rc = tool.main(static_cast<int>(args.size()), args.data());
  }

  if (rc == OpenMS::TOPPBase::EXECUTION_OK && !help && descriptor == DescriptorRequest::Ini)
  {
    const std::string written = valueOf(argc, argv, "-write_ini");
    try { if (!written.empty()) { finishWrittenIni(written, ini, instance); } }
    catch (const std::exception& e)
    {
      // No half-finished INI: it would carry -threads 1.
      std::error_code ec;
      fs::remove(written, ec);
      std::cerr << "Error: -write_ini " << written << ": " << e.what() << std::endl;
      rc = OpenMS::TOPPBase::ILLEGAL_PARAMETERS;
    }
  }
  return rc;
}
