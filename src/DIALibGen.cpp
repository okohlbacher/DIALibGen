// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

#include "DIALibGen.h"
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/FORMAT/ParamXMLFile.h>
#include <OpenMS/APPLICATIONS/ParameterInformation.h>
#include <odia/AtomicFile.h>
#include <odia/DIANNLibraryFile.h>
#include <odia/Trace.h>
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
namespace fs = std::filesystem;

DIALibGen::DIALibGen()
  : TOPPBase("DIALibGen", "Generate, refine and tune DIA spectral libraries.", false,
      {{"Zeng WF et al.", "AlphaPeptDeep: a modular deep learning framework to predict peptide properties for proteomics",
        "Nat Commun 2022; 13: 7238", "10.1038/s41467-022-34904-3"}})
{
  version_ = DIALIBGEN_VERSION;
  verboseVersion_ = OpenMS::String(DIALIBGEN_VERSION) + " (OpenMS " + OpenMS::VersionInfo::getVersion() + ")";
}

void DIALibGen::registerOptionsAndFlags_()
{
  registerStringOption_("mode", "<mode>", "generate", "generate: FASTA to library; refine: apply observed values; tune: train RT/CCS models and re-predict the whole library; convert: rewrite a library in the format of -out, content unchanged.", false);
  setValidStrings_("mode", {"generate", "refine", "tune", "convert"});
  registerInputFile_("in", "<file>", "", "Protein FASTA (generate) or spectral library (refine/tune/convert).", false);
  setValidFormats_("in", {"fasta", "parquet", "tsv"}, false);
  registerOutputFile_("out", "<file>", "", "Output library. .tsv: DIA-NN's TSV library format, read by DIA-NN 1.7 and 2.x. .parquet: DIALibGen's own layout, which embeds provenance; DIA-NN does not read it.", false);
  setValidFormats_("out", {"parquet", "tsv"}, false);
  registerInputFile_("config", "<file>", "", "Optional JSON configuration for the selected mode. Explicit CLI/INI generation settings override JSON.", false);
  setValidFormats_("config", {"json"}, false);
  registerOutputFile_("write_config", "<file>", "", "Write the effective mode configuration and exit.", false);
  setValidFormats_("write_config", {"json"}, false);
  registerGenerationOptions_();
  registerRefinementOptions_();
}

OpenMS::TOPPBase::ExitCodes DIALibGen::main_(int argc, const char** argv)
{
  // ODIA_TRACE=<dir>: phase events and /proc samples for this run (docs/tracing.md).
  const ODIA::trace::Session trace_session(argc, argv);
  std::string ini, instance = "1";
  for (int i = 1; i < argc; ++i)
  {
    const std::string token = argv[i];
    if (token.starts_with("-")) { supplied_.insert(token.substr(1)); }
    if (i + 1 < argc && token == "-ini") { ini = argv[i + 1]; }
    if (i + 1 < argc && token == "-instance") { instance = argv[i + 1]; }
  }
  if (!ini.empty())
  {
    OpenMS::Param parameters;
    OpenMS::ParamXMLFile().load(ini, parameters);
    for (auto it = parameters.begin(); it != parameters.end(); ++it)
    {
      const std::string name = it.getName();
      for (const auto& prefix : {"DIALibGen:" + instance + ":", std::string("common:DIALibGen:"), std::string("common:")})
      { if (name.starts_with(prefix)) { supplied_.insert(name.substr(prefix.size())); break; } }
    }
  }
  const std::string mode = getStringOption_("mode");
  const auto& parameters = getParam_();
  // Arrow's pools follow -threads, not the hardware-thread count.
  try { ODIA::DIANNLibraryFile::useThreads(static_cast<unsigned>(std::max(0, getIntOption_("threads")))); }
  catch (const std::exception& e) { writeLogWarn_(std::string("cannot size Arrow's thread pools: ") + e.what()); }
  // Inspect the raw path before TOPP's output getter can probe a dangling symlink.
  const std::string write_config = parameters.getValue("write_config").toString();
  if (!write_config.empty() && (fs::exists(write_config) || fs::is_symlink(write_config)))
  { writeLogError_("refusing to overwrite config output: " + write_config); return CANNOT_WRITE_OUTPUT_FILE; }
  for (const auto& option : supplied_)
  {
    // The GPU options steer prediction in generate as well as in tune.
    const bool prediction_device = option == "machine:require_gpu" || option == "machine:tf32" ||
                                   option == "machine:gpu_sessions" || option == "machine:gpu_devices";
    if (prediction_device && mode == "generate") { continue; }
    const bool training = option.starts_with("tune") || option.starts_with("filter:") ||
      option.starts_with("cohort:") || option.starts_with("train:") || option.starts_with("stop:") || option.starts_with("machine:");
    const bool refinement = refinement_options_.count(option) || training;
    const bool generation = option.starts_with("generation:") || option == "irt_standards";
    bool wrong_mode = mode == "generate" ? refinement : generation;
    if (mode == "convert") { wrong_mode = refinement || generation || option == "config" || option == "write_config"; }
    if (mode == "refine" && training && !getFlag_("tune")) { wrong_mode = true; }
    if (mode == "tune" && refinement && !training && option != "ids" && option != "out_report" &&
        option != "no_filter" && option != "no_write_rt") { wrong_mode = true; }
    // A generated INI contains every mode's defaults. Only reject an explicit
    // value that would change behavior if it belonged to the selected mode.
    if (!wrong_mode) { continue; }
    const auto& info = findEntry_(option);
    const auto default_value = info.type == OpenMS::ParameterInformation::FLAG ? OpenMS::ParamValue("false") : info.default_value;
    if (parameters.getValue(option) != default_value)
    {
      writeLogError_("-" + option + " has no effect in -mode " + mode +
                     (mode == "refine" && training ? "; enable -tune or remove the training option" :
                                                    "; select the matching mode or remove the option"));
      return ILLEGAL_PARAMETERS;
    }
  }
  if (mode == "convert") { return convert_(); }
  return mode == "generate" ? generate_() : refine_(mode == "tune");
}

OpenMS::TOPPBase::ExitCodes DIALibGen::convert_()
{
  // Read the raw path, as refine does: TOPP's writable-file probe would follow
  // and then remove a dangling symlink before it could be refused.
  const std::string in = getStringOption_("in");
  const std::string out = getParam_().getValue("out").toString();
  if (in.empty() || out.empty()) { writeLogError_("-mode convert needs -in and -out"); return ILLEGAL_PARAMETERS; }
  if (fs::exists(out) || fs::is_symlink(out))
  { writeLogError_("refusing to overwrite existing output: " + out); return CANNOT_WRITE_OUTPUT_FILE; }
  ODIA::Library library;
  try
  {
    ODIA::DIANNLibraryFile::load(in, library);
    // Parquet output is the compact layout every other mode writes; it carries
    // no recipe, because conversion changes none.
    if (out.ends_with(".parquet"))
    { ODIA::DIANNLibraryFile::storeParquetCompact(out, library, ODIA::DIANNLibraryFile::Fingerprint{}); }
    else { ODIA::DIANNLibraryFile::store(out, library); }
  }
  catch (const std::exception& e) { writeLogError_(e.what()); return INPUT_FILE_CORRUPT; }
  writeLogInfo_("converted " + std::to_string(library.precursorCount()) + " precursors (" +
                std::to_string(library.decoyCount()) + " decoys) and " +
                std::to_string(library.transitionCount()) + " transitions: " + in + " -> " + out);
  return EXECUTION_OK;
}

ODIA::GpuOptions DIALibGen::gpuOptions_(std::vector<std::string>* warnings)
{
  ODIA::GpuOptions gpu;
  gpu.require = getFlag_("machine:require_gpu");
  gpu.tf32 = getStringOption_("machine:tf32") == "true";
  gpu.sessions_per_device = getIntOption_("machine:gpu_sessions");
  for (const int device : getIntList_("machine:gpu_devices")) { gpu.devices.push_back(device); }
  gpu.warnings = warnings;
  return gpu;
}

namespace
{
  enum class DescriptorRequest { None, Ini, Ctd, CwlOrJson };

  DescriptorRequest requestedDescriptor(int argc, const char** argv)
  {
    auto given = [argc, argv](const char* opt) {
      for (int i = 1; i < argc; ++i)
      { if (std::strcmp(argv[i], opt) == 0) { return true; } }
      return false;
    };
    if (given("-write_ini")) { return DescriptorRequest::Ini; }
    if (given("-write_ctd")) { return DescriptorRequest::Ctd; }
    if (given("-write_cwl") || given("-write_nested_cwl") ||
        given("-write_json") || given("-write_nested_json"))
    { return DescriptorRequest::CwlOrJson; }
    return DescriptorRequest::None;
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

  bool help = false;
  for (int i = 1; i < argc; ++i)
  { help = help || std::strcmp(argv[i], "--help") == 0 || std::strcmp(argv[i], "--helphelp") == 0; }
  const DescriptorRequest descriptor = requestedDescriptor(argc, argv);
  const bool want_descriptor =
    !help && (descriptor == DescriptorRequest::Ctd || descriptor == DescriptorRequest::CwlOrJson);

#ifndef ENABLE_TDL
  if (want_descriptor && descriptor == DescriptorRequest::CwlOrJson)
  {
    std::cerr << "DIALibGen: -write_cwl and -write_json need an OpenMS built with\n"
                 "                     ENABLE_TDL=ON; the OpenMS this binary links has it off.\n"
                 "                     Use -write_ctd instead." << std::endl;
    return OpenMS::TOPPBase::UNKNOWN_ERROR;
  }
#endif

  std::unique_ptr<ToolHandlerRegistration> ttd;
  if (want_descriptor) { ttd = std::make_unique<ToolHandlerRegistration>(); }

  DIALibGen tool;
  OpenMS::TOPPBase::ExitCodes rc;
  if (ttd)
  {
    try { rc = tool.main(argc, argv); }
    catch (const std::exception& e)
    {
      std::cerr << "DIALibGen: tool description export failed: " << e.what() << std::endl;
      rc = OpenMS::TOPPBase::UNKNOWN_ERROR;
    }
  }
  else
  {
    rc = tool.main(argc, argv);
  }

  return rc;
}
