// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

#include "DIALibGen.h"
#include <odia/DIANNLibraryFile.h>
#include <odia/AtomicFile.h>
#include <odia/Library.h>
#include <odia/LibraryGenerator.h>
#include <odia/PeptDeepEncoder.h>
#include <OpenMS/SYSTEM/File.h>
#include <filesystem>
#include <fstream>
#include <cstdlib>
#include <cstdint>
#include <iostream>
#include <thread>
#include <type_traits>
#include <algorithm>
using json = nlohmann::json;
namespace fs = std::filesystem;
namespace
{
  constexpr const char* kRtModelFile = "peptdeep_rt_dynamic.onnx";
  constexpr const char* kMs2ModelFile = "peptdeep_ms2_dynamic.onnx";
  constexpr const char* kCcsModelFile = "peptdeep_ccs_dynamic.onnx";

  constexpr const char* kDecoyMethods[] = {"none", "mutate", "pseudo_reverse",
                                           "reverse", "shuffle"};

  constexpr int kSchemaVersion = 1;

  fs::path executableDir()
  {
    const std::string p = OpenMS::File::getExecutablePath();
    return fs::path(p).lexically_normal();
  }

  std::vector<fs::path> dataDirCandidates()
  {
    std::vector<fs::path> out;
    if (const char* env = std::getenv("DIALIBGEN_DATA_DIR"); env != nullptr && *env != '\0')
    {
      out.emplace_back(env);
    }
    const fs::path exe = executableDir();
    out.push_back(exe / ".." / "share" / "DIALibGen");  // installed, and the bundles
    out.push_back(exe / ".." / "data");                           // an uninstalled build tree
#ifdef ODIA_DATA_DIR
    out.emplace_back(ODIA_DATA_DIR);
#endif
    for (auto& p : out) { p = p.lexically_normal(); }
    return out;
  }

  std::string findDataFile(const char* name)
  {
    std::error_code ec;
    for (const auto& dir : dataDirCandidates())
    {
      const fs::path c = dir / name;
      if (fs::exists(c, ec)) { return c.string(); }
    }
    return {};
  }

  std::vector<fs::path> modelDirCandidates()
  {
    std::vector<fs::path> out;
    if (const char* env = std::getenv("DIALIBGEN_MODEL_DIR"); env != nullptr && *env != '\0')
    {
      out.emplace_back(env);
    }
    const fs::path exe = executableDir();
    out.push_back(exe / ".." / "share" / "DIALibGen" / "models");
    try
    {
      out.push_back(fs::path(std::string(OpenMS::File::getOpenMSDataPath())) / "models");
    }
    catch (...)
    {
    }
    for (auto& p : out) { p = p.lexically_normal(); }
    return out;
  }

  std::string findModel(const char* name)
  {
    std::error_code ec;
    for (const auto& dir : modelDirCandidates())
    {
      const fs::path c = dir / name;
      if (fs::exists(c, ec)) { return c.string(); }
    }
    return {};
  }

  json effectiveConfig(const ODIA::DigestParams& p, const std::string& decoys,
                       const std::string& rt_model, const std::string& ms2_model,
                       const std::string& ccs_model, double nce,
                       const std::string& instrument, bool irt_rescale,
                       bool recompute_decoy_mz)
  {
    return json{
      {"schema_version", kSchemaVersion},
      {"enzyme", p.enzyme},
      {"missed_cleavages", p.missed_cleavages},
      {"peptide_length", {p.min_length, p.max_length}},
      {"precursor_charges", p.charges},
      {"precursor_mz", {p.precursor_mz_min, p.precursor_mz_max}},
      {"fragment_mz", {p.fragment_mz_min, p.fragment_mz_max}},
      {"max_fragment_charge", p.max_fragment_charge},
      {"fragments", {p.min_fragments, p.max_fragments}},
      {"fixed_modifications", p.fixed_modifications},
      {"variable_modifications", p.variable_modifications},
      {"max_variable_modifications", p.max_variable_modifications},
      {"n_terminal_methionine_excision", p.n_terminal_methionine_excision},
      {"free_cysteine_rt_correction", p.free_cysteine_rt_correction},
      {"derive_ion_mobility", p.derive_ion_mobility},
      {"min_relative_intensity", p.min_relative_intensity},
      {"reserved_doubly_charged", p.reserved_doubly_charged},
      {"decoys", decoys},
      {"rt_model", rt_model}, {"ms2_model", ms2_model}, {"ccs_model", ccs_model},
      {"instrument", instrument}, {"nce", nce},
      {"irt_rescale", irt_rescale},
      {"recompute_decoy_mz", recompute_decoy_mz}};
  }
}


namespace
{
  json generationDefaults()
  {
    json j = effectiveConfig(ODIA::DigestParams{}, "none", "", "", "", 30, "QE", false, false);
    j.erase("schema_version");
    j["nce"] = -1.0; // TOPP sentinel: use the instrument default.
    return j;
  }
  void apply_(const json& j, ODIA::DigestParams& p, std::string& decoys,
              std::string& rt_model, std::string& ms2_model, std::string& ccs_model,
              double& nce, std::string& instrument, bool& irt_rescale,
              bool& recompute_decoy_mz, bool& nce_was_set, std::string& instrument_alias_of)
  {
    if (!j.is_object()) { throw std::runtime_error("config must be a JSON object"); }
    const json ref = effectiveConfig(p, decoys, rt_model, ms2_model, ccs_model,
                                     nce, instrument, irt_rescale,
                                     recompute_decoy_mz);
    for (const auto& [k, v] : j.items())
    {
      if (k == "nce_source" || k == "instrument_named") { continue; }
      if (!ref.contains(k)) { throw std::runtime_error("unknown config key: " + k); }
      (void)v;
    }
    if (j.contains("schema_version"))
    {
      if (!j["schema_version"].is_number_integer())
      { throw std::runtime_error("schema_version must be an integer"); }
      if (const int v = j["schema_version"]; v != kSchemaVersion)
      {
        throw std::runtime_error("schema_version " + std::to_string(v) +
                                 " is not supported by this build (expected " +
                                 std::to_string(kSchemaVersion) + ")");
      }
    }
    auto count = [&](const char* k, std::size_t& out, std::size_t max) {
      if (!j.contains(k)) { return; }
      if (!j[k].is_number_unsigned() || j[k].get<std::uint64_t>() > max)
      {
        throw std::runtime_error(std::string(k) + " must be a whole number in 0.." +
                                 std::to_string(max));
      }
      out = j[k].get<std::size_t>();
    };
    auto range = [&](const char* k, auto& lo, auto& hi, double max) {
      if (!j.contains(k)) { return; }
      if (!j[k].is_array() || j[k].size() != 2)
      { throw std::runtime_error(std::string(k) + " must be [min, max]"); }
      for (const auto& v : j[k])
      {
        if (!v.is_number() || v.get<double>() < 0.0 || v.get<double>() > max)
        {
          throw std::runtime_error(std::string(k) + " values must be between 0 and " +
                                   std::to_string(static_cast<long long>(max)));
        }
      }
      if constexpr (std::is_integral_v<std::remove_reference_t<decltype(lo)>>)
      {
        for (const auto& v : j[k])
        { if (!v.is_number_integer()) { throw std::runtime_error(std::string(k) + " values must be integers"); } }
      }
      lo = j[k][0]; hi = j[k][1];
      if (!(lo < hi)) { throw std::runtime_error(std::string(k) + " needs min < max"); }
    };
    if (j.contains("enzyme")) { p.enzyme = j["enzyme"]; }
    count("missed_cleavages", p.missed_cleavages, 10);
    range("peptide_length", p.min_length, p.max_length, 200);
    range("precursor_mz", p.precursor_mz_min, p.precursor_mz_max, 100000);
    range("fragment_mz", p.fragment_mz_min, p.fragment_mz_max, 100000);
    range("fragments", p.min_fragments, p.max_fragments, 1000);
    if (j.contains("precursor_charges"))
    {
      if (!j["precursor_charges"].is_array())
      { throw std::runtime_error("precursor_charges must be an integer list"); }
      for (const auto& value : j["precursor_charges"])
      { if (!value.is_number_integer()) { throw std::runtime_error("precursor_charges must be integers"); } }
      p.charges = j["precursor_charges"].get<std::vector<int>>();
      for (const int z : p.charges)
      {
        if (z < 1 || z > 8)
        { throw std::runtime_error("precursor_charges must each be between 1 and 8"); }
      }
    }
    if (j.contains("max_fragment_charge"))
    {
      if (!j["max_fragment_charge"].is_number_integer() ||
          j["max_fragment_charge"] < 1 || j["max_fragment_charge"] > 2)
      { throw std::runtime_error("max_fragment_charge must be between 1 and 2"); }
      p.max_fragment_charge = j["max_fragment_charge"];
    }
    if (j.contains("fixed_modifications"))
    { p.fixed_modifications = j["fixed_modifications"].get<std::vector<std::string>>(); }
    if (j.contains("variable_modifications"))
    { p.variable_modifications = j["variable_modifications"].get<std::vector<std::string>>(); }
    count("max_variable_modifications", p.max_variable_modifications, 10);
    if (j.contains("n_terminal_methionine_excision"))
    { p.n_terminal_methionine_excision = j["n_terminal_methionine_excision"]; }
    if (j.contains("min_relative_intensity"))
    {
      if (!j["min_relative_intensity"].is_number() ||
          j["min_relative_intensity"] < 0.0 || j["min_relative_intensity"] > 1.0)
      { throw std::runtime_error("min_relative_intensity must be between 0 and 1"); }
      p.min_relative_intensity = j["min_relative_intensity"];
    }
    if (j.contains("derive_ion_mobility"))
    { p.derive_ion_mobility = j["derive_ion_mobility"]; }
    if (j.contains("free_cysteine_rt_correction"))
    { p.free_cysteine_rt_correction = j["free_cysteine_rt_correction"]; }
    count("reserved_doubly_charged", p.reserved_doubly_charged, 100);
    if (j.contains("decoys")) { decoys = j["decoys"]; }
    if (j.contains("rt_model")) { rt_model = j["rt_model"]; }
    if (j.contains("ms2_model")) { ms2_model = j["ms2_model"]; }
    if (j.contains("ccs_model")) { ccs_model = j["ccs_model"]; }
    if (j.contains("instrument")) { instrument = j["instrument"]; }
    if (j.contains("nce"))
    {
      if (!j["nce"].is_number() || j["nce"] <= 0 || j["nce"] > 100)
      { throw std::runtime_error("nce must be greater than 0 and at most 100"); }
      nce = j["nce"]; nce_was_set = true;
    }
    if (j.contains("irt_rescale")) { irt_rescale = j["irt_rescale"]; }
    if (j.contains("recompute_decoy_mz"))
    { recompute_decoy_mz = j["recompute_decoy_mz"]; }
    if (p.charges.empty()) { throw std::runtime_error("precursor_charges must not be empty"); }
    auto charges = p.charges;
    std::sort(charges.begin(), charges.end());
    if (std::adjacent_find(charges.begin(), charges.end()) != charges.end())
    { throw std::runtime_error("precursor_charges must not contain duplicates"); }

    bool known = false;
    for (const char* m : kDecoyMethods) { known = known || decoys == m; }
    if (!known)
    {
      std::string all;
      for (const char* m : kDecoyMethods) { all += (all.empty() ? "" : ", "); all += m; }
      throw std::runtime_error("unknown decoys method: '" + decoys + "' (known: " + all + ")");
    }

    const std::string canonical = ODIA::PeptDeepEncoder::canonicalInstrument(instrument);
    if (canonical.empty())
    {
      throw std::runtime_error(
        "unknown instrument: '" + instrument + "'. Known: QE, Lumos, timsTOF, SciexTOF, ThermoTOF, and the "
        "aliases upstream groups onto them -- Astral, Fusion, Eclipse, Velos, Elite and the Tribrids are "
        "Lumos; QE+, QEHF, QEHFX, Q Exactive and Exploris are QE; timsTOF Pro/SCP/HT/Ultra/flex are "
        "timsTOF; TripleTOF and ZenoTOF are SciexTOF. A leading 'Orbitrap' and a trailing model number "
        "are ignored, so 'Orbitrap Exploris 480' is accepted. An unrecognised name would index a slot no "
        "training addressed, which is close to Lumos in practice but would be recorded as though you had "
        "chosen an instrument. Use Lumos to ask for no instrument correction deliberately.");
    }
    if (canonical != instrument) { instrument_alias_of = instrument; }
    instrument = canonical;
  }

}

std::string DIALibGen::bundledModelDir()
{
  const std::string model = findModel(kRtModelFile);
  return model.empty() ? std::string() : fs::path(model).parent_path().string();
}

void DIALibGen::registerGenerationOptions_()
{
  registerTOPPSubsection_("generation", "FASTA digestion, prediction and decoys.");
  const auto defaults = generationDefaults();
  for (const auto& [key, value] : defaults.items())
  {
    const std::string name = "generation:" + key;
    std::string description = key;
    std::replace(description.begin(), description.end(), '_', ' ');
    if (key == "nce") { description = "Normalized collision energy (>0 and <=100). -1 selects the instrument default."; }
    if (value.is_boolean())
    {
      registerStringOption_(name, "<true/false>", value.get<bool>() ? "true" : "false", description, false);
      setValidStrings_(name, {"true", "false"});
    }
    else if (value.is_string())
    { registerStringOption_(name, "<value>", value.get<std::string>(), description, false); }
    else if (value.is_number_integer())
    { registerIntOption_(name, "<number>", value.get<int>(), description, false); setMinInt_(name, 0); }
    else if (value.is_number())
    { registerDoubleOption_(name, "<number>", value.get<double>(), description, false); setMinFloat_(name, key == "nce" ? -1 : 0); }
    else if (value.is_array())
    {
      if (key == "fixed_modifications" || key == "variable_modifications")
      {
        OpenMS::StringList values;
        for (const auto& item : value) { values.emplace_back(item.get<std::string>()); }
        registerStringList_(name, "<names>", values, description, false);
      }
      else if (key == "precursor_mz" || key == "fragment_mz")
      { registerDoubleList_(name, "<min max>", value.get<std::vector<double>>(), description, false); }
      else
      { registerIntList_(name, "<numbers>", value.get<std::vector<int>>(), description, false); }
    }
  }
  setValidStrings_("generation:decoys", {"none", "mutate", "pseudo_reverse", "reverse", "shuffle"});
  setMinInt_("generation:max_fragment_charge", 1);
  setMaxInt_("generation:max_fragment_charge", 2);
  setMinInt_("generation:precursor_charges", 1);
  setMaxInt_("generation:precursor_charges", 8);
  setMaxFloat_("generation:nce", 100);
  registerInputFile_("irt_standards", "<file>", "", "iRT calibration standards; defaults to the bundled table.", false);
  setValidFormats_("irt_standards", {"tsv"}, false);
}

json DIALibGen::generationValue_(const std::string& key, const json& value)
{
  const std::string name = "generation:" + key;
  if (value.is_boolean()) { return getStringOption_(name) == "true"; }
  if (value.is_string()) { return std::string(getStringOption_(name)); }
  if (value.is_number_unsigned()) { return static_cast<std::uint64_t>(getIntOption_(name)); }
  if (value.is_number_integer()) { return getIntOption_(name); }
  if (value.is_number()) { return getDoubleOption_(name); }
  if (key == "fixed_modifications" || key == "variable_modifications") { return getStringList_(name); }
  if (key == "precursor_mz" || key == "fragment_mz") { return getDoubleList_(name); }
  return getIntList_(name);
}
OpenMS::TOPPBase::ExitCodes DIALibGen::generate_()
  {
    ODIA::DigestParams p;                       // defaults live in DigestParams
    std::string decoys = "none", rt_model, ms2_model, ccs_model, instrument = "QE";
    bool recompute_decoy_mz = false;
    double nce = 30.0;
    bool irt_rescale = false;
    bool nce_was_set = false;
    std::string instrument_alias_of;   // the name the caller wrote, when it was an alias
    std::string nce_source = "config";

    const std::string cfg = getStringOption_("config");
    if (!cfg.empty())
    {
      std::ifstream in(cfg);
      if (!in) { writeLogError_("cannot read config: " + cfg); return INPUT_FILE_NOT_FOUND; }
      try { apply_(json::parse(in, nullptr, true, true), p, decoys, rt_model,
                   ms2_model, ccs_model, nce, instrument, irt_rescale,
                   recompute_decoy_mz, nce_was_set, instrument_alias_of); }
      catch (const std::exception& e)
      { writeLogError_(std::string("config: ") + e.what()); return ILLEGAL_PARAMETERS; }

      const auto base = fs::path(cfg).parent_path();
      for (std::string* m : {&rt_model, &ms2_model, &ccs_model})
      {
        if (!m->empty() && fs::path(*m).is_relative())
        { *m = (base / *m).lexically_normal().string(); }
      }
    }

    try
    {
      json overrides = json::object();
      const auto defaults = generationDefaults();
      for (const auto& [key, value] : defaults.items())
      {
        const std::string option = "generation:" + key;
        if (cfg.empty() || supplied_.count(option)) { overrides[key] = generationValue_(key, value); }
      }
      const bool automatic_nce = overrides.contains("nce") && overrides["nce"] == -1.0;
      if (automatic_nce) { overrides.erase("nce"); }
      else if (overrides.contains("nce")) { nce_source = "TOPP"; }
      if (overrides.contains("instrument")) { instrument_alias_of.clear(); }
      apply_(overrides, p, decoys, rt_model, ms2_model, ccs_model, nce, instrument,
             irt_rescale, recompute_decoy_mz, nce_was_set, instrument_alias_of);
      if (automatic_nce) { nce_was_set = false; }
    }
    catch (const std::exception& e)
    { writeLogError_(std::string("generation: ") + e.what()); return ILLEGAL_PARAMETERS; }

    if (!nce_was_set)
    {
      nce = ODIA::PeptDeepEncoder::defaultNce(instrument);
      nce_source = "instrument-default:" + instrument;
    }
    if (!instrument_alias_of.empty())
    { writeLogInfo_("instrument '" + instrument_alias_of + "' is '" + instrument + "' to this model (upstream's instrument_group)"); }
    if (!nce_was_set)
    {
      writeLogInfo_("nce not set; using " + std::to_string(nce) + " for " + instrument +
                    " (recorded as " + nce_source + ")");

    }
    if (instrument == "SciexTOF" || instrument == "ThermoTOF")
    { writeLogWarn_(instrument + " is in upstream's instrument list but carries no trained weights in the shipped "
                    "MS2 checkpoint (whose own constants name four instruments, and whose SciexTOF column sits at "
                    "its initialisation). Predictions will be close to Lumos, the no-correction baseline; name "
                    "Lumos if that is what you want."); }

    struct { const char* name; const char* file; std::string* path; } models[] = {
      {"rt_model", kRtModelFile, &rt_model},
      {"ms2_model", kMs2ModelFile, &ms2_model},
      {"ccs_model", kCcsModelFile, &ccs_model}};
    for (auto& m : models)
    {
      if (m.path->empty()) { *m.path = findModel(m.file); }
    }

    if (const std::string wc = getStringOption_("write_config"); !wc.empty())
    {
      json eff = effectiveConfig(p, decoys, rt_model, ms2_model, ccs_model,
                                 nce, instrument, irt_rescale,
                                 recompute_decoy_mz);
      eff["nce_source"] = nce_source;
      try
      {
        ODIA::AtomicFile staged(wc);
        std::ofstream os(staged.temporaryPath());
        os << eff.dump(2) << '\n'; os.close();
        if (!os) { throw std::runtime_error("cannot write config to " + wc); }
        staged.commit();
      }
      catch (const std::exception& e) { writeLogError_(e.what()); return CANNOT_WRITE_OUTPUT_FILE; }
      writeLogInfo_("wrote effective config to " + wc);
      return EXECUTION_OK;
    }

    const std::string fasta = getStringOption_("in");
    const std::string out = getParam_().getValue("out").toString();
    if (!out.empty() && (fs::exists(out) || fs::is_symlink(out)))
    { writeLogError_("refusing to overwrite output: " + out); return CANNOT_WRITE_OUTPUT_FILE; }
    (void)getStringOption_("out");
    if (fasta.empty() || out.empty())
    { writeLogError_("-in and -out are required"); return ILLEGAL_PARAMETERS; }
    const bool parquet = out.ends_with(".parquet");
    if (!parquet && !out.ends_with(".tsv"))
    {
      writeLogError_("-out must end in .parquet or .tsv (got '" + out + "')");
      return ILLEGAL_PARAMETERS;
    }

    json eff = effectiveConfig(p, decoys, rt_model, ms2_model, ccs_model,
                               nce, instrument, irt_rescale,
                               recompute_decoy_mz);
    eff["nce_source"] = nce_source;
    if (!instrument_alias_of.empty()) { eff["instrument_named"] = instrument_alias_of; }
    std::error_code ec;
    for (const auto& m : models)
    {
      if (m.path->empty())
      {
        std::ostringstream e;
        e << "no " << m.name << " configured and " << m.file << " was not found. "
          << "No tagged OpenMS release ships the PeptDeep models; set \"" << m.name
          << "\" in the config, or put all three .onnx files in one directory and "
          << "point DIALIBGEN_MODEL_DIR at it. Searched:";
        for (const auto& d : modelDirCandidates()) { e << "\n  " << d.string(); }
        writeLogError_(e.str());
        return INPUT_FILE_NOT_FOUND;
      }
      if (!fs::exists(*m.path, ec))
      {
        writeLogError_(std::string(m.name) + " does not exist: " + *m.path);
        return INPUT_FILE_NOT_FOUND;
      }
    }

    const unsigned threads = static_cast<unsigned>(std::max(0, getIntOption_("threads")));
    ODIA::Library library;
    std::string standards;
    try
    {
      const auto st = ODIA::LibraryGenerator::generate(fasta, p, library);
      writeLogInfo_("digest: " + std::to_string(st.proteins) + " proteins, " +
                    std::to_string(st.peptides) + " peptides, " +
                    std::to_string(st.precursors) + " precursors, " +
                    std::to_string(st.transitions) + " transitions");
      if (st.dropped_proteins || st.dropped_ambiguous_peptides)
      { writeLogWarn_("digest skipped " + std::to_string(st.dropped_proteins) +
          " invalid proteins and " + std::to_string(st.dropped_ambiguous_peptides) +
          " peptides with ambiguous or undefined residue masses"); }

      if (const auto miss = ODIA::LibraryGenerator::predictRetentionTimes(
            library, rt_model, true, threads, p.free_cysteine_rt_correction); miss)
      {
        writeLogError_(std::to_string(miss) + " precursors have no predicted RT. "
                       "Refusing rather than shipping a library with NaN retention "
                       "times, which extracts from nowhere and looks like a search bug.");
        return UNEXPECTED_RESULT;
      }
      if (irt_rescale)
      {
        standards = getStringOption_("irt_standards");
        if (standards.empty()) { standards = findDataFile("irt_standards.tsv"); }
        if (standards.empty())
        {
          std::ostringstream e;
          e << "no iRT standards found; pass -irt_standards. Searched:";
          for (const auto& d : dataDirCandidates()) { e << "\n  " << (d / "irt_standards.tsv").string(); }
          writeLogError_(e.str());
          return INPUT_FILE_NOT_FOUND;
        }
        const auto cal = ODIA::LibraryGenerator::fitIrtCalibration(rt_model, standards, true);
        ODIA::LibraryGenerator::applyIrtCalibration(library, cal);
        std::ostringstream m;
        m.setf(std::ios::fixed); m.precision(4);
        m << "rescaled to iRT: " << cal.slope << " * raw + " << cal.intercept
          << " from " << cal.peptides << " standards (worst off by "
          << cal.max_abs_error << " iRT)";
        writeLogInfo_(m.str());
      }
      else
      {
        writeLogInfo_("RT column is the model's RAW 0..1 output, NOT iRT "
                      "(irt_rescale=false). The per-run map is scale-invariant, so "
                      "this is the pipeline domain; it is NOT interchangeable with an "
                      "iRT library in another tool.");
      }

      if (const auto miss = ODIA::LibraryGenerator::predictFragmentIntensities(
            library, ms2_model, p, static_cast<float>(nce), instrument, true, threads); miss)
      { throw std::runtime_error(std::to_string(miss) + " precursors have no predicted MS2 spectrum"); }
      writeLogInfo_("after fragment selection: " + std::to_string(library.precursorCount()) + " precursors");
      if (library.precursorCount() == 0)
      { throw std::runtime_error("no precursors remain after digestion and fragment selection"); }
      if (const auto miss = ODIA::LibraryGenerator::predictCollisionCrossSections(
            library, ccs_model, true, threads, p.derive_ion_mobility); miss)
      { throw std::runtime_error(std::to_string(miss) + " precursors have no predicted CCS"); }
      std::size_t skipped = 0;
      const auto method = ODIA::parseDecoyMethod(decoys);
      const auto made = ODIA::LibraryGenerator::appendDecoys(library, method, &skipped,
                                                             p.min_fragments,
                                                             recompute_decoy_mz);
      writeLogInfo_("decoys: " + std::to_string(made) + " (" + decoys + "), " +
                    std::to_string(skipped) + " skipped");
    }
    catch (const std::exception& e)
    { writeLogError_(std::string("library generation failed: ") + e.what()); return UNEXPECTED_RESULT; }

    ODIA::DIANNLibraryFile::Fingerprint fp;
    try { fp = ODIA::DIANNLibraryFile::fingerprintFasta(fasta); }
    catch (const std::exception& e)
    { writeLogError_(std::string("cannot fingerprint the FASTA: ") + e.what()); return INPUT_FILE_NOT_FOUND; }
    fp.target_params = ODIA::LibraryGenerator::fingerprintParams(
      p, ODIA::DIANNLibraryFile::hashFile(rt_model),
      ODIA::DIANNLibraryFile::hashFile(ms2_model),
      ODIA::DIANNLibraryFile::hashFile(ccs_model), nce, instrument, irt_rescale);
    if (irt_rescale)
    { fp.target_params += ";irt_standards=" + ODIA::DIANNLibraryFile::hashFile(standards); }
    fp.decoy_method = decoys;
    fp.params = fp.target_params + ";decoy=" + fp.decoy_method +
                ";recompute_decoy_mz=" + (recompute_decoy_mz ? "1" : "0");

    try
    {
      if (parquet)
      { ODIA::DIANNLibraryFile::storeParquetCompact(out, library, fp, eff.dump()); }
      else
      {
        ODIA::DIANNLibraryFile::storeTSV(out, library);
        writeLogWarn_("TSV cannot carry the recipe: the config is NOT embedded and this "
                      "file is not cache-eligible. Use .parquet to keep them together.");
      }
    }
    catch (const std::exception& e)
    { writeLogError_(std::string("cannot write ") + out + ": " + e.what()); return CANNOT_WRITE_OUTPUT_FILE; }
    writeLogInfo_("wrote " + out);
    return EXECUTION_OK;
  }
