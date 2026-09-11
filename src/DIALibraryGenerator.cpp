// Copyright (c) 2026, Oliver Kohlbacher and the DIALibraryGenerator authors.
// SPDX-License-Identifier: BSD-3-Clause

/// DIALibraryGenerator: FASTA + JSON config -> in-silico DIA library.
///
/// A library is a CROSS-RUN artefact built once and reused; a search is per-run.
/// They were welded together behind OpenDIAlyzer's -out_lib, which is why the
/// library's own parameters were scattered across ~15 CLI flags and why the
/// alkylation state -- the defect that cost a quarter of the search space --
/// was hard-coded and unreachable.
///
/// Every algorithm here already exists in odia_library. This file is a config
/// parser, a call sequence, and the standalone-tool scaffolding that OpenMS
/// does not supply to a tool living outside its own tree (version, update
/// check, -threads default, tool-description export, data and model lookup).
/// If it grows a second algorithm it is in the wrong file.

#include <odia/DIANNLibraryFile.h>
#include <odia/Library.h>
#include <odia/LibraryGenerator.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/SYSTEM/File.h>

#include <nlohmann/json.hpp>

#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <thread>
#include <vector>

using json = nlohmann::json;
namespace fs = std::filesystem;

namespace
{
  /// The three PeptDeep exports, by the file names OpenMS publishes them under.
  constexpr const char* kRtModelFile = "peptdeep_rt_dynamic.onnx";
  constexpr const char* kMs2ModelFile = "peptdeep_ms2_dynamic.onnx";
  constexpr const char* kCcsModelFile = "peptdeep_ccs_dynamic.onnx";

  /// Every decoy method the generator implements. Kept here, not inferred from
  /// parseDecoyMethod, because that function maps every name it does not know
  /// to Mutate: an unvalidated "reverese" produced mutation decoys and then
  /// recorded "reverese" in the provenance as though it were real.
  constexpr const char* kDecoyMethods[] = {"none", "mutate", "pseudo_reverse",
                                           "reverse", "shuffle"};

  /// The only schema this build can read. Accepted-and-ignored was worse than
  /// unsupported: a future config would have been silently misinterpreted.
  constexpr int kSchemaVersion = 1;

  /// Where the executable lives, without a trailing separator. OpenMS returns
  /// it with one and, on a path that is not resolvable, returns ".".
  fs::path executableDir()
  {
    const std::string p = OpenMS::File::getExecutablePath();
    return fs::path(p).lexically_normal();
  }

  /// Directories that may hold the shipped data tables (irt_standards.tsv), in
  /// priority order.
  ///
  /// The compile-time ODIA_DATA_DIR used to be the ONLY answer, which made an
  /// installed binary read its standards out of the source tree it was built
  /// from -- and made a relocated binary, a package, or a macOS .app read
  /// nothing at all. It is still consulted, last, because it is what makes an
  /// uninstalled build tree work.
  std::vector<fs::path> dataDirCandidates()
  {
    std::vector<fs::path> out;
    if (const char* env = std::getenv("DIALIBGEN_DATA_DIR"); env != nullptr && *env != '\0')
    {
      out.emplace_back(env);
    }
    const fs::path exe = executableDir();
    out.push_back(exe / ".." / "share" / "DIALibraryGenerator");  // installed, and the bundles
    out.push_back(exe / ".." / "data");                           // an uninstalled build tree
#ifdef ODIA_DATA_DIR
    out.emplace_back(ODIA_DATA_DIR);
#endif
    for (auto& p : out) { p = p.lexically_normal(); }
    return out;
  }

  /// First existing @p name across the data directories, or "".
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

  /// Directories that may hold the PeptDeep ONNX exports, in priority order.
  ///
  /// No tagged OpenMS release ships them (WITH_ONNX exists only on develop and
  /// defaults OFF), so this search usually comes up empty and the tool has to
  /// say so in words rather than dying inside an ONNX session constructor with
  /// "Load model from  failed".
  std::vector<fs::path> modelDirCandidates()
  {
    std::vector<fs::path> out;
    if (const char* env = std::getenv("DIALIBGEN_MODEL_DIR"); env != nullptr && *env != '\0')
    {
      out.emplace_back(env);
    }
    const fs::path exe = executableDir();
    out.push_back(exe / ".." / "share" / "DIALibraryGenerator" / "models");
    // Where OpenMS puts them when built WITH_ONNX=ON.
    try
    {
      out.push_back(fs::path(OpenMS::File::getOpenMSDataPath()) / "models");
    }
    catch (...)
    {
      // getOpenMSDataPath throws when share/OpenMS is unreachable. Not fatal
      // here: the other candidates may still hold the models.
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

  /// Defaults come from DigestParams so the two cannot drift. The JSON is
  /// EFFECTIVE: every default materialised, so the embedded recipe is complete
  /// and a reader never has to know our defaults to reproduce the library.
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
      // Recorded in the embedded recipe so a library states how its decoys were
      // built. NOTE, because an earlier comment here claimed more than this: it
      // does NOT protect cache reuse. The reuse decision is made by
      // fingerprintParams in LibraryGenerator.cpp, which does not include this
      // flag, and the embedded config is never read back to decide reuse.
      {"recompute_decoy_mz", recompute_decoy_mz}};
  }
}

class DIALibraryGenerator final : public OpenMS::TOPPBase
{
public:
  // official = false: this tool lives outside the OpenMS tree and so is not in
  // ToolHandler's list. Passing true makes the TOPPBase constructor throw.
  DIALibraryGenerator()
    : TOPPBase("DIALibraryGenerator",
               "Build an in-silico DIA spectral library from a FASTA.",
               false,
               {{"Zeng WF, Zhou XX, Willems S, Ammar C, Wahle M, Bludau I, Voytik E, "
                 "Strauss MT, Mann M",
                 "AlphaPeptDeep: a modular deep learning framework to predict peptide "
                 "properties for proteomics",
                 "Nat Commun 2022; 13: 7238", "10.1038/s41467-022-34904-3"}})
  {
    // TOPPBase prints OpenMS's version when version_ is empty. For a tool
    // shipped inside OpenMS that is right; for this one it means --help
    // announces the version of whatever OpenMS tree it happened to be built
    // against, and a user cannot report a bug against a version that does not
    // exist. DIALIBGEN_VERSION comes from project(DIALibraryGenerator VERSION).
#ifdef DIALIBGEN_VERSION
    version_ = DIALIBGEN_VERSION;
    verboseVersion_ = OpenMS::String(DIALIBGEN_VERSION) + " (OpenMS " +
                      OpenMS::VersionInfo::getVersion() + ")";
#endif
  }

protected:
  void registerOptionsAndFlags_() override
  {
    // Not required: -write_config and the tool-description exports are useful
    // on their own, and TOPPBase would otherwise refuse them for want of an
    // input file nothing reads. Enforced by hand in main_ instead.
    registerInputFile_("in", "<file>", "", "Protein FASTA.", false);
    setValidFormats_("in", {"fasta"}, false);
    registerInputFile_("config", "<file>", "",
                       "JSON library configuration. Every content-affecting parameter lives "
                       "here, so a library cannot be built with settings nobody wrote down. "
                       "Omit it for the defaults; -write_config dumps the effective config "
                       "including everything defaulted.", false);
    setValidFormats_("config", {"json"}, false);
    registerOutputFile_("out", "<file>", "",
                        "Library. .parquet carries the recipe in its schema metadata "
                        "(odia.config_json); .tsv is the DIA-NN dialect and CANNOT carry it "
                        "-- a TSV is not reproducible from itself. No other extension is "
                        "accepted.", false);
    setValidFormats_("out", {"parquet", "tsv"}, false);
    registerOutputFile_("write_config", "<file>", "",
                        "Write the effective config here and exit without building. The "
                        "quickest way to see what the defaults actually are.", false);
    setValidFormats_("write_config", {"json"}, false);
    registerInputFile_("irt_standards", "<file>", "",
                       "Standards defining the iRT scale. The RT model emits a raw 0..1 "
                       "value; without this the library ships THAT, and every consumer "
                       "then calibrates against a scale nobody declared. Defaults to the "
                       "shipped data/irt_standards.tsv.", false);
    setValidFormats_("irt_standards", {"tsv"}, false);
    // NOTE: -threads is NOT registered here. TOPPBase registers it itself,
    // after this function runs and with a default of 1, and registration only
    // appends -- a second entry of the same name is dead weight and TOPPBase's
    // wins at parse time. This tool's default of 0 (= all cores) is applied to
    // argv in main() instead. See the note printed after --help.
  }

  /// Strict: an unknown key is a typo that would otherwise silently become a
  /// default, and a library built from a typo looks exactly like one built
  /// correctly.
  void apply_(const json& j, ODIA::DigestParams& p, std::string& decoys,
              std::string& rt_model, std::string& ms2_model, std::string& ccs_model,
              double& nce, std::string& instrument, bool& irt_rescale,
              bool& recompute_decoy_mz)
  {
    const json ref = effectiveConfig(p, decoys, rt_model, ms2_model, ccs_model,
                                     nce, instrument, irt_rescale,
                                     recompute_decoy_mz);
    for (const auto& [k, v] : j.items())
    {
      if (!ref.contains(k)) { throw std::runtime_error("unknown config key: " + k); }
      (void)v;
    }
    // Refused, not ignored. A config written for a schema this build does not
    // know would otherwise be read with THIS build's meaning for every key.
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
    auto pair = [&](const char* k, auto& lo, auto& hi) {
      if (!j.contains(k)) { return; }
      if (!j[k].is_array() || j[k].size() != 2)
      { throw std::runtime_error(std::string(k) + " must be [min, max]"); }
      lo = j[k][0]; hi = j[k][1];
      if (!(lo < hi)) { throw std::runtime_error(std::string(k) + " needs min < max"); }
    };
    if (j.contains("enzyme")) { p.enzyme = j["enzyme"]; }
    if (j.contains("missed_cleavages")) { p.missed_cleavages = j["missed_cleavages"]; }
    pair("peptide_length", p.min_length, p.max_length);
    pair("precursor_mz", p.precursor_mz_min, p.precursor_mz_max);
    pair("fragment_mz", p.fragment_mz_min, p.fragment_mz_max);
    pair("fragments", p.min_fragments, p.max_fragments);
    if (j.contains("precursor_charges")) { p.charges = j["precursor_charges"].get<std::vector<int>>(); }
    if (j.contains("max_fragment_charge")) { p.max_fragment_charge = j["max_fragment_charge"]; }
    if (j.contains("fixed_modifications"))
    { p.fixed_modifications = j["fixed_modifications"].get<std::vector<std::string>>(); }
    if (j.contains("variable_modifications"))
    { p.variable_modifications = j["variable_modifications"].get<std::vector<std::string>>(); }
    if (j.contains("max_variable_modifications"))
    { p.max_variable_modifications = j["max_variable_modifications"]; }
    if (j.contains("n_terminal_methionine_excision"))
    { p.n_terminal_methionine_excision = j["n_terminal_methionine_excision"]; }
    if (j.contains("min_relative_intensity"))
    { p.min_relative_intensity = j["min_relative_intensity"]; }
    if (j.contains("derive_ion_mobility"))
    { p.derive_ion_mobility = j["derive_ion_mobility"]; }
    if (j.contains("free_cysteine_rt_correction"))
    { p.free_cysteine_rt_correction = j["free_cysteine_rt_correction"]; }
    if (j.contains("reserved_doubly_charged"))
    { p.reserved_doubly_charged = j["reserved_doubly_charged"]; }
    if (j.contains("decoys")) { decoys = j["decoys"]; }
    if (j.contains("rt_model")) { rt_model = j["rt_model"]; }
    if (j.contains("ms2_model")) { ms2_model = j["ms2_model"]; }
    if (j.contains("ccs_model")) { ccs_model = j["ccs_model"]; }
    if (j.contains("instrument")) { instrument = j["instrument"]; }
    if (j.contains("nce")) { nce = j["nce"]; }
    if (j.contains("irt_rescale")) { irt_rescale = j["irt_rescale"]; }
    if (j.contains("recompute_decoy_mz"))
    { recompute_decoy_mz = j["recompute_decoy_mz"]; }
    if (p.charges.empty()) { throw std::runtime_error("precursor_charges must not be empty"); }

    // parseDecoyMethod maps every name it does not know to Mutate, so an
    // unchecked typo produced mutation decoys AND wrote the typo into the
    // provenance as the method used.
    bool known = false;
    for (const char* m : kDecoyMethods) { known = known || decoys == m; }
    if (!known)
    {
      std::string all;
      for (const char* m : kDecoyMethods) { all += (all.empty() ? "" : ", "); all += m; }
      throw std::runtime_error("unknown decoys method: '" + decoys + "' (known: " + all + ")");
    }
  }

  ExitCodes main_(int, const char**) override
  {
    ODIA::DigestParams p;                       // defaults live in DigestParams
    // DEFAULT: no decoys. A generated library is an interchange artefact, and the
    // consumer decides its own null -- DIA-NN's README is explicit that it
    // "will search these decoys in addition to the regular decoys it generates",
    // so shipping ours DOUBLED its decoy population and made its FDR far more
    // conservative. ODIA appends its own on load when a library has none.
    std::string decoys = "none", rt_model, ms2_model, ccs_model, instrument = "QE";
    // Off pending its first measured arm, not because of a hazard: the decoy's
    // 1/K0 is re-derived from its CCS at the new mass, so it stays a physically
    // consistent ion inside the two-dimensional diaPASEF window.
    bool recompute_decoy_mz = false;
    double nce = 30.0;
    // RAW MODEL UNITS are the pipeline domain: the library carries the RT
    // model's own 0..1 output and the per-run map takes it to seconds. The map
    // is scale-invariant -- the constant cancels -- so the iRT rescale buys
    // nothing internally, and its 11 spiked standards had a worst-case error of
    // 8.76 iRT. Set true only to write a library another tool must read as iRT.
    bool irt_rescale = false;

    const std::string cfg = getStringOption_("config");
    if (!cfg.empty())
    {
      std::ifstream in(cfg);
      if (!in) { writeLogError_("cannot read config: " + cfg); return INPUT_FILE_NOT_FOUND; }
      try { apply_(json::parse(in, nullptr, true, true), p, decoys, rt_model,
                   ms2_model, ccs_model, nce, instrument, irt_rescale,
                   recompute_decoy_mz); }
      catch (const std::exception& e)
      { writeLogError_(std::string("config: ") + e.what()); return ILLEGAL_PARAMETERS; }

      // Relative model paths are relative to the CONFIG, so a config is
      // portable with its models and does not depend on where the tool was
      // invoked from.
      const auto base = fs::path(cfg).parent_path();
      for (std::string* m : {&rt_model, &ms2_model, &ccs_model})
      {
        if (!m->empty() && fs::path(*m).is_relative())
        { *m = (base / *m).lexically_normal().string(); }
      }
    }

    if (const std::string wc = getStringOption_("write_config"); !wc.empty())
    {
      const json eff = effectiveConfig(p, decoys, rt_model, ms2_model, ccs_model,
                                       nce, instrument, irt_rescale,
                                       recompute_decoy_mz);
      std::ofstream os(wc);
      if (!os) { writeLogError_("cannot write config to " + wc); return CANNOT_WRITE_OUTPUT_FILE; }
      os << eff.dump(2) << '\n';
      writeLogInfo_("wrote effective config to " + wc);
      return EXECUTION_OK;
    }

    const std::string fasta = getStringOption_("in"), out = getStringOption_("out");
    if (fasta.empty() || out.empty())
    { writeLogError_("-in and -out are required"); return ILLEGAL_PARAMETERS; }
    // setValidFormats_ only warns on an unknown output extension, and the
    // writer below treated "anything that is not .parquet" as TSV -- so
    // `-out library.parqet` silently produced a DIA-NN TSV under a Parquet
    // name, which no consumer would open and every consumer would misreport.
    const bool parquet = out.ends_with(".parquet");
    if (!parquet && !out.ends_with(".tsv"))
    {
      writeLogError_("-out must end in .parquet or .tsv (got '" + out + "')");
      return ILLEGAL_PARAMETERS;
    }

    // Fill in unset models from the shipped/installed locations before the ONNX
    // session constructor gets a chance to fail with "Load model from  failed".
    struct { const char* name; const char* file; std::string* path; } models[] = {
      {"rt_model", kRtModelFile, &rt_model},
      {"ms2_model", kMs2ModelFile, &ms2_model},
      {"ccs_model", kCcsModelFile, &ccs_model}};
    for (auto& m : models)
    {
      if (m.path->empty()) { *m.path = findModel(m.file); }
    }
    // The recipe is built HERE, after resolution, and not before: a model
    // supplied through DIALIBGEN_MODEL_DIR would otherwise be embedded as
    // "rt_model": "", so a library produced that way would not state what
    // produced it -- which is the one thing the embedded recipe is for. The
    // cache fingerprint was unaffected (it hashes model CONTENT), so this was
    // invisible to every check except reading the recipe back.
    const json eff = effectiveConfig(p, decoys, rt_model, ms2_model, ccs_model,
                                     nce, instrument, irt_rescale,
                                     recompute_decoy_mz);
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
    try
    {
      const auto st = ODIA::LibraryGenerator::generate(fasta, p, library);
      writeLogInfo_("digest: " + std::to_string(st.proteins) + " proteins, " +
                    std::to_string(st.peptides) + " peptides, " +
                    std::to_string(st.precursors) + " precursors, " +
                    std::to_string(st.transitions) + " transitions");

      // Predictions BEFORE decoys: a decoy carries its target's predictions, so
      // appending first would leave them unpredicted.
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
        std::string standards = getStringOption_("irt_standards");
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

      ODIA::LibraryGenerator::predictFragmentIntensities(
        library, ms2_model, p, static_cast<float>(nce), instrument, true, threads);
      ODIA::LibraryGenerator::predictCollisionCrossSections(library, ccs_model, true, threads,
                                                    p.derive_ion_mobility);
      std::size_t skipped = 0;
      // parseDecoyMethod, NOT a second hand-rolled ternary chain. This used to
      // enumerate the methods itself and silently mapped every name it did not
      // know to Mutate, so "reverse" and "shuffle" produced mutation decoys and
      // the three libraries came out byte-identical. The name itself is checked
      // in apply_, which is where an unknown one is now refused.
      const auto method = ODIA::parseDecoyMethod(decoys);
      // Same fragment bar as the targets: applying it to one class only is an
      // anti-conservative FDR.
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
    // Models by CONTENT, not path (see DIANNLibraryFile::hashFile).
    fp.target_params = ODIA::LibraryGenerator::fingerprintParams(
      p, ODIA::DIANNLibraryFile::hashFile(rt_model),
      ODIA::DIANNLibraryFile::hashFile(ms2_model),
      ODIA::DIANNLibraryFile::hashFile(ccs_model), nce, instrument, irt_rescale);
    fp.decoy_method = decoys;
    fp.params = fp.target_params + ";decoy=" + fp.decoy_method;

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
};

// ---------------------------------------------------------------------------
// -write_ctd / -write_cwl / -write_json for a tool that is not part of OpenMS
// ---------------------------------------------------------------------------
// All three route through TOPPBase::writeToolDescription_(), whose first act is
//     StringList type_list = ToolHandler::getTypes(tool_name_);
// getTypes() looks the name up in ToolHandler's hard-coded registry of official
// TOPP tools and THROWS on a miss, so any tool built outside the OpenMS tree
// dies with "Requested tool 'DIALibraryGenerator' does not exist!" and exit 8.
// `official=false` does not help: inside writeToolDescription_ that flag guards
// only the CATEGORY lookup, whose getCategory() is the non-throwing sibling.
// -write_ini is unaffected because its branch never calls getTypes at all,
// which is why the GUI's parameter manifest works while CTD export did not.
//
// The escape hatch stock OpenMS offers is ToolHandler's INTERNAL tool registry:
// it merges tools parsed from *.ttd files found in, among other places,
// whatever OPENMS_TTD_INTERNAL_PATH points at. Registering this tool there
// makes getTypes() return an empty type list and the export proceeds. Every
// byte of the descriptor is then written by OpenMS's own ParamCTDFile over
// OpenMS's own getDefaultParameters_(); nothing here composes XML. That is
// deliberate -- a descriptor we generated ourselves would be a second,
// drifting source of truth for the parameter contract.
namespace
{
  /// Which of TOPPBase's descriptor branches this command line will take,
  /// mirroring TOPPBase::main's own precedence: write_ini wins over everything,
  /// then ctd, then the cwl/json spellings.
  enum class DescriptorRequest { None, Ini, Ctd, CwlOrJson };

  DescriptorRequest requestedDescriptor(int argc, const char** argv)
  {
    // TOPPBase::parseCommandLine_ matches option tokens by exact string ("-" +
    // registered name), so "--write_ctd" and "-write_ctd=DIR" are unknown
    // options to it too; a literal comparison here recognises exactly what it
    // recognises.
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

  /// Registers the tool with OpenMS's ToolHandler for the lifetime of the
  /// object, by writing a one-element .ttd into a private temporary directory
  /// and pointing OPENMS_TTD_INTERNAL_PATH at it.
  ///
  /// Every failure path is a silent no-op that leaves the environment
  /// untouched; the caller then gets OpenMS's stock behaviour, which is what
  /// happens today anyway. The one thing this must never do is leave a
  /// half-written .ttd where ToolHandler can find it: a truncated file raises a
  /// FATAL Xerces error from inside getTOPPToolList(), which runs before
  /// TOPPBase::main's try block. Hence write to a scratch name, size-check,
  /// then rename into place.
  class ToolHandlerRegistration
  {
  public:
    ToolHandlerRegistration()
    {
      // Never override a value the user set: OpenMS reads exactly one directory
      // from this variable, so there is nothing to append our entry to. An
      // empty value counts as unset -- it reaches QDir(""), which means the
      // current working directory, not "disabled".
      const char* existing = std::getenv("OPENMS_TTD_INTERNAL_PATH");
      if (existing != nullptr && *existing != '\0') { return; }

      std::error_code ec;
      const fs::path base = fs::temp_directory_path(ec);
      if (ec) { return; }
      const fs::path dir =
        base / ("dialibgen-ttd-" + std::string(OpenMS::File::getUniqueName(false)));
      fs::create_directories(dir, ec);
      if (ec || !fs::is_directory(dir, ec)) { return; }

      static const char kTtd[] = "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n"
                                 "<ttd><tool status=\"internal\">"
                                 "<name>DIALibraryGenerator</name>"
                                 "<category/><type/></tool></ttd>\n";
      constexpr std::size_t kTtdLen = sizeof(kTtd) - 1;
      const fs::path scratch = dir / "DIALibraryGenerator.ttd.part";
      const fs::path target = dir / "DIALibraryGenerator.ttd";
      {
        std::ofstream os(scratch, std::ios::binary | std::ios::trunc);
        os.write(kTtd, static_cast<std::streamsize>(kTtdLen));
        os.close();
        if (!os) { fs::remove_all(dir, ec); return; }
      }
      fs::rename(scratch, target, ec);
      if (ec || fs::file_size(target, ec) != kTtdLen || ec)
      { fs::remove_all(dir, ec); return; }

      const std::string dir_str = dir.string();
#ifdef _WIN32
      const bool set_ok = (_putenv_s("OPENMS_TTD_INTERNAL_PATH", dir_str.c_str()) == 0);
#else
      const bool set_ok = (::setenv("OPENMS_TTD_INTERNAL_PATH", dir_str.c_str(), 1) == 0);
#endif
      if (!set_ok) { fs::remove_all(dir, ec); return; }
      dir_ = dir;
      active_ = true;
    }

    ~ToolHandlerRegistration()
    {
      if (!active_) { return; }
      std::error_code ec;
      fs::remove_all(dir_, ec);   // best effort, never throws
    }

    ToolHandlerRegistration(const ToolHandlerRegistration&) = delete;
    ToolHandlerRegistration& operator=(const ToolHandlerRegistration&) = delete;

  private:
    fs::path dir_;
    bool active_ = false;
  };
}

int main(int argc, const char** argv)
{
  // OpenMS asks its REST server whether a newer OPENMS exists. For a tool that
  // is not OpenMS the answer is not actionable, and when the request fails --
  // offline, behind a proxy, or sandboxed inside an .app -- Qt prints
  // "QIODevice::read (QNetworkReplyHttpImpl): device not open" to stderr, which
  // reads like an error from this tool and gets reported as one.
  //
  // Set only if the user has NOT: exporting it yourself still wins, in both
  // directions, so the check can be turned back on.
#ifdef _WIN32
  size_t sz = 0;
  if (getenv_s(&sz, nullptr, 0, "OPENMS_DISABLE_UPDATE_CHECK") != 0 || sz == 0)
  { _putenv_s("OPENMS_DISABLE_UPDATE_CHECK", "ON"); }
#else
  ::setenv("OPENMS_DISABLE_UPDATE_CHECK", "ON", 0);
#endif

  // -threads defaults to 0 = all cores (capped at 16 sessions; see
  // inferenceSessions). TOPPBase registers that option itself, with default 1,
  // AFTER registerOptionsAndFlags_ and offers no hook to change it -- so the
  // default is applied here, to argv, before TOPPBase parses it. Left alone
  // when the user passes -threads, and when they pass -ini: the command line
  // outranks the INI, so an injected value would silently override a threads=
  // line in the file. Also left alone for a bare invocation, where TOPPBase
  // prints the usage and an injected pair would turn that into a parameter
  // error instead.
  std::vector<const char*> args(argv, argv + argc);
  bool explicit_threads = false, help = false;
  for (int i = 1; i < argc; ++i)
  {
    if (std::strcmp(argv[i], "-threads") == 0 || std::strcmp(argv[i], "-ini") == 0)
    { explicit_threads = true; }
    if (std::strcmp(argv[i], "--help") == 0 || std::strcmp(argv[i], "--helphelp") == 0)
    { help = true; }
  }
  if (!explicit_threads && argc > 1)
  {
    args.push_back("-threads");
    args.push_back("0");
  }

  // Tool description export (see ToolHandlerRegistration). Gated on the usage
  // banner: registering with ToolHandler flips exactly one --help line,
  // "Common UTIL options:" -> "Common TOPP options:", chosen before any option
  // is parsed. Keeping the shim off whenever --help was asked for keeps that
  // output byte-identical.
  const DescriptorRequest descriptor = requestedDescriptor(argc, argv);
  const bool want_descriptor =
    !help && (descriptor == DescriptorRequest::Ctd || descriptor == DescriptorRequest::CwlOrJson);

#ifndef ENABLE_TDL
  // ENABLE_TDL is a PUBLIC compile definition on OpenMS's exported CMake
  // target, so it reaches this file exactly when the OpenMS being linked can
  // write CWL/JSON at all. Refuse here rather than letting ParamCWLFile::store()
  // open (and therefore TRUNCATE) the target before throwing from a depth
  // TOPPBase does not catch: that would replace a good descriptor with an empty
  // one.
  if (want_descriptor && descriptor == DescriptorRequest::CwlOrJson)
  {
    std::cerr << "DIALibraryGenerator: -write_cwl and -write_json need an OpenMS built with\n"
                 "                     ENABLE_TDL=ON; the OpenMS this binary links has it off.\n"
                 "                     Use -write_ctd instead." << std::endl;
    return OpenMS::TOPPBase::UNKNOWN_ERROR;
  }
#endif

  std::unique_ptr<ToolHandlerRegistration> ttd;
  if (want_descriptor) { ttd = std::make_unique<ToolHandlerRegistration>(); }

  DIALibraryGenerator tool;
  OpenMS::TOPPBase::ExitCodes rc;
  if (ttd)
  {
    // Only descriptor runs get this net, so every other invocation keeps
    // today's exact failure behaviour. It is needed because a malformed or
    // duplicated .ttd raises from getTOPPToolList(), outside TOPPBase::main's
    // own try block, so nothing else would catch it.
    try { rc = tool.main(static_cast<int>(args.size()), args.data()); }
    catch (const std::exception& e)
    {
      std::cerr << "DIALibraryGenerator: tool description export failed: " << e.what() << std::endl;
      rc = OpenMS::TOPPBase::UNKNOWN_ERROR;
    }
  }
  else
  {
    rc = tool.main(static_cast<int>(args.size()), args.data());
  }

  // The same missing hook leaves --help printing TOPPBase's own -threads line,
  // "(default: '1')" -- wrong for this tool and not editable: printUsage_ is
  // not virtual and wraps at the console width, so patching the captured text
  // would be fragile. TOPPBase prints help to stderr and returns EXECUTION_OK,
  // so the correction follows on the same stream, and only when help was really
  // printed. Nothing parses this: the GUI manifest comes from -write_ini, and
  // CI looks for "Version:" only.
  if (help && rc == OpenMS::TOPPBase::EXECUTION_OK)
  {
    std::cerr << "Note: in DIALibraryGenerator, -threads defaults to 0 = all available cores;\n"
                 "      the -threads line above is OpenMS's and shows OpenMS's own default.\n"
                 "      With -ini, the file's threads value applies (1 if it has none) unless\n"
                 "      -threads is also given on the command line.\n" << std::endl;
  }
  return rc;
}
