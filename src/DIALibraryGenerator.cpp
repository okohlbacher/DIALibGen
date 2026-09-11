// Copyright (c) 2026, Oliver Kohlbacher and the DIALibraryGenerator authors.
// SPDX-License-Identifier: BSD-3-Clause

/// DIALibraryGenerator: FASTA + JSON config -> in-silico DIA library.
///
/// A library is a CROSS-RUN artefact built once and reused; a search is per-run.
/// They were welded together behind OpenDIAlyzer's -out_lib, which is why the
/// library's own parameters were scattered across ~15 CLI flags and why the
/// alkylation state -- the defect that cost a quarter of the search space
/// (doc/27) -- was hard-coded and unreachable.
///
/// Every algorithm here already exists in odia_core. This file is a config
/// parser and a call sequence; if it grows a second algorithm it is in the
/// wrong file.

#include <odia/DIANNLibraryFile.h>
#include <odia/Library.h>
#include <odia/LibraryGenerator.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <nlohmann/json.hpp>

#include <filesystem>
#include <fstream>
#include <sstream>

using json = nlohmann::json;

namespace
{
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
      {"schema_version", 1},
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
      // flag, and the embedded config is never read back to decide reuse. Two
      // OpenDIAlyzer paths also call appendDecoys without it. Until those are
      // closed, a cached library built with this ON can be silently rebuilt with
      // inherited-m/z decoys.
      {"recompute_decoy_mz", recompute_decoy_mz}};
  }
}

class DIALibraryGenerator final : public OpenMS::TOPPBase
{
public:
  DIALibraryGenerator()
    : TOPPBase("DIALibraryGenerator",
               "Build an in-silico DIA spectral library from a FASTA.",
               false) {}   // not an official TOPP tool; lives in this repo

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Protein FASTA.");
    setValidFormats_("in", {"fasta"});
    registerInputFile_("config", "<file>", "",
                       "JSON library configuration. Every content-affecting parameter lives "
                       "here, so a library cannot be built with settings nobody wrote down. "
                       "Omit it for the defaults; -write_config dumps the effective config "
                       "including everything defaulted.", false);
    setValidFormats_("config", {"json"}, false);
    registerOutputFile_("out", "<file>", "",
                        "Library. .parquet is ODIA-compatible and carries the recipe in its "
                        "schema metadata (odia.config_json); .tsv is the DIA-NN dialect and "
                        "CANNOT carry it -- a TSV is not reproducible from itself.");
    setValidFormats_("out", {"parquet", "tsv"}, false);
    registerOutputFile_("write_config", "<file>", "",
                        "Write the effective config here and exit without building. The "
                        "quickest way to see what the defaults actually are.", false);
    setValidFormats_("write_config", {"json"}, false);
    registerInputFile_("irt_standards", "<file>", "",
                       "Standards defining the iRT scale. The RT model emits a raw 0..1 "
                       "value; without this the library ships THAT, and every consumer "
                       "then calibrates against a scale nobody declared. Defaults to "
                       "data/irt_standards.tsv.", false);
    setValidFormats_("irt_standards", {"tsv"}, false);
    registerIntOption_("threads", "<n>", 0, "Threads. 0 = all.", false);
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
    // RAW MODEL UNITS are the pipeline domain (doc/28): the library carries the
    // RT model's own 0..1 output and the per-run map takes it to seconds. The
    // map is scale-invariant -- the constant cancels -- so the iRT rescale buys
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
    }
    // Relative model paths are relative to the CONFIG, so a config is portable
    // with its models and does not depend on where the tool was invoked from.
    if (!cfg.empty())
    {
      const auto base = std::filesystem::path(cfg).parent_path();
      for (std::string* m : {&rt_model, &ms2_model, &ccs_model})
      {
        if (!m->empty() && std::filesystem::path(*m).is_relative())
        { *m = (base / *m).lexically_normal().string(); }
      }
    }
    const json eff = effectiveConfig(p, decoys, rt_model, ms2_model, ccs_model,
                                     nce, instrument, irt_rescale,
                                     recompute_decoy_mz);

    if (const std::string wc = getStringOption_("write_config"); !wc.empty())
    {
      std::ofstream(wc) << eff.dump(2) << '\n';
      writeLogInfo_("wrote effective config to " + wc);
      return EXECUTION_OK;
    }

    const std::string fasta = getStringOption_("in"), out = getStringOption_("out");
    if (fasta.empty() || out.empty())
    { writeLogError_("-in and -out are required"); return ILLEGAL_PARAMETERS; }

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
        if (standards.empty())
        {
          for (const auto& c : {std::filesystem::path("data/irt_standards.tsv"),
                                std::filesystem::path(ODIA_DATA_DIR) / "irt_standards.tsv"})
          { if (std::filesystem::exists(c)) { standards = c.string(); break; } }
        }
        if (standards.empty())
        { writeLogError_("no iRT standards found; pass -irt_standards"); return INPUT_FILE_NOT_FOUND; }
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
      // the three libraries came out byte-identical.
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
      ODIA::DIANNLibraryFile::hashFile(ccs_model), nce, instrument);
    fp.decoy_method = decoys;
    fp.params = fp.target_params + ";decoy=" + fp.decoy_method;

    if (out.ends_with(".parquet"))
    { ODIA::DIANNLibraryFile::storeParquetCompact(out, library, fp, eff.dump()); }
    else
    {
      ODIA::DIANNLibraryFile::storeTSV(out, library);
      writeLogWarn_("TSV cannot carry the recipe: the config is NOT embedded and this "
                    "file is not cache-eligible. Use .parquet to keep them together.");
    }
    writeLogInfo_("wrote " + out);
    return EXECUTION_OK;
  }
};

int main(int argc, const char** argv)
{
  DIALibraryGenerator tool;
  return tool.main(argc, argv);
}
