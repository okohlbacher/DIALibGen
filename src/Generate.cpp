// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

#include "DIALibGen.h"
#include <odia/DIANNLibraryFile.h>
#include <odia/AtomicFile.h>
#include <odia/Library.h>
#include <odia/LibraryExport.h>
#include <odia/LibraryGenerator.h>
#include <odia/PeptDeepPredictor.h>
#include <odia/Hardware.h>
#include <odia/PeptDeepEncoder.h>
#include <odia/Trace.h>
#include <OpenMS/SYSTEM/File.h>
#include "library/Parallel.h"
#include <atomic>
#include <filesystem>
#include <future>
#include <fstream>
#include <condition_variable>
#include <cstdlib>
#include <deque>
#include <exception>
#include <functional>
#include <memory>
#include <mutex>
#include <string_view>
#include <cstdint>
#include <iostream>
#include <system_error>
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

  /// Hands MS2 results to the library writer while the prediction goes on.
  ///
  /// generate predicts RT and CCS before MS2, or on a thread of their own
  /// while a GPU predicts MS2: a precursor is complete once its MS2 spectrum
  /// is in and the RT and CCS phases are over. A consumer thread takes each
  /// block's kept precursors (a queue bounded at about one row group of
  /// precursors, at least two blocks: MS2 runs on while a row group is
  /// encoded, and the copies held stay near 200 MB) and holds them until
  /// the CCS guard is there (setGuard: the RT and CCS phases are done),
  /// settles their CCS batching (CcsBatchGuard, once the old order's CCS
  /// call they belonged to is complete) and passes them to the StreamWriter,
  /// which writes every complete row group. Strictly in precursor order;
  /// the output is the one-pass writer's, byte for byte.
  class Ms2Stream final : public ODIA::LibraryGenerator::Ms2Sink
  {
  public:
    using Open = std::function<std::unique_ptr<ODIA::DIANNLibraryFile::StreamWriter>()>;

    /// Starts the consumer thread at once, so that a process without a
    /// thread to spare finds out before CCS is predicted in the new order
    /// (std::system_error) and can keep 0.11.1's order instead.
    Ms2Stream(ODIA::Library& library, Open open)
      : library_(library), open_(std::move(open)), thread_([this] { run(); }) {}

    /// The RT and CCS values are final: from now on the consumer settles
    /// and writes. Any time, from any thread; nothing is written before.
    void setGuard(ODIA::LibraryGenerator::CcsBatchGuard* guard)
    {
      {
        const std::lock_guard<std::mutex> lock(mutex_);
        guard_ = guard;
      }
      changed_.notify_all();
    }

    /// The RT or CCS phase failed: the next rows() and drain() rethrow it,
    /// and the consumer stops.
    void fail(std::exception_ptr error)
    {
      {
        const std::lock_guard<std::mutex> lock(mutex_);
        if (!error_) { error_ = std::move(error); }
      }
      changed_.notify_all();
    }

    ~Ms2Stream() override
    {
      {
        const std::lock_guard<std::mutex> lock(mutex_);
        done_ = true;
        abandon_ = !drained_;
      }
      changed_.notify_all();
      if (thread_.joinable()) { thread_.join(); }
    }

    void rows(ODIA::Library::PrecursorRows rows) override
    {
      std::unique_lock<std::mutex> lock(mutex_);
      // On one H100 a two-block bound made MS2 wait while each row group
      // was encoded (+7 s for the proteome); a row group's worth does not.
      changed_.wait(lock, [&] { return queue_.size() < 2 || queued_rows_ < QUEUE_ROWS || error_; });
      if (error_) { std::rethrow_exception(error_); }
      queued_rows_ += rows.precursor.size();
      queue_.push_back(std::move(rows));
      changed_.notify_all();
    }

    void drain() override
    {
      {
        const std::lock_guard<std::mutex> lock(mutex_);
        done_ = true;
        drained_ = true;
      }
      changed_.notify_all();
      thread_.join();
      if (error_) { std::rethrow_exception(error_); }
    }

    std::unique_ptr<ODIA::DIANNLibraryFile::StreamWriter> takeWriter() { return std::move(writer_); }
    std::size_t resettled() const { return resettled_; }

  private:
    /// A whole old-order CCS block of kept precursors is not settled yet.
    /// kept_ and settled_ belong to the consumer thread; read here only by
    /// it (inside the wait predicate).
    bool settleable() const
    {
      const std::size_t block = ODIA::LibraryGenerator::ccsBlock();
      return kept_.size() / block * block > settled_;
    }

    void run()
    {
      try
      {
        for (;;)
        {
          ODIA::Library::PrecursorRows rows;
          bool got = false, finished = false;
          ODIA::LibraryGenerator::CcsBatchGuard* guard = nullptr;
          {
            std::unique_lock<std::mutex> lock(mutex_);
            // Rows are taken as they come, guard or no guard (the queue must
            // not fill while RT and CCS are still running); settling and
            // writing wait for the guard.
            changed_.wait(lock, [&] {
              return !queue_.empty() || error_ || abandon_ || (guard_ != nullptr && (done_ || settleable()));
            });
            if (abandon_ || error_) { return; }
            guard = guard_;
            if (!queue_.empty())
            {
              rows = std::move(queue_.front());
              queue_.pop_front();
              queued_rows_ -= rows.precursor.size();
              got = true;
            }
            else if (done_) { finished = true; }   // and the guard is there
          }
          if (got)
          {
            changed_.notify_all();
            kept_.insert(kept_.end(), rows.precursor.begin(), rows.precursor.end());
            held_end_.push_back(kept_.size());
            held_.push_back(std::move(rows));
          }
          if (guard != nullptr)
          {
            const std::size_t block = ODIA::LibraryGenerator::ccsBlock();
            const std::size_t complete = kept_.size() / block * block;
            if (complete > settled_)
            {
              resettled_ += guard->settle(library_, kept_, settled_, complete);
              settled_ = complete;
            }
            release();
          }
          if (finished) { break; }
        }
        // Every kept precursor is known: the last old-order CCS call is complete.
        resettled_ += guard_->settle(library_, kept_, settled_, kept_.size());
        settled_ = kept_.size();
        release();
      }
      catch (...)
      {
        const std::lock_guard<std::mutex> lock(mutex_);
        if (!error_) { error_ = std::current_exception(); }
        changed_.notify_all();
      }
    }

    void release()
    {
      while (!held_.empty() && held_end_.front() <= settled_)
      {
        if (!writer_) { writer_ = open_(); }
        writer_->add(std::move(held_.front()));
        held_.pop_front();
        held_end_.pop_front();
      }
    }

    ODIA::Library& library_;
    ODIA::LibraryGenerator::CcsBatchGuard* guard_ = nullptr;   ///< guarded by mutex_
    Open open_;
    std::unique_ptr<ODIA::DIANNLibraryFile::StreamWriter> writer_;
    std::vector<std::uint32_t> kept_;          ///< kept precursors so far, in order
    std::deque<ODIA::Library::PrecursorRows> held_;
    std::deque<std::size_t> held_end_;         ///< kept_.size() after each held block
    std::size_t settled_ = 0;
    std::size_t resettled_ = 0;
    std::mutex mutex_;
    std::condition_variable changed_;
    std::deque<ODIA::Library::PrecursorRows> queue_;
    static constexpr std::size_t QUEUE_ROWS = std::size_t(1) << 20;
    std::size_t queued_rows_ = 0;
    bool done_ = false, drained_ = false, abandon_ = false;
    std::exception_ptr error_;
    std::thread thread_;   ///< last: starts once everything above exists
  };

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

OpenMS::TOPPBase::ExitCodes DIALibGen::prepare_()
{
  // The MS2 model's CoreML shapes for the length window, compiled into the
  // cache now rather than during the first -machine:device coreml run (a
  // cold cache cost that run 66 s on a whole proteome; 42 shapes take
  // about 105 s here, Core ML compiling one at a time). RT and CCS do not
  // run on CoreML (no LSTM), so
  // there is nothing to prepare for them. Elsewhere there is nothing to do.
  if (!ODIA::PeptDeepPredictor::coremlProviderAvailable())
  {
    writeLogInfo_("nothing to prepare: this ONNX Runtime has no CoreML provider (macOS only)");
    return EXECUTION_OK;
  }
  std::string ms2_model = getStringOption_("generation:ms2_model");
  if (ms2_model.empty()) { ms2_model = findModel(kMs2ModelFile); }
  std::error_code ec;
  if (ms2_model.empty() || !fs::exists(ms2_model, ec))
  { writeLogError_("no MS2 model found; set DIALIBGEN_MODEL_DIR or -generation:ms2_model"); return INPUT_FILE_NOT_FOUND; }
  const auto window = getIntList_("generation:peptide_length");
  if (window.size() != 2 || window[0] < 1 || window[1] < window[0])
  { writeLogError_("-generation:peptide_length must be [min, max]"); return ILLEGAL_PARAMETERS; }
  std::vector<std::string> warnings;
  ODIA::GpuOptions gpu = gpuOptions_(&warnings);
  gpu.coreml = true;
  gpu.prebuild_min_residues = window[0];
  gpu.prebuild_max_residues = window[1];
  try
  {
    ODIA::trace::Scope phase("prepare", 0);
    ODIA::PeptDeepPredictor predictor(ms2_model, true, 1, 1, -1, gpu);
    if (predictor.provider() != ODIA::PeptDeepPredictor::Provider::CoreML)
    {
      writeLogError_("the MS2 model cannot run on CoreML: " + predictor.fallbackReason());
      return UNEXPECTED_RESULT;
    }
    const auto [shapes, seconds] = predictor.waitForShapes();
    std::ostringstream m;
    m.setf(std::ios::fixed); m.precision(1);
    m << "prepared " << shapes << " CoreML shape(s) of " << fs::path(ms2_model).filename().string() << " for "
      << window[0] << "-" << window[1] << " residues in " << seconds << " s (summed over the builders)";
    const auto cache = predictor.coremlCacheDir();
    if (!cache.empty()) { m << "; cache " << cache.string(); } else { m << "; the cache is off (ODIA_COREML_CACHE=0), nothing was kept"; }
    writeLogInfo_(m.str());
  }
  catch (const std::exception& e) { writeLogError_(std::string("prepare failed: ") + e.what()); return UNEXPECTED_RESULT; }
  return EXECUTION_OK;
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
    // Configuration, model lookup and option checks, up to the FASTA read.
    ODIA::trace::Scope setup_phase("setup", 0);
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
    // Neither instrument has trained weights in the shipped MS2 checkpoint
    // (its own constants name four instruments). Measured against Lumos at
    // the same NCE (median spectral cosine per precursor, synthetic
    // peptides): ThermoTOF 0.993 -- effectively Lumos -- but SciexTOF only
    // 0.914, with a tenth of the precursors below 0.64.
    if (instrument == "ThermoTOF")
    { writeLogWarn_("ThermoTOF carries no trained weights in the shipped MS2 checkpoint; its predictions are "
                    "effectively Lumos's (median spectral cosine 0.993 at the same NCE). Name Lumos to say so."); }
    if (instrument == "SciexTOF")
    { writeLogWarn_("SciexTOF (also TripleTOF and ZenoTOF) carries no trained weights in the shipped MS2 "
                    "checkpoint, and its predictions are NOT those of Lumos: median spectral cosine 0.914 against "
                    "Lumos at the same NCE, below 0.64 for a tenth of precursors. Neither has been validated on "
                    "SCIEX data; name Lumos to get the no-correction baseline instead."); }

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
    ODIA::LibraryFormat format;
    try { format = ODIA::LibraryExport::resolve(out, getStringOption_("out_format")); }
    catch (const std::invalid_argument& e) { writeLogError_(e.what()); return ILLEGAL_PARAMETERS; }
    const bool parquet = format == ODIA::LibraryFormat::DIALibGenParquet;

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
    if (threads == 0)
    { writeLogInfo_("-threads 0: " + std::to_string(ODIA::physicalCores()) + " threads (" +
                    ODIA::physicalCoresReason() + ")"); }
    {
      // Interleaving the run's memory over the NUMA nodes (what numactl
      // --interleave=all does) was the default above 32 threads in 0.11.2:
      // the threads started per call crowded onto one node. Pool workers
      // alternate between the nodes and touch their own memory first, and
      // interleaving now costs time: 131.2 s against 123.9 s without, at 96
      // threads on 2 x 96 cores. ODIA_NUMA_INTERLEAVE=1 still turns it on.
      const unsigned workers = threads == 0 ? ODIA::physicalCores() : threads;
      const char* forced = std::getenv("ODIA_NUMA_INTERLEAVE");
      const bool want = forced != nullptr && std::string_view(forced) == "1";
      const unsigned nodes = ODIA::numaNodes();
      const bool interleaved = want && nodes > 1 && ODIA::interleaveMemory();
      if (interleaved)
      { writeLogInfo_("memory interleaved over " + std::to_string(nodes) + " NUMA nodes (" +
                      std::to_string(workers) + " threads; ODIA_NUMA_INTERLEAVE=1)"); }
      ODIA::trace::note("numa_interleave", interleaved ? "on" : "off");
    }
    writeLogInfo_("MS2 model: " + ms2_model);
    std::vector<std::string> device_warnings;
    ODIA::GpuOptions gpu = gpuOptions_(&device_warnings);
    // CoreML compiles one MS2 shape per residue count: the window's shapes
    // are built in the background from the predictor's construction on.
    gpu.prebuild_min_residues = static_cast<int>(p.min_length);
    gpu.prebuild_max_residues = static_cast<int>(p.max_length);
    ODIA::Library library;
    std::string standards;
    setup_phase.fields().rows = threads;
    setup_phase.end();

    // The library's identity, known before generation so that the writer can
    // start while MS2 runs.
    ODIA::trace::Scope fingerprint_phase("fingerprint", 0);
    ODIA::DIANNLibraryFile::Fingerprint fp;
    try { fp = ODIA::DIANNLibraryFile::fingerprintFasta(fasta); }
    catch (const std::exception& e)
    { writeLogError_(std::string("cannot fingerprint the FASTA: ") + e.what()); return INPUT_FILE_NOT_FOUND; }
    fp.target_params = ODIA::LibraryGenerator::fingerprintParams(
      p, ODIA::DIANNLibraryFile::hashFile(rt_model),
      ODIA::DIANNLibraryFile::hashFile(ms2_model),
      ODIA::DIANNLibraryFile::hashFile(ccs_model), nce, instrument, irt_rescale);
    fp.decoy_method = decoys;
    fingerprint_phase.end();
    const auto finishFingerprint = [&]() {
      if (irt_rescale && fp.target_params.find(";irt_standards=") == std::string::npos)
      { fp.target_params += ";irt_standards=" + ODIA::DIANNLibraryFile::hashFile(standards); }
      fp.params = fp.target_params + ";decoy=" + fp.decoy_method +
                  ";recompute_decoy_mz=" + (recompute_decoy_mz ? "1" : "0");
    };
    // Recorded only when there is something to say, so a library built where
    // nothing fell back is byte-identical to one from a build without this.
    const auto extraMetadata = [&]() {
      std::vector<std::pair<std::string, std::string>> extra;
      if (!device_warnings.empty()) { extra.emplace_back("odia.warnings", json(device_warnings).dump()); }
      return extra;
    };
    // RT and CCS first and MS2 last, with the writer streaming behind MS2.
    // ODIA_STREAM=0 restores 0.11.1's order (MS2 before CCS) and the one-pass
    // writer: the same file, for comparison.
    const char* stream_env = std::getenv("ODIA_STREAM");
    // The other layouts are written in one pass, from the final library.
    const bool stream = (stream_env == nullptr || std::string_view(stream_env) != "0") &&
                        (parquet || format == ODIA::LibraryFormat::DIANNTSV);
    // When MS2 runs on a GPU, RT and CCS run on the CPU meanwhile, on a
    // thread of their own: MS2 does not need them, the writer waits for
    // both, and the CPU's RT and CCS calls are the same calls as before, so
    // the library is the same. On a whole human library on an M4 Max they
    // took 77 s of 16 idle cores before 236 s of GPU time. ODIA_OVERLAP=0
    // keeps the phases in order, for comparison.
    const char* overlap_env = std::getenv("ODIA_OVERLAP");
    const bool overlap_wanted = overlap_env == nullptr || std::string_view(overlap_env) != "0";
    std::unique_ptr<ODIA::DIANNLibraryFile::StreamWriter> writer;
    // What a failing phase on the CPU thread reports: the message and the
    // exit code the phase returned when it ran on the main thread.
    struct PhaseFailure : std::runtime_error
    {
      PhaseFailure(ExitCodes c, const std::string& what) : std::runtime_error(what), code(c) {}
      ExitCodes code;
    };
    try
    {
      const auto st = ODIA::LibraryGenerator::generate(fasta, p, library, threads);
      writeLogInfo_("digest: " + std::to_string(st.proteins) + " proteins, " +
                    std::to_string(st.peptides) + " peptides, " +
                    std::to_string(st.precursors) + " precursors, " +
                    std::to_string(st.transitions) + " transitions");
      if (st.dropped_proteins || st.dropped_ambiguous_peptides)
      { writeLogWarn_("digest skipped " + std::to_string(st.dropped_proteins) +
          " invalid proteins and " + std::to_string(st.dropped_ambiguous_peptides) +
          " peptides with ambiguous or undefined residue masses"); }
      for (const auto& [name, count] : st.variable_modification_precursors)
      { writeLogInfo_("variable modification '" + name + "': " + std::to_string(count) + " precursors"); }

      // Sessions of the CPU phases: all of -threads, or fewer while they
      // overlap the GPU's MS2 phase (below).
      unsigned cpu_threads = threads;
      // The RT phase and the iRT calibration; @p info takes the log lines
      // (written at once on the main thread, collected on the CPU thread).
      const auto rtPhase = [&](const std::function<void(const std::string&)>& info) {
        if (const auto miss = ODIA::LibraryGenerator::predictRetentionTimes(
              library, rt_model, true, cpu_threads, p.free_cysteine_rt_correction, gpu); miss)
        {
          throw PhaseFailure(UNEXPECTED_RESULT, std::to_string(miss) + " precursors have no predicted RT. "
                             "Refusing rather than shipping a library with NaN retention "
                             "times, which extracts from nowhere and looks like a search bug.");
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
            throw PhaseFailure(INPUT_FILE_NOT_FOUND, e.str());
          }
          const auto cal = ODIA::LibraryGenerator::fitIrtCalibration(rt_model, standards, true, gpu);
          ODIA::LibraryGenerator::applyIrtCalibration(library, cal);
          std::ostringstream m;
          m.setf(std::ios::fixed); m.precision(4);
          m << "rescaled to iRT: " << cal.slope << " * raw + " << cal.intercept
            << " from " << cal.peptides << " standards (worst off by "
            << cal.max_abs_error << " iRT)";
          info(m.str());
        }
        else
        {
          info("RT column is the model's RAW 0..1 output, NOT iRT "
               "(irt_rescale=false). The per-run map is scale-invariant, so "
               "this is the pipeline domain; it is NOT interchangeable with an "
               "iRT library in another tool.");
        }
        finishFingerprint();
      };
      const auto logInfo = [&](const std::string& line) { writeLogInfo_(line); };

      // The MS2 predictor first: its provider decides whether RT and CCS
      // overlap it, and on CoreML its shape builds start now.
      std::unique_ptr<ODIA::PeptDeepPredictor> ms2;
      if (stream) { ms2 = ODIA::LibraryGenerator::openMs2Predictor(ms2_model, true, threads, gpu); }
      const bool overlap = overlap_wanted && ms2 && ms2->provider() != ODIA::PeptDeepPredictor::Provider::CPU;
      // While the GPU predicts MS2, its host side (Core ML's own threads,
      // the encoding workers) needs cores: with all 16 CPU sessions busy on
      // CCS, the Astral MS2 phase on an M4 Max took 229 s instead of 153 s
      // and the overlap gained 6 s of 238 (with 12 sessions 208 s, with 8
      // 189 s under a noisier load). The CPU phases therefore leave half
      // the cores free on CoreML (ODIA_OVERLAP_RESERVE overrides the count);
      // their chunking does not depend on the session count, so the library
      // does not either. CUDA's host side is lighter: nothing reserved.
      if (overlap)
      {
        const unsigned all = threads == 0 ? ODIA::physicalCores() : threads;
        unsigned reserve = ms2->provider() == ODIA::PeptDeepPredictor::Provider::CoreML ? all / 2 : 0;
        if (const char* v = std::getenv("ODIA_OVERLAP_RESERVE"); v != nullptr && *v != '\0')
        {
          char* end = nullptr;
          const long long x = std::strtoll(v, &end, 10);
          if (end != v && x >= 0) { reserve = static_cast<unsigned>(x); }
        }
        cpu_threads = std::max(1u, all > reserve ? all - reserve : 1u);
      }
      if (!overlap) { rtPhase(logInfo); }

      // Guard before sink: the sink's thread may use it until the sink is gone.
      std::unique_ptr<ODIA::LibraryGenerator::CcsBatchGuard> guard;
      std::unique_ptr<Ms2Stream> sink;
      const std::string config_json = parquet ? eff.dump() : std::string();
      if (stream)
      {
        try
        {
          sink = std::make_unique<Ms2Stream>(library, [&]() {
            return std::make_unique<ODIA::DIANNLibraryFile::StreamWriter>(out, library, fp, config_json,
                                                                           extraMetadata(), threads);
          });
        }
        catch (const std::system_error& e)
        {
          writeLogWarn_(std::string("no thread for the streaming writer (") + e.what() +
                        "); predicting in 0.11.1's order and writing after MS2");
        }
      }
      if (sink)
      {
        const auto ccsPhase = [&]() {
          std::vector<std::uint16_t> lengths;
          bool ccs_on_cpu = false;
          if (const auto miss = ODIA::LibraryGenerator::predictCollisionCrossSections(
                library, ccs_model, true, cpu_threads, p.derive_ion_mobility, gpu, &lengths, &ccs_on_cpu); miss)
          { throw std::runtime_error(std::to_string(miss) + " precursors have no predicted CCS"); }
          guard = std::make_unique<ODIA::LibraryGenerator::CcsBatchGuard>(ccs_model, p.derive_ion_mobility,
                                                                           std::move(lengths), ccs_on_cpu);
        };
        // The CPU phases on a thread of their own (after the sink: waited for
        // before the sink goes; their log lines follow once they are done).
        std::vector<std::string> cpu_log;
        std::mutex cpu_log_mutex;
        std::atomic<bool> cancel{false};
        std::future<void> cpu_phases;
        struct WaitFor
        {
          std::future<void>& f;
          ~WaitFor() { if (f.valid()) { f.wait(); } }
        } wait_for_cpu_phases{cpu_phases};
        if (overlap)
        {
          writeLogInfo_("RT and CCS predicted on a thread of their own, " + std::to_string(cpu_threads) +
                        " sessions, while the GPU predicts MS2 (ODIA_OVERLAP=0: one after the other)");
          cpu_phases = ODIA::parallel::async([&]() {
            try
            {
              rtPhase([&](const std::string& line) {
                const std::lock_guard<std::mutex> lock(cpu_log_mutex);
                cpu_log.push_back(line);
              });
              if (cancel.load()) { return; }
              ccsPhase();
              sink->setGuard(guard.get());
            }
            catch (...)
            {
              sink->fail(std::current_exception());
              throw;
            }
          });
        }
        else
        {
          ccsPhase();
          sink->setGuard(guard.get());
        }
        try
        {
          if (const auto miss = ODIA::LibraryGenerator::predictFragmentIntensities(
                library, *ms2, p, static_cast<float>(nce), instrument, threads, sink.get()); miss)
          { throw std::runtime_error(std::to_string(miss) + " precursors have no predicted MS2 spectrum"); }
        }
        catch (...)
        {
          cancel = true;
          throw;
        }
        if (cpu_phases.valid()) { cpu_phases.get(); }
        for (const auto& line : cpu_log) { writeLogInfo_(line); }
        ms2.reset();   // its summary line, and on CoreML the cache's pruning, now
        writer = sink->takeWriter();
        if (sink->resettled())
        { writeLogInfo_("CCS of " + std::to_string(sink->resettled()) + " precursors predicted again in their 0.11.1 batching"); }
        writeLogInfo_("after fragment selection: " + std::to_string(library.precursorCount()) + " precursors");
        if (library.precursorCount() == 0)
        { throw std::runtime_error("no precursors remain after digestion and fragment selection"); }
      }
      else
      {
        if (!ms2) { ms2 = ODIA::LibraryGenerator::openMs2Predictor(ms2_model, true, threads, gpu); }
        if (const auto miss = ODIA::LibraryGenerator::predictFragmentIntensities(
              library, *ms2, p, static_cast<float>(nce), instrument, threads); miss)
        { throw std::runtime_error(std::to_string(miss) + " precursors have no predicted MS2 spectrum"); }
        ms2.reset();
        writeLogInfo_("after fragment selection: " + std::to_string(library.precursorCount()) + " precursors");
        if (library.precursorCount() == 0)
        { throw std::runtime_error("no precursors remain after digestion and fragment selection"); }
        if (const auto miss = ODIA::LibraryGenerator::predictCollisionCrossSections(
              library, ccs_model, true, threads, p.derive_ion_mobility, gpu); miss)
        { throw std::runtime_error(std::to_string(miss) + " precursors have no predicted CCS"); }
      }
      std::size_t skipped = 0;
      const auto method = ODIA::parseDecoyMethod(decoys);
      const auto made = ODIA::LibraryGenerator::appendDecoys(library, method, &skipped,
                                                             p.min_fragments,
                                                             recompute_decoy_mz);
      writeLogInfo_("decoys: " + std::to_string(made) + " (" + decoys + "), " +
                    std::to_string(skipped) + " skipped");
    }
    catch (const PhaseFailure& e)
    {
      writeLogError_(e.what());
      return e.code;
    }
    catch (const std::invalid_argument& e)
    {
      // A request the library cannot honour (a modification that modifies
      // nothing, an unsupported specificity, an impossible charge range).
      writeLogError_(std::string("library generation refused: ") + e.what());
      return ILLEGAL_PARAMETERS;
    }
    catch (const std::exception& e)
    { writeLogError_(std::string("library generation failed: ") + e.what()); return UNEXPECTED_RESULT; }

    if (!device_warnings.empty())
    {
      writeLogWarn_(std::to_string(device_warnings.size()) + " model(s) requested on the GPU ran on the CPU; "
                    "recorded as odia.warnings" + std::string(parquet ? "" : " (TSV: log only)"));
    }

    ODIA::trace::Scope write_phase("write", 0);
    write_phase.fields().rows = static_cast<std::int64_t>(library.transitionCount());
    try
    {
      if (writer)
      {
        write_phase.fields().idx = static_cast<std::int64_t>(writer->written());
        writer->finish(library);
      }
      else if (parquet)
      { ODIA::DIANNLibraryFile::storeParquetCompact(out, library, fp, eff.dump(), extraMetadata(), threads); }
      else if (format == ODIA::LibraryFormat::DIANNTSV) { ODIA::DIANNLibraryFile::storeTSV(out, library, threads); }
      else
      {
        // DIA-NN ignores the file's metadata, so the DIA-NN 2.0 layout
        // carries the recipe as DIALibGen's own does.
        auto metadata = extraMetadata();
        if (ODIA::LibraryExport::carriesRecipe(format)) { metadata.insert(metadata.begin(), {"odia.config_json", eff.dump()}); }
        ODIA::LibraryExport::store(out, library, format, metadata, threads);
      }
      if (format == ODIA::LibraryFormat::DIANNTSV)
      {
        writeLogWarn_("TSV cannot carry the recipe: the config is NOT embedded and this "
                      "file is not cache-eligible. Use .parquet to keep them together.");
      }
      else if (!ODIA::LibraryExport::carriesRecipe(format))
      {
        writeLogWarn_(std::string(ODIA::LibraryExport::name(format)) + " cannot carry the recipe: the config "
                      "is NOT embedded. Use .parquet to keep them together.");
      }
    }
    catch (const std::exception& e)
    { writeLogError_(std::string("cannot write ") + out + ": " + e.what()); return CANNOT_WRITE_OUTPUT_FILE; }
    write_phase.end();
    writeLogInfo_("wrote " + out + (format == ODIA::LibraryFormat::DIALibGenParquet || format == ODIA::LibraryFormat::DIANNTSV
                                     ? std::string() : std::string(" (") + ODIA::LibraryExport::name(format) + ")"));
    {
      // Releasing the library: what the destructor would do on return, done
      // here so that the trace can time it.
      ODIA::trace::Scope teardown_phase("teardown", 0);
      library = ODIA::Library();
    }
    return EXECUTION_OK;
  }
