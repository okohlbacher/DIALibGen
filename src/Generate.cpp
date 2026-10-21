// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

#include "DIALibGen.h"
#include "InstrumentDefaults.h"
#include <optional>
#include <odia/DIANNLibraryFile.h>
#include <odia/AtomicFile.h>
#include <odia/DeviceRoute.h>
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
#include <map>
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

  /// Where the models are looked for: @p first (-models), $DIALIBGEN_MODEL_DIR,
  /// the bundled directory, OpenMS's data path.
  std::vector<fs::path> modelDirCandidates(const std::string& first = {})
  {
    std::vector<fs::path> out;
    if (!first.empty()) { out.emplace_back(first); }
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

  std::string findModel(const char* name, const std::string& first = {})
  {
    std::error_code ec;
    for (const auto& dir : modelDirCandidates(first))
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
                       bool recompute_decoy_mz, const std::string& ion_mobility,
                       const std::string& irt_standards,
                       const std::string& mobility_correction,
                       const std::map<int, double>& nce_by_charge)
  {
    json by_charge = nullptr;
    if (!nce_by_charge.empty())
    {
      by_charge = json::object();
      for (const auto& [z, v] : nce_by_charge) { by_charge[std::to_string(z)] = v; }
    }
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
      {"predict_ion_mobility", ion_mobility},
      {"mobility_correction", mobility_correction},
      {"min_relative_intensity", p.min_relative_intensity},
      {"reserved_doubly_charged", p.reserved_doubly_charged},
      {"decoys", decoys},
      {"rt_model", rt_model}, {"ms2_model", ms2_model}, {"ccs_model", ccs_model},
      {"instrument", instrument}, {"nce", nce}, {"nce_by_charge", by_charge},
      {"irt_rescale", irt_rescale},
      {"irt_standards", irt_standards},
      {"recompute_decoy_mz", recompute_decoy_mz}};
  }
}


namespace
{
  json generationDefaults()
  {
    json j = effectiveConfig(ODIA::DigestParams{}, "none", "", "", "", 30, "QE", false, false, "auto", "", "auto", {});
    j.erase("schema_version");
    j.erase("nce_by_charge");   // the instrument's, or a -config JSON's: not a TOPP option
    j.erase("irt_standards");   // the top-level -irt_standards
    j["nce"] = -1.0; // TOPP sentinel: use the instrument default.
    // Empty: the instrument's windows, else DigestParams' (resolved in generate_).
    j["precursor_mz"] = json::array();
    j["fragment_mz"] = json::array();
    return j;
  }
  void apply_(const json& j, ODIA::DigestParams& p, std::string& decoys,
              std::string& rt_model, std::string& ms2_model, std::string& ccs_model,
              double& nce, std::string& instrument, bool& irt_rescale,
              bool& recompute_decoy_mz, bool& nce_was_set, std::string& instrument_alias_of,
              std::string& ion_mobility, std::string& irt_standards,
              std::string& mobility_correction, std::optional<std::map<int, double>>& nce_by_charge)
  {
    if (!j.is_object()) { throw std::runtime_error("config must be a JSON object"); }
    const json ref = effectiveConfig(p, decoys, rt_model, ms2_model, ccs_model,
                                     nce, instrument, irt_rescale,
                                     recompute_decoy_mz, ion_mobility, irt_standards, mobility_correction, {});
    for (const auto& [k, v] : j.items())
    {
      // Recorded for the reader, not settings.
      if (k == "nce_source" || k == "instrument_named" || k == "model_instrument" || k == "precursor_mz_source" || k == "fragment_mz_source" ||
          k == "precursor_charges_source" || k == "instrument_defaults" ||
          k == "mobility_correction_applied" || k == "nce_by_charge_source") { continue; }
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
    if (j.contains("predict_ion_mobility"))
    {
      const json& v = j["predict_ion_mobility"];
      if (v.is_boolean()) { ion_mobility = v.get<bool>() ? "true" : "false"; }
      else if (v.is_string() && (v == "auto" || v == "true" || v == "false")) { ion_mobility = v.get<std::string>(); }
      else { throw std::runtime_error("predict_ion_mobility must be auto, true or false"); }
    }
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
    if (j.contains("nce_by_charge"))
    {
      // null: one NCE for every charge; an object: "charge": NCE over nce.
      const json& v = j["nce_by_charge"];
      std::map<int, double> by;
      if (!v.is_null())
      {
        if (!v.is_object()) { throw std::runtime_error("nce_by_charge must be null or an object of \"charge\": NCE"); }
        for (const auto& [key, value] : v.items())
        {
          std::size_t used = 0;
          int z = 0;
          try { z = std::stoi(key, &used); } catch (const std::exception&) { used = 0; }
          if (used != key.size() || z < 1 || z > 8) { throw std::runtime_error("nce_by_charge: '" + key + "' is not a charge from 1 to 8"); }
          if (!value.is_number() || value.get<double>() <= 0 || value.get<double>() > 100)
          { throw std::runtime_error("nce_by_charge: the NCE of charge " + key + " must be greater than 0 and at most 100"); }
          by[z] = value.get<double>();
        }
      }
      nce_by_charge = by;
    }
    if (j.contains("irt_rescale")) { irt_rescale = j["irt_rescale"]; }
    if (j.contains("mobility_correction"))
    {
      const json& v = j["mobility_correction"];
      if (!v.is_string() || (v != "auto" && v != "off"))
      { throw std::runtime_error("mobility_correction must be auto or off"); }
      mobility_correction = v.get<std::string>();
    }
    if (j.contains("irt_standards"))
    {
      if (!j["irt_standards"].is_string()) { throw std::runtime_error("irt_standards must be a path (\"\" = the bundled table)"); }
      irt_standards = j["irt_standards"];
    }
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

    // The name as written: generate_ resolves it through the instrument
    // files and the model's own names (instrument_alias_of is set there).
    (void)instrument_alias_of;
    if (instrument.empty()) { throw std::runtime_error("instrument must not be empty"); }
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
  if (ms2_model.empty()) { ms2_model = findModel(kMs2ModelFile, getParam_().getValue("models").toString()); }
  std::error_code ec;
  if (ms2_model.empty() || !fs::exists(ms2_model, ec))
  { writeLogError_("no MS2 model found; give -models or -generation:ms2_model, or set DIALIBGEN_MODEL_DIR"); return INPUT_FILE_NOT_FOUND; }
  const auto window = getIntList_("generation:peptide_length");
  if (window.size() != 2 || window[0] < 1 || window[1] < window[0])
  { writeLogError_("-generation:peptide_length must be [min, max]"); return ILLEGAL_PARAMETERS; }
  std::vector<std::string> warnings;
  ODIA::GpuOptions gpu = gpuOptions_(&warnings);
  gpu.coreml = true;
  gpu.prebuild_min_residues = window[0];
  gpu.prebuild_max_residues = window[1];
  gpu.prebuild_discard = true;   // the cache is the product, not the sessions
  try
  {
    ODIA::trace::Scope phase("prepare", 0);
    ODIA::PeptDeepPredictor predictor(ms2_model, true, 1, 1, -1, gpu);
    if (predictor.provider() != ODIA::PeptDeepPredictor::Provider::CoreML)
    {
      writeLogError_("nothing prepared: " + predictor.fallbackReason());
      return UNEXPECTED_RESULT;
    }
    const auto [shapes, seconds] = predictor.waitForShapes();
    std::ostringstream m;
    m.setf(std::ios::fixed); m.precision(1);
    m << "prepared " << shapes << " CoreML shape(s) for peptide lengths " << window[0] << "-" << window[1];
    if (ODIA::diagnosticsLevel() > 0)
    { m << " (" << fs::path(ms2_model).filename().string() << ", " << seconds << " s summed over the builders)"; }
    const auto cache = predictor.coremlCacheDir();
    if (!cache.empty()) { m << "; cache: " << cache.string(); } else { m << "; nothing kept: the cache is off (ODIA_COREML_CACHE=0)"; }
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
  registerTOPPSubsection_("generation", "Generate (FASTA to library)");
  // Registration order is the order --help prints: the settings a method
  // needs first, then the rest (advanced: --helphelp).
  struct Text { const char* key; const char* metavar; const char* text; };
  static const Text texts[] = {
    {"instrument", "<name>", "Instrument of the run: Astral, Exploris, Lumos, QE, timsTOF, or a vendor name such as "
                             "'Orbitrap Astral' or 'timsTOF Pro'. Its defaults (instrument files, see --helphelp "
                             "-instrument_file) set the NCE, ion mobility (timsTOF) and m/z windows and charges "
                             "(Astral, timsTOF, Exploris 480, Q Exactive HF/HF-X, Fusion Lumos/Eclipse, "
                             "TripleTOF/ZenoTOF). Default: default.yaml's (QE)"},
    {"nce", "<number>", "Collision energy; -1 = the instrument's (its file's: Astral 25, Q Exactive HF-X and "
                        "TripleTOF/ZenoTOF 27, Fusion Lumos/Eclipse 30; else Lumos 25, others 30)"},
    {"precursor_mz", "<min max>", "Precursor m/z window; match it to the acquisition method. Default: the "
                                  "instrument's (Astral 380-980), else 350-1200"},
    {"fragment_mz", "<min max>", "Fragment m/z window. Default: the instrument's (Astral 150-2000), else 200-1800"},
    {"precursor_charges", "<z>", "Precursor charges. Default: the instrument's (timsTOF 2-4), else 1-4"},
    {"fixed_modifications", "<names>", "OpenMS names, e.g. 'Carbamidomethyl (C)'; give the option without a value for none"},
    {"variable_modifications", "<names>", "OpenMS names, e.g. 'Oxidation (M)'; at most "
                                          "-generation:max_variable_modifications (1) per peptide"},
    {"enzyme", "<name>", "Protease; Trypsin/P also cuts before P"},
    {"missed_cleavages", "<n>", "Missed cleavages per peptide"},
    {"peptide_length", "<min max>", "Peptide length window, in residues"},
    {"max_variable_modifications", "<n>", "Most variable modifications on one peptide"},
    {"n_terminal_methionine_excision", "<true/false>", "Also digest each protein without its initiator methionine"},
    {"fragments", "<min max>", "Fragments per precursor: fewer than min drops the precursor, the most intense up to "
                               "max are kept"},
    {"max_fragment_charge", "<z>", "Highest fragment charge (1 or 2)"},
    {"min_relative_intensity", "<fraction>", "Drop fragments below this fraction of the precursor's most intense one"},
    {"reserved_doubly_charged", "<n>", "Fragment slots reserved for doubly charged ions; 0 ranks all fragments by "
                                       "predicted intensity"},
    {"predict_ion_mobility", "<auto/true/false>", "Predict CCS (and 1/K0, see derive_ion_mobility). auto: only for an "
                                                  "instrument with ion mobility (timsTOF); true: always; false: never, "
                                                  "the library then carries neither"},
    {"derive_ion_mobility", "<true/false>", "Also write 1/K0, derived from the predicted CCS (when CCS is predicted)"},
    {"mobility_correction", "<auto/off>", "auto: correct the derived 1/K0 with the instrument file's mobility_correction "
                                          "(timsTOF; only for the CCS model it was fitted for); off: the plain Mason-Schamp "
                                          "1/K0 of the predicted CCS"},
    {"irt_rescale", "<true/false>", "Rescale RT to iRT with the standards of -irt_standards; false: RT is the "
                                    "model's raw 0..1 scale (DIA-NN fits its own RT map)"},
    {"free_cysteine_rt_correction", "<true/false>", "Add the measured RT offset of unmodified cysteines (nothing to do "
                                                    "when cysteines are alkylated)"},
    {"decoys", "<method>", "Decoys in the library: none (DIA-NN makes its own), mutate, pseudo_reverse, reverse, shuffle"},
    {"recompute_decoy_mz", "<true/false>", "Recompute a decoy's precursor m/z from its decoy sequence (with "
                                           "-generation:decoys)"},
    {"rt_model", "<file>", "RT model (ONNX); default: peptdeep_rt_dynamic.onnx in -models"},
    {"ms2_model", "<file>", "MS2 model (ONNX); default: peptdeep_ms2_dynamic.onnx in -models"},
    {"ccs_model", "<file>", "CCS model (ONNX); default: peptdeep_ccs_dynamic.onnx in -models"}};
  const auto defaults = generationDefaults();
  if (defaults.size() != std::size(texts)) { throw std::logic_error("registerGenerationOptions_: a key has no text"); }
  for (const auto& [key, metavar, text] : texts)
  {
    const json& value = defaults.at(key);
    const std::string name = std::string("generation:") + key;
    const bool advanced = advanced_(name);
    if (value.is_boolean())
    {
      registerStringOption_(name, metavar, value.get<bool>() ? "true" : "false", text, false, advanced);
      setValidStrings_(name, {"true", "false"});
    }
    else if (value.is_string())
    { registerStringOption_(name, metavar, value.get<std::string>(), text, false, advanced); }
    else if (value.is_number_integer())
    { registerIntOption_(name, metavar, value.get<int>(), text, false, advanced); setMinInt_(name, 0); }
    else if (value.is_number())
    { registerDoubleOption_(name, metavar, value.get<double>(), text, false, advanced); setMinFloat_(name, std::string(key) == "nce" ? -1 : 0); }
    else if (std::string(key) == "fixed_modifications" || std::string(key) == "variable_modifications")
    {
      OpenMS::StringList values;
      for (const auto& item : value) { values.emplace_back(item.get<std::string>()); }
      registerStringList_(name, metavar, values, text, false, advanced);
    }
    else if (std::string(key) == "precursor_mz" || std::string(key) == "fragment_mz")
    { registerDoubleList_(name, metavar, value.get<std::vector<double>>(), text, false, advanced); }
    else
    { registerIntList_(name, metavar, value.get<std::vector<int>>(), text, false, advanced); }
  }
  setValidStrings_("generation:decoys", {"none", "mutate", "pseudo_reverse", "reverse", "shuffle"});
  setValidStrings_("generation:predict_ion_mobility", {"auto", "true", "false"});
  setValidStrings_("generation:mobility_correction", {"auto", "off"});
  setMinInt_("generation:max_fragment_charge", 1);
  setMaxInt_("generation:max_fragment_charge", 2);
  setMinInt_("generation:precursor_charges", 1);
  setMaxInt_("generation:precursor_charges", 8);
  setMaxFloat_("generation:nce", 100);
  registerInputFile_("irt_standards", "<file>", "", "Generate: iRT standards for -generation:irt_rescale true; "
                     "default: the bundled table", false, advanced_("irt_standards"));
  setValidFormats_("irt_standards", {"tsv"}, false);
  registerInputFile_("instrument_file", "<file>", "", "Generate: the instrument's defaults (NCE, ion mobility, m/z "
                     "windows, charges) from this YAML file. Otherwise one file per instrument from "
                     "$DIALIBGEN_INSTRUMENT_DIR, then the user's DIALibGen/instruments directory (~/.config, "
                     "~/Library/Application Support or %APPDATA%), then the shipped share/DIALibGen/instruments; "
                     "default.yaml for the rest (docs/usage.md#instrument-defaults)", false, advanced_("instrument_file"));
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
    std::string ion_mobility = "auto";   // -generation:predict_ion_mobility, resolved below
    std::string irt_standards;           // "" = the bundled table
    std::string mobility_correction = "auto";   // -generation:mobility_correction, resolved below
    std::optional<std::map<int, double>> nce_by_charge;   // set by a -config JSON; else the instrument's (below)
    std::string nce_by_charge_source = "none";
    // Where the m/z windows came from: the recipe, the command line or INI,
    // else the instrument's method (or DigestParams' generic default).
    std::string precursor_mz_source = "default", fragment_mz_source = "default", charges_source = "default";
    bool instrument_in_config = false;

    const std::string cfg = getStringOption_("config");
    if (!cfg.empty())
    {
      std::ifstream in(cfg);
      if (!in) { writeLogError_("cannot read config: " + cfg); return INPUT_FILE_NOT_FOUND; }
      try
      {
        const json recipe = json::parse(in, nullptr, true, true);
        apply_(recipe, p, decoys, rt_model, ms2_model, ccs_model, nce, instrument, irt_rescale,
               recompute_decoy_mz, nce_was_set, instrument_alias_of, ion_mobility, irt_standards, mobility_correction, nce_by_charge);
        if (nce_by_charge) { nce_by_charge_source = "config"; }
        // A recipe written before 0.13 (it records nce_source, but has no
        // mobility_correction) replays with the 1/K0 it was made with.
        if (recipe.contains("nce_source") && !recipe.contains("mobility_correction"))
        {
          mobility_correction = "off";
          writeLogInfo_("the recipe predates the 1/K0 correction (no mobility_correction): replayed without it "
                        "(-generation:mobility_correction auto to apply it)");
        }
        instrument_in_config = recipe.contains("instrument");
        if (recipe.contains("precursor_mz")) { precursor_mz_source = "config"; }
        if (recipe.contains("fragment_mz")) { fragment_mz_source = "config"; }
        if (recipe.contains("precursor_charges")) { charges_source = "config"; }
      }
      catch (const std::exception& e)
      { writeLogError_(std::string("config: ") + e.what()); return ILLEGAL_PARAMETERS; }

      const auto base = fs::path(cfg).parent_path();
      for (std::string* m : {&rt_model, &ms2_model, &ccs_model, &irt_standards})
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
      // An empty window is the default: the instrument's, resolved below.
      for (const char* key : {"precursor_mz", "fragment_mz"})
      { if (overrides.contains(key) && overrides[key].empty()) { overrides.erase(key); } }
      if (overrides.contains("precursor_mz")) { precursor_mz_source = "TOPP"; }
      if (overrides.contains("fragment_mz")) { fragment_mz_source = "TOPP"; }
      if (given_("generation:precursor_charges")) { charges_source = "TOPP"; }
      // The raw path: TOPP's input-file getter refuses a missing file, which
      // matters only when the standards are used (irt_rescale).
      if (cfg.empty() || supplied_.count("irt_standards")) { overrides["irt_standards"] = getParam_().getValue("irt_standards").toString(); }
      const bool automatic_nce = overrides.contains("nce") && overrides["nce"] == -1.0;
      if (automatic_nce) { overrides.erase("nce"); }
      else if (overrides.contains("nce")) { nce_source = "TOPP"; }
      apply_(overrides, p, decoys, rt_model, ms2_model, ccs_model, nce, instrument,
             irt_rescale, recompute_decoy_mz, nce_was_set, instrument_alias_of, ion_mobility, irt_standards, mobility_correction, nce_by_charge);
      if (automatic_nce) { nce_was_set = false; }
      // An explicit -generation:nce is one NCE for every charge, also over a
      // recipe's nce_by_charge.
      if (nce_source == "TOPP") { nce_by_charge = std::map<int, double>{}; nce_by_charge_source = "TOPP"; }
    }
    catch (const std::exception& e)
    { writeLogError_(std::string("generation: ") + e.what()); return ILLEGAL_PARAMETERS; }

    // The instrument's defaults: an explicit option > its YAML file > default.yaml
    // (docs/usage.md#instrument-defaults).
    namespace inst = ODIA::instruments;
    std::unique_ptr<inst::Catalog> catalog;
    std::optional<inst::InstrumentFile> instrument_file;
    bool ion_mobility_default = false;
    std::string named;   // the instrument as the caller (or default.yaml) named it
    try
    {
      catalog = std::make_unique<inst::Catalog>(inst::searchPath(dataDirCandidates()));
      const auto& generic = catalog->defaults();
      const bool explicit_instrument = given_("generation:instrument") || instrument_in_config;
      named = explicit_instrument ? instrument : generic.instrument;
      if (const std::string file = getParam_().getValue("instrument_file").toString(); !file.empty())
      {
        if (!fs::exists(file)) { throw std::runtime_error("-instrument_file " + file + " does not exist"); }
        instrument_file = inst::readInstrumentFile(file);
        if (!explicit_instrument) { named = instrument_file->name; }
        else if (!instrument_file->names(named))
        {
          throw std::runtime_error("-generation:instrument " + named + " is not the instrument of " + file +
                                   " (" + instrument_file->name + "); give one of them");
        }
      }
      else if (const inst::InstrumentFile* f = catalog->find(named)) { instrument_file = *f; }
      const std::string model = instrument_file ? instrument_file->model_instrument
                                                : ODIA::PeptDeepEncoder::canonicalInstrument(named);
      if (model.empty())
      {
        std::string files;
        for (const auto& f : catalog->files()) { files += (files.empty() ? "" : ", ") + f.name; }
        throw std::runtime_error(
          "unknown instrument: '" + named + "'. Use QE, Lumos, timsTOF, SciexTOF or ThermoTOF, a vendor name such as "
          "Astral, Exploris 480 or timsTOF Pro, or an instrument file (" + files + "; -instrument_file, "
          "docs/usage.md#instrument-defaults); Lumos means no instrument correction.");
      }
      instrument = model;
      instrument_alias_of = named != model ? named : std::string();
      const auto slot = generic.models.find(model);
      const auto pick = [&](const auto& from_file, const auto& generic_value) {
        return from_file ? *from_file : generic_value;
      };
      const std::string label = "instrument-default:" + (instrument_file ? instrument_file->name : model);
      if (precursor_mz_source == "default")
      {
        const inst::Window w = pick(instrument_file ? instrument_file->precursor_mz : std::nullopt, generic.precursor_mz);
        p.precursor_mz_min = w[0]; p.precursor_mz_max = w[1];
        if (instrument_file && instrument_file->precursor_mz) { precursor_mz_source = label; }
      }
      if (fragment_mz_source == "default")
      {
        const inst::Window w = pick(instrument_file ? instrument_file->fragment_mz : std::nullopt, generic.fragment_mz);
        p.fragment_mz_min = w[0]; p.fragment_mz_max = w[1];
        if (instrument_file && instrument_file->fragment_mz) { fragment_mz_source = label; }
      }
      if (charges_source == "default")
      {
        p.charges = pick(instrument_file ? instrument_file->precursor_charges : std::nullopt, generic.precursor_charges);
        if (instrument_file && instrument_file->precursor_charges) { charges_source = label; }
      }
      if (!nce_was_set)
      {
        nce = instrument_file && instrument_file->nce ? *instrument_file->nce
            : slot != generic.models.end() && slot->second.nce ? *slot->second.nce
            : ODIA::PeptDeepEncoder::defaultNce(model);
        nce_source = label;
        // The instrument's per-charge NCE goes with its NCE: an explicit
        // -generation:nce (one value) or a -config nce_by_charge replaces it.
        if (!nce_by_charge && instrument_file && instrument_file->nce_by_charge)
        {
          nce_by_charge = *instrument_file->nce_by_charge;
          nce_by_charge_source = label;
        }
      }
      ion_mobility_default = instrument_file && instrument_file->ion_mobility ? *instrument_file->ion_mobility
                           : slot != generic.models.end() && slot->second.ion_mobility ? *slot->second.ion_mobility
                           : ODIA::PeptDeepEncoder::hasIonMobility(model);
    }
    catch (const std::exception& e)
    { writeLogError_(std::string("generation: ") + e.what()); return ILLEGAL_PARAMETERS; }
    {
      std::ostringstream m;
      m << std::defaultfloat << "instrument defaults: "
        << (instrument_file ? instrument_file->path.filename().string() + " (" + instrument_file->name + ") over " : std::string())
        << catalog->defaults().path.filename().string();
      if (ODIA::diagnosticsLevel() > 0)
      {
        m << " [" << (instrument_file ? instrument_file->path.string() + ", " : std::string())
          << catalog->defaults().path.string() << "; searched:";
        for (const auto& d : catalog->dirs()) { m << " " << d.string(); }
        m << "]";
      }
      writeLogInfo_(m.str());
      std::vector<std::string> implied;
      std::ostringstream w;
      w << std::defaultfloat;
      if (precursor_mz_source.starts_with("instrument-default:"))
      { w << "precursor m/z " << p.precursor_mz_min << "-" << p.precursor_mz_max; implied.push_back("-generation:precursor_mz"); }
      if (fragment_mz_source.starts_with("instrument-default:"))
      {
        w << (implied.empty() ? "" : ", ") << "fragment m/z " << p.fragment_mz_min << "-" << p.fragment_mz_max;
        implied.push_back("-generation:fragment_mz");
      }
      if (charges_source.starts_with("instrument-default:"))
      {
        w << (implied.empty() ? "" : ", ") << "charges";
        for (const int z : p.charges) { w << " " << z; }
        implied.push_back("-generation:precursor_charges");
      }
      if (!implied.empty())
      {
        std::string options;
        for (const auto& o : implied) { options += (options.empty() ? "" : ", ") + o; }
        writeLogInfo_(instrument_file->name + " method defaults: " + w.str() + " (" + options + " to change them" +
                      (implied.front() == "-generation:precursor_mz" && (p.precursor_mz_min > 350 || p.precursor_mz_max < 1200)
                         ? ", e.g. -generation:precursor_mz 350 1200 for a "
                                                                        "wider window)" : std::string(")")));
      }
    }
    if (!instrument_alias_of.empty())
    { writeLogInfo_("instrument '" + instrument_alias_of + "' is predicted as " + instrument); }
    if (!nce_was_set)
    {
      std::ostringstream m;
      m << "NCE not set; using " << named << "'s default, " << nce;
      if (nce_by_charge && !nce_by_charge->empty())
      {
        m << std::defaultfloat << " (by precursor charge:";
        for (const auto& [z, v] : *nce_by_charge) { m << " " << z << "+ " << v; }
        m << "; other charges " << nce << ")";
      }
      writeLogInfo_(m.str());
    }
    else if (nce_by_charge && !nce_by_charge->empty())
    {
      std::ostringstream m;
      m << std::defaultfloat << "NCE by precursor charge (" << nce_by_charge_source << "):";
      for (const auto& [z, v] : *nce_by_charge) { m << " " << z << "+ " << v; }
      m << "; other charges " << nce;
      writeLogInfo_(m.str());
    }
    const std::map<int, double> nce_charge = nce_by_charge ? *nce_by_charge : std::map<int, double>{};
    std::map<int, float> nce_charge_f;
    for (const auto& [z, v] : nce_charge) { nce_charge_f[z] = static_cast<float>(v); }
    // CCS and 1/K0 only where the instrument can use them, unless the
    // caller says otherwise. One line either way: the columns' absence (or
    // presence) is otherwise a surprise in the consumer.
    p.predict_ion_mobility = ion_mobility == "auto" ? ion_mobility_default : ion_mobility == "true";
    {
      const bool instrument_named = given_("generation:instrument") || instrument_in_config;
      if (ion_mobility != "auto")
      {
        writeLogInfo_(std::string(p.predict_ion_mobility ? "CCS and 1/K0 predicted" : "CCS and 1/K0 not predicted") +
                      " (-generation:predict_ion_mobility " + ion_mobility + ")");
      }
      else if (!p.predict_ion_mobility && !instrument_named && getParam_().getValue("instrument_file").toString().empty())
      {
        writeLogInfo_("no -generation:instrument: " + named + " assumed, so CCS and 1/K0 are not predicted "
                      "(name the instrument, e.g. -generation:instrument timsTOF for ion mobility)");
      }
      else if (!p.predict_ion_mobility)
      {
        writeLogInfo_(named + " has no ion mobility: CCS and 1/K0 not predicted "
                      "(-generation:predict_ion_mobility true to predict them)");
      }
      else
      { writeLogInfo_(named + " has ion mobility: CCS and 1/K0 predicted"); }
      if (p.predict_ion_mobility && !p.derive_ion_mobility)
      { writeLogInfo_("1/K0 not derived from CCS (-generation:derive_ion_mobility false)"); }
    }
    // Settings without an effect here are refused (a WARNING in 0.12);
    // -irt_standards only warns: alone it is harmless.
    {
      const ODIA::DigestParams defaults;
      const auto inert = [this](const std::string& option, const std::string& because, const std::string& fix)
      { writeLogError_(option + " has no effect " + because + ": " + fix); };
      if (!irt_standards.empty() && !irt_rescale)
      { warn_("-irt_standards has no effect without -generation:irt_rescale true"); }
      if (recompute_decoy_mz && decoys == "none")
      {
        inert("-generation:recompute_decoy_mz true", "with -generation:decoys none",
              "add -generation:decoys or remove -generation:recompute_decoy_mz");
        return ILLEGAL_PARAMETERS;
      }
      if (p.max_variable_modifications != defaults.max_variable_modifications && p.variable_modifications.empty())
      {
        inert("-generation:max_variable_modifications", "without -generation:variable_modifications",
              "add -generation:variable_modifications or remove -generation:max_variable_modifications");
        return ILLEGAL_PARAMETERS;
      }
      if (!p.predict_ion_mobility && !p.derive_ion_mobility)
      {
        inert("-generation:derive_ion_mobility false", "when CCS is not predicted",
              "remove it, or predict CCS (-generation:predict_ion_mobility true, or an ion-mobility instrument)");
        return ILLEGAL_PARAMETERS;
      }
    }
    // Neither instrument has trained weights in the shipped MS2 checkpoint
    // (its own constants name four instruments). Measured against Lumos at
    // the same NCE (median spectral cosine per precursor, synthetic
    // peptides): ThermoTOF 0.993 -- effectively Lumos -- but SciexTOF only
    // 0.914, with a tenth of the precursors below 0.64. On SCIEX data (five
    // public TripleTOF and ZenoTOF runs, all-charge b/y cosine from the raw
    // spectra) the SciexTOF slot scores 0.74 against Lumos 27's 0.83, which is
    // why sciextof.yaml maps the TripleTOF and ZenoTOF names to Lumos.
    if (instrument == "ThermoTOF")
    { warn_("ThermoTOF has no trained MS2 weights; its spectra are effectively Lumos's. Use "
            "-generation:instrument Lumos to say so."); }
    if (instrument == "SciexTOF")
    { warn_("SciexTOF has no trained MS2 weights; its spectra are NOT those of Lumos and match SCIEX data less "
            "well. -generation:instrument TripleTOF or ZenoTOF uses sciextof.yaml (Lumos, NCE 27)."); }

    const std::string models_dir = getParam_().getValue("models").toString();
    struct { const char* name; const char* file; std::string* path; } models[] = {
      {"rt_model", kRtModelFile, &rt_model},
      {"ms2_model", kMs2ModelFile, &ms2_model},
      {"ccs_model", kCcsModelFile, &ccs_model}};
    for (auto& m : models)
    {
      if (m.path->empty()) { *m.path = findModel(m.file, models_dir); }
    }

    // The instrument's 1/K0 correction (timstof.yaml, mobility_correction):
    // for the CCS model it was fitted for, so a tuned CCS model (trained on a
    // run's own 1/K0) is never corrected twice.
    ODIA::MobilityCorrection mobility;
    std::string mobility_note;
    if (mobility_correction == "auto" && p.predict_ion_mobility && p.derive_ion_mobility && instrument_file &&
        instrument_file->mobility_correction)
    {
      std::error_code exists_ec;
      if (!ccs_model.empty() && fs::exists(ccs_model, exists_ec))
      {
        try
        {
          if (ODIA::DIANNLibraryFile::hashFile(ccs_model) == instrument_file->mobility_correction->ccs_model_fnv1a64)
          { mobility = instrument_file->mobility_correction->correction; }
          else { mobility_note = "the CCS model is not the one " + instrument_file->path.filename().string() +
                                 "'s mobility_correction was fitted for (a tuned model needs none)"; }
        }
        catch (const std::exception& e) { writeLogError_(e.what()); return INPUT_FILE_NOT_FOUND; }
      }
    }
    if (mobility.active)
    {
      std::ostringstream m;
      m << std::defaultfloat << "1/K0 corrected with " << instrument_file->path.filename().string()
        << "'s mobility_correction (intercept + slope x the Mason-Schamp 1/K0 of the predicted CCS:";
      for (const auto& [z, l] : mobility.by_charge) { m << " " << z << "+ " << l.intercept << " + " << l.slope << "x;"; }
      m << " other charges " << mobility.other.intercept << " + " << mobility.other.slope
        << "x); -generation:mobility_correction off for the uncorrected value";
      writeLogInfo_(m.str());
    }
    else if (!mobility_note.empty()) { writeLogInfo_("1/K0 not corrected: " + mobility_note); }

    const auto recipe = [&]() {
      json eff = effectiveConfig(p, decoys, rt_model, ms2_model, ccs_model, nce, instrument, irt_rescale,
                                 recompute_decoy_mz, ion_mobility, irt_standards, mobility_correction,
                                 nce_charge);
      eff["nce_by_charge_source"] = nce_by_charge_source;
      // The 1/K0 correction applied, with the file it came from.
      if (mobility.active)
      {
        json lines = json::object();
        for (const auto& [z, l] : mobility.by_charge) { lines[std::to_string(z)] = {l.intercept, l.slope}; }
        lines["other"] = {mobility.other.intercept, mobility.other.slope};
        eff["mobility_correction_applied"] = {{"file", instrument_file->path.filename().string()},
                                              {"intercept_slope", lines}};
      }
      else { eff["mobility_correction_applied"] = nullptr; }
      // The instrument as named, so that a replay resolves it through the same
      // instrument file (a custom instrument keeps its ion mobility); the MS2
      // model's slot is recorded beside it.
      eff["instrument"] = named;
      eff["model_instrument"] = instrument;
      eff["nce_source"] = nce_source;
      eff["precursor_mz_source"] = precursor_mz_source;
      eff["fragment_mz_source"] = fragment_mz_source;
      eff["precursor_charges_source"] = charges_source;
      // Which files gave the defaults (name and FNV-1a64): the values are in
      // the recipe and the fingerprint already.
      eff["instrument_defaults"] = {
        {"defaults", catalog->defaults().path.filename().string()}, {"defaults_fnv1a64", catalog->defaults().hash},
        {"instrument", instrument_file ? json(instrument_file->path.filename().string()) : json(nullptr)},
        {"instrument_fnv1a64", instrument_file ? json(instrument_file->hash) : json(nullptr)}};
      return eff;
    };
    if (const std::string wc = getStringOption_("write_config"); !wc.empty())
    {
      const json eff = recipe();
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
    if ((format == ODIA::LibraryFormat::OpenSwathPQP || format == ODIA::LibraryFormat::OpenSwathTSV) && decoys == "none")
    {
      writeLogInfo_("OpenSWATH scores against the library's decoys and this one has none: add -generation:decoys "
                    "shuffle (or mutate, pseudo_reverse, reverse)");
    }

    json eff = recipe();
    std::error_code ec;
    for (const auto& m : models)
    {
      if (m.path == &ccs_model && !p.predict_ion_mobility) { continue; }   // not used
      if (m.path->empty())
      {
        std::ostringstream e;
        e << m.file << " not found; give -models (or set DIALIBGEN_MODEL_DIR) with the directory holding the three "
          << "PeptDeep models, or -generation:" << m.name;
        if (ODIA::diagnosticsLevel() > 0)
        { e << ". Searched:"; for (const auto& d : modelDirCandidates(models_dir)) { e << "\n  " << d.string(); } }
        else { e << " (-debug 1 lists the directories searched)"; }
        writeLogError_(e.str());
        return INPUT_FILE_NOT_FOUND;
      }
      if (!fs::exists(*m.path, ec))
      {
        writeLogError_("-generation:" + std::string(m.name) + " " + *m.path + " does not exist");
        return INPUT_FILE_NOT_FOUND;
      }
    }

    if (const std::string device = getStringOption_("machine:device"); !ODIA::validDeviceName(device))
    { writeLogError_("-machine:device '" + device + "': use auto, cpu, cuda, cuda:N, coreml or mps"); return ILLEGAL_PARAMETERS; }
    else if (device == "cpu" && getFlag_("machine:require_gpu"))
    { writeLogError_("-machine:require_gpu contradicts -machine:device cpu; drop one of them"); return ILLEGAL_PARAMETERS; }
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
    // CoreML compiles one MS2 shape per residue count and batch size: the
    // shapes the library needs are built in the background from the
    // predictor's construction on (filled in after the digest, below).
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
      p.predict_ion_mobility ? ODIA::DIANNLibraryFile::hashFile(ccs_model) : std::string("none"), nce, instrument,
      irt_rescale, mobility, nce_charge_f);
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
    // Where RT and CCS run on the CPU while MS2 runs on a GPU, they can run
    // meanwhile, on a thread of their own: MS2 does not need them, the
    // writer waits for both, and the CPU's RT and CCS calls are the same
    // calls as before, so the library is the same. That is the default on
    // CUDA without cuDNN; on CoreML only with ODIA_OVERLAP=1, since the
    // GPU's host side then lost more than the overlap gained (0.13.0, whole
    // human library on an M4 Max: 343.9 s against 297.1 s in order). With
    // RT and CCS on the GPU they stay in order
    // (LibraryGenerator::overlapsCpuPhases). ODIA_OVERLAP=0 keeps the phases
    // in order everywhere, for comparison.
    const char* overlap_env = std::getenv("ODIA_OVERLAP");
    const bool overlap_wanted = overlap_env == nullptr || std::string_view(overlap_env) != "0";
    const bool coreml_overlap = overlap_env != nullptr && std::string_view(overlap_env) == "1";
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

      // Where each model runs: -machine:device, and for auto what this
      // machine offers and the library's size. One line, and a WARNING with
      // the fix when a GPU is there but RT and CCS cannot use it.
      const ODIA::DeviceCaps caps = ODIA::probeDeviceCaps();
      if (ODIA::diagnosticsLevel() > 0) { for (const auto& line : caps.diagnostics) { writeLogInfo_("device check: " + line); } }
      const ODIA::DeviceRoute route = ODIA::chooseDeviceRoute(getStringOption_("machine:device"), gpu.require, caps,
                                                              library.precursorCount(),
                                                              ODIA::coremlMinPrecursors(threads == 0 ? ODIA::physicalCores()
                                                                                         : std::min(threads, ODIA::physicalCores())),
                                                              p.predict_ion_mobility);
      if (!route.error.empty()) { throw std::runtime_error(route.error); }
      writeLogInfo_("device: " + route.summary);
      if (!route.warning.empty())
      {
        warn_(route.warning);
        device_warnings.push_back(route.warning);
      }
      ODIA::GpuOptions lstm_gpu = gpu;
      lstm_gpu.coreml = route.lstm_coreml;

      // Sessions of the CPU phases: all of -threads, or fewer while they
      // overlap the GPU's MS2 phase (below).
      unsigned cpu_threads = threads;
      // The RT phase and the iRT calibration; @p info takes the log lines
      // (written at once on the main thread, collected on the CPU thread).
      const auto rtPhase = [&](const std::function<void(const std::string&)>& info) {
        if (const auto miss = ODIA::LibraryGenerator::predictRetentionTimes(
              library, rt_model, route.lstm_gpu, cpu_threads, p.free_cysteine_rt_correction, lstm_gpu); miss)
        {
          throw PhaseFailure(UNEXPECTED_RESULT, std::to_string(miss) + " precursors have no predicted RT (a "
                             "modification the RT model cannot encode?); no library written");
        }
        if (irt_rescale)
        {
          standards = irt_standards;
          if (!standards.empty() && !fs::exists(standards))
          { throw PhaseFailure(INPUT_FILE_NOT_FOUND, "-irt_standards " + standards + " does not exist"); }
          if (standards.empty()) { standards = findDataFile("irt_standards.tsv"); }
          if (standards.empty())
          {
            std::ostringstream e;
            e << "no iRT standards file found; give -irt_standards";
            if (ODIA::diagnosticsLevel() > 0)
            { e << ". Searched:"; for (const auto& d : dataDirCandidates()) { e << "\n  " << (d / "irt_standards.tsv").string(); } }
            throw PhaseFailure(INPUT_FILE_NOT_FOUND, e.str());
          }
          const auto cal = ODIA::LibraryGenerator::fitIrtCalibration(rt_model, standards, route.lstm_gpu, lstm_gpu);
          ODIA::LibraryGenerator::applyIrtCalibration(library, cal);
          std::ostringstream m;
          m.setf(std::ios::fixed); m.precision(1);
          m << "rescaled to iRT with " << cal.peptides << " standards (largest error " << cal.max_abs_error << " iRT)";
          if (ODIA::diagnosticsLevel() > 0)
          { m.precision(4); m << ": " << cal.slope << " * raw + " << cal.intercept; }
          info(m.str());
        }
        else
        {
          info("RT is the model's raw 0..1 scale, not iRT; -generation:irt_rescale true writes iRT for "
               "a library another tool reads as iRT");
        }
        finishFingerprint();
      };
      const auto logInfo = [&](const std::string& line) { writeLogInfo_(line); };

      // The MS2 predictor first: its provider decides whether RT and CCS
      // overlap it, and on CoreML its shape builds start now, for the
      // residue counts the library holds (counted from the modified
      // sequences: the residues are the letters outside the modification
      // brackets) and only the batch sizes their chunks can need.
      {
        const auto& precursors = library.precursors();
        std::map<int, std::uint64_t> per_length;
        for (std::size_t i = 0; i < library.precursorCount(); ++i)
        {
          int residues = 0, depth = 0;
          for (const char c : library.strings().get(precursors.modified_sequence[i]))
          {
            if (c == '(' || c == '[') { ++depth; }
            else if (c == ')' || c == ']') { if (depth > 0) { --depth; } }
            else if (depth == 0 && c >= 'A' && c <= 'Z') { ++residues; }
          }
          if (residues > 0) { ++per_length[residues]; }
        }
        gpu.prebuild_lengths.assign(per_length.begin(), per_length.end());
      }
      ODIA::GpuOptions ms2_gpu = gpu;
      ms2_gpu.coreml = route.ms2_coreml;
      std::unique_ptr<ODIA::PeptDeepPredictor> ms2;
      if (stream) { ms2 = ODIA::LibraryGenerator::openMs2Predictor(ms2_model, route.ms2_gpu, threads, ms2_gpu); }
      const bool overlap = overlap_wanted && ms2 &&
                           ODIA::LibraryGenerator::overlapsCpuPhases(ms2->provider(), !route.lstm_gpu, coreml_overlap);
      // While the GPU predicts MS2, its host side (Core ML's own threads,
      // the encoding workers) needs cores: with all 16 CPU sessions busy on
      // CCS, the Astral MS2 phase on an M4 Max took 229 s instead of 153 s
      // and the overlap gained 6 s of 238 (with 12 sessions 208 s, with 8
      // 189 s under a noisier load). The CPU phases therefore leave half
      // the cores free (ODIA_OVERLAP_RESERVE overrides the count); their
      // chunking does not depend on the session count, so the library does
      // not either.
      if (overlap)
      {
        const unsigned all = threads == 0 ? ODIA::physicalCores() : threads;
        unsigned reserve = all / 2;
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
          warn_(std::string("cannot start the streaming writer (") + e.what() +
                "); the library is written after prediction instead (same output, more memory)");
        }
      }
      if (sink)
      {
        const auto ccsPhase = [&]() {
          if (!p.predict_ion_mobility)
          {
            // Nothing to settle: the guard does nothing.
            guard = std::make_unique<ODIA::LibraryGenerator::CcsBatchGuard>(ccs_model, false,
                                                                             std::vector<std::uint16_t>(), false);
            return;
          }
          std::vector<std::uint16_t> lengths;
          bool ccs_on_cpu = false;
          if (const auto miss = ODIA::LibraryGenerator::predictCollisionCrossSections(
                library, ccs_model, route.lstm_gpu, cpu_threads, p.derive_ion_mobility, lstm_gpu, &lengths, &ccs_on_cpu,
                mobility); miss)
          { throw std::runtime_error(std::to_string(miss) + " precursors have no predicted CCS"); }
          guard = std::make_unique<ODIA::LibraryGenerator::CcsBatchGuard>(ccs_model, p.derive_ion_mobility,
                                                                           std::move(lengths), ccs_on_cpu, mobility);
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
          writeLogInfo_(std::string(p.predict_ion_mobility ? "RT and CCS" : "RT") + " predicted on the CPU (" + std::to_string(cpu_threads) +
                        " sessions) while the GPU predicts MS2");
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
                library, *ms2, p, static_cast<float>(nce), instrument, threads, sink.get(), nce_charge_f); miss)
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
        {
          if (ODIA::diagnosticsLevel() > 0)
          { writeLogInfo_("CCS of " + std::to_string(sink->resettled()) + " precursors re-predicted in their one-pass batching"); }
        }
        writeLogInfo_("after fragment selection: " + std::to_string(library.precursorCount()) + " precursors");
        if (library.precursorCount() == 0)
        { throw std::runtime_error("no precursors remain after digestion and fragment selection"); }
      }
      else
      {
        // Overlap chosen but no streaming writer: RT has not run yet, and
        // runs on all of -threads, as without the overlap.
        if (overlap)
        {
          cpu_threads = threads;
          rtPhase(logInfo);
        }
        if (!ms2) { ms2 = ODIA::LibraryGenerator::openMs2Predictor(ms2_model, route.ms2_gpu, threads, ms2_gpu); }
        if (const auto miss = ODIA::LibraryGenerator::predictFragmentIntensities(
              library, *ms2, p, static_cast<float>(nce), instrument, threads, nullptr, nce_charge_f); miss)
        { throw std::runtime_error(std::to_string(miss) + " precursors have no predicted MS2 spectrum"); }
        ms2.reset();
        writeLogInfo_("after fragment selection: " + std::to_string(library.precursorCount()) + " precursors");
        if (library.precursorCount() == 0)
        { throw std::runtime_error("no precursors remain after digestion and fragment selection"); }
        if (p.predict_ion_mobility)
        {
          if (const auto miss = ODIA::LibraryGenerator::predictCollisionCrossSections(
                library, ccs_model, route.lstm_gpu, threads, p.derive_ion_mobility, lstm_gpu, nullptr, nullptr,
                mobility); miss)
          { throw std::runtime_error(std::to_string(miss) + " precursors have no predicted CCS"); }
        }
      }
      std::size_t skipped = 0;
      const auto method = ODIA::parseDecoyMethod(decoys);
      const auto made = ODIA::LibraryGenerator::appendDecoys(library, method, &skipped,
                                                             p.min_fragments,
                                                             recompute_decoy_mz, mobility);
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
    {
      const std::string what = e.what();
      writeLogError_(what.rfind("the GPU is required", 0) == 0 ? what : "library generation failed: " + what);
      return UNEXPECTED_RESULT;
    }

    if (!device_warnings.empty())
    {
      warn_("part of the prediction ran on the CPU instead of the GPU (see the WARNINGs above)" +
            std::string(parquet ? "; recorded in the library's metadata (odia.warnings)" : ""));
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
        writeLogInfo_("the TSV does not embed the generation settings; .parquet keeps them with the library");
      }
      else if (!ODIA::LibraryExport::carriesRecipe(format))
      {
        writeLogInfo_(std::string(ODIA::LibraryExport::name(format)) + " does not embed the generation settings; "
                      ".parquet keeps them with the library");
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
