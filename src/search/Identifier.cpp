// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <odia/search/Identifier.h>

#include <odia/search/EvidencePrefilter.h>

#include <nlohmann/json.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace ODIA::search
{
  namespace
  {
    using json = nlohmann::json;
    using Clock = std::chrono::steady_clock;

    double since(Clock::time_point t) { return std::chrono::duration<double>(Clock::now() - t).count(); }

    std::string seconds(double s)
    {
      std::ostringstream o;
      o.setf(std::ios::fixed);
      o.precision(1);
      o << s << " s";
      return o.str();
    }

    json num(double v) { return std::isfinite(v) ? json(v) : json(nullptr); }

    std::string fixed2(double v)
    {
      std::ostringstream o;
      o.setf(std::ios::fixed);
      o.precision(2);
      o << v;
      return o.str();
    }

    json object(const std::string& text) { return text.empty() ? json::object() : json::parse(text); }

    json selectionJson(const SelectionStats& s)
    {
      return {{"library_precursors", s.library_precursors}, {"library_decoys_ignored", s.library_decoys_ignored},
              {"targets", s.targets}, {"ineligible_mz", s.ineligible_mz}, {"ineligible_charge", s.ineligible_charge},
              {"ineligible_rt", s.ineligible_rt}, {"ineligible_fragments", s.ineligible_fragments},
              {"ineligible_window", s.ineligible_window}, {"isolation_windows", s.windows},
              {"duplicate_key", s.duplicate_key}, {"eligible", s.eligible}, {"drawn", s.drawn},
              {"no_decoy", s.no_decoy},
              {"no_decoy_reasons", {{"unparsable", s.decoy_unparsable}, {"unshufflable", s.decoy_unshufflable},
                                    {"out_of_range", s.decoy_out_of_range}, {"copy", s.decoy_copy},
                                    {"too_few_fragments", s.decoy_too_few_fragments}}},
              {"decoys_redrawn", s.decoy_redrawn}, {"fragment_slots_dropped", s.fragment_slots_dropped},
              {"decoy_fragment_mz_range", {s.fragment_mz_min, s.fragment_mz_max}},
              {"capped", s.capped}, {"pairs", s.pairs}};
    }

    json scoringJson(const ScoringOutcome& o)
    {
      const auto& d = o.scored.diagnostics;
      json j = {
        {"peak_groups", d.rows}, {"target_precursors", d.target_groups}, {"decoy_precursors", d.decoy_groups},
        {"pairs_complete", d.pairs_complete}, {"targets_unpaired", d.targets_unpaired}, {"decoys_unpaired", d.decoys_unpaired},
        {"features_used", d.features_used}, {"features_dropped", d.features_dropped},
        {"features_excluded", d.features_excluded}, {"exclusions_unmatched", d.exclusions_unmatched},
        {"cells_imputed", d.cells_imputed}, {"folds", d.n_folds}, {"iterations_trained", d.iterations_trained},
        {"iterations_skipped", d.iterations_skipped}, {"folds_unscaled", d.folds_unscaled},
        {"target_winners", d.target_winners}, {"decoy_winners", d.decoy_winners},
        {"estimator", "concatenated pair competition, q = (D + 1) / T; peptides and protein groups by picked competition"},
        {"null_balance", {{"fraction", SearchParams::null_balance_fraction}, {"targets", o.null_targets},
                          {"decoys", o.null_decoys}, {"ratio", num(o.null_ratio)}, {"z", num(o.null_z)},
                          {"warn_z", SearchParams::null_balance_warn_z}}}};
      return j;
    }

    json identificationsJson(const ScoringOutcome& o)
    {
      const auto& d = o.scored.diagnostics;
      return {{"q", SearchParams::identification_q}, {"precursors", d.targets_at_q}, {"decoy_precursors", d.decoys_at_q},
              {"pooled_precursors", d.pooled_targets_at_q}, {"pooled_vs_paired", num(d.pooled_vs_paired)},
              {"peptides", o.peptides_at_q}, {"peptide_pairs", o.peptide_pairs},
              {"protein_groups", o.proteins_at_q}, {"protein_group_pairs", o.protein_pairs}};
    }

    json entrapmentJson(const ScoringOutcome& o)
    {
      if (!o.entrapment) { return nullptr; }
      const auto& e = o.entrapment_estimate;
      return {{"estimator", "combined"}, {"valid", e.valid}, {"identified", e.n_reported}, {"entrapment", e.n_entrapment},
              {"shared_left_out", o.entrapment_shared}, {"db_target", e.db_target}, {"db_entrapment", e.db_entrapment},
              {"ratio", num(e.ratio)}, {"fdp", num(e.fdp)}};
    }

    json selftestJson(const ScoringOutcome& o)
    {
      if (!o.selftest) { return nullptr; }
      return {{"label_swap_ids", o.selftest_label_swap_ids}, {"random_label_ids", o.selftest_random_label_ids},
              {"limit", o.selftest_limit}};
    }
  }

  Identifier::Identifier(SearchParams params, Log info, Log warn)
    : params_(std::move(params)), info_(std::move(info)), warn_(std::move(warn))
  {
  }

  Identifier::~Identifier() = default;

  void Identifier::info(const std::string& message) const
  {
    if (info_) { info_(message); } else { std::cout << message << std::endl; }
  }

  void Identifier::warn(const std::string& message) const
  {
    if (warn_) { warn_(message); } else { std::cerr << "Warning: " << message << std::endl; }
  }

  std::filesystem::path Identifier::scratchParent() const
  {
    return output_dir_.empty() ? std::filesystem::temp_directory_path() : output_dir_;
  }

  std::vector<ReportRow> Identifier::reportRows(const SearchSet& set, const PeakGroups& groups,
                                                const ScoringOutcome& outcome, const SearchParams& params)
  {
    const auto& p = set.library.precursors();
    const auto& result = outcome.scored.groups;
    std::vector<ReportRow> rows;
    for (std::size_t g = 0; g < result.size(); ++g)
    {
      const auto& gr = result[g];
      if (!gr.winner || gr.qvalue > params.report_max_q) { continue; }
      const auto i = static_cast<std::size_t>(gr.group);
      const std::size_t r = gr.best_row;
      ReportRow row;
      row.precursor_id = set.precursorId(i);
      row.modified_sequence = std::string(set.modifiedSequence(i));
      row.charge = set.charge(i);
      row.precursor_mz = fromFixed(p.mz[i]);
      row.protein_group = std::string(set.proteinGroup(i));
      row.decoy = set.isDecoy(i);
      row.rt = groups.apex_rt[r] / 60.0;
      row.rt_start = groups.rt_start[r] / 60.0;
      row.rt_stop = groups.rt_stop[r] / 60.0;
      row.irt = p.irt[i];
      row.im = groups.im[r];
      row.q = gr.qvalue;
      row.global_q = outcome.peptide_q[g];
      row.pg_q = outcome.protein_q[g];
      row.pep = gr.pep;
      row.evidence = gr.score;
      rows.push_back(std::move(row));
    }
    return rows;
  }

  std::vector<IsolationWindow> Identifier::isolationWindows(const RunData& run)
  {
    std::vector<IsolationWindow> out;
    for (const auto& m : run.maps)
    {
      if (!m.ms1) { out.push_back({m.lower, m.upper}); }
    }
    std::sort(out.begin(), out.end(), [](const IsolationWindow& a, const IsolationWindow& b) {
      return a.lower != b.lower ? a.lower < b.lower : a.upper < b.upper;
    });
    return out;
  }

  IdentificationResult Identifier::identify(const Library& library, const std::string& run_path,
                                            const std::string& out_ids, const std::string& tool_version)
  {
    params_.validate();
    if (std::filesystem::exists(out_ids) || std::filesystem::is_symlink(out_ids))
    { throw std::runtime_error("refusing to overwrite existing output: " + out_ids); }
    output_dir_ = std::filesystem::absolute(out_ids).parent_path();
    const auto started = Clock::now();
    json timing = json::object();
    json warnings = json::array();
    auto warning = [&](const std::string& message) { warn("search: " + message); warnings.push_back(message); };

    // The scratch directory of a cache-mode load goes when this scope ends,
    // whether the search succeeded or not. It is declared BEFORE the run, so
    // the run -- and with it every open cache file -- is destroyed first:
    // Windows cannot delete a file that is still open, and a failure is said.
    struct Scratch
    {
      const Identifier* owner = nullptr;
      std::filesystem::path path;
      void remove()
      {
        if (path.empty()) { return; }
        std::error_code ec;
        std::uintmax_t bytes = 0;
        for (auto it = std::filesystem::recursive_directory_iterator(path, ec);
             !ec && it != std::filesystem::recursive_directory_iterator(); it.increment(ec))
        {
          std::error_code size_ec;
          if (it->is_regular_file(size_ec)) { bytes += it->file_size(size_ec); }
        }
        ec.clear();
        std::filesystem::remove_all(path, ec);
        std::error_code exists_ec;
        if (ec || std::filesystem::exists(path, exists_ec))
        {
          owner->warn("search: could not remove the run's cache directory " + path.string() + " (" +
                      std::to_string(bytes / 1000000) + " MB" + (ec ? ": " + ec.message() : std::string()) +
                      "); remove it by hand");
        }
        path.clear();
      }
      ~Scratch()
      {
        try { remove(); } catch (...) {}
      }
    } scratch;
    scratch.owner = this;
    RunData run;
    PeakGroups groups;
    Calibration calibration;
    json run_json, calibration_json;

    // 1. The run, first: candidate selection needs its isolation windows.
    auto t = Clock::now();
    run = loadRun(run_path);
    scratch.path = run.scratch;
    run.meta.reset();   // SwathFile's per-spectrum metadata: nothing reads it (0.56 GB on a 6 GB run)
    if (run.name.empty()) { run.name = std::filesystem::path(run_path).stem().string(); }
    if (run.path.empty()) { run.path = run_path; }
    timing["load"] = since(t);
    const std::vector<IsolationWindow> windows = isolationWindows(run);
    {
      std::size_t ms1 = 0, ms2 = 0;
      for (const auto& m : run.maps) { (m.ms1 ? ms1 : ms2)++; }
      run_json = {{"name", run.name}, {"ion_mobility", run.ion_mobility}, {"read_mode", run.read_mode},
                  {"ms1_maps", ms1}, {"ms2_windows", ms2}, {"loader", object(run.provenance_json)}};
      if (!windows.empty())
      { run_json["isolation_mz_range"] = {windows.front().lower, std::max_element(windows.begin(), windows.end(),
                                           [](const IsolationWindow& a, const IsolationWindow& b) { return a.upper < b.upper; })->upper}; }
      info("search run: " + run.name + ", " + std::to_string(ms2) + " isolation windows, " + std::to_string(ms1) + " MS1 maps, " +
           (run.ion_mobility ? "diaPASEF (ion mobility)" : "no ion mobility") + ", read " +
           (run.read_mode.empty() ? std::string("?") : run.read_mode) + " (" + seconds(timing["load"].get<double>()) + ")");
    }

    // 2. Candidates and their in-memory decoys.
    t = Clock::now();
    const bool evidence = params_.candidates == "evidence";
    std::string prefilter_seconds;
    const SearchSet set = evidence ? EvidencePrefilter::select(library, params_, windows, run.maps,
                                                               [this](const std::string& m) { info(m); }, &prefilter_seconds)
                                   : CandidateSelector::select(library, params_, windows);
    timing["candidates"] = since(t);
    if (!prefilter_seconds.empty()) { timing["prefilter"] = json::parse(prefilter_seconds); }
    const auto& st = set.stats;
    info("search candidates: " + std::to_string(st.pairs) + " target-decoy pairs (" + searchDecoyName(params_.decoys) +
         " decoys, " + (evidence ? "chosen by fragment evidence" : "drawn with seed " + std::to_string(params_.seed)) +
         ") from " + std::to_string(st.eligible) + " eligible of " + std::to_string(st.targets) + " library targets; " +
         (evidence ? "kept by the prefilter " : "drawn ") + std::to_string(st.drawn) + ", no decoy " +
         std::to_string(st.no_decoy) + " (" + std::to_string(st.decoy_unshufflable) + " unshufflable, " +
         std::to_string(st.decoy_out_of_range) + " fragments out of range, " + std::to_string(st.decoy_copy) +
         " fragment copies, " + std::to_string(st.decoy_unparsable + st.decoy_too_few_fragments) + " other), capped " +
         std::to_string(st.capped) + "; ineligible: " +
         std::to_string(st.ineligible_mz) + " m/z, " + std::to_string(st.ineligible_charge) + " charge, " +
         std::to_string(st.ineligible_rt) + " RT, " + std::to_string(st.ineligible_fragments) + " < " +
         std::to_string(SearchParams::min_assay_fragments) + " fragments, " + std::to_string(st.ineligible_window) +
         " outside the " + std::to_string(st.windows) + " isolation windows, " + std::to_string(st.duplicate_key) +
         " duplicate keys; " + std::to_string(st.library_decoys_ignored) + " library decoys not searched (" +
         seconds(timing["candidates"].get<double>()) + ")");
    if (st.windows == 0) { warning("the run reported no isolation windows; candidates were not checked against them"); }
    if (evidence && params_.subset != 0)
    { warning("search:subset " + std::to_string(params_.subset) + " applies to search:candidates random only and was ignored"); }
    if (st.pairs == 0)
    {
      throw SearchAbort("search: no target-decoy pair to search (" + std::to_string(st.eligible) + " eligible targets, " +
                        std::to_string(st.ineligible_window) + " outside every isolation window, " +
                        std::to_string(st.no_decoy) + " without a decoy" +
                        (evidence ? ", none with fragment evidence at search:prefilter_depth " + std::to_string(params_.prefilter_depth)
                                  : std::string()) + ")");
    }

    // 3-4. Calibration and extraction.
    {
      t = Clock::now();
      calibration = calibrate(library, set, run);
      timing["calibrate"] = since(t);
      calibration_json = {{"seeds", calibration.seeds}, {"points", calibration.points}, {"rsq", num(calibration.rsq)},
                          {"points_per_seed", num(calibration.coverage)}, {"bootstrap", calibration.bootstrap},
                          {"rt_window_s", num(calibration.rt_window)}, {"mz_ppm", num(calibration.mz_ppm)},
                          {"ms1_mz_ppm", num(calibration.ms1_mz_ppm)}, {"im_window", num(calibration.im_window)},
                          {"detail", object(calibration.provenance_json)}};
      const json& cd = calibration_json["detail"];
      std::string model = cd.contains("model") ? cd["model"].get<std::string>() : std::string("linear");
      if (cd.contains("nonlinear") && cd["nonlinear"].contains("window_lowess_s") && cd["nonlinear"]["window_lowess_s"].is_number() &&
          cd["nonlinear"].contains("window_linear_s") && cd["nonlinear"]["window_linear_s"].is_number())
      {
        model += " (LOWESS window " + seconds(cd["nonlinear"]["window_lowess_s"].get<double>()) + ", linear " +
                 seconds(cd["nonlinear"]["window_linear_s"].get<double>()) + ")";
      }
      info("search calibration: " + std::to_string(calibration.points) + " points from " + std::to_string(calibration.seeds) +
           " seeds, " + model + ", r^2 " + (std::isfinite(calibration.rsq) ? std::to_string(calibration.rsq) : std::string("n/a")) +
           ", windows RT " + std::to_string(calibration.rt_window) + " s, m/z " + std::to_string(calibration.mz_ppm) +
           " ppm, 1/K0 " + (calibration.im_window > 0 ? std::to_string(calibration.im_window) : std::string("off")) +
           (calibration.bootstrap ? " (BOOTSTRAP)" : "") + " (" +
           seconds(timing["calibrate"].get<double>()) + ")");
      if (calibration.bootstrap) { warning("the RT calibration FAILED and search:allow_bootstrap replaced it with a linear map"); }

      t = Clock::now();
      extract(set, run, calibration, groups);
      timing["extract"] = since(t);
    }
    // The run's memory and its cache files are not needed for scoring.
    run.maps.clear();
    scratch.remove();

    groups.validate(set);
    checkExtraction(set, groups, params_);
    info("search extraction: " + std::to_string(groups.rows()) + " peak groups for " + std::to_string(set.size()) +
         " precursors, " + std::to_string(groups.scores.width()) + " sub-scores (" + seconds(timing["extract"].get<double>()) + ")");

    // 5. Classifier, competition, roll-ups, guards.
    t = Clock::now();
    const ScoringOutcome outcome = scorePeakGroups(set, groups, params_);
    timing["score"] = since(t);
    const auto& d = outcome.scored.diagnostics;
    for (const auto& w : outcome.warnings) { warning(w); }
    info("search scoring: " + std::to_string(d.targets_at_q) + " target precursors at q <= 0.01 (" +
         std::to_string(d.decoys_at_q) + " decoys; pooled estimate " + std::to_string(d.pooled_targets_at_q) + "), " +
         std::to_string(outcome.peptides_at_q) + " peptides, " + std::to_string(outcome.proteins_at_q) + " protein groups; " +
         std::to_string(d.features_used.size()) + " sub-scores, " + std::to_string(d.iterations_trained) +
         " iterations trained; null-pair balance " + std::to_string(outcome.null_targets) + " : " +
         std::to_string(outcome.null_decoys) + " (" + seconds(timing["score"].get<double>()) + ")");
    if (outcome.selftest)
    {
      info("search self-check: label swap " + std::to_string(outcome.selftest_label_swap_ids) + ", random pair labels " +
           std::to_string(outcome.selftest_random_label_ids) + " identifications (limit " + std::to_string(outcome.selftest_limit) + ")");
    }
    // Every guard but search:min_ids aborts here, before a report exists: a
    // search they stop is not one to keep. Too few identifications is an
    // honest result, and its report is written for inspection first.
    checkGuards(outcome, params_, false);
    const bool too_few = d.targets_at_q < params_.min_ids;

    // 6. The report.
    t = Clock::now();
    const std::vector<ReportRow> rows = reportRows(set, groups, outcome, params_);
    IdentificationResult result;
    result.report = out_ids;
    result.run_name = run.name;
    result.identified = d.targets_at_q;
    result.report_rows = rows.size();
    for (const auto& r : rows) { (r.decoy ? result.report_decoys : result.report_targets)++; }

    json candidates = selectionJson(st);
    candidates["method"] = params_.candidates;
    candidates["prefilter"] = set.prefilter_json.empty() ? json(nullptr) : json::parse(set.prefilter_json);
    json search = {
      {"experimental", true},
      {"settings", json::parse(params_.toJson(false))},
      {"candidates", candidates},
      {"run", run_json},
      {"calibration", calibration_json},
      {"scoring", scoringJson(outcome)},
      {"identifications", identificationsJson(outcome)},
      {"entrapment", entrapmentJson(outcome)},
      {"selftest", selftestJson(outcome)}};
    if (too_few)
    {
      search["aborted"] = {{"guard", "search:min_ids"}, {"min_ids", params_.min_ids}, {"identified", d.targets_at_q}};
    }
    // The report carries everything that is a function of the inputs and
    // settings, and nothing that is not (timings, threads, chunking, paths of
    // this machine): the same search at any -threads or search:chunk writes
    // the same bytes. The provenance sidecar gets the full settings.
    const json embedded = {{"tool", "DIALibGen"}, {"tool_version", tool_version},
                           {"producer", "DIALibGen built-in identification (-run); DIA-NN column names as a compatibility "
                                        "dialect, not a DIA-NN result"},
                           {"search", search}};
    ReportWriter::write(out_ids, run.name, rows, {{"odia.identifier", embedded.dump()}});
    if (too_few)
    {
      // The guard fires only now, with the report on disk for inspection. This
      // invocation never uses it as -ids; the same command is refused while
      // it exists.
      try { checkGuards(outcome, params_, true); }
      catch (const SearchAbort& e)
      {
        throw SearchAbort(std::string(e.what()) + ". The report " + out_ids + " was written for inspection (" +
                          std::to_string(result.report_targets) + " targets and " + std::to_string(result.report_decoys) +
                          " decoys at q <= " + fixed2(params_.report_max_q) + "); remove it before running the same command again");
      }
    }
    search["settings"] = json::parse(params_.toJson(true));
    timing["report"] = since(t);
    timing["total"] = since(started);
    info("search report: " + std::to_string(rows.size()) + " precursors (" + std::to_string(result.report_targets) + " targets, " +
         std::to_string(result.report_decoys) + " decoys) at q <= " + std::to_string(params_.report_max_q) + " written to " +
         out_ids + " (" + seconds(timing["total"].get<double>()) + " in total)");

    search["report"] = {{"path", std::filesystem::absolute(out_ids).string()}, {"rows", rows.size()},
                        {"targets", result.report_targets}, {"decoys", result.report_decoys}, {"max_q", params_.report_max_q}};
    search["warnings"] = warnings;
    search["resources"] = {{"threads", params_.threads},
                           {"openswath_threads", std::min(params_.threads, SearchParams::openswath_max_threads)},
                           {"seconds", timing}};
    result.provenance_json = search.dump();
    return result;
  }
}
