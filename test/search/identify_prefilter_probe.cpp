// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

// A measurement tool, not a test: the evidence prefilter of search:candidates
// evidence on a real library and run, at several settings, without extraction.
// It reads the run once (as -run does), builds the pair universe once, and
// sweeps once per (top peaks, ppm); every depth is evaluated on that sweep.
//
//   identify_prefilter_probe <library> <run.mzML> <out prefix> [key=value ...]
//     threads=16  readoptions=auto|normal|cache  cache_dir=<dir>  decoys=shuffle
//     top=1000[,...]  ppm=10[,...]  depth=3,4,5,6  max_pairs=200000
//     intensities=predicted|library  model=<peptdeep_ms2_dynamic.onnx>  instrument=<name>  nce=<nce>
//     (predicted needs model; instrument empty = timsTOF on an ion-mobility run, QE otherwise)
//
// Writes <prefix>.json (counts, timings, memory) and, per sweep,
// <prefix>.top<N>_ppm<P>.tsv: every pair whose better member reached depth 3,
// with both members' evidence and, per depth, whether the capped selection
// keeps it ("kept_d<D>"). Precursor key: modified sequence and charge.

#include <odia/DIANNLibraryFile.h>
#include <odia/search/EvidencePrefilter.h>
#include <odia/search/Identifier.h>

#include <nlohmann/json.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

using namespace ODIA::search;
using json = nlohmann::json;

namespace
{
  class Loader : public Identifier
  {
  public:
    explicit Loader(SearchParams p)
      : Identifier(std::move(p), [](const std::string& m) { std::cout << m << std::endl; },
                   [](const std::string& m) { std::cout << "warning: " << m << std::endl; }) {}
    using Identifier::loadRun;
  };

  double since(std::chrono::steady_clock::time_point t)
  { return std::chrono::duration<double>(std::chrono::steady_clock::now() - t).count(); }

  /// Peak resident memory in MB where the OS tells (Linux /proc), else -1.
  double peakMB()
  {
    std::ifstream in("/proc/self/status");
    std::string line;
    while (std::getline(in, line))
    {
      if (line.rfind("VmHWM:", 0) == 0) { return std::stod(line.substr(6)) / 1024.0; }
    }
    return -1.0;
  }

  template <typename T>
  std::vector<T> list(const std::string& s)
  {
    std::vector<T> out;
    std::stringstream in(s);
    std::string item;
    while (std::getline(in, item, ',')) { std::stringstream v(item); T x{}; v >> x; out.push_back(x); }
    return out;
  }
}

int main(int argc, char** argv)
{
  if (argc < 4)
  {
    std::cerr << "usage: identify_prefilter_probe <library> <run.mzML> <out prefix> [key=value ...]\n";
    return 2;
  }
  std::map<std::string, std::string> opt = {{"threads", "16"}, {"readoptions", "auto"}, {"cache_dir", ""},
                                            {"decoys", "shuffle"}, {"top", "1000"}, {"ppm", "10"},
                                            {"depth", "3,4,5,6"}, {"max_pairs", "200000"}, {"intensities", "predicted"},
                                            {"model", ""}, {"instrument", ""}, {"nce", "-1"}};
  for (int a = 4; a < argc; ++a)
  {
    const std::string kv = argv[a];
    const auto eq = kv.find('=');
    if (eq == std::string::npos || !opt.count(kv.substr(0, eq))) { std::cerr << "bad option " << kv << "\n"; return 2; }
    opt[kv.substr(0, eq)] = kv.substr(eq + 1);
  }
  const std::string prefix = argv[3];
  SearchParams params;
  params.threads = std::stoi(opt["threads"]);
  params.readoptions = parseReadMode(opt["readoptions"]);
  params.cache_dir = opt["cache_dir"];
  params.decoys = parseSearchDecoyMethod(opt["decoys"]);
  params.max_pairs = std::stoul(opt["max_pairs"]);
  params.intensities = parseIntensities(opt["intensities"]);
#ifdef _OPENMP
  omp_set_num_threads(params.threads);
#endif
  json record = {{"library", argv[1]}, {"run", argv[2]}, {"options", opt}};

  auto t = std::chrono::steady_clock::now();
  ODIA::Library library;
  ODIA::DIANNLibraryFile::load(argv[1], library);
  record["library_precursors"] = library.precursorCount();
  record["seconds_library"] = since(t);
  std::cout << "library: " << library.precursorCount() << " precursors (" << since(t) << " s)" << std::endl;

  t = std::chrono::steady_clock::now();
  Loader loader(params);
  RunData run = loader.loadRun(argv[2]);
  record["seconds_run"] = since(t);
  const std::vector<IsolationWindow> windows = Identifier::isolationWindows(run);
  record["windows"] = windows.size();
  record["mb_after_load"] = peakMB();

  std::unique_ptr<PeptDeepFragmentModel> model;
  if (params.intensities == Intensities::Predicted)
  {
    if (opt["model"].empty()) { std::cerr << "intensities=predicted needs model=<peptdeep_ms2_dynamic.onnx>\n"; return 2; }
    const std::string instrument = !opt["instrument"].empty() ? opt["instrument"] : (run.ion_mobility ? "timsTOF" : "QE");
    model = std::make_unique<PeptDeepFragmentModel>(opt["model"], instrument, std::stod(opt["nce"]), params.threads);
    record["model"] = model->describe();
    std::cout << "model: " << model->describe() << std::endl;
  }
  t = std::chrono::steady_clock::now();
  const PrefilterPairs universe = EvidencePrefilter::pairs(library, params, windows, model.get());
  record["seconds_pairs"] = since(t);
  record["pairs"] = universe.size();
  record["eligible"] = universe.stats.eligible;
  record["no_decoy"] = universe.stats.no_decoy;
  record["ineligible_window"] = universe.stats.ineligible_window;
  record["mb_after_pairs"] = peakMB();
  std::cout << "pairs: " << universe.size() << " of " << universe.stats.eligible << " eligible (" << since(t) << " s)" << std::endl;

  const auto& pre = library.precursors();
  json sweeps = json::array();
  for (const std::size_t top : list<std::size_t>(opt["top"]))
  {
    for (const double ppm : list<double>(opt["ppm"]))
    {
      SearchParams p = params;
      p.prefilter_top_peaks = top;
      p.prefilter_ppm = ppm;
      const std::vector<int> depths = list<int>(opt["depth"]);
      // Spectra counts are for the configured depth: sweep at the smallest.
      p.prefilter_depth = *std::min_element(depths.begin(), depths.end());
      SweepStats sw;
      t = std::chrono::steady_clock::now();
      const std::vector<MemberEvidence> evidence = EvidencePrefilter::sweep(universe, run.maps, p, &sw);
      json s = {{"top_peaks", top}, {"ppm", ppm}, {"spectra", sw.spectra}, {"maps", sw.maps}, {"entries", sw.entries},
                {"index_mb", sw.index_bytes / 1e6}, {"seconds_index", sw.index_seconds}, {"seconds_sweep", sw.sweep_seconds},
                {"peaks", sw.peaks}, {"mb_peak", peakMB()}, {"spectra_counted_at_depth", p.prefilter_depth}};
      std::vector<std::vector<char>> kept(depths.size(), std::vector<char>(universe.size(), 0));
      json per_depth = json::array();
      for (std::size_t j = 0; j < depths.size(); ++j)
      {
        SearchParams q = p;
        q.prefilter_depth = depths[j];
        const auto c = std::chrono::steady_clock::now();
        const PrefilterSelection sel = EvidencePrefilter::choose(universe, evidence, q, windows);
        for (const std::size_t k : sel.kept) { kept[j][k] = 1; }
        per_depth.push_back({{"depth", depths[j]}, {"targets_passing", sel.targets_passing}, {"decoys_passing", sel.decoys_passing},
                             {"both_passing", sel.both_passing}, {"union_pairs", sel.union_pairs}, {"capped", sel.capped},
                             {"kept", sel.kept.size()}, {"kept_targets_passing", sel.kept_targets_passing},
                             {"kept_decoys_passing", sel.kept_decoys_passing}, {"strata", sel.strata},
                             {"depth_targets", sel.depth_targets}, {"depth_decoys", sel.depth_decoys},
                             {"seeds", EvidencePrefilter::seeds(library, universe, evidence, q, CandidateSelector::robustRtRange(library)).size()},
                             {"seconds_choose", since(c)}});
      }
      s["depths"] = per_depth;
      std::ostringstream name;
      name << prefix << ".top" << top << "_ppm" << ppm << ".tsv";
      std::ofstream out(name.str());
      out << "Modified.Sequence\tPrecursor.Charge\tlibrary_rt\ttarget_depth\tdecoy_depth\ttarget_spectra\tdecoy_spectra\ttarget_rt\tdecoy_rt";
      for (const int d : depths) { out << "\tkept_d" << d; }
      out << "\n";
      std::size_t rows = 0;
      for (std::size_t k = 0; k < universe.size(); ++k)
      {
        const MemberEvidence& a = evidence[2 * k];
        const MemberEvidence& b = evidence[2 * k + 1];
        if (std::max(a.depth, b.depth) < 3) { continue; }
        const std::size_t i = universe.targets[k].second;
        out << library.strings().get(pre.modified_sequence[i]) << '\t' << static_cast<int>(pre.charge[i]) << '\t'
            << pre.irt[i] << '\t' << static_cast<int>(a.depth) << '\t' << static_cast<int>(b.depth) << '\t' << a.spectra << '\t'
            << b.spectra << '\t' << a.rt << '\t' << b.rt;
        for (std::size_t j = 0; j < depths.size(); ++j) { out << '\t' << static_cast<int>(kept[j][k]); }
        out << '\n';
        ++rows;
      }
      s["tsv"] = name.str();
      s["tsv_rows"] = rows;
      sweeps.push_back(s);
      std::cout << "sweep top " << top << " ppm " << ppm << ": " << sw.sweep_seconds << " s, " << per_depth.dump() << std::endl;
      std::ofstream(prefix + ".json") << json(record).dump(1) << "\n";
    }
  }
  record["sweeps"] = sweeps;
  record["mb_peak"] = peakMB();
  std::ofstream(prefix + ".json") << record.dump(1) << "\n";
  run.maps.clear();
  if (!run.scratch.empty()) { std::error_code ec; std::filesystem::remove_all(run.scratch, ec); }
  std::cout << "done" << std::endl;
  return 0;
}
