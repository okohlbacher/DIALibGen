// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <odia/search/SearchParams.h>

#include <odia/Hardware.h>

#include <OpenMS/config.h>

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include <string>

namespace ODIA::search
{
  const char* toString(ReadMode m)
  {
    switch (m)
    {
      case ReadMode::Normal: return "normal";
      case ReadMode::Cache: return "cache";
      case ReadMode::Auto: break;
    }
    return "auto";
  }

  ReadMode parseReadMode(const std::string& s)
  {
    if (s == "auto") { return ReadMode::Auto; }
    if (s == "normal") { return ReadMode::Normal; }
    if (s == "cache") { return ReadMode::Cache; }
    throw std::invalid_argument("search:readoptions must be auto, normal or cache, not '" + s + "'");
  }

  const char* toString(Ms1InMemory m)
  {
    switch (m)
    {
      case Ms1InMemory::On: return "on";
      case Ms1InMemory::Off: return "off";
      case Ms1InMemory::Auto: break;
    }
    return "auto";
  }

  Ms1InMemory parseMs1InMemory(const std::string& s)
  {
    if (s == "auto") { return Ms1InMemory::Auto; }
    if (s == "on") { return Ms1InMemory::On; }
    if (s == "off") { return Ms1InMemory::Off; }
    throw std::invalid_argument("search:ms1_in_memory must be auto, on or off, not '" + s + "'");
  }

  const char* toString(Intensities i)
  {
    return i == Intensities::Library ? "library" : "predicted";
  }

  Intensities parseIntensities(const std::string& s)
  {
    if (s == "predicted") { return Intensities::Predicted; }
    if (s == "library") { return Intensities::Library; }
    throw std::invalid_argument("search:intensities must be predicted or library, not '" + s + "'");
  }

  DecoyMethod parseSearchDecoyMethod(const std::string& s)
  {
    // Spelled out rather than routed through parseDecoyMethod, which maps every
    // unknown name to Mutate -- the one method the search must never use.
    if (s == "shuffle") { return DecoyMethod::Shuffle; }
    if (s == "pseudo_reverse") { return DecoyMethod::PseudoReverse; }
    if (s == "mutate")
    { throw std::invalid_argument("search:decoys mutate is not offered: its substitution table is DIA-NN's; use shuffle or pseudo_reverse"); }
    // Whole-sequence reversal moves the C-terminal K/R of a tryptic peptide to
    // the N-terminus: the classifier then learns how decoys were made rather
    // than whether a peptide is present (measured on an Astral run: null pairs
    // won 1.66 : 1 by the target instead of 1.04 : 1 with shuffle).
    if (s == "reverse")
    { throw std::invalid_argument("search:decoys reverse is not offered: reversing the termini makes decoys separable from "
                                  "tryptic targets by construction; use shuffle or pseudo_reverse, which keep both termini"); }
    throw std::invalid_argument("search:decoys must be shuffle or pseudo_reverse, not '" + s + "'");
  }

  const char* searchDecoyName(DecoyMethod m)
  {
    switch (m)
    {
      case DecoyMethod::Shuffle: return "shuffle";
      case DecoyMethod::PseudoReverse: return "pseudo_reverse";
      case DecoyMethod::Reverse: return "reverse";
      case DecoyMethod::Mutate: return "mutate";
      case DecoyMethod::None: break;
    }
    return "none";
  }

  int SearchParams::resolveThreads(int requested)
  {
    return requested > 0 ? requested : static_cast<int>(std::max(1u, physicalCores()));
  }

  namespace
  {
    /// An integer environment switch in [lo, hi], or @p fallback when it is
    /// unset, empty, malformed or out of range.
    int envInt(const char* name, int fallback, long lo, long hi)
    {
      const char* e = std::getenv(name);
      if (!e || !*e) { return fallback; }
      char* end = nullptr;
      const long v = std::strtol(e, &end, 10);
      return (end && *end == '\0' && v >= lo && v <= hi) ? static_cast<int>(v) : fallback;
    }

    /// ODIA_OPENSWATH_OUTER, read once: -2 unset, else its value (<= 0: off).
    int outerSwitch()
    {
      static const int v = envInt("ODIA_OPENSWATH_OUTER", -2, -1, 4096);
      return v;
    }
  }

  bool SearchParams::openswathNestedLoop()
  {
#if defined(_OPENMP) && defined(MT_ENABLE_NESTED_OPENMP)
    return true;
#else
    return false;
#endif
  }

  int SearchParams::openswathThreads()
  {
    static const int cap = [] {
      const int nested = openswathNestedLoop() && outerSwitch() != 0 && outerSwitch() != -1 ? openswath_max_threads
                                                                                             : openswath_max_threads_flat;
      return envInt("ODIA_OPENSWATH_THREADS", nested, 1, 4096);
    }();
    return cap;
  }

  int SearchParams::openswathOuterLoop(int team)
  {
    if (!openswathNestedLoop() || team < 1) { return -1; }
    const int k = outerSwitch();
    if (k == 0 || k == -1) { return -1; }
    if (k >= 1) { return std::min(k, team); }
    // From a team of 8: at 4 the flat loop keeps all four threads busy and
    // two windows at once were slower (extraction 840 -> 946 s on the timsTOF
    // cell); at 8 the nested loop is level (512 -> 495 s), at 16 it wins.
    return team >= 8 ? team / 2 : -1;
  }

  int SearchParams::calibrationThreads(int team)
  {
    static const int k = envInt("ODIA_CALIB_PAR", -1, 0, 4096);
    const int wanted = k >= 0 ? k : calibration_max_threads;
    return std::max(1, std::min(wanted, team));
  }

  int SearchParams::loadThreads()
  {
    static const int k = envInt("ODIA_LOAD_THREADS", load_max_threads, 0, 4096);
    return k;
  }

  bool SearchParams::asyncIo()
  {
    static const bool on = [] {
      const char* e = std::getenv("ODIA_ASYNC_IO");
      return !(e && std::string(e) == "0");
    }();
    return on;
  }

  int SearchParams::predictionSessions(int threads)
  {
    // One per thread, never more than the physical cores, as generate runs
    // them since 0.11.2 (LibraryGenerator inferenceSessions): with 64-row
    // batches (PeptDeepPredictor::CPU_CHUNK_ROWS) a session's MS2 arena is
    // about 53 MB, so the cap of 16 sessions that kept 2048-row arenas in
    // memory is gone; sessions on SMT siblings share a core and are slower.
    return std::clamp(threads, 1, static_cast<int>(std::max(1u, physicalCores())));
  }

  void SearchParams::validate() const
  {
    auto fail = [](const std::string& what) { throw std::invalid_argument(what); };
    auto finite = [](double v) { return std::isfinite(v); };
    if (candidates != "evidence" && candidates != "random")
    { fail("search:candidates must be evidence or random, not '" + candidates + "'"); }
    if (prefilter_depth < 1 || prefilter_depth > static_cast<int>(prefilter_fragments))
    { fail("search:prefilter_depth must be in [1, " + std::to_string(prefilter_fragments) + "]"); }
    if (prefilter_top_peaks == 0) { fail("search:prefilter_top_peaks must be at least 1"); }
    if (!finite(prefilter_ppm) || prefilter_ppm <= 0 || prefilter_ppm > 100)
    { fail("search:prefilter_ppm must be in (0, 100]"); }
    if (decoys != DecoyMethod::Shuffle && decoys != DecoyMethod::PseudoReverse)
    { fail(std::string("search:decoys ") + searchDecoyName(decoys) + " is not offered; use shuffle or pseudo_reverse"); }
    if (passes != 1) { fail("search:passes must be 1 in this version, not " + std::to_string(passes)); }
    if (!finite(rt_window) || rt_window < 0) { fail("search:rt_window must be >= 0 seconds (0 = from the calibration)"); }
    if (!finite(mz_ppm) || mz_ppm < 0) { fail("search:mz_ppm must be >= 0 (0 = automatic)"); }
    if (!finite(im_window) || (im_window < 0 && im_window != -1.0)) { fail("search:im_window must be >= 0, or -1 for off"); }
    if (im_window > 1.0) { fail("search:im_window is a full 1/K0 width; " + std::to_string(im_window) + " exceeds the whole mobility range"); }
    if (batch_size == 0) { fail("search:batch_size must be at least 1"); }
    if (chunk < 2) { fail("search:chunk must be at least 2 (one target-decoy pair)"); }
    if (!finite(calibration_min_rsq) || calibration_min_rsq < 0 || calibration_min_rsq > 1)
    { fail("search:calibration_min_rsq must be in [0, 1]"); }
    if (!finite(calibration_min_coverage) || calibration_min_coverage < 0 || calibration_min_coverage > 1)
    { fail("search:calibration_min_coverage must be in [0, 1]"); }
    if (!finite(max_target_fraction) || max_target_fraction <= 0 || max_target_fraction > 1)
    { fail("search:max_target_fraction must be in (0, 1]"); }
    if (!finite(report_max_q) || report_max_q < identification_q || report_max_q > 1)
    { fail("search:report_max_q must be in [0.01, 1]: the report must carry every identification the gates can pass"); }
    if (threads < 1) { fail("threads must be at least 1"); }
    if (!finite(nce) || nce > 100) { fail("search:nce must be <= 100 (<= 0: the instrument's default)"); }
  }

  std::string SearchParams::toJson(bool execution) const
  {
    nlohmann::json j = {
      {"candidates", candidates},
      {"subset", subset},
      {"max_pairs", max_pairs},
      {"prefilter_depth", prefilter_depth},
      {"prefilter_top_peaks", prefilter_top_peaks},
      {"prefilter_ppm", prefilter_ppm},
      {"prefilter_fragments", prefilter_fragments},
      {"calibration_seeds", calibration_seeds},
      {"decoys", searchDecoyName(decoys)},
      {"intensities", toString(intensities)},
      {"ms2_model", ms2_model},
      {"ms2_model_hash", ms2_model_hash},
      {"instrument", instrument},
      {"nce", nce},
      {"instrument_source", instrument_source},
      {"seed", seed},
      {"passes", passes},
      {"rt_window", rt_window},
      {"mz_ppm", mz_ppm},
      {"im_window", im_window},
      {"ms1", ms1},
      {"rt_im_scores", rt_im_scores},
      {"batch_size", batch_size},
      {"chunk", chunk},
      {"calibration_min_rsq", calibration_min_rsq},
      {"calibration_min_coverage", calibration_min_coverage},
      {"allow_bootstrap", allow_bootstrap},
      {"readoptions", toString(readoptions)},
      {"cache_dir", cache_dir},
      {"ms1_in_memory", toString(ms1_in_memory)},
      {"min_ids", min_ids},
      {"max_target_fraction", max_target_fraction},
      {"report_max_q", report_max_q},
      {"entrapment_tag", entrapment_tag},
      {"selftest", selftest},
      // threads deliberately absent: it must not change the result, and the
      // report embeds these settings, so the report must not depend on it.
      {"identification_q", identification_q},
      {"decoy_ratio_band", {decoy_ratio_low, decoy_ratio_high}},
      {"min_assay_fragments", min_assay_fragments}};
    if (!execution)
    {
      for (const char* key : {"chunk", "batch_size", "readoptions", "cache_dir", "ms1_in_memory", "ms2_model"}) { j.erase(key); }
    }
    return j.dump();
  }
}
