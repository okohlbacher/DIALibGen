// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <odia/search/SearchParams.h>

#include <nlohmann/json.hpp>

#include <cmath>
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
      for (const char* key : {"chunk", "batch_size", "readoptions", "cache_dir", "ms2_model"}) { j.erase(key); }
    }
    return j.dump();
  }
}
