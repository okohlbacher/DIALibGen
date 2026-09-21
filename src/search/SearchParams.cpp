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

  DecoyMethod parseSearchDecoyMethod(const std::string& s)
  {
    // Spelled out rather than routed through parseDecoyMethod, which maps every
    // unknown name to Mutate -- the one method the search must never use.
    if (s == "shuffle") { return DecoyMethod::Shuffle; }
    if (s == "pseudo_reverse") { return DecoyMethod::PseudoReverse; }
    if (s == "reverse") { return DecoyMethod::Reverse; }
    if (s == "mutate")
    { throw std::invalid_argument("search:decoys mutate is not offered: its substitution table is DIA-NN's; use shuffle, pseudo_reverse or reverse"); }
    throw std::invalid_argument("search:decoys must be shuffle, pseudo_reverse or reverse, not '" + s + "'");
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
    if (candidates != "random") { fail("search:candidates must be random in this version, not '" + candidates + "'"); }
    if (decoys != DecoyMethod::Shuffle && decoys != DecoyMethod::PseudoReverse && decoys != DecoyMethod::Reverse)
    { fail(std::string("search:decoys ") + searchDecoyName(decoys) + " is not offered; use shuffle, pseudo_reverse or reverse"); }
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
  }

  std::string SearchParams::toJson() const
  {
    const nlohmann::json j = {
      {"candidates", candidates},
      {"subset", subset},
      {"max_pairs", max_pairs},
      {"decoys", searchDecoyName(decoys)},
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
    return j.dump();
  }
}
