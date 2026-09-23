// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <odia/search/EvidencePrefilter.h>

#include <odia/search/PredictedAssays.h>
#include <odia/search/Scoring.h>
#include <odia/search/SearchDecoys.h>

#include <OpenMS/OPENSWATHALGO/DATAACCESS/DataStructures.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/ISpectrumAccess.h>

#include <nlohmann/json.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>

namespace ODIA::search
{
  namespace
  {
    using json = nlohmann::json;
    using Clock = std::chrono::steady_clock;
    constexpr std::size_t F = SearchParams::prefilter_fragments;
    static_assert(F <= 8, "a member's matched fragments are an 8-bit mask");

    double since(Clock::time_point t) { return std::chrono::duration<double>(Clock::now() - t).count(); }

    std::string fixed(double v, int digits)
    {
      std::ostringstream s;
      s.setf(std::ios::fixed);
      s.precision(digits);
      s << v;
      return s.str();
    }

    int bits(std::uint8_t v)
    {
      int n = 0;
      for (; v; v = static_cast<std::uint8_t>(v & (v - 1))) { ++n; }
      return n;
    }

    /// Uppercase residues outside (...) and [...]: "AC(UniMod:4)K" -> "ACK".
    std::string stripped(std::string_view modified)
    {
      std::string out;
      int depth = 0;
      for (const char c : modified)
      {
        if (c == '(' || c == '[') { ++depth; continue; }
        if (c == ')' || c == ']') { depth = std::max(0, depth - 1); continue; }
        if (depth == 0 && c >= 'A' && c <= 'Z') { out.push_back(c); }
      }
      return out;
    }

    /// A pair member's evidence as one comparable number: depth, then spectra.
    std::uint64_t strength(const MemberEvidence& e)
    { return (static_cast<std::uint64_t>(e.depth) << 32) | e.spectra; }

    /// One isolation window: the pairs whose precursor it holds, and every
    /// indexed fragment of both their members, sorted by m/z. An entry's code
    /// is (local member << 3) | fragment slot, local member = 2 * local pair
    /// + class.
    struct WindowIndex
    {
      double lower = 0.0, upper = 0.0;
      std::vector<std::uint32_t> pairs;
      std::vector<float> mz;
      std::vector<std::uint32_t> code;
    };

    /// One map's sweep result, per local member of its window.
    struct MapResult
    {
      std::size_t window = 0;
      std::vector<std::uint8_t> depth;
      std::vector<std::uint32_t> spectra;
      std::vector<float> rt;
      std::vector<float> intensity;
      std::size_t n_spectra = 0;
      std::size_t n_peaks = 0;
    };
  }

  // ---- the pair universe ----------------------------------------------------------

  PrefilterPairs EvidencePrefilter::pairs(const Library& library, const SearchParams& params,
                                          const std::vector<IsolationWindow>& windows, FragmentModel* model, const Log& progress)
  {
    params.validate();
    SelectionStats st;
    std::vector<std::pair<std::uint64_t, std::size_t>> eligible = CandidateSelector::eligible(library, params, windows, st);
    std::sort(eligible.begin(), eligible.end(),
              [&](const auto& a, const auto& b) { return CandidateSelector::drawLess(library, a, b); });
    return pairsOf(library, params, eligible, st, model, 0, progress);
  }

  PrefilterPairs EvidencePrefilter::pairsOf(const Library& library, const SearchParams& params,
                                            const std::vector<std::pair<std::uint64_t, std::size_t>>& ordered,
                                            const SelectionStats& stats, FragmentModel* model, std::size_t stop_after,
                                            const Log& progress)
  {
    const auto started = Clock::now();
    PrefilterPairs out;
    out.stats = stats;
    SelectionStats& st = out.stats;
    const DecoyRules rules = CandidateSelector::decoyRules(library, params);
    st.fragment_mz_min = fromFixed(rules.fragment_min);
    st.fragment_mz_max = fromFixed(rules.fragment_max);
    const bool predicted = params.intensities == Intensities::Predicted;
    if (predicted && !model)
    { throw std::invalid_argument("search:intensities predicted needs a fragment model (the PeptDeep MS2 model)"); }

    // Each target's decoy, exactly as the search will build it, and each
    // member's top F fragments:
    //   predicted: both members predicted by the model from their own
    //     sequences, each its own most intense (PredictedAssays.h);
    //   library: the top F slots by the TARGET's library intensity (ties to
    //     the earlier slot), and the decoy's m/z in the same slots.
    const std::size_t n = ordered.size();
    const auto& pre = library.precursors();
    const auto& tr = library.transitions();
    std::vector<DecoyOutcome> outcome(n, DecoyOutcome::Made);
    std::vector<std::uint8_t> count(n, 0), assay(n, 0);
    std::vector<float> mz(n * 2 * F, std::numeric_limits<float>::quiet_NaN());
    std::string failure;
    const std::size_t block = predicted ? PredictedAssays::block_pairs : std::max<std::size_t>(1, n);
    std::size_t processed = 0, made_so_far = 0;
    auto reported = Clock::now();
    for (std::size_t base = 0; base < n; base += block)
    {
      const std::size_t last = std::min(n, base + block);
      std::vector<std::string> sequences(predicted ? last - base : 0);
      const auto total = static_cast<std::ptrdiff_t>(last - base);
#pragma omp parallel for schedule(dynamic, 256)
      for (std::ptrdiff_t s = 0; s < total; ++s)
      {
        const std::size_t k = base + static_cast<std::size_t>(s);
        try
        {
          DecoyAssay a = searchDecoy(library, ordered[k].second, rules);
          outcome[k] = a.outcome;
          if (a.outcome != DecoyOutcome::Made) { continue; }
          if (predicted) { sequences[k - base] = std::move(a.sequence); continue; }
          std::vector<std::size_t> order(a.slots.size());
          std::iota(order.begin(), order.end(), std::size_t{0});
          std::stable_sort(order.begin(), order.end(), [&](std::size_t x, std::size_t y) {
            return tr.library_intensity[a.slots[x]] > tr.library_intensity[a.slots[y]];
          });
          const std::size_t m = std::min(F, order.size());
          for (std::size_t j = 0; j < m; ++j)
          {
            mz[(2 * k) * F + j] = static_cast<float>(fromFixed(tr.product_mz[a.slots[order[j]]]));
            mz[(2 * k + 1) * F + j] = static_cast<float>(fromFixed(a.mz[order[j]]));
          }
          count[k] = static_cast<std::uint8_t>(m);
          assay[k] = static_cast<std::uint8_t>(std::min<std::size_t>(255, a.slots.size()));
        }
        catch (const std::exception& e)
        {
#pragma omp critical(odia_prefilter_failure)
          { if (failure.empty()) { failure = e.what(); } }
        }
      }
      if (!failure.empty()) { throw std::runtime_error("search prefilter: " + failure); }
      if (predicted)
      {
        std::vector<std::size_t> targets;
        targets.reserve(last - base);
        for (std::size_t k = base; k < last; ++k) { targets.push_back(ordered[k].second); }
        std::vector<PairPrediction> predictions;
        PredictedAssays::predict(library, targets, sequences, rules, *model, predictions);
        for (std::size_t k = base; k < last; ++k)
        {
          if (sequences[k - base].empty()) { continue; }   // no decoy: its outcome stands
          const PairPrediction& p = predictions[k - base];
          outcome[k] = p.outcome;
          if (p.outcome != DecoyOutcome::Made) { continue; }
          const std::size_t m = std::min(F, p.count());
          for (std::size_t j = 0; j < m; ++j)
          {
            mz[(2 * k) * F + j] = static_cast<float>(p.target[j].mz);
            mz[(2 * k + 1) * F + j] = static_cast<float>(p.decoy[j].mz);
          }
          count[k] = static_cast<std::uint8_t>(m);
          assay[k] = static_cast<std::uint8_t>(std::min<std::size_t>(255, p.count()));
        }
      }
      processed = last;
      for (std::size_t k = base; k < last; ++k) { made_so_far += outcome[k] == DecoyOutcome::Made ? 1 : 0; }
      if (predicted && progress && (since(reported) >= 30.0 || last == n))
      {
        const double s = since(started);
        progress("search prefilter: " + std::to_string(last) + " of " + std::to_string(n) + " targets and their decoys predicted (" +
                 fixed(s, 0) + " s" + (last < n ? ", about " + fixed(s * static_cast<double>(n - last) / static_cast<double>(last), 0) +
                                                      " s to go" : std::string()) + ")");
        reported = Clock::now();
      }
      if (stop_after != 0 && made_so_far >= stop_after) { break; }
    }

    out.targets.reserve(made_so_far);
    out.fragments.reserve(made_so_far);
    out.assay_fragments.reserve(made_so_far);
    out.fragment_mz.reserve(made_so_far * 2 * F);
    out.precursor_mz.reserve(made_so_far);
    out.library_rt.reserve(made_so_far);
    out.charge.reserve(made_so_far);
    out.processed = processed;
    for (std::size_t k = 0; k < processed; ++k)
    {
      switch (outcome[k])
      {
        case DecoyOutcome::Made: break;
        case DecoyOutcome::Unparsable: ++st.decoy_unparsable; break;
        case DecoyOutcome::Unshufflable: ++st.decoy_unshufflable; break;
        case DecoyOutcome::OutOfRange: ++st.decoy_out_of_range; break;
        case DecoyOutcome::Copy: ++st.decoy_copy; break;
        case DecoyOutcome::TooFewFragments: ++st.decoy_too_few_fragments; break;
        case DecoyOutcome::Unpredictable: ++st.decoy_unpredictable; break;
      }
      if (outcome[k] != DecoyOutcome::Made) { ++st.no_decoy; continue; }
      const std::size_t i = ordered[k].second;
      out.targets.push_back(ordered[k]);
      out.fragments.push_back(count[k]);
      out.assay_fragments.push_back(assay[k]);
      out.fragment_mz.insert(out.fragment_mz.end(), mz.begin() + static_cast<std::ptrdiff_t>(2 * k * F),
                             mz.begin() + static_cast<std::ptrdiff_t>(2 * (k + 1) * F));
      out.precursor_mz.push_back(fromFixed(pre.mz[i]));
      out.library_rt.push_back(pre.irt[i]);
      out.charge.push_back(pre.charge[i]);
    }
    out.decoy_seconds = since(started);
    return out;
  }

  // ---- the sweep ---------------------------------------------------------------------

  std::vector<MemberEvidence> EvidencePrefilter::sweep(const PrefilterPairs& pairs, const std::vector<OpenSwath::SwathMap>& maps,
                                                       const SearchParams& params, SweepStats* stats, bool swap_classes)
  {
    params.validate();
    SweepStats local;
    SweepStats& ss = stats ? *stats : local;
    ss = SweepStats();
    const std::size_t P = pairs.size();
    std::vector<MemberEvidence> evidence(2 * P);

    // Distinct isolation windows of the MS2 maps (IM-split maps share one).
    std::vector<std::pair<double, double>> bounds;
    std::vector<std::size_t> ms2;
    for (std::size_t i = 0; i < maps.size(); ++i)
    {
      if (maps[i].ms1 || !maps[i].sptr) { continue; }
      ms2.push_back(i);
      bounds.emplace_back(maps[i].lower, maps[i].upper);
    }
    std::sort(bounds.begin(), bounds.end());
    bounds.erase(std::unique(bounds.begin(), bounds.end()), bounds.end());
    ss.windows = bounds.size();
    ss.maps = ms2.size();
    if (P == 0 || ms2.empty()) { return evidence; }

    // 1. One index per window, both classes, sorted by (m/z, code).
    auto t = Clock::now();
    std::vector<WindowIndex> index(bounds.size());
    const auto W = static_cast<std::ptrdiff_t>(bounds.size());
    std::string failure;
#pragma omp parallel for schedule(dynamic, 1)
    for (std::ptrdiff_t w = 0; w < W; ++w)
    {
      // An exception must not escape the parallel region (std::terminate): an
      // allocation failure on a very large library is reported like the
      // other loops' failures.
      try
      {
        WindowIndex& x = index[static_cast<std::size_t>(w)];
        x.lower = bounds[static_cast<std::size_t>(w)].first;
        x.upper = bounds[static_cast<std::size_t>(w)].second;
        std::size_t entries = 0;
        for (std::size_t k = 0; k < P; ++k)
        {
          const double mz = pairs.precursor_mz[k];
          if (x.lower < mz && mz < x.upper)
          {
            x.pairs.push_back(static_cast<std::uint32_t>(k));
            entries += 2u * pairs.fragments[k];
          }
        }
        if (2 * x.pairs.size() >= (std::size_t{1} << 29))
        {
#pragma omp critical(odia_prefilter_failure)
          { failure = "an isolation window holds " + std::to_string(x.pairs.size()) + " pairs, more than the index can address"; }
          continue;
        }
        std::vector<std::pair<float, std::uint32_t>> e;
        e.reserve(entries);
        for (std::size_t l = 0; l < x.pairs.size(); ++l)
        {
          const std::size_t k = x.pairs[l];
          for (std::uint32_t cls = 0; cls < 2; ++cls)
          {
            const std::uint32_t member = static_cast<std::uint32_t>(2 * l) + (swap_classes ? 1u - cls : cls);
            const float* f = pairs.mz(2 * k + cls);
            for (std::uint32_t j = 0; j < pairs.fragments[k]; ++j) { e.emplace_back(f[j], (member << 3) | j); }
          }
        }
        std::sort(e.begin(), e.end());
        x.mz.resize(e.size());
        x.code.resize(e.size());
        for (std::size_t j = 0; j < e.size(); ++j) { x.mz[j] = e[j].first; x.code[j] = e[j].second; }
      }
      catch (const std::exception& e)
      {
#pragma omp critical(odia_prefilter_failure)
        { if (failure.empty()) { failure = std::string("index: ") + e.what(); } }
      }
    }
    if (!failure.empty()) { throw std::runtime_error("search prefilter: " + failure); }
    for (const auto& x : index)
    {
      ss.entries += x.mz.size();
      ss.index_bytes += x.mz.size() * (sizeof(float) + sizeof(std::uint32_t)) + x.pairs.size() * sizeof(std::uint32_t);
    }
    ss.index_seconds = since(t);

    // 2. The sweep, one map at a time per thread, then merged in map order.
    t = Clock::now();
    const int depth_min = params.prefilter_depth;
    const std::size_t top = params.prefilter_top_peaks;
    const double tol = params.prefilter_ppm * 1e-6;
    std::vector<MapResult> results(ms2.size());
    const auto Q = static_cast<std::ptrdiff_t>(ms2.size());
#pragma omp parallel for schedule(dynamic, 1)
    for (std::ptrdiff_t q = 0; q < Q; ++q)
    {
      try
      {
        const OpenSwath::SwathMap& m = maps[ms2[static_cast<std::size_t>(q)]];
        MapResult& r = results[static_cast<std::size_t>(q)];
        r.window = static_cast<std::size_t>(std::lower_bound(bounds.begin(), bounds.end(), std::make_pair(m.lower, m.upper)) -
                                            bounds.begin());
        const WindowIndex& x = index[r.window];
        const std::size_t M = 2 * x.pairs.size();
        if (M == 0) { continue; }
        r.depth.assign(M, 0);
        r.spectra.assign(M, 0);
        r.rt.assign(M, std::numeric_limits<float>::quiet_NaN());
        r.intensity.assign(M, 0.0f);
        std::vector<std::uint32_t> stamp(M, std::numeric_limits<std::uint32_t>::max());
        std::vector<std::uint8_t> mask(M, 0);
        std::vector<float> summed(M, 0.0f);
        std::vector<std::uint32_t> touched;
        std::vector<std::uint32_t> order;
        std::vector<std::pair<double, double>> peaks;

        // A private handle: a cached map reads through one file stream.
        OpenSwath::SpectrumAccessPtr access = m.sptr->lightClone();
        if (!access) { access = m.sptr; }
        const std::size_t spectra = access->getNrSpectra();
        for (std::size_t s = 0; s < spectra; ++s)
        {
          const OpenSwath::SpectrumPtr spectrum = access->getSpectrumById(static_cast<int>(s));
          const double rt = access->getSpectrumMetaById(static_cast<int>(s)).RT;
          ++r.n_spectra;
          if (!spectrum || !spectrum->getMZArray() || !spectrum->getIntensityArray()) { continue; }
          const auto& mzs = spectrum->getMZArray()->data;
          const auto& ints = spectrum->getIntensityArray()->data;
          const std::size_t n = std::min(mzs.size(), ints.size());
          if (n == 0) { continue; }
          // The `top` most intense peaks (ties to the lower m/z, then the
          // earlier peak), matched in m/z order.
          order.resize(n);
          std::iota(order.begin(), order.end(), 0u);
          auto stronger = [&](std::uint32_t a, std::uint32_t b) {
            if (ints[a] != ints[b]) { return ints[a] > ints[b]; }
            if (mzs[a] != mzs[b]) { return mzs[a] < mzs[b]; }
            return a < b;
          };
          if (n > top)
          {
            std::nth_element(order.begin(), order.begin() + static_cast<std::ptrdiff_t>(top), order.end(), stronger);
            order.resize(top);
          }
          peaks.clear();
          for (const std::uint32_t p : order) { peaks.emplace_back(mzs[p], ints[p]); }
          std::sort(peaks.begin(), peaks.end());
          r.n_peaks += peaks.size();

          const auto stamp_now = static_cast<std::uint32_t>(s);
          auto from = x.mz.begin();
          for (const auto& [mz, intensity] : peaks)
          {
            const double lo = mz * (1.0 - tol), hi = mz * (1.0 + tol);
            from = std::lower_bound(from, x.mz.end(), lo, [](float a, double v) { return static_cast<double>(a) < v; });
            for (auto it = from; it != x.mz.end() && static_cast<double>(*it) <= hi; ++it)
            {
              const std::uint32_t code = x.code[static_cast<std::size_t>(it - x.mz.begin())];
              const std::uint32_t member = code >> 3;
              const auto bit = static_cast<std::uint8_t>(1u << (code & 7u));
              if (stamp[member] != stamp_now)
              {
                stamp[member] = stamp_now;
                mask[member] = 0;
                summed[member] = 0.0f;
                touched.push_back(member);
              }
              if (!(mask[member] & bit))
              {
                mask[member] = static_cast<std::uint8_t>(mask[member] | bit);
                summed[member] += static_cast<float>(intensity);
              }
            }
          }
          for (const std::uint32_t member : touched)
          {
            const int d = bits(mask[member]);
            if (d >= depth_min) { ++r.spectra[member]; }
            if (d > r.depth[member] || (d == r.depth[member] && summed[member] > r.intensity[member]))
            {
              r.depth[member] = static_cast<std::uint8_t>(d);
              r.intensity[member] = summed[member];
              r.rt[member] = static_cast<float>(rt);
            }
          }
          touched.clear();
        }
      }
      catch (const std::exception& e)
      {
#pragma omp critical(odia_prefilter_failure)
        { if (failure.empty()) { failure = std::string("sweep: ") + e.what(); } }
      }
    }
    if (!failure.empty()) { throw std::runtime_error("search prefilter: " + failure); }

    for (const MapResult& r : results)
    {
      ss.spectra += r.n_spectra;
      ss.peaks += r.n_peaks;
      if (r.depth.empty()) { continue; }
      const WindowIndex& x = index[r.window];
      for (std::size_t l = 0; l < r.depth.size(); ++l)
      {
        MemberEvidence& e = evidence[2 * static_cast<std::size_t>(x.pairs[l >> 1]) + (l & 1)];
        if (r.depth[l] > e.depth || (r.depth[l] == e.depth && r.intensity[l] > e.intensity))
        {
          e.depth = r.depth[l];
          e.rt = r.rt[l];
          e.intensity = r.intensity[l];
        }
        e.spectra += r.spectra[l];
      }
    }
    ss.sweep_seconds = since(t);
    return evidence;
  }

  // ---- the choice --------------------------------------------------------------------

  void EvidencePrefilter::checkRatio(std::size_t targets, std::size_t decoys)
  {
    const double ratio = targets > 0 ? static_cast<double>(decoys) / static_cast<double>(targets) : 0.0;
    if (targets == 0 && decoys == 0) { return; }
    if (!(ratio >= SearchParams::candidate_ratio_low && ratio <= SearchParams::candidate_ratio_high))
    {
      throw SearchAbort("search: the candidate set holds " + std::to_string(targets) + " targets and " + std::to_string(decoys) +
                        " decoys, a decoy:target ratio of " + fixed(ratio, 3) + " outside [" +
                        fixed(SearchParams::candidate_ratio_low, 2) + ", " + fixed(SearchParams::candidate_ratio_high, 2) +
                        "]: a selection that keeps members instead of whole pairs makes the FDR estimate meaningless");
    }
  }

  PrefilterSelection EvidencePrefilter::choose(const PrefilterPairs& pairs, const std::vector<MemberEvidence>& evidence,
                                               const SearchParams& params, const std::vector<IsolationWindow>& windows)
  {
    params.validate();
    const std::size_t P = pairs.size();
    if (evidence.size() != 2 * P)
    { throw std::invalid_argument("search prefilter: " + std::to_string(evidence.size()) + " members of evidence for " + std::to_string(P) + " pairs"); }
    const int D = params.prefilter_depth;
    PrefilterSelection sel;
    sel.depth_targets.assign(F + 1, 0);
    sel.depth_decoys.assign(F + 1, 0);
    std::vector<std::size_t> pass;
    for (std::size_t k = 0; k < P; ++k)
    {
      const int dt = evidence[2 * k].depth, dd = evidence[2 * k + 1].depth;
      ++sel.depth_targets[std::min<std::size_t>(F, static_cast<std::size_t>(dt))];
      ++sel.depth_decoys[std::min<std::size_t>(F, static_cast<std::size_t>(dd))];
      const bool t = dt >= D, d = dd >= D;
      sel.targets_passing += t;
      sel.decoys_passing += d;
      sel.both_passing += t && d;
      if (t || d) { pass.push_back(k); }   // pair-union: the label is never read
    }
    sel.union_pairs = pass.size();
    // Whole pairs only: as many decoys as targets, by construction.
    checkRatio(pass.size(), pass.size());

    const std::size_t cap = params.max_pairs;
    if (cap != 0 && pass.size() > cap)
    {
      // Strata: window x library-RT decile x charge, all read from the
      // target, which the decoy shares.
      std::vector<float> rts;
      rts.reserve(pass.size());
      for (const std::size_t k : pass) { rts.push_back(pairs.library_rt[k]); }
      std::sort(rts.begin(), rts.end());
      std::vector<float> edges;
      for (std::size_t j = 1; j < 10; ++j) { edges.push_back(rts[j * rts.size() / 10]); }
      auto stratum = [&](std::size_t k) {
        std::uint64_t window = 0;
        for (std::size_t w = 0; w < windows.size(); ++w)
        {
          if (windows[w].lower < pairs.precursor_mz[k] && pairs.precursor_mz[k] < windows[w].upper) { window = w; break; }
        }
        const auto decile = static_cast<std::uint64_t>(std::upper_bound(edges.begin(), edges.end(), pairs.library_rt[k]) - edges.begin());
        return (window * 10 + decile) * 16 + std::min<std::uint64_t>(pairs.charge[k], 15);
      };
      std::vector<std::uint64_t> key(P, 0), score(P, 0);
      for (const std::size_t k : pass)
      {
        key[k] = stratum(k);
        score[k] = std::max(strength(evidence[2 * k]), strength(evidence[2 * k + 1]));   // symmetric in the members
      }
      // Within a stratum: stronger evidence first, then the draw key (pairs
      // are in draw order, so the pair index is the draw rank).
      std::sort(pass.begin(), pass.end(), [&](std::size_t a, std::size_t b) {
        if (key[a] != key[b]) { return key[a] < key[b]; }
        if (score[a] != score[b]) { return score[a] > score[b]; }
        return a < b;
      });
      // Largest-remainder shares of the cap; remainder ties to the lower stratum key.
      struct Stratum { std::size_t begin, size, quota; std::uint64_t remainder; };
      std::vector<Stratum> strata;
      for (std::size_t a = 0; a < pass.size();)
      {
        std::size_t b = a + 1;
        while (b < pass.size() && key[pass[b]] == key[pass[a]]) { ++b; }
        strata.push_back({a, b - a, 0, 0});
        a = b;
      }
      const auto U = static_cast<std::uint64_t>(pass.size());
      std::uint64_t given = 0;
      for (auto& s : strata)
      {
        const std::uint64_t share = static_cast<std::uint64_t>(cap) * s.size;
        s.quota = static_cast<std::size_t>(share / U);
        s.remainder = share % U;
        given += s.quota;
      }
      std::vector<std::size_t> by_remainder(strata.size());
      std::iota(by_remainder.begin(), by_remainder.end(), std::size_t{0});
      std::stable_sort(by_remainder.begin(), by_remainder.end(),
                       [&](std::size_t a, std::size_t b) { return strata[a].remainder > strata[b].remainder; });
      for (std::size_t j = 0; given < cap && j < by_remainder.size(); ++j)
      {
        Stratum& s = strata[by_remainder[j]];
        if (s.quota < s.size) { ++s.quota; ++given; }
      }
      std::vector<std::size_t> kept;
      kept.reserve(cap);
      for (const auto& s : strata)
      {
        for (std::size_t j = 0; j < s.quota; ++j) { kept.push_back(pass[s.begin + j]); }
      }
      sel.strata = strata.size();
      sel.capped = pass.size() - kept.size();
      pass.swap(kept);
    }
    std::sort(pass.begin(), pass.end());
    for (const std::size_t k : pass)
    {
      sel.kept_targets_passing += evidence[2 * k].depth >= D;
      sel.kept_decoys_passing += evidence[2 * k + 1].depth >= D;
    }
    sel.kept = std::move(pass);
    return sel;
  }

  // ---- calibration seeds -------------------------------------------------------------

  std::vector<SeedHint> EvidencePrefilter::seeds(const Library& library, const PrefilterPairs& pairs,
                                                 const std::vector<MemberEvidence>& evidence, const SearchParams& params,
                                                 const RtScale& scale)
  {
    const int D = params.prefilter_depth;
    const std::size_t B = SearchParams::calibration_seed_bins;
    const std::size_t per_bin = std::max<std::size_t>(1, SearchParams::calibration_seeds / B);
    // Rank: deeper, then more spectra, then the draw (the pair index).
    auto better = [&](std::size_t a, std::size_t b) {
      const std::uint64_t sa = strength(evidence[2 * a]), sb = strength(evidence[2 * b]);
      if (sa != sb) { return sa > sb; }
      return a < b;
    };
    // One seed per peptide: its best precursor.
    std::unordered_map<std::string, std::size_t> best;
    const auto& pre = library.precursors();
    for (std::size_t k = 0; k < pairs.size(); ++k)
    {
      const MemberEvidence& t = evidence[2 * k];
      const MemberEvidence& d = evidence[2 * k + 1];
      if (t.depth < D || !std::isfinite(t.rt)) { continue; }
      // A seed must beat its own decoy clearly: a target that only matches as
      // well as a random rearrangement of itself is a chance match, and on
      // an entrapment library a sparse RT bin was otherwise filled with them.
      if (t.spectra < SearchParams::seed_min_spectra ||
          static_cast<std::uint64_t>(t.spectra) < std::uint64_t{SearchParams::seed_decoy_factor} * d.spectra)
      { continue; }
      // Only the central library-RT range: a library row with an absurd RT
      // must neither squeeze the bins nor stretch the seeds' RT span.
      if (!(pairs.library_rt[k] >= scale.min && pairs.library_rt[k] <= scale.max)) { continue; }
      const std::string peptide = stripped(library.strings().get(pre.modified_sequence[pairs.targets[k].second]));
      const auto it = best.find(peptide);
      if (it == best.end()) { best.emplace(peptide, k); }
      else if (better(k, it->second)) { it->second = k; }
    }
    std::vector<std::vector<std::size_t>> bin(B);
    for (const auto& [peptide, k] : best)
    {
      const double assay = scale.toAssay(pairs.library_rt[k]);
      const auto b = static_cast<std::size_t>(std::min<double>(static_cast<double>(B - 1), std::max(0.0, std::floor(assay / 100.0 * static_cast<double>(B)))));
      bin[b].push_back(k);
    }
    std::vector<std::size_t> chosen;
    for (auto& v : bin)
    {
      std::sort(v.begin(), v.end(), better);
      if (v.size() > per_bin) { v.resize(per_bin); }
      chosen.insert(chosen.end(), v.begin(), v.end());
    }
    std::sort(chosen.begin(), chosen.end(), better);
    std::vector<SeedHint> out;
    out.reserve(chosen.size());
    for (const std::size_t k : chosen)
    {
      const MemberEvidence& e = evidence[2 * k];
      out.push_back({pairs.targets[k].second, static_cast<double>(e.rt), e.depth, e.spectra});
    }
    return out;
  }

  // ---- everything --------------------------------------------------------------------

  namespace
  {
    /// The searched set's predicted assays: pair k takes the count its
    /// target got in @p universe.
    void predictSet(SearchSet& set, const Library& library, const SearchParams& params, const PrefilterPairs& universe,
                    FragmentModel& model)
    {
      std::unordered_map<std::size_t, std::uint8_t> count_of;
      count_of.reserve(universe.size());
      for (std::size_t k = 0; k < universe.size(); ++k) { count_of.emplace(universe.targets[k].second, universe.assay_fragments[k]); }
      std::vector<std::uint8_t> counts(set.pairs(), 0);
      for (std::size_t k = 0; k < set.pairs(); ++k)
      {
        const auto it = count_of.find(set.source[k]);
        if (it == count_of.end()) { throw std::logic_error("search: a searched pair is not in the prefilter's universe"); }
        counts[k] = it->second;
      }
      PredictedAssays::apply(set, library, CandidateSelector::decoyRules(library, params), model, counts);
      set.stats.fragment_slots_dropped = 0;   // a library-slot notion; predicted assays drop none
    }
  }

  SearchSet EvidencePrefilter::selectRandom(const Library& library, const SearchParams& params, const std::vector<IsolationWindow>& windows,
                                            FragmentModel& model)
  {
    params.validate();
    SelectionStats st;
    std::vector<std::pair<std::uint64_t, std::size_t>> draws = CandidateSelector::eligible(library, params, windows, st);
    auto less = [&](const auto& a, const auto& b) { return CandidateSelector::drawLess(library, a, b); };
    std::sort(draws.begin(), draws.end(), less);
    if (params.subset != 0 && params.subset < draws.size()) { draws.resize(params.subset); }
    st.drawn = draws.size();
    const std::size_t cap = params.max_pairs != 0 ? params.max_pairs : std::numeric_limits<std::size_t>::max();
    const PrefilterPairs universe = pairsOf(library, params, draws, st, &model, cap == std::numeric_limits<std::size_t>::max() ? 0 : cap);
    std::vector<std::pair<std::uint64_t, std::size_t>> chosen(universe.targets.begin(),
                                                              universe.targets.begin() + static_cast<std::ptrdiff_t>(std::min(cap, universe.size())));
    SearchSet set = CandidateSelector::fromTargets(library, params, chosen, universe.stats);
    set.stats.drawn = st.drawn;
    set.stats.capped = st.drawn - set.stats.no_decoy - set.stats.pairs;
    predictSet(set, library, params, universe, model);
    return set;
  }

  SearchSet EvidencePrefilter::select(const Library& library, const SearchParams& params, const std::vector<IsolationWindow>& windows,
                                      const std::vector<OpenSwath::SwathMap>& maps, const Log& info, std::string* seconds,
                                      FragmentModel* model)
  {
    auto say = [&](const std::string& m) { if (info) { info(m); } };
    const auto started = Clock::now();
    const bool predicted = params.intensities == Intensities::Predicted;
    const PrefilterPairs universe = pairs(library, params, windows, model, info);
    const SelectionStats& ps = universe.stats;
    say("search prefilter: " + std::to_string(universe.size()) + " target-decoy pairs indexed from " + std::to_string(ps.eligible) +
        " eligible targets (" + std::to_string(ps.no_decoy) + " without a decoy; " + std::to_string(ps.ineligible_window) +
        " library targets outside the isolation windows); " +
        (predicted ? "both members predicted by " + model->describe()
                   : std::string("library intensities, decoys in their target's slots")) +
        " (" + fixed(universe.decoy_seconds, 1) + " s)");

    SweepStats sw;
    const std::vector<MemberEvidence> evidence = sweep(universe, maps, params, &sw);
    say("search prefilter: " + std::to_string(sw.spectra) + " MS2 spectra of " + std::to_string(sw.maps) + " maps swept (" +
        std::to_string(sw.windows) + " isolation windows, top " + std::to_string(params.prefilter_top_peaks) + " peaks, +-" +
        fixed(params.prefilter_ppm, 1) + " ppm) against " + std::to_string(sw.entries) + " fragments (" +
        std::to_string(sw.index_bytes / 1000000) + " MB index; " + fixed(sw.index_seconds, 1) + " s index, " +
        fixed(sw.sweep_seconds, 1) + " s sweep)");

    auto t = Clock::now();
    const PrefilterSelection sel = choose(universe, evidence, params, windows);
    const double choose_s = since(t);
    say("search prefilter: at depth >= " + std::to_string(params.prefilter_depth) + " of " + std::to_string(F) + ": " +
        std::to_string(sel.targets_passing) + " targets, " + std::to_string(sel.decoys_passing) + " decoys (" +
        std::to_string(sel.both_passing) + " pairs both); " + std::to_string(sel.union_pairs) + " pairs kept by either member, " +
        std::to_string(sel.capped) + " capped by search:max_pairs " + std::to_string(params.max_pairs) + " over " +
        std::to_string(sel.strata) + " strata -> " + std::to_string(sel.kept.size()) + " pairs");

    std::vector<std::pair<std::uint64_t, std::size_t>> chosen;
    chosen.reserve(sel.kept.size());
    for (const std::size_t k : sel.kept) { chosen.push_back(universe.targets[k]); }
    t = Clock::now();
    SearchSet set = CandidateSelector::fromTargets(library, params, std::move(chosen), ps);
    if (predicted) { predictSet(set, library, params, universe, *model); }
    const double build_s = since(t);
    set.stats.capped = sel.capped;
    checkRatio(set.pairs(), set.size() - set.pairs());
    set.seeds = seeds(library, universe, evidence, params, set.rt_robust);
    if (!params.entrapment_tag.empty())
    {
      set.entrapment_db_universe = true;
      for (const auto& [draw, i] : universe.targets)
      {
        const Entrapment c = entrapmentClass(library.strings().get(library.precursors().protein_group[i]), params.entrapment_tag);
        if (c == Entrapment::Real) { ++set.entrapment_db_real; }
        else if (c == Entrapment::Trap) { ++set.entrapment_db_trap; }
      }
    }

    auto ratio = [](std::size_t a, std::size_t b) { return b > 0 ? json(static_cast<double>(a) / static_cast<double>(b)) : json(nullptr); };
    std::size_t seed_bins = 0, seeds_trap = 0;
    {
      std::vector<char> filled(SearchParams::calibration_seed_bins, 0);
      for (const auto& s : set.seeds)
      {
        if (!params.entrapment_tag.empty() &&
            entrapmentClass(library.strings().get(library.precursors().protein_group[s.index]), params.entrapment_tag) == Entrapment::Trap)
        { ++seeds_trap; }
        const double a = set.rt_robust.toAssay(library.precursors().irt[s.index]);
        filled[static_cast<std::size_t>(std::min<double>(SearchParams::calibration_seed_bins - 1,
                                                         std::max(0.0, std::floor(a / 100.0 * SearchParams::calibration_seed_bins))))] = 1;
      }
      seed_bins = static_cast<std::size_t>(std::count(filled.begin(), filled.end(), 1));
    }
    const json record = {
      {"rule", "a pair is kept when its target or its decoy has depth distinct fragments of its top fragments among one "
               "spectrum's top peaks; the same rule for both classes"},
      {"depth", params.prefilter_depth}, {"fragments", F}, {"top_peaks", params.prefilter_top_peaks},
      {"ppm", params.prefilter_ppm},
      {"pairs_indexed", universe.size()},
      {"isolation_windows", sw.windows}, {"ms2_maps", sw.maps}, {"spectra", sw.spectra}, {"peaks_matched", sw.peaks},
      {"index_entries", sw.entries},
      {"members_by_depth", {{"targets", sel.depth_targets}, {"decoys", sel.depth_decoys}}},
      {"passing", {{"targets", sel.targets_passing}, {"decoys", sel.decoys_passing}, {"both", sel.both_passing},
                   {"pairs_union", sel.union_pairs}, {"target_decoy_ratio", ratio(sel.targets_passing, sel.decoys_passing)}}},
      {"cap", {{"max_pairs", params.max_pairs}, {"strata", sel.strata}, {"capped", sel.capped}, {"kept", sel.kept.size()},
               {"rank", "window x library-RT decile x charge; within a stratum the better member's (depth, spectra), then the draw key"},
               {"kept_passing", {{"targets", sel.kept_targets_passing}, {"decoys", sel.kept_decoys_passing},
                                 {"target_decoy_ratio", ratio(sel.kept_targets_passing, sel.kept_decoys_passing)}}}}},
      {"candidate_ratio_band", {SearchParams::candidate_ratio_low, SearchParams::candidate_ratio_high}},
      {"intensities", {{"source", toString(params.intensities)}, {"model", predicted ? json(model->describe()) : json(nullptr)},
                       {"rule", predicted ? "both members predicted from their own sequences by one model; each takes its own most "
                                            "intense b/y fragments, the same number for both"
                                          : "the target's library assay; the decoy in its target's slots with its target's intensities"}}},
      {"seeds", {{"count", set.seeds.size()}, {"bins", SearchParams::calibration_seed_bins}, {"bins_filled", seed_bins},
                 {"library_rt_range", {set.rt_robust.min, set.rt_robust.max}},
                 {"library_rt_full", {set.rt_scale.min, set.rt_scale.max}},
                 {"entrapment", params.entrapment_tag.empty() ? json(nullptr) : json(seeds_trap)},
                 {"rule", "targets whose own depth passes, with >= " + std::to_string(SearchParams::seed_min_spectra) +
                          " spectra and >= " + std::to_string(SearchParams::seed_decoy_factor) +
                          "x their decoy's; one per peptide, best (depth, spectra) per bin of the central library-RT range"}}}};
    set.prefilter_json = record.dump();
    if (seconds)
    {
      *seconds = json({{"decoys", universe.decoy_seconds}, {"index", sw.index_seconds}, {"sweep", sw.sweep_seconds},
                       {"choose", choose_s}, {"build", build_s}, {"total", since(started)}, {"index_bytes", sw.index_bytes}})
                   .dump();
    }
    say("search prefilter: " + std::to_string(set.pairs()) + " pairs searched, " + std::to_string(set.seeds.size()) +
        " calibration seeds over " + std::to_string(seed_bins) + " of " + std::to_string(SearchParams::calibration_seed_bins) +
        " library-RT bins (" + fixed(since(started), 1) + " s in all)");
    return set;
  }
}
