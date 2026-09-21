// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause
//
// scoreAndControl, the pipeline's entry point: honesty and determinism.
//
//   S1 planted signal: most true hits identified, few null targets among them, diagnostics add up
//   S2 pure null: exactly 0 identifications at q <= 0.01, over several seeds
//   S3 label swap (decoys carry the signal): about 0 identifications
//   S4 row permutation: bit-identical result
//   S5 thread count: bit-identical result at 1, 2, 3 and 8 threads
//   S6 all-missing and constant columns are dropped, label-blind, with no effect on the result
//   S7 missing cells are imputed and counted
//   S8 input that cannot be scored honestly is refused
//   S9 peptide (pooled) and protein (picked) q-values from the precursor result

#include "synthetic.h"

#include <odia-core/odia_core.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

using synth::check;
using odia::core::ScoreTable;
using odia::core::ScoredResult;

namespace
{
ScoreTable toTable(const synth::Data& d)
{
  ScoreTable t;
  std::vector<std::string> names;
  for (std::size_t j = 0; j < d.m; ++j) { names.push_back("var_s" + std::to_string(j)); }
  t.setColumns(names);
  t.reserve(d.rows());
  for (std::size_t i = 0; i < d.rows(); ++i)
  {
    t.append(d.group[i], d.pair[i], d.label[i] == 0, d.feature[i], d.x.data() + i * d.m);
  }
  return t;
}

// Identifications: winning targets at q <= 0.01; and how many of them are null.
struct Ids { std::size_t ids = 0, null_ids = 0; };
Ids count(const synth::Data& d, const ScoredResult& r)
{
  Ids out;
  for (const auto& g : r.groups)
  {
    if (g.is_decoy || !g.winner || g.qvalue > 0.01) { continue; }
    ++out.ids;
    if (!d.group_true[static_cast<std::size_t>(g.group)]) { ++out.null_ids; }
  }
  return out;
}

bool identical(const ScoredResult& a, const ScoredResult& b)
{
  if (a.dscore != b.dscore || a.groups.size() != b.groups.size()) { return false; }
  for (std::size_t g = 0; g < a.groups.size(); ++g)
  {
    const auto& x = a.groups[g];
    const auto& y = b.groups[g];
    if (x.group != y.group || x.best_row != y.best_row || x.score != y.score || x.winner != y.winner ||
        x.qvalue != y.qvalue || x.pep != y.pep || x.pooled_qvalue != y.pooled_qvalue) { return false; }
  }
  return a.diagnostics.targets_at_q == b.diagnostics.targets_at_q &&
         a.diagnostics.iterations_trained == b.diagnostics.iterations_trained;
}
} // namespace

int main()
{
  const synth::Data d = synth::make({});
  const ScoreTable table = toTable(d);
  const ScoredResult base = odia::core::scoreAndControl(table);
  const auto& diag = base.diagnostics;

  // ---- S1 --------------------------------------------------------------------------------
  {
    const Ids ids = count(d, base);
    std::fprintf(stderr, "[S1] %zu IDs (%zu null) of %d true; decoys at q: %zu; pooled %zu; "
                         "trained %d, skipped %d\n",
                 ids.ids, ids.null_ids, d.n_true, diag.decoys_at_q, diag.pooled_targets_at_q,
                 diag.iterations_trained, diag.iterations_skipped);
    check(ids.ids >= static_cast<std::size_t>(0.8 * d.n_true), "S1 most true hits are identified");
    check(static_cast<double>(ids.null_ids) <= 0.02 * static_cast<double>(ids.ids), "S1 at most 2 % null among IDs");
    check(ids.ids == diag.targets_at_q, "S1 targets_at_q counts the identifications");
    check(diag.rows == d.rows() && diag.groups == d.groups() && diag.target_groups == 2000 &&
          diag.decoy_groups == 2000 && diag.pairs_complete == 2000 && diag.targets_unpaired == 0 &&
          diag.decoys_unpaired == 0, "S1 row, precursor and pair counts");
    check(diag.target_winners + diag.decoy_winners == 2000, "S1 one winner per pair");
    check(diag.n_folds == 3 && diag.iterations_trained > 0, "S1 three folds, trained");
    check(diag.features_used.size() == 8 && diag.features_dropped.empty(), "S1 all eight sub-scores used");
    check(diag.pooled_vs_paired > 0.5 && diag.pooled_vs_paired < 2.0, "S1 pooled and paired agree within 2x");
    bool rows_ok = true;
    for (const auto& g : base.groups)
    {
      rows_ok = rows_ok && g.rows == 4 && table.group[g.best_row] == g.group &&
                base.dscore[g.best_row] == g.score && (g.winner || g.qvalue == 1.0);
    }
    check(rows_ok, "S1 best row belongs to its precursor; losers have q = 1");
  }

  // ---- S2 pure null ----------------------------------------------------------------------
  {
    std::size_t worst = 0;
    for (unsigned seed = 101; seed < 106; ++seed)
    {
      synth::Spec spec;
      spec.seed = seed;
      spec.true_frac = 0.0;
      const synth::Data null_data = synth::make(spec);
      const auto r = odia::core::scoreAndControl(toTable(null_data));
      worst = std::max(worst, r.diagnostics.targets_at_q);
    }
    std::fprintf(stderr, "[S2] pure null, 5 seeds: at most %zu IDs\n", worst);
    check(worst == 0, "S2 pure null gives 0 identifications at q <= 0.01");
  }

  // ---- S3 label swap -----------------------------------------------------------------------
  {
    ScoreTable swapped = table;
    for (auto& v : swapped.is_decoy) { v = v ? 0 : 1; }
    const auto r = odia::core::scoreAndControl(swapped);
    std::fprintf(stderr, "[S3] label swap: %zu IDs (unswapped %zu)\n", r.diagnostics.targets_at_q,
                 diag.targets_at_q);
    check(r.diagnostics.targets_at_q <= std::max<std::size_t>(2, diag.targets_at_q / 100),
          "S3 swapping labels leaves about 0 identifications");
  }

  // ---- S4 permutation ----------------------------------------------------------------------
  {
    std::vector<std::size_t> perm(table.rows());
    std::iota(perm.begin(), perm.end(), std::size_t{0});
    std::mt19937 rng(77);
    std::shuffle(perm.begin(), perm.end(), rng);
    ScoreTable shuffled;
    shuffled.setColumns(table.feature_names);
    for (const std::size_t i : perm)
    {
      shuffled.append(table.group[i], table.pair[i], table.is_decoy[i] != 0, table.feature_id[i], table.row(i));
    }
    const auto r = odia::core::scoreAndControl(shuffled);
    ScoredResult mapped = r;   // back to the original row numbering
    for (std::size_t k = 0; k < perm.size(); ++k) { mapped.dscore[perm[k]] = r.dscore[k]; }
    for (auto& g : mapped.groups) { g.best_row = perm[g.best_row]; }
    check(identical(mapped, base), "S4 row order does not change a single bit");
  }

  // ---- S5 thread count -----------------------------------------------------------------------
  {
#ifdef _OPENMP
    std::fprintf(stderr, "[S5] OpenMP on\n");
#else
    std::fprintf(stderr, "[S5] OpenMP off: the thread count is ignored\n");
#endif
    bool same = true;
    for (const int threads : {1, 2, 3, 8})
    {
      odia::core::Options o;
      o.lda.threads = threads;
      same = same && identical(odia::core::scoreAndControl(table, o), base);
    }
    check(same, "S5 1, 2, 3 and 8 threads give bit-identical results");
  }

  // ---- S6 uninformative columns ----------------------------------------------------------------
  {
    ScoreTable wide;
    std::vector<std::string> names = table.feature_names;
    names.insert(names.begin() + 3, "var_all_missing");
    names.push_back("var_constant");
    wide.setColumns(names);
    const float nan = std::numeric_limits<float>::quiet_NaN();
    for (std::size_t i = 0; i < table.rows(); ++i)
    {
      std::vector<float> row(table.row(i), table.row(i) + table.width());
      row.insert(row.begin() + 3, nan);
      row.push_back(i % 2 ? 4.0f : nan);   // constant where present
      wide.append(table.group[i], table.pair[i], table.is_decoy[i] != 0, table.feature_id[i], row.data());
    }
    const auto r = odia::core::scoreAndControl(wide);
    check(r.diagnostics.features_dropped == std::vector<std::string>({"var_all_missing", "var_constant"}),
          "S6 the all-missing and the constant column are dropped by name");
    check(r.diagnostics.features_used == table.feature_names, "S6 the informative columns stay, in order");
    check(identical(r, base), "S6 dropping them changes nothing else");
  }

  // ---- S7 missing cells ------------------------------------------------------------------------
  {
    ScoreTable holes = table;
    std::mt19937 rng(3);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    std::size_t made = 0;
    for (auto& v : holes.values)
    {
      if (u(rng) < 0.05) { v = std::numeric_limits<float>::quiet_NaN(); ++made; }
    }
    const auto r = odia::core::scoreAndControl(holes);
    const Ids ids = count(d, r);
    std::fprintf(stderr, "[S7] %zu cells missing: %zu IDs (%zu null)\n", made, ids.ids, ids.null_ids);
    check(r.diagnostics.cells_imputed == made, "S7 every missing cell is counted");
    check(ids.ids > static_cast<std::size_t>(0.6 * d.n_true), "S7 still identifies with 5 % missing");
    check(static_cast<double>(ids.null_ids) <= 0.02 * static_cast<double>(ids.ids), "S7 still honest");
  }

  // ---- S8 refusals ---------------------------------------------------------------------------------
  {
    const auto throws = [](const ScoreTable& t) {
      try { odia::core::scoreAndControl(t); }
      catch (const std::invalid_argument& e) { std::fprintf(stderr, "[S8] refused: %s\n", e.what()); return true; }
      return false;
    };
    check(throws(ScoreTable()), "S8 empty table");
    ScoreTable dup = table;
    dup.feature_id[1] = dup.feature_id[0];
    check(throws(dup), "S8 duplicate (group, feature id)");
    ScoreTable no_decoys = table;
    for (auto& v : no_decoys.is_decoy) { v = 0; }
    check(throws(no_decoys), "S8 no decoy precursor");
    ScoreTable mixed = table;
    mixed.is_decoy[1] = mixed.is_decoy[0] ? 0 : 1;
    check(throws(mixed), "S8 a precursor with target and decoy rows");
    ScoreTable two_targets = table;
    for (std::size_t i = 0; i < two_targets.rows(); ++i)
    {
      if (two_targets.group[i] == 3) { two_targets.is_decoy[i] = 0; }   // pair 1's decoy becomes a target
    }
    check(throws(two_targets), "S8 a pair with two targets");
    ScoreTable ragged = table;
    ragged.pair.pop_back();
    check(throws(ragged), "S8 vectors of different lengths");
    ScoreTable flat = table;
    std::fill(flat.values.begin(), flat.values.end(), 1.0f);
    check(throws(flat), "S8 no column varies");
  }

  // ---- S9 peptide and protein level -------------------------------------------------------------
  {
    // Ten pairs per protein; a decoy precursor carries its target's keys.
    std::vector<std::string> peptide(base.groups.size()), protein(base.groups.size());
    for (std::size_t g = 0; g < base.groups.size(); ++g)
    {
      const std::int64_t p = base.groups[g].pair;
      peptide[g] = "PEPTIDE" + std::to_string(p);
      protein[g] = "PROT" + std::to_string(p / 10);
    }
    const auto pep = odia::core::entityQValues(base, peptide, false);
    const auto prot = odia::core::entityQValues(base, protein, true);
    std::size_t pep_ids = 0, prot_ids = 0;
    for (const auto& e : pep.entities) { if (e.label == 1 && e.qvalue <= 0.01) { ++pep_ids; } }
    for (const auto& e : prot.entities) { if (e.label == 1 && e.qvalue <= 0.01) { ++prot_ids; } }
    std::fprintf(stderr, "[S9] peptides at 1 %%: %zu (%zu pairs); proteins at 1 %%: %zu of 200 (%zu pairs)\n",
                 pep_ids, pep.pairs, prot_ids, prot.pairs);
    check(pep.entities.size() == 4000 && pep.pairs == 2000, "S9 one peptide entity per precursor here, paired");
    check(prot.entities.size() == 400 && prot.pairs == 200, "S9 200 target and 200 decoy proteins, paired");
    check(pep_ids > 1000 && prot_ids > 150, "S9 entities are identified");
    bool mapped = true;
    for (std::size_t g = 0; g < base.groups.size(); ++g)
    {
      for (const auto& e : prot.entities)
      {
        if (e.id == protein[g] && (e.label == 1) == !base.groups[g].is_decoy)
        {
          mapped = mapped && prot.group_qvalue[g] == e.qvalue;
        }
      }
    }
    check(mapped, "S9 each precursor gets its own entity's q");

    synth::Spec spec;
    spec.true_frac = 0.0;
    spec.seed = 9;
    const auto null_result = odia::core::scoreAndControl(toTable(synth::make(spec)));
    const auto null_prot = odia::core::entityQValues(null_result, protein, true);
    std::size_t null_ids = 0;
    for (const auto& e : null_prot.entities) { if (e.label == 1 && e.qvalue <= 0.01) { ++null_ids; } }
    check(null_ids == 0, "S9 pure null gives 0 proteins at q <= 0.01");
  }

  return synth::finish("odia_score_and_control_test");
}
