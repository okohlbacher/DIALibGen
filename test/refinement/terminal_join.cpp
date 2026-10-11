// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

// Terminal modifications join a DIA-NN report to the library whichever way
// either side spells them. A library from generation stores OpenMS's
// `.(Acetyl)PEPTIDEK` and `PEPTIDEK.(Amidated)`; one read from a DIA-NN TSV or
// Parquet stores `(UniMod:1)PEPTIDEK` and `PEPTIDEK(UniMod:2)`; DIA-NN reports
// the latter. `Methyl` can sit on the C-terminus and on K or R, so its two
// spellings keep their sites and are joined by refine itself when the library
// holds only one of them.

#include <odia/LibraryRefiner.h>

#include <iostream>
#include <string>
#include <vector>

using namespace ODIA;

namespace
{
  int failures = 0;
  void check(bool ok, const std::string& what)
  { if (!ok) { ++failures; std::cerr << "FAIL: " << what << "\n"; } }

  Library library(const std::vector<std::string>& sequences)
  {
    Library lib;
    auto& p = lib.precursors();
    auto& t = lib.transitions();
    for (std::size_t i = 0; i < sequences.size(); ++i)
    {
      p.mz.push_back(toFixed(500.25 + i)); p.irt.push_back(0.5f); p.im.push_back(1.0f); p.ccs.push_back(400.0f);
      p.charge.push_back(2); p.decoy.push_back(0);
      p.modified_sequence.push_back(lib.strings().intern(sequences[i]));
      p.protein_group.push_back(lib.strings().intern("sp|P1|TEST"));
      p.transition_begin.push_back(static_cast<std::uint32_t>(t.product_mz.size()));
      p.transition_count.push_back(3);
      for (int k = 0; k < 3; ++k)
      {
        t.product_mz.push_back(toFixed(300.0 + 100.0 * k)); t.library_intensity.push_back(1.0f);
        t.type.push_back(FragmentType::Y); t.ordinal.push_back(static_cast<std::uint8_t>(3 + k));
        t.charge.push_back(1); t.loss.push_back(LossType::None);
      }
    }
    return lib;
  }

  /// Refine @p library_sequences with observations named @p report_sequences
  /// (RT 10, 20, ... in order); returns each kept precursor's sequence and RT.
  std::vector<std::pair<std::string, float>> refine(const std::vector<std::string>& library_sequences,
                                                    const std::vector<std::string>& report_sequences,
                                                    RefineStats& stats)
  {
    Library lib = library(library_sequences);
    LibraryRefiner::ObsMap obs;
    for (std::size_t i = 0; i < report_sequences.size(); ++i)
    {
      Observation o;
      o.rt = 10.0f * static_cast<float>(i + 1); o.q = 0.001f;
      obs.emplace(LibraryRefiner::key(report_sequences[i], 2), o);
    }
    RefineParams p;
    LibraryRefiner::refine(lib, obs, p, stats);
    std::vector<std::pair<std::string, float>> out;
    for (std::size_t i = 0; i < lib.precursorCount(); ++i)
    { out.emplace_back(std::string(lib.strings().get(lib.precursors().modified_sequence[i])), lib.precursors().irt[i]); }
    return out;
  }
}

int main()
{
  try
  {
    const std::vector<std::string> report{"(UniMod:1)PEPTIDEK", "PEPTIDEK(UniMod:2)", "LLPEPTIDER(UniMod:34)", "PEPTIDER"};
    // Generation's spelling, and the spelling of a library read from DIA-NN's.
    for (const auto& lib : {std::vector<std::string>{".(Acetyl)PEPTIDEK", "PEPTIDEK.(Amidated)", "LLPEPTIDER.(Methyl)", "PEPTIDER"},
                            std::vector<std::string>{"(UniMod:1)PEPTIDEK", "PEPTIDEK(UniMod:2)", "LLPEPTIDER.(UniMod:34)", "PEPTIDER"}})
    {
      RefineStats stats;
      const auto kept = refine(lib, report, stats);
      check(stats.matched == 4 && stats.ids_unmatched == 0 && kept.size() == 4,
            lib[1] + ": all four precursors joined (" + std::to_string(stats.matched) + " matched, " +
            std::to_string(stats.ids_unmatched) + " unmatched)");
      for (std::size_t i = 0; i < kept.size() && i < 4; ++i)
      {
        check(kept[i].first == lib[i] && kept[i].second == 10.0f * static_cast<float>(i + 1),
              kept[i].first + " carries its own observation (RT " + std::to_string(kept[i].second) + ")");
      }
    }
    // The other way round: an OpenMS-spelled reference against a DIA-NN-spelled library.
    {
      RefineStats stats;
      const auto kept = refine({"PEPTIDEK(UniMod:2)", "LLPEPTIDER(UniMod:34)"}, {"PEPTIDEK.(Amidated)", "LLPEPTIDER.(Methyl)"}, stats);
      check(stats.matched == 2 && kept.size() == 2, "OpenMS-spelled observations join a DIA-NN-spelled library");
    }
    // A library holding BOTH Methyl forms keeps each to its own observation and
    // does not hand the residue form's observation to the C-terminal one.
    {
      RefineStats stats;
      const auto kept = refine({"LLPEPTIDER(Methyl)", "LLPEPTIDER.(Methyl)"}, {"LLPEPTIDER(UniMod:34)"}, stats);
      check(stats.matched == 1 && kept.size() == 1 && kept[0].first == "LLPEPTIDER(Methyl)",
            "with both Methyl forms in the library, only the residue form takes the residue-spelled observation");
    }
  }
  catch (const std::exception& e) { std::cerr << "FAIL: " << e.what() << "\n"; return 1; }
  if (failures == 0) { std::cout << "terminal_join: N- and C-terminal modifications join in either spelling\n"; }
  return failures == 0 ? 0 : 1;
}
