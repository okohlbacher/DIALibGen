// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <odia/PeptDeepEncoder.h>

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/Element.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>

#include <cctype>
#include <map>
#include <stdexcept>

using namespace OpenMS;

namespace ODIA
{
  namespace
  {
    /// OpenMS spells isotopes "(13)C"; AlphaPeptDeep's element list spells them
    /// "13C". These four renames are the only translation needed -- every one of
    /// the 40 element symbols appearing in the shipped unimod.xml is otherwise
    /// present in the 109-element list verbatim. Without them, isotope-labelled
    /// modifications land silently in the '?' bucket.
    std::string toPeptDeepSymbol(const std::string& openms_symbol)
    {
      static const std::map<std::string, std::string> renames{
        // Heavy isotopes have their own slots in AlphaPeptDeep's list.
        {"(2)H", "2H"}, {"(13)C", "13C"}, {"(15)N", "15N"}, {"(18)O", "18O"},
        // Light isotopes are the base element there, and they do occur:
        // ModificationsDB loads PSI-MOD and XLMOD as well as UniMod (3610
        // modifications, not 2859), and Label:13C(8)15N(2) is spelled with
        // explicit (12)C and (14)N. Without these, its C and N counts were
        // dropped into the '?' bucket instead of cancelling the heavy ones.
        {"(1)H", "H"}, {"(12)C", "C"}, {"(14)N", "N"}, {"(16)O", "O"},
      };
      const auto it = renames.find(openms_symbol);
      return it == renames.end() ? openms_symbol : it->second;
    }
  } // namespace

  std::vector<float> PeptDeepEncoder::modificationVector(const std::string& mod_id)
  {
    std::vector<float> feature(PEPTDEEP_MOD_ELEMENTS.size(), 0.0f);

    const ResidueModification* mod = nullptr;
    try
    {
      mod = ModificationsDB::getInstance()->getModification(mod_id);
    }
    catch (const std::exception&)
    {
      throw std::runtime_error("unknown modification: " + mod_id);
    }
    if (mod == nullptr) { throw std::runtime_error("unknown modification: " + mod_id); }

    // A modification with no elemental composition would encode to an all-zero
    // vector -- indistinguishable from no modification at all. OpenMS accepts
    // bare mass shifts such as C[999] and registers them with an empty diff
    // formula, so this is reachable from ordinary input, and zero is not a safe
    // default for something the models are supposed to see.
    if (mod->getDiffFormula().isEmpty())
    {
      throw std::runtime_error(
        "modification '" + mod_id + "' has no elemental composition, so it "
        "cannot be encoded; PeptDeep needs a composition, not a mass shift");
    }

    for (const auto& [element, count] : mod->getDiffFormula())
    {
      const auto symbol = toPeptDeepSymbol(element->getSymbol());
      const auto index = peptDeepElementIndex(symbol);
      if (index)
      {
        // Assigned, not added -- upstream writes `feature[idx] = num`.
        feature[*index] = static_cast<float>(count);
      }
      else
      {
        // Accumulated -- upstream writes `feature[-1] += num`, so two unknown
        // symbols in one formula sum rather than the last one winning.
        feature[PEPTDEEP_UNKNOWN_ELEMENT] += static_cast<float>(count);
      }
    }
    return feature;
  }

  void PeptDeepEncoder::addModification_(Batch& batch, std::size_t row,
                                         std::size_t site, const std::string& mod_id)
  {
    const auto feature = modificationVector(mod_id);
    const std::size_t base =
      (row * batch.sequence_length + site) * PEPTDEEP_MOD_ELEMENTS.size();
    for (std::size_t i = 0; i < feature.size(); ++i)
    {
      // Accumulates: two modifications on one site add, they do not overwrite.
      batch.mod_x[base + i] += feature[i];
    }
  }

  PeptDeepEncoder::Batch PeptDeepEncoder::encode(const AASequence& peptide)
  {
    return encode(std::vector<AASequence>{peptide});
  }

  PeptDeepEncoder::Batch PeptDeepEncoder::encode(const std::vector<AASequence>& peptides)
  {
    if (peptides.empty()) { throw std::invalid_argument("empty peptide batch"); }

    const std::size_t n = peptides.front().size();
    if (n == 0) { throw std::invalid_argument("empty peptide sequence"); }
    for (const auto& p : peptides)
    {
      if (p.size() != n)
      {
        throw std::invalid_argument(
          "a batch must be length-homogeneous: padding is not inert, so mixing "
          "lengths changes the prediction");
      }
    }

    Batch batch;
    batch.rows = peptides.size();
    batch.sequence_length = n + 2;      // one terminal token at each end
    batch.aa_indices.assign(batch.rows * batch.sequence_length, 0);
    batch.mod_x.assign(batch.rows * batch.sequence_length * PEPTDEEP_MOD_ELEMENTS.size(),
                       0.0f);

    for (std::size_t row = 0; row < peptides.size(); ++row)
    {
      const auto& peptide = peptides[row];

      for (std::size_t k = 0; k < n; ++k)
      {
        const auto letter = peptide[k].getOneLetterCode();
        if (letter.size() != 1 || letter[0] < 'A' || letter[0] > 'Z')
        {
          throw std::invalid_argument("residue outside A-Z in " + peptide.toString());
        }
        // A->1 .. Z->26; 0 is the pad and terminal token. Residue k sits at
        // index k+1, so index 0 stays the N-terminal slot.
        batch.aa_indices[row * batch.sequence_length + k + 1] =
          static_cast<std::int64_t>(letter[0] - 'A' + 1);

        if (const auto* mod = peptide[k].getModification(); mod != nullptr)
        {
          addModification_(batch, row, k + 1, mod->getFullId());
        }
      }

      // Terminal modifications belong on the terminal rows, not on the adjacent
      // residue. Placing a C-terminal modification on the last residue instead
      // of row nAA+1 shifted a predicted iRT by more than the modification's
      // own effect.
      if (const auto* mod = peptide.getNTerminalModification(); mod != nullptr)
      {
        addModification_(batch, row, 0, mod->getFullId());
      }
      if (const auto* mod = peptide.getCTerminalModification(); mod != nullptr)
      {
        addModification_(batch, row, n + 1, mod->getFullId());
      }
    }
    return batch;
  }

  namespace
  {
    /// Upper-case, and drop '-', '_' and ' ', so "timsTOF Pro" and "TIMSTOF-PRO"
    /// are one name. Instrument names are written a dozen ways in the wild and a
    /// hyphen is not a different mass spectrometer.
    std::string foldInstrument(const std::string& name)
    {
      std::string f;
      for (const char c : name)
      {
        if (c == '-' || c == '_' || c == ' ') { continue; }
        f.push_back(static_cast<char>(std::toupper(static_cast<unsigned char>(c))));
      }
      return f;
    }

    /// Upstream's `instrument_group` (peptdeep/constants/default_settings.yaml),
    /// all 17 entries, PLUS spellings of our own -- so this is a superset, not a
    /// copy. Pinned in data/peptdeep_meta_inputs.txt.
    ///
    /// An Astral is a Lumos here because that is the group upstream trained it
    /// with. They tried the other reading: commit 74a9816 mapped Astral onto
    /// ThermoTOF, and 4b0aaf1 reverted it, "FIX use lumos model for Astral".
    struct InstrumentAlias { const char* folded; const char* canonical; };
    constexpr InstrumentAlias INSTRUMENT_ALIASES[] = {
      // the five the model indexes, each mapping to itself
      {"QE", "QE"}, {"LUMOS", "Lumos"}, {"TIMSTOF", "timsTOF"},
      {"SCIEXTOF", "SciexTOF"}, {"THERMOTOF", "ThermoTOF"},
      // upstream instrument_group
      {"ASTRAL", "Lumos"}, {"FUSION", "Lumos"}, {"ECLIPSE", "Lumos"},
      {"VELOS", "Lumos"}, {"ELITE", "Lumos"},
      {"ORBITRAPTRIBRID", "Lumos"}, {"THERMOTRIBRID", "Lumos"},
      {"QE+", "QE"}, {"QEHF", "QE"}, {"QEHFX", "QE"},
      {"EXPLORIS", "QE"}, {"EXPLORIS480", "QE"},
      // ours: the spellings PSI-MS, SDRF and PRIDE actually carry
      {"QEPLUS", "QE"}, {"QEXACTIVE", "QE"}, {"QEXACTIVEPLUS", "QE"},
      {"QEXACTIVEHF", "QE"}, {"QEXACTIVEHFX", "QE"},
      {"TIMSTOFPRO", "timsTOF"}, {"TIMSTOFPRO2", "timsTOF"}, {"TIMSTOFSCP", "timsTOF"},
      {"TIMSTOFHT", "timsTOF"}, {"TIMSTOFULTRA", "timsTOF"}, {"TIMSTOFULTRA2", "timsTOF"},
      {"TIMSTOFFLEX", "timsTOF"},
      {"FUSIONLUMOS", "Lumos"}, {"ASTRALZOOM", "Lumos"},
      {"TRIPLETOF", "SciexTOF"}, {"ZENOTOF", "SciexTOF"},
    };
  }


  std::string PeptDeepEncoder::canonicalInstrument(const std::string& name)
  {
    std::string folded = foldInstrument(name);
    for (const InstrumentAlias& a : INSTRUMENT_ALIASES)
    { if (folded == a.folded) { return a.canonical; } }

    // Then the two things a vendor name carries that do not change which of the
    // five it is: a leading "Orbitrap", and a trailing model number. That turns
    // "Orbitrap Exploris 480" and "ZenoTOF 7600" -- the names a user copies out
    // of their own metadata -- into an answer rather than a refusal.
    //
    // The result must match the table EXACTLY. A prefix match here would be a
    // hole big enough to drive the whole feature through: "Astrall" begins with
    // "Astral", and an earlier version of this function duly accepted it, which
    // is precisely the typo the refusal exists to catch.
    if (folded.rfind("ORBITRAP", 0) == 0) { folded.erase(0, 8); }
    while (!folded.empty() && std::isdigit(static_cast<unsigned char>(folded.back()))) { folded.pop_back(); }
    for (const InstrumentAlias& a : INSTRUMENT_ALIASES)
    { if (folded == a.folded) { return a.canonical; } }
    return {};
  }

  float PeptDeepEncoder::defaultNce(const std::string& canonical_instrument)
  {
    // Read out of the shipped checkpoint: meta_nn is Linear(9,7) over
    // one_hot(instrument,8) + nce, and only the QE, timsTOF and nce columns
    // carry weights outside the layer's initialisation bound. Lumos, SciexTOF
    // and ThermoTOF sit at init -- Lumos is effectively the BASELINE, the
    // no-correction case that QE and timsTOF are deltas from. So the label is
    // worth choosing for those two, and the NCE is worth choosing for everyone.
    //
    // timsTOF 30 is upstream's value, kept deliberately over the 40 we measured.
    // On one timsTOF diaPASEF method, normalised spectral angle against the run's own observed
    // fragment areas was 0.8939 +/- 0.0011 at NCE 30 and 0.9041 +/- 0.0004 at 40
    // over three replicates -- a real difference, ten times the replicate
    // spread, and worth +1,028 precursors end to end. It is still not a default:
    // the curve falls about four times more steeply above its peak than below
    // it, so defaulting AT the peak puts every method with a cooler ramp on the
    // steep side, while defaulting below costs our own ramp 0.010. One ramp was
    // measured and the caller's is unknown, so the default sits below the peak.
    // A timsTOF method like ours should set nce 40 explicitly; the README says so.
    //
    // QE 30 is peptdeep's default and Lumos 25 is AlphaDIA's. NEITHER is
    // measured here: every row of our sweep was scored on timsTOF spectra, so it
    // says which LABEL suits timsTOF data, not what NCE a real QE run wants.
    if (canonical_instrument == "timsTOF") { return 30.0f; }
    if (canonical_instrument == "Lumos") { return 25.0f; }
    // QE, SciexTOF and ThermoTOF: upstream's generic default. Named explicitly
    // rather than left to a sentinel, so that no instrument can reach the caller
    // with "no default" and have the fallback reported as its own.
    return 30.0f;
  }

  std::int64_t PeptDeepEncoder::instrumentIndex(const std::string& name)
  {
    // AlphaPeptDeep's list, in its order, pinned in data/peptdeep_meta_inputs.txt
    // against the upstream repo constants -- do not "correct" it from a comment.
    //
    // ThermoTOF=4 was briefly removed here on the grounds that the SHIPPED
    // checkpoint's model_const.yaml names only QE, Lumos, timsTOF and SciexTOF.
    // That was wrong twice over: the pin's source is the upstream constants,
    // which are ahead of this checkpoint, and peptdeep_meta_inputs caught the
    // change immediately -- which is what it is for. Restored.
    //
    // What IS true, and worth knowing before selecting one: max_instrument_num
    // is 8, so slots 4-7 exist in the embedding, but the checkpoint we load
    // lists four instruments, so index 4 is untrained FOR THIS MODEL.
    static const std::map<std::string, std::int64_t> known{
      {"QE", 0}, {"LUMOS", 1}, {"TIMSTOF", 2}, {"SCIEXTOF", 3}, {"THERMOTOF", 4},
    };
    std::string upper;
    for (const char c : name) { upper.push_back(static_cast<char>(std::toupper(static_cast<unsigned char>(c)))); }
    const auto it = known.find(upper);
    // max_instrument_num - 1, the "unknown" slot.
    return it == known.end() ? 7 : it->second;
  }

  PeptDeepEncoder::Batch
  PeptDeepEncoder::encode(const std::vector<AASequence>& peptides,
                          const std::vector<int>& charges,
                          float nce, const std::string& instrument)
  {
    if (charges.size() != peptides.size())
    {
      throw std::invalid_argument("one charge per peptide is required");
    }
    Batch batch = encode(peptides);

    batch.charges.reserve(batch.rows);
    batch.nces.assign(batch.rows, nce * NCE_SCALE);
    batch.instrument_indices.assign(batch.rows, instrumentIndex(instrument));
    for (const int z : charges)
    {
      // A charge below 1 is not a charge, and the model does not say so: it
      // returns a plausible, base-peak-normalised spectrum for 0 and for -3
      // alike. The DIA-NN reader has a live path to charge 0, because a null
      // Precursor.Charge parses to it, so this is reachable from a file rather
      // than only from a caller's mistake. Refusing here rather than in
      // predictMS2 means the per-peptide retry records the offending row and
      // keeps the rest of the chunk.
      if (z < 1 || z > MAX_PRECURSOR_CHARGE)
      {
        throw std::invalid_argument("precursor charge " + std::to_string(z) +
                                    " is outside 1.." +
                                    std::to_string(MAX_PRECURSOR_CHARGE));
      }
      batch.charges.push_back(static_cast<float>(z) * CHARGE_SCALE);
    }
    return batch;
  }

  std::vector<std::vector<std::size_t>>
  PeptDeepEncoder::groupByLength(const std::vector<AASequence>& peptides)
  {
    std::map<std::size_t, std::vector<std::size_t>> by_length;
    for (std::size_t i = 0; i < peptides.size(); ++i)
    {
      by_length[peptides[i].size()].push_back(i);
    }
    std::vector<std::vector<std::size_t>> groups;
    groups.reserve(by_length.size());
    for (auto& [length, indices] : by_length) { groups.push_back(std::move(indices)); }
    return groups;
  }

} // namespace ODIA
