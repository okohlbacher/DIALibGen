// Copyright (c) 2026, Oliver Kohlbacher and the DIALibraryGenerator authors.
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <odia/Library.h>

#include <cstddef>
#include <string>
#include <vector>

namespace ODIA
{

  /// How decoy precursors are constructed.
  ///
  /// The families attested in the DIA literature, selectable because the choice
  /// is a modelling decision with no settled answer and the tools disagree:
  /// DIA-NN 2.x defaults to SHUFFLING (its Generic mode; the precursor mass is
  /// deliberately left unchanged), DIA-NN 1.7.12-1.8 to two-residue MUTATION,
  /// the OpenSWATH lineage to shuffling with an anti-similarity filter
  /// (Schubert et al., Nat Protoc 2015, citing Rost et al. 2014: decoys "need
  /// to represent the targets well but at the same time they have to be
  /// different from the target assays"), and mProphet to REVERSAL.
  ///
  /// Every method here leaves the precursor m/z, RT, ion mobility and fragment
  /// intensities equal to the target's and changes only fragment masses, which
  /// is what DIA-NN does in both its 1.x mutation and its 2.x Generic mode.
  enum class DecoyMethod
  {
    None,
    Mutate,         ///< DIA-NN 1.7.12-1.8: substitute one residue near each terminus
    PseudoReverse,  ///< reverse the interior, both termini fixed
    Reverse,        ///< reverse the whole sequence
    Shuffle         ///< permute the interior; DIA-NN 2.x's default family
  };

  /// Residues held fixed at each terminus by PseudoReverse and Shuffle.
  ///
  /// One each, matching DIA-NN's --dg-keep-nterm / --dg-keep-cterm defaults.
  /// It matters when the targets share terminal residues by construction --
  /// synthetic peptides, or a protease that fixes the C-terminus.
  inline constexpr std::size_t DECOY_KEEP_NTERM = 1;
  inline constexpr std::size_t DECOY_KEEP_CTERM = 1;

  DecoyMethod parseDecoyMethod(const std::string& s);

  struct DigestParams
  {
    /// "Trypsin/P" rather than "Trypsin": cleave after K and R even when the
    /// next residue is proline.
    ///
    /// The no-cut-before-proline rule is a partial one, and search engines
    /// have largely abandoned it -- DIA-NN's default is `--cut K*,R*`, which
    /// cuts regardless. Measured against DIA-NN's library on the human
    /// proteome, keeping "Trypsin" cost 169,044 peptides, and 48.6% of those
    /// began with a proline while none of ours did. That was the single
    /// largest source of the coverage gap, well ahead of methionine excision.
    std::string enzyme = "Trypsin/P";
    std::size_t missed_cleavages = 1;
    std::size_t min_length = 7;
    std::size_t max_length = 30;

    /// MEASURED, not chosen: against DIA-NN's 12,308 confident precursors on
    /// Astral, a 2,3 library covers 92.76% and every one of the 891 it misses
    /// is charge 1 (193) or charge 4 (698); 1,2,3,4 covers 100.00%. A precursor
    /// the library cannot express is a CEILING on identifications, not a tuning
    /// knob, which is why the default is the one that reaches parity.
    std::vector<int> charges{1, 2, 3, 4};

    /// UniMod-style names as OpenMS knows them, e.g. "Carbamidomethyl (C)".
    std::vector<std::string> fixed_modifications{"Carbamidomethyl (C)"};
    std::vector<std::string> variable_modifications{};
    std::size_t max_variable_modifications = 1;

    /// Add a measured offset to the predicted RT of peptides carrying an
    /// UNMODIFIED cysteine (doc/30).
    ///
    /// AlphaPeptDeep and DIA-NN's predictor were both trained on corpora in
    /// which essentially every cysteine was carbamidomethylated, so the
    /// modification is confounded with the residue and free cysteine cannot be
    /// represented. Free cysteine is more hydrophobic than CAM-cysteine, so it
    /// elutes later than either model predicts -- by +0.90 min per cysteine on
    /// Astral and +0.74 min on S08, and DIA-NN's own predictor misses it by
    /// +0.90 and +1.03 min on the same peptides. Not an ODIA defect, but ODIA
    /// pays for it because it windows on the library RT.
    ///
    /// Self-gating, and therefore safe to leave on: the count is of cysteines
    /// that carry NO modification, so a carbamidomethylated library gets no
    /// correction at all and a library with variable alkylation gets it only
    /// where it applies.
    bool free_cysteine_rt_correction = true;

    /// Also emit 1/K0, derived from the predicted CCS (doc/32).
    ///
    /// CCS remains the authoritative value; this adds the nominal mobility a
    /// timsTOF would report for it. Off means the library carries angstroms
    /// only, which is what ODIA did before and which costs a diaPASEF consumer
    /// the entire mobility dimension -- DIA-NN saw `iIM` = 0 for 100% of its
    /// identifications when searching our library against S08.
    bool derive_ion_mobility = true;

    /// Digest the sequence again with the initiator methionine removed.
    ///
    /// Not an option in practice, though it is one here. The initiator Met is
    /// cleaved co-translationally from a large fraction of proteins, so the
    /// observed N-terminal peptide usually begins one residue in. Leaving this
    /// off cost 272,565 precursors against DIA-NN on the human proteome --
    /// 7.1% of the library, every one of them inside the m/z window and evenly
    /// split across charge, which is what identified it: a filter would not be
    /// charge-neutral.
    bool n_terminal_methionine_excision = true;

    double precursor_mz_min = 350.0;
    double precursor_mz_max = 1200.0;

    double fragment_mz_min = 200.0;
    double fragment_mz_max = 1800.0;
    int max_fragment_charge = 2;

    /// A precursor with fewer usable fragments than this is not searchable.
    ///
    /// Three, not four, to match DIA-NN: its predicted library from the same
    /// FASTA carries precursors with as few as 3 transitions, and the extra
    /// requirement was most of our 0.6% precursor surplus (4,984,739 against
    /// 4,954,236). A precursor the reference engine will search and we will not
    /// is a ceiling on identifications, not a quality filter.
    std::size_t min_fragments = 3;
    /// Cap per precursor, taken in descending fragment m/z.
    std::size_t max_fragments = 12;

    /// Smallest predicted intensity, relative to the precursor's base peak,
    /// that may enter the library.
    ///
    /// The guard used to be `intensity > 0`, which admits denormals. Measured
    /// on the shipped S08 library: 1,980 of 5,000 precursors' transitions had a
    /// stored intensity of median 2.0e-07 against a library median of 0.255 --
    /// six orders of magnitude down, i.e. numerical noise. Only 13.6% of them
    /// were ever observed, against 81.8% for ordinary picks, and 1,081 of them
    /// were b2+, which was most of our b2+ over-representation.
    ///
    /// It also made the library non-reproducible: on the H100 those values come
    /// out as tiny positives and on CPU as exact zeros, so the same config
    /// produced different libraries on different hardware.
    double min_relative_intensity = 1e-4;

    /// Reserve this many of the cap for doubly-charged fragments, if the
    /// precursor has any worth keeping.
    ///
    /// Zero means pure intensity ranking, which is the default and is faithful
    /// to the model. The reason to consider otherwise is measured: DIA-NN's
    /// library carries 21.7% doubly-charged fragments for a 2+ precursor and
    /// 35.7% for a 3+, where ours carries 17.5% and 25.8% -- and the precursors
    /// we fail to identify are enriched in charge 3 at 37.0% against 24.5%.
    ///
    /// But the cause is a disagreement between the two INTENSITY models, not a
    /// difference in the rule: for a 2+ precursor our model puts a doubly
    /// charged fragment at median 0.098 where DIA-NN's puts it at 0.155, so
    /// ranking by intensity demotes them exactly as it should given what we
    /// believe. Setting this above zero therefore overrides the model rather
    /// than correcting it, and is only justified if a search says so.
    std::size_t reserved_doubly_charged = 0;

    DecoyMethod decoy_method = DecoyMethod::Mutate;
  };

  /// Builds an assay library from protein sequences.
  ///
  /// Produces the peptide, precursor and fragment skeleton with placeholder
  /// intensities and retention times; those are filled in by prediction.
  class LibraryGenerator
  {
  public:
    struct Stats
    {
      std::size_t proteins = 0;
      std::size_t peptides = 0;
      std::size_t precursors = 0;
      std::size_t transitions = 0;
      std::size_t dropped_too_few_fragments = 0;
      std::size_t dropped_precursor_mz = 0;
    };

    static Stats generate(const std::string& fasta_file,
                          const DigestParams& params,
                          Library& library);

    /// Fill in predicted retention times, replacing the placeholders.
    ///
    /// Predicts once per distinct modified sequence rather than once per
    /// precursor: the RT model takes no charge input, so the charge states of
    /// one peptide would otherwise be predicted identically several times over.
    /// On the human proteome that is 1.25 M predictions instead of 4.0 M.
    ///
    /// @returns the number of precursors whose iRT could not be predicted;
    ///          theirs are left NaN rather than given a made-up value.
    /// The canonical content-affecting parameter string for the library cache.
    ///
    /// ONE definition, called by every tool that builds a library. It used to
    /// be assembled inline in OpenDIAlyzer; a second tool assembling its own
    /// would silently disagree about what "the same library" is, and a
    /// fingerprint miss regenerates for hours without saying why.
    /// @returns the TARGET key: everything content-affecting EXCEPT the decoy
    /// method, which callers append as ";decoy=<method>" for the full key. No
    /// part of the inference depends on the decoy method, so changing it must
    /// not discard a library.
    static std::string fingerprintParams(const DigestParams& p,
                                         const std::string& rt_model,
                                         const std::string& ms2_model,
                                         const std::string& ccs_model,
                                         double nce, const std::string& instrument,
                                         bool irt_rescale);

    static std::size_t predictRetentionTimes(Library& library,
                                             const std::string& rt_model_path,
                                             bool prefer_gpu = true,
                                             unsigned sessions = 0,
                                             bool free_cysteine_rt_correction = true);

    /// The measured free-cysteine offset, in the RT model's own normalised
    /// units, for a peptide with @p free_cysteines unmodified cysteines.
    ///
    /// Bucketed rather than linear because the effect saturates, and stated in
    /// rt_norm rather than minutes because minutes are gradient-dependent.
    /// Fitted on Astral (Orbitrap, 8.5-37.9 min) and S08 (timsTOF, 7-30 min)
    /// and cross-validated by holding each out: fitting on S08 and applying to
    /// Astral takes the cysteine-peptide residual p95 from 3.193 to 2.416 min,
    /// and the reverse takes S08 from 3.327 to 2.678 min. See doc/30.
    static float freeCysteineRtOffset(std::size_t free_cysteines);

    /// Replace placeholder intensities with predicted ones, and re-choose the
    /// fragments now that there is a basis for ranking.
    ///
    /// This necessarily rebuilds the transition arrays rather than editing
    /// them: generate() caps fragments by descending m/z because it has no
    /// intensities, and the top twelve by m/z are not the top twelve by
    /// intensity, so re-ranking what survived that cap would not give the same
    /// answer as ranking before it.
    ///
    /// Predicts in blocks. A whole-proteome call would materialise ~4 M
    /// spectra -- about 2.2 GiB of payload in as many allocations -- all live
    /// at once, on top of the library being built.
    ///
    /// Call before appendDecoys: a decoy copies its target's intensity
    /// pattern, so predicting afterwards would leave every decoy with the
    /// placeholder.
    ///
    /// @returns the number of precursors left with placeholder intensities
    ///          because prediction failed for them.
    static std::size_t predictFragmentIntensities(Library& library,
                                                  const std::string& ms2_model_path,
                                                  const DigestParams& params,
                                                  float nce = 30.0f,
                                                  const std::string& instrument = "QE",
                                                  bool prefer_gpu = true,
                                                  unsigned sessions = 0);

    /// The line mapping the RT model's raw output onto the iRT scale.
    struct IrtCalibration
    {
      double slope = 1.0;
      double intercept = 0.0;
      std::size_t peptides = 0;      ///< standards that could be predicted
      double max_abs_error = 0.0;    ///< worst standard, in iRT units

      double apply(double raw) const { return slope * raw + intercept; }
    };

    /// Fit the raw-to-iRT line from the Biognosys standards, using @p model.
    ///
    /// Recomputed rather than pinned: the line depends on the model checkpoint,
    /// and a hardcoded one would rot silently on a model swap -- the outputs
    /// would still look like iRT and simply be wrong.
    static IrtCalibration fitIrtCalibration(const std::string& rt_model_path,
                                            const std::string& standards_file,
                                            bool prefer_gpu = true);

    /// Rescale a library's retention times onto the iRT scale in place.
    ///
    /// This does NOT make the predictions more accurate, and must not be
    /// described as though it does. It is monotone, so any consumer that fits
    /// its own retention-time calibration -- DIA-NN does -- is unaffected:
    /// measured on the human proteome, applying it left DIA-NN's search window
    /// unchanged to the last digit (2.18905 both ways). It is applied because a
    /// column named iRT holding a 0..1 training-gradient coordinate is a
    /// silent-units error waiting for a consumer that applies a tolerance in
    /// iRT units without calibrating first.
    static void applyIrtCalibration(Library& library, const IrtCalibration& calibration);

    /// Fill in predicted collision cross-sections, in square angstroms.
    ///
    /// Per precursor, not per peptide: the CCS model takes charge, and a
    /// peptide's charge states have genuinely different cross-sections.
    ///
    /// This does NOT populate the ion-mobility column. CCS and 1/K0 are
    /// different quantities, related through the drift gas and the
    /// instrument's calibration, and that conversion belongs downstream where
    /// the instrument is known.
    ///
    /// @returns the number of precursors left without a value; theirs are NaN.
    static std::size_t predictCollisionCrossSections(Library& library,
                                                    const std::string& ccs_model_path,
                                                    bool prefer_gpu = true,
                                                    unsigned sessions = 0,
                                                    bool derive_mobility = true);

    /// Append a decoy for every target currently in @p library.
    ///
    /// Decoys keep the target's precursor m/z, iRT and intensity pattern; only
    /// fragment m/z values move. See DecoyMethod.
    /// @param skipped_out receives the number of targets for which no decoy
    ///        could be built. Reported rather than dropped silently: 1.3% of the
    ///        DIA-NN fixture is affected, all of them N-terminally modified.
    /// @param min_fragments the same bar the targets had to clear. Applying it
    ///        to one class only is an anti-conservative FDR (D7 rule 2).
    /// Idempotent: targets that already have a decoy are skipped.
    /// @param recompute_decoy_mz recompute the decoy's PRECURSOR m/z from its
    ///        own sequence instead of inheriting the target's, and re-derive its
    ///        1/K0 from its collision cross-section at that new mass so it stays
    ///        a consistent ion in the two-dimensional diaPASEF window. Off by
    ///        default pending its first measured arm.
    static std::size_t appendDecoys(Library& library, DecoyMethod method,
                                    std::size_t* skipped_out = nullptr,
                                    std::size_t min_fragments = 0,
                                    bool recompute_decoy_mz = false);
  };

} // namespace ODIA
