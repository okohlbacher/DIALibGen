// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <odia/PeptDeepElements.h>

#include <cstdint>
#include <string>
#include <vector>

namespace OpenMS { class AASequence; }

namespace ODIA
{

  /// Builds the input tensors AlphaPeptDeep's PeptDeep models expect.
  ///
  /// OpenMS ships these models and a binding for them, but the binding rejects
  /// modified peptides outright and its headers are not installed, so ODIA
  /// encodes for itself. The encoding is specified in doc/04-peptdeep-encoding.md
  /// and cross-checked against test/peptdeep_reference.py, an independent
  /// implementation written from that spec.
  ///
  /// Input is an OpenMS AASequence rather than a string: the terminal-versus-
  /// residue distinction is exactly where a string parser goes wrong, and
  /// AASequence already carries it.
  class PeptDeepEncoder
  {
  public:
    /// Tensors for one batch. All peptides in a batch must have the same
    /// residue count -- padding is not inert, because index 0 is one-hot
    /// encoded and no model applies a padding mask.
    struct Batch
    {
      std::size_t rows = 0;             ///< peptides
      std::size_t sequence_length = 0;  ///< nAA + 2, including both terminal tokens
      std::vector<std::int64_t> aa_indices;   ///< [rows, sequence_length]
      std::vector<float> mod_x;              ///< [rows, sequence_length, 109]

      // MS2 and CCS only, empty otherwise. The scale factors are applied here
      // rather than left to the caller: passing raw NCE gives a spectrum with
      // cosine 0.0028 against the correct one, and raw charge 0.6377 -- both
      // unrelated output with no error, so the scaling must not be optional.
      std::vector<float> charges;                  ///< [rows, 1], already x0.1
      std::vector<float> nces;                     ///< [rows, 1], already x0.01
      std::vector<std::int64_t> instrument_indices; ///< [rows], rank 1, not [rows, 1]
    };

    static constexpr float CHARGE_SCALE = 0.1f;

    /// Upper bound on a precursor charge the model is asked to predict for.
    /// Not a model limit -- there is none -- but a bound past which the input
    /// is certainly wrong rather than unusual. AlphaPeptDeep's own training
    /// data tops out well below this.
    static constexpr int MAX_PRECURSOR_CHARGE = 10;
    static constexpr float NCE_SCALE = 0.01f;

    /// Instrument index. Anything unrecognised maps to max_instrument_num - 1,
    /// not to 0: mapping an unknown instrument onto QE would quietly predict
    /// for the wrong one.
    static std::int64_t instrumentIndex(const std::string& name);

    /// Encode one peptide. Throws if it is empty or carries a residue outside
    /// A-Z, which would one-hot to an all-off row and so be indistinguishable
    /// from padding.
    static Batch encode(const OpenMS::AASequence& peptide);

    /// Encode several peptides of equal length into one batch.
    /// Throws if the lengths differ.
    static Batch encode(const std::vector<OpenMS::AASequence>& peptides);

    /// As above, plus the meta inputs the MS2 model needs.
    /// @param charges one per peptide, unscaled; scaling is applied here.
    static Batch encode(const std::vector<OpenMS::AASequence>& peptides,
                        const std::vector<int>& charges,
                        float nce, const std::string& instrument);

    /// Group by residue count, so each batch is length-homogeneous. Returns the
    /// index groups; the caller encodes and runs each in turn.
    static std::vector<std::vector<std::size_t>>
    groupByLength(const std::vector<OpenMS::AASequence>& peptides);

    /// Known limitation: OpenMS's AASequence keeps only one modification per
    /// residue, so a doubly-modified residue loses all but the last one before
    /// the encoder ever sees it. AlphaPeptDeep would accumulate them. The
    /// accumulation in addModification_ is therefore correct but currently
    /// unreachable, and a doubly-modified residue is silently under-encoded.
    /// See doc/BACKLOG.md.
    ///
    /// The 109-wide feature vector for a modification's elemental composition.
    ///
    /// Counts are signed -- roughly a fifth of UniMod has negative entries, and
    /// Deamidated is H(-1)N(-1)O(1). Known elements are assigned, unknown ones
    /// accumulate into the '?' slot, matching upstream exactly; the two differ,
    /// so the distinction has to be preserved.
    static std::vector<float> modificationVector(const std::string& mod_id);

  private:
    static void addModification_(Batch& batch, std::size_t row, std::size_t site,
                                 const std::string& mod_id);
  };

} // namespace ODIA
