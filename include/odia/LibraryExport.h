// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <odia/Library.h>

#include <string>
#include <utility>
#include <vector>

namespace OpenMS { class TargetedExperiment; }

namespace ODIA
{

  /// The formats a library is written in, selected by `-out_format` and the
  /// extension of `-out`.
  ///
  /// The first two are what `-out x.parquet` and `-out x.tsv` have always
  /// written; the others are opt-in. Each opt-in name carries the version of
  /// the layout it follows, so that a later layout is a new name rather than a
  /// silent change to an existing one.
  enum class LibraryFormat
  {
    DIALibGenParquet,  ///< `dialibgen_parquet`: DIALibGen's compact layout, with provenance
    DIANNTSV,          ///< `diann_tsv`: DIA-NN's TSV library, read by DIA-NN 1.7 and 2.x
    DIANN2Parquet,     ///< `diann2_parquet`: DIA-NN 2.0's own Parquet library layout
    OpenSwathPQP,      ///< `openswath_pqp`: OpenSWATH's PQP (SQLite), written by OpenMS
    OpenSwathTSV       ///< `openswath_tsv`: OpenSWATH's assay TSV, written by OpenMS
  };

  class LibraryExport
  {
  public:
    /// Valid `-out_format` values: `auto` and the name of every format.
    static std::vector<std::string> formatNames();

    /// The `-out_format` name of @p format.
    static const char* name(LibraryFormat format);

    /// The format for @p out and the `-out_format` value @p requested.
    ///
    /// `auto` picks by extension: .parquet DIALibGen's layout, .tsv DIA-NN's
    /// TSV, .pqp OpenSWATH's PQP. An explicit format must match the
    /// extension (.parquet for diann2_parquet, .tsv for openswath_tsv).
    /// @throws std::invalid_argument with a message for the user otherwise.
    static LibraryFormat resolve(const std::string& out, const std::string& requested);

    /// Whether @p format can carry DIALibGen's recipe (the embedded config).
    static bool carriesRecipe(LibraryFormat format);

    /// Write @p library in DIA-NN 2.0's Parquet library layout.
    ///
    /// The columns, their order and their types are those of a library DIA-NN
    /// 2.0 writes itself (`--out-lib x.parquet`), one row per transition:
    ///
    ///   Precursor.Id, Modified.Sequence, Stripped.Sequence: string
    ///   Precursor.Charge, Proteotypic, Decoy, N.Term, C.Term: int64
    ///   RT, IM, Q.Value, Peptidoform.Q.Value, PTM.Site.Confidence, PG.Q.Value,
    ///   Precursor.Mz, Product.Mz, Relative.Intensity: float
    ///   Fragment.Type: string, Fragment.Charge, Fragment.Series.Number: int64
    ///   Fragment.Loss.Type: string, Exclude.From.Quant: int64
    ///   Protein.Ids, Protein.Group, Protein.Names, Genes: string, Flags: int64
    ///
    /// all non-null. Sequences in DIA-NN's UniMod notation (diannModifiedSequence),
    /// proteins as accessions and entry names (`sp|P02768|ALBU_HUMAN` gives
    /// P02768 and ALBU_HUMAN). What a predicted library cannot know is written
    /// as DIA-NN writes it for its own predicted library: the q-values 0, the
    /// site confidence 1, Exclude.From.Quant and Flags 0, an absent 1/K0 0.
    /// N.Term and C.Term (protein-terminal peptide) are 0 and Genes is empty:
    /// the library does not keep them. DIA-NN does not need either to search.
    ///
    /// @param metadata key/value pairs for the file's metadata (DIA-NN ignores
    ///        them), for example `odia.config_json`.
    /// @param threads host threads for building row groups; the file does not
    ///        depend on it.
    static void storeDIANN2Parquet(const std::string& filename, const Library& library,
                                   const std::vector<std::pair<std::string, std::string>>& metadata = {},
                                   unsigned threads = 1);

    /// The library as an OpenMS TargetedExperiment, the input of OpenMS's
    /// OpenSWATH writers. One peptide per precursor (decoys prefixed
    /// `DECOY_`), one transition per fragment, the precursor's proteins with
    /// their accessions, RT as the peptide's normalised retention time and
    /// 1/K0 as its drift time.
    static void toTargetedExperiment(const Library& library, OpenMS::TargetedExperiment& experiment);

    /// Write @p library as OpenSWATH's PQP (@p pqp) or assay TSV, through
    /// OpenMS's TransitionPQPFile and TransitionTSVFile.
    ///
    /// Exports what the library holds, decoys included. A library without
    /// decoys needs OpenSwathDecoyGenerator before OpenSwathWorkflow.
    static void storeOpenSwath(const std::string& filename, const Library& library, bool pqp);

    /// Dispatch on @p format.
    static void store(const std::string& filename, const Library& library, LibraryFormat format,
                      const std::vector<std::pair<std::string, std::string>>& metadata = {},
                      unsigned threads = 1);
  };

} // namespace ODIA
