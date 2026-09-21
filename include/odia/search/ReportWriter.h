// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

/// The identification report of a built-in search (-out_ids).
///
/// Written in the column names the existing report readers expect
/// (LibraryRefiner::readObservations for refine, the trainer's report loader
/// for tune), so that after it is written it simply BECOMES -ids for the rest
/// of the invocation and no reader, gate or provenance rule changes. The
/// column names are a compatibility dialect; the Parquet metadata key
/// "odia.identifier" says which tool and settings produced the file.
///
/// One row per reported precursor, targets and decoys (Decoy 1). Units: RT,
/// RT.Start and RT.Stop in MINUTES of run time; iRT is the library's own RT
/// value; IM is 1/K0 (null when the run has none). Missing values are nulls.
#pragma once

#include <string>
#include <utility>
#include <vector>

namespace ODIA::search
{
  struct ReportRow
  {
    std::string precursor_id;        ///< Precursor.Id
    std::string modified_sequence;   ///< Modified.Sequence: the library string verbatim
    int charge = 0;                  ///< Precursor.Charge
    double precursor_mz = 0.0;       ///< Precursor.Mz
    std::string protein_group;       ///< Protein.Group: the library string verbatim
    bool decoy = false;              ///< Decoy
    double rt = 0.0;                 ///< RT, minutes: apex of the chosen peak group
    double rt_start = 0.0;           ///< RT.Start, minutes
    double rt_stop = 0.0;            ///< RT.Stop, minutes
    double irt = 0.0;                ///< iRT: the library's RT value
    double im = 0.0;                 ///< IM, 1/K0; NaN = none (written as null)
    double q = 1.0;                  ///< Q.Value: precursor
    double global_q = 1.0;           ///< Global.Q.Value: peptide
    double pg_q = 1.0;               ///< PG.Q.Value: protein group
    double pep = 1.0;                ///< PEP
    double evidence = 0.0;           ///< Evidence: the classifier's d-score
  };

  class ReportWriter
  {
  public:
    /// Column names, in file order.
    static const std::vector<std::string>& columns();

    /// Write @p rows as Parquet (ZSTD) to @p path, with Run = @p run on every
    /// row and @p metadata as Arrow schema key-value metadata. Staged and
    /// committed atomically; refuses a @p path that already exists (the
    /// report is never overwritten). Throws std::runtime_error.
    static void write(const std::string& path, const std::string& run, const std::vector<ReportRow>& rows,
                      const std::vector<std::pair<std::string, std::string>>& metadata);
  };
}
