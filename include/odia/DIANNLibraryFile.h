// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <odia/Library.h>
#include <arrow/type_fwd.h>

#include <string>
#include <utility>
#include <vector>

namespace ODIA
{

  /// Reads and writes assay libraries: DIA-NN's TSV dialect, and DIALibGen's
  /// own Parquet layout.
  ///
  /// The TSV is one row per transition, with every precursor field repeated --
  /// exactly the redundancy the CSR layout in Library collapses on load. Its
  /// column names are `TSVColumns`, chosen so that DIA-NN 1.7 and 2.x both read
  /// the whole file; the Parquet layouts use `Columns`.
  class DIANNLibraryFile
  {
  public:
    /// Column names, as DIA-NN writes them.
    struct Columns
    {
      static constexpr const char* PRECURSOR_ID = "Precursor.Id";
      static constexpr const char* MODIFIED_SEQUENCE = "Modified.Sequence";
      static constexpr const char* STRIPPED_SEQUENCE = "Stripped.Sequence";
      static constexpr const char* PRECURSOR_CHARGE = "Precursor.Charge";
      static constexpr const char* DECOY = "Decoy";
      static constexpr const char* RT = "RT";
      static constexpr const char* IM = "IM";
      /// Predicted collision cross-section, in square angstroms.
      ///
      /// Its own column, and NOT interchangeable with IM: `IM` is 1/K0 in
      /// Vs/cm^2 (Bruker's convention, ~0.6-1.6 for tryptic peptides), CCS is
      /// A^2 (~300-700). The two are related by Mason-Schamp through the drift
      /// gas and the instrument's calibration.
      ///
      /// Both are now written and both are read. A library carrying only one
      /// gets the other filled by `completeMobility` on load, so a consumer
      /// never has to know which the producer chose. Verified against 37,193
      /// MEASURED 1/K0 values on IH1: our derived mobility sits at 2.8%
      /// relative error, and the fitted coefficient is within 2.0% of the
      /// textbook constant, which is what establishes the units agree with
      /// DIA-NN's.
      static constexpr const char* CCS = "CCS";
      static constexpr const char* PRECURSOR_MZ = "Precursor.Mz";
      static constexpr const char* PRODUCT_MZ = "Product.Mz";
      static constexpr const char* RELATIVE_INTENSITY = "Relative.Intensity";
      static constexpr const char* FRAGMENT_TYPE = "Fragment.Type";
      static constexpr const char* FRAGMENT_CHARGE = "Fragment.Charge";
      static constexpr const char* FRAGMENT_SERIES_NUMBER = "Fragment.Series.Number";
      static constexpr const char* FRAGMENT_LOSS_TYPE = "Fragment.Loss.Type";
      static constexpr const char* PROTEIN_GROUP = "Protein.Group";
    };

    /// Column names of the TSV that `storeTSV` writes.
    ///
    /// DIA-NN matches a library TSV's header against fixed lists of accepted
    /// names. DIA-NN 2.x accepts both its own dotted names (`Precursor.Mz`) and
    /// the older undotted ones; DIA-NN 1.7 accepts ONLY the undotted ones and
    /// refuses a file that lacks any of sequence, charge, precursor m/z, RT,
    /// product m/z or intensity. So the undotted names are the only spelling
    /// both read. Verified by loading libraries in DIA-NN 1.7.12 and 2.0.
    ///
    /// Protein identity: neither version reads `Protein.Group` from a library
    /// TSV -- with it alone both report "0 protein groups". Both read the
    /// accessions from `UniprotID` (2.x: `Protein.Ids`), entry names from
    /// `ProteinName` and the flag from `Proteotypic`. `Precursor.Id`,
    /// `Protein.Group` and `CCS` are not DIA-NN columns; DIA-NN ignores them,
    /// and they carry DIALibGen's own id, verbatim protein group and
    /// cross-section so that `loadTSV` reads back exactly what was written.
    struct TSVColumns
    {
      static constexpr const char* PRECURSOR_ID = "Precursor.Id";
      static constexpr const char* MODIFIED_PEPTIDE = "ModifiedPeptide";
      static constexpr const char* PRECURSOR_CHARGE = "PrecursorCharge";
      static constexpr const char* DECOY = "Decoy";
      /// DIA-NN treats library RT as a reference scale and fits its own map
      /// to the run, so the unit is whatever the library holds: a generated
      /// library's normalised 0..1 prediction, or a tuned/refined library's
      /// run minutes. Written as stored, never rescaled.
      static constexpr const char* RT = "iRT";
      static constexpr const char* IM = "IonMobility";
      static constexpr const char* PRECURSOR_MZ = "PrecursorMz";
      static constexpr const char* PRODUCT_MZ = "ProductMz";
      static constexpr const char* RELATIVE_INTENSITY = "LibraryIntensity";
      static constexpr const char* FRAGMENT_TYPE = "FragmentType";
      static constexpr const char* FRAGMENT_CHARGE = "FragmentCharge";
      static constexpr const char* FRAGMENT_SERIES_NUMBER = "FragmentSeriesNumber";
      static constexpr const char* FRAGMENT_LOSS_TYPE = "FragmentLossType";
      /// Protein.Group's members' accessions, ';'-separated: `P02768` from
      /// `sp|P02768|ALBU_HUMAN`, the whole member when it has no `db|acc|name` form.
      static constexpr const char* PROTEIN_IDS = "UniprotID";
      /// The members' entry names, `ALBU_HUMAN`, by the same rule.
      static constexpr const char* PROTEIN_NAMES = "ProteinName";
      /// 1 when the precursor maps to exactly one protein, else 0 -- the rule
      /// DIA-NN itself applies when the column is absent.
      static constexpr const char* PROTEOTYPIC = "Proteotypic";
      static constexpr const char* PROTEIN_GROUP = "Protein.Group";
      static constexpr const char* CCS = "CCS";
    };

    /// Load a library. Accepts .parquet or .tsv (case insensitive).
    ///
    /// Rows are grouped into precursors by a change in Precursor.Id, so the
    /// input must keep a precursor's transitions together -- which DIA-NN does.
    ///
    /// A TSV may use either spelling of a column: the `TSVColumns` names that
    /// `storeTSV` writes now, the dotted `Columns` names that 0.11.1 and earlier wrote,
    /// or DIA-NN's other common aliases (`Tr_recalibrated`, `FragmentMz`,
    /// `RelativeIntensity`, ...). Without a Precursor.Id column the id is
    /// built from sequence, charge and decoy flag; without Protein.Group the
    /// protein group is taken from Protein.Ids/UniprotID.
    static void load(const std::string& filename, Library& library);

    /// Fill in whichever of 1/K0 and CCS the file omitted, from the other.
    ///
    /// Idempotent, and a no-op where both or neither are present. Public
    /// because a caller that builds a library in memory rather than reading
    /// one wants the same completion.
    static void completeMobility(Library& library);

    static void loadTSV(const std::string& filename, Library& library);
    static void loadParquet(const std::string& filename, Library& library);

    /// Write in DIA-NN's TSV dialect, with the `TSVColumns` header.
    ///
    /// The sequence is written in UniMod notation (`C(UniMod:4)`), which is
    /// what DIA-NN parses; a modification without a UniMod accession is
    /// written verbatim and counted in a warning. Precursor.Id is synthesised
    /// as <sequence><charge>, with `_decoy` appended for a decoy.
    static void storeTSV(const std::string& filename, const Library& library);

    /// Write the library as flat Parquet, one row per transition, with the
    /// `Columns` names, so `loadParquet` reads back what this writes.
    ///
    /// Worth having because the TSV is one row per TRANSITION with the sequence
    /// and protein group repeated 12 times: a proteome library is 12.35 GiB of
    /// text that takes 67-86 s to re-parse on every single search. Parquet
    /// dictionary-encodes exactly those repeated strings.
    static void storeParquet(const std::string& filename, const Library& library);

    /// Dispatches on the extension, mirroring `load`.
    static void store(const std::string& filename, const Library& library);

    /// Identity of a generated library: the FASTA it came from and every
    /// parameter that changes its contents.
    ///
    /// Predicting a proteome library costs ~23 minutes of GPU-less inference,
    /// and the inputs rarely change between runs -- so a library is worth
    /// reusing, but ONLY when it was built from the same FASTA under the same
    /// rules. Reusing one built with different charges or a different m/z
    /// window would silently answer a different question, which is the failure
    /// this whole benchmark keeps tripping over.
    ///
    /// The FASTA is identified by CONTENT, not path or mtime: the same file
    /// copied to /scratch must hit the cache, and an edited file with the same
    /// size and timestamp must miss it.
    struct Fingerprint
    {
      std::string fasta_hash;      ///< 64-bit FNV-1a over the file's bytes, hex
      std::uint64_t fasta_bytes = 0;
      std::string params;          ///< every content-affecting parameter, canonical

      /// Everything in `params` EXCEPT the decoy method.
      ///
      /// The expensive half of building a library is inference -- ~19 minutes
      /// of retention time, fragment intensity and CCS prediction on the human
      /// proteome -- and none of it depends on how decoys are made. Keying the
      /// cache on the full parameter set therefore threw all of it away when
      /// only the decoy method changed, which is exactly what happens when
      /// comparing Mutate against PseudoReverse.
      std::string target_params;
      std::string decoy_method;

      std::string key() const
      { return fasta_hash + ":" + std::to_string(fasta_bytes) + ":" + params; }

      /// Identifies the PREDICTIONS, independent of the decoy layer.
      std::string targetKey() const
      { return fasta_hash + ":" + std::to_string(fasta_bytes) + ":" + target_params; }


    };

    /// Hash a FASTA by content. Throws if it cannot be read.
    static Fingerprint fingerprintFasta(const std::string& fasta);

    /// FNV-1a over a file's CONTENT, as 16 hex digits. Empty path -> "bundled".
    ///
    /// Used to identify PREDICTOR MODELS in the cache key. A path is not an
    /// identity: a changed model at the same path silently collides, and the
    /// same model staged elsewhere silently misses.
    static std::string hashFile(const std::string& path);

    /// Write with a fingerprint recorded in the file's metadata. Parquet only:
    /// TSV has nowhere to put it, so a TSV library is never cache-eligible.
    static void storeParquet(const std::string& filename, const Library& library,
                             const Fingerprint& fp);

    /// The fingerprint recorded in a Parquet library, or empty if the file is
    /// missing, unreadable, or carries none. Reads ONLY the metadata -- it must
    /// not cost a full table read to decide whether to use the table.
    static std::string readFingerprint(const std::string& filename,
                                       const char* key = "odia.fingerprint");

  private:
    /// Read the compact (one row per precursor, transitions nested) layout.
    static void loadParquetCompact(const std::shared_ptr<arrow::Table>& table,
                                   Library& library);
  public:

    /// Write the library in ODIA's COMPACT Parquet layout: one row per
    /// PRECURSOR, with the transitions as Arrow list columns.
    ///
    /// The flat layout writes one row per TRANSITION, so every precursor-level
    /// value -- id, sequence, protein, precursor m/z, RT, IM, CCS, charge,
    /// decoy -- is repeated for each of that precursor's ~12 fragments. On the
    /// parity library that is ~44 of 56 bytes per row describing only
    /// 9,983,789 precursors across 119,088,506 rows.
    ///
    /// Storing them once takes the raw form from ~6.2 GiB to ~1.7 GiB, which is
    /// about what the in-memory library costs (1,845.7 MiB) -- the right
    /// target, since that representation is already the well-designed one.
    ///
    /// `loadParquet` detects the layout from the schema (Product.Mz being a
    /// list) and reads either, so old files keep working.
    /// @param config_json when non-empty, embedded verbatim as the schema
    /// metadata key `odia.config_json`: the recipe that produced the library,
    /// travelling with it. Readers ignore unknown keys, so this is additive.
    /// @param extra_metadata further key/value pairs appended after it (for
    /// example `odia.warnings`); nothing is written for an empty list.
    static void storeParquetCompact(const std::string& filename, const Library& library,
                                    const Fingerprint& fp,
                                    const std::string& config_json = "",
                                    const std::vector<std::pair<std::string, std::string>>& extra_metadata = {});

  };

} // namespace ODIA
