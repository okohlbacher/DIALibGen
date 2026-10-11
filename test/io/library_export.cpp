// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

// The opt-in output layouts, read back:
//
//   diann2 <dir>     DIA-NN 2.0's Parquet layout: the schema column by column
//                    against the one DIA-NN 2.0 writes itself, the values, and
//                    DIALibGen reading the file back.

#include <odia/DIANNLibraryFile.h>
#include <odia/LibraryExport.h>

#include <arrow/api.h>
#include <arrow/io/file.h>
#include <parquet/arrow/reader.h>

#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace ODIA;
namespace fs = std::filesystem;

namespace
{
  void check(bool condition, const std::string& message)
  {
    if (!condition) { throw std::runtime_error(message); }
  }

  /// Five precursors: N- and C-terminal, residue and no modifications, a
  /// decoy, a two-protein group, one precursor without 1/K0. Sequences as
  /// generation stores them (OpenMS names).
  Library fixture()
  {
    Library lib;
    auto& p = lib.precursors();
    auto& t = lib.transitions();
    const float absent = std::numeric_limits<float>::quiet_NaN();
    const char* sequences[] = {".(Acetyl)PEPTIDEK", "AC(Carbamidomethyl)DEFGHK", "LLM(Oxidation)EPTIDER",
                               "PEPTIDEK.(Amidated)", "AC(Carbamidomethyl)DEFGHK"};
    const char* groups[] = {"sp|P10001|ONE_HUMAN;sp|P10002|TWO_HUMAN", "sp|P10003|THREE_HUMAN",
                            "sp|P10003|THREE_HUMAN", "plain_id", "sp|P10003|THREE_HUMAN"};
    for (std::size_t i = 0; i < 5; ++i)
    {
      p.modified_sequence.push_back(lib.strings().intern(sequences[i]));
      p.protein_group.push_back(lib.strings().intern(groups[i]));
      p.mz.push_back(toFixed(500.25 + 10.0 * i));
      p.irt.push_back(0.1f * static_cast<float>(i + 1));
      p.im.push_back(i == 2 ? absent : 0.8f + 0.05f * static_cast<float>(i));
      p.ccs.push_back(i == 2 ? absent : 400.0f + static_cast<float>(i));
      p.charge.push_back(static_cast<std::uint8_t>(2 + i % 2));
      p.decoy.push_back(i == 4 ? 1 : 0);
      p.transition_begin.push_back(static_cast<std::uint32_t>(t.product_mz.size()));
      p.transition_count.push_back(3);
      for (std::uint8_t k = 0; k < 3; ++k)
      {
        t.product_mz.push_back(toFixed(300.125 + 100.0 * k + i));
        t.library_intensity.push_back(k == 0 ? 1.0f : 0.5f / static_cast<float>(k));
        t.type.push_back(k == 1 ? FragmentType::B : FragmentType::Y);
        t.ordinal.push_back(static_cast<std::uint8_t>(2 + k));
        t.charge.push_back(static_cast<std::int8_t>(k == 2 ? 2 : 1));
        t.loss.push_back(i == 1 && k == 2 ? LossType::Water : LossType::None);
      }
    }
    return lib;
  }

  // ------------------------------------------------------------------ diann2

  std::shared_ptr<arrow::Table> readParquet(const std::string& file, std::shared_ptr<const arrow::KeyValueMetadata>* kv)
  {
    auto in = arrow::io::ReadableFile::Open(file);
    check(in.ok(), "cannot open " + file);
    auto reader = parquet::arrow::OpenFile(*in, arrow::default_memory_pool());
    check(reader.ok(), "not Parquet: " + file);
    std::shared_ptr<arrow::Table> table;
    check((*reader)->ReadTable(&table).ok(), "cannot read " + file);
    *kv = (*reader)->parquet_reader()->metadata()->key_value_metadata();
    auto combined = table->CombineChunks();
    check(combined.ok(), "cannot combine chunks");
    return *combined;
  }

  void diann2(const fs::path& dir)
  {
    fs::create_directories(dir);
    const auto file = (dir / "diann2.parquet").string();
    fs::remove(file);
    const Library lib = fixture();
    LibraryExport::storeDIANN2Parquet(file, lib, {{"odia.config_json", "{}"}}, 2);

    // The schema of a library DIA-NN 2.0 writes (--out-lib x.parquet), in order.
    const std::vector<std::pair<std::string, std::string>> expected = {
      {"Precursor.Id", "string"}, {"Modified.Sequence", "string"}, {"Stripped.Sequence", "string"},
      {"Precursor.Charge", "int64"}, {"Proteotypic", "int64"}, {"Decoy", "int64"}, {"N.Term", "int64"},
      {"C.Term", "int64"}, {"RT", "float"}, {"IM", "float"}, {"Q.Value", "float"},
      {"Peptidoform.Q.Value", "float"}, {"PTM.Site.Confidence", "float"}, {"PG.Q.Value", "float"},
      {"Precursor.Mz", "float"}, {"Product.Mz", "float"}, {"Relative.Intensity", "float"},
      {"Fragment.Type", "string"}, {"Fragment.Charge", "int64"}, {"Fragment.Series.Number", "int64"},
      {"Fragment.Loss.Type", "string"}, {"Exclude.From.Quant", "int64"}, {"Protein.Ids", "string"},
      {"Protein.Group", "string"}, {"Protein.Names", "string"}, {"Genes", "string"}, {"Flags", "int64"}};
    std::shared_ptr<const arrow::KeyValueMetadata> kv;
    const auto table = readParquet(file, &kv);
    const auto& schema = *table->schema();
    check(schema.num_fields() == static_cast<int>(expected.size()), "DIA-NN 2.0 has 27 columns, file has " +
          std::to_string(schema.num_fields()));
    for (std::size_t c = 0; c < expected.size(); ++c)
    {
      const auto& f = schema.field(static_cast<int>(c));
      check(f->name() == expected[c].first && f->type()->ToString() == expected[c].second && !f->nullable(),
            "column " + std::to_string(c) + " is " + f->name() + ": " + f->type()->ToString() +
            (f->nullable() ? " (nullable)" : "") + ", expected " + expected[c].first + ": " + expected[c].second);
    }
    check(table->num_rows() == 15, "one row per transition");
    check(kv && kv->Contains("odia.config_json"), "the metadata is written");

    auto str = [&](const char* name, std::int64_t row) {
      return std::string(std::static_pointer_cast<arrow::StringArray>(table->GetColumnByName(name)->chunk(0))->GetView(row));
    };
    auto i64 = [&](const char* name, std::int64_t row) {
      return std::static_pointer_cast<arrow::Int64Array>(table->GetColumnByName(name)->chunk(0))->Value(row);
    };
    auto f32 = [&](const char* name, std::int64_t row) {
      return std::static_pointer_cast<arrow::FloatArray>(table->GetColumnByName(name)->chunk(0))->Value(row);
    };
    // DIA-NN's UniMod notation, terminal modifications as DIA-NN writes them.
    const char* modified[] = {"(UniMod:1)PEPTIDEK", "AC(UniMod:4)DEFGHK", "LLM(UniMod:35)EPTIDER",
                              "PEPTIDEK(UniMod:2)", "AC(UniMod:4)DEFGHK"};
    const char* stripped[] = {"PEPTIDEK", "ACDEFGHK", "LLMEPTIDER", "PEPTIDEK", "ACDEFGHK"};
    for (std::int64_t i = 0; i < 5; ++i)
    {
      const std::int64_t row = 3 * i;
      check(str("Modified.Sequence", row) == modified[i], "Modified.Sequence " + str("Modified.Sequence", row));
      check(str("Stripped.Sequence", row) == stripped[i], "Stripped.Sequence " + str("Stripped.Sequence", row));
      check(str("Precursor.Id", row) == std::string(modified[i]) + std::to_string(2 + i % 2) + (i == 4 ? "_decoy" : ""),
            "Precursor.Id " + str("Precursor.Id", row));
      check(i64("Decoy", row) == (i == 4 ? 1 : 0), "Decoy");
      check(i64("Precursor.Charge", row) == 2 + i % 2, "Precursor.Charge");
      check(std::abs(f32("Precursor.Mz", row) - (500.25 + 10.0 * i)) < 1e-3, "Precursor.Mz");
      check(std::abs(f32("RT", row) - 0.1f * static_cast<float>(i + 1)) < 1e-6, "RT");
      check(f32("IM", row) == (i == 2 ? 0.0f : 0.8f + 0.05f * static_cast<float>(i)), "IM, 0 when absent");
      check(f32("PTM.Site.Confidence", row) == 1.0f && f32("Q.Value", row) == 0.0f && i64("Flags", row) == 0,
            "constant columns as DIA-NN writes them");
    }
    check(str("Protein.Ids", 0) == "P10001;P10002" && str("Protein.Group", 0) == "P10001;P10002" &&
          str("Protein.Names", 0) == "ONE_HUMAN;TWO_HUMAN" && i64("Proteotypic", 0) == 0,
          "a two-protein group: accessions, entry names, not proteotypic");
    check(str("Protein.Ids", 3) == "P10003" && str("Protein.Names", 3) == "THREE_HUMAN" && i64("Proteotypic", 3) == 1,
          "a one-protein group is proteotypic");
    check(str("Protein.Ids", 9) == "plain_id" && str("Protein.Names", 9) == "plain_id", "a plain identifier is both");
    check(str("Fragment.Type", 1) == "b" && i64("Fragment.Series.Number", 1) == 3 && i64("Fragment.Charge", 2) == 2 &&
          str("Fragment.Loss.Type", 5) == "H2O" && str("Fragment.Loss.Type", 0) == "noloss",
          "fragment annotation");
    check(std::abs(f32("Product.Mz", 4) - 401.125f) < 1e-3 && f32("Relative.Intensity", 4) == 0.5f, "fragments");

    // DIALibGen reads its DIA-NN layout back: the same precursors and transitions.
    Library back;
    DIANNLibraryFile::load(file, back);
    check(back.precursorCount() == 5 && back.transitionCount() == 15 && back.decoyCount() == 1,
          "read back: " + std::to_string(back.precursorCount()) + " precursors, " +
          std::to_string(back.transitionCount()) + " transitions");
    std::cout << "diann2_parquet: 27 columns as DIA-NN 2.0 writes them, 15 rows, read back whole\n";
  }

}

int main(int argc, char** argv)
{
  if (argc != 3) { std::cerr << "usage: library_export diann2 <dir>\n"; return 2; }
  try
  {
    const std::string mode = argv[1];
    if (mode == "diann2") { diann2(argv[2]); }
    else { std::cerr << "unknown mode " << mode << '\n'; return 2; }
  }
  catch (const std::exception& e)
  {
    std::cerr << "FAIL: " << e.what() << '\n';
    return 1;
  }
  return 0;
}
