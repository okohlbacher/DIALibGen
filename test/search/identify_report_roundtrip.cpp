// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

// The hand-off: a report written by ReportWriter must read back through the
// UNCHANGED readers that consume -ids -- LibraryRefiner::readObservations for
// refine and, in builds with the fine-tuning stage, the trainer's report
// loader (through validateTrainingReport) for tune -- with the counts, units
// and gates they apply to a DIA-NN report.

#include "synthetic_library.h"

#include <odia/LibraryRefiner.h>
#include <odia/search/ReportWriter.h>
#ifdef HAS_TUNE
#include <odia/tune/Trainer.h>
#endif

#include <arrow/api.h>
#include <arrow/io/file.h>
#include <parquet/arrow/reader.h>

#include <nlohmann/json.hpp>

#include <cmath>
#include <filesystem>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

using namespace ODIA::search;
namespace fs = std::filesystem;

int main(int argc, char** argv)
{
  const fs::path dir = argc > 1 ? fs::path(argv[1]) : fs::temp_directory_path() / "dialibgen-report-roundtrip";
  fs::remove_all(dir);
  fs::create_directories(dir);
  const std::string path = (dir / "run.ids.parquet").string();

  // 1500 target precursors of distinct protein groups (held-out tuning cohorts
  // need >= 100 units each) and 300 decoys; q-values spread across the gates.
  const auto peptides = synth::peptides(1500, 11, 1);
  std::vector<ReportRow> rows;
  std::mt19937_64 rng(5);
  std::size_t pass_all = 0, pass_tune_rt = 0, pass_tune_ccs = 0;
  for (std::size_t k = 0; k < peptides.size(); k += 2)   // charge 2 of every peptide
  {
    const auto& p = peptides[k];
    ReportRow r;
    r.modified_sequence = p.sequence;
    r.charge = p.charge;
    r.precursor_id = p.sequence + std::to_string(p.charge);
    r.precursor_mz = OpenMS::AASequence::fromString(p.sequence).getMZ(p.charge);
    r.protein_group = p.protein;
    r.rt = 5.0 + p.rt * 0.3;                // minutes
    r.rt_start = r.rt - 0.1;
    r.rt_stop = r.rt + 0.1;
    r.irt = p.rt;
    r.im = (k % 10 == 0) ? std::nan("") : 0.7 + 0.001 * static_cast<double>(k % 300);
    r.q = (k % 7 == 0) ? 0.05 : 0.001 + 0.00001 * static_cast<double>(k % 100);
    r.global_q = (k % 11 == 0) ? 0.02 : 0.004;
    r.pg_q = 0.003;
    r.pep = 0.01;
    r.evidence = 3.0 - 0.001 * static_cast<double>(k);
    if (r.q <= 0.01 && r.global_q <= 0.01 && r.pg_q <= 0.01) { ++pass_all; }
    if (r.q <= 0.01) { ++pass_tune_rt; if (std::isfinite(r.im)) { ++pass_tune_ccs; } }
    rows.push_back(r);
  }
  const std::size_t targets = rows.size();
  for (std::size_t k = 0; k < 300; ++k)
  {
    ReportRow d = rows[k];
    d.decoy = true;
    d.precursor_id += "_decoy";
    d.q = 0.008;
    d.evidence = -1.0;
    rows.push_back(d);
  }
  const std::string meta = R"({"tool":"DIALibGen","search":{"settings":{"seed":42}}})";
  ReportWriter::write(path, "synthetic_run", rows, {{"odia.identifier", meta}});
  CHECK(fs::exists(path));

  // Never overwritten.
  {
    bool refused = false;
    try { ReportWriter::write(path, "synthetic_run", rows, {}); }
    catch (const std::runtime_error& e) { refused = std::string(e.what()).find("refusing to overwrite") != std::string::npos; }
    CHECK(refused);
  }

  // The file itself: column order, nulls for absent 1/K0, the metadata.
  {
    auto infile = arrow::io::ReadableFile::Open(path);
    CHECK(infile.ok());
    auto reader = parquet::arrow::OpenFile(*infile, arrow::default_memory_pool());
    CHECK(reader.ok());
    std::shared_ptr<arrow::Table> table;
    CHECK((*reader)->ReadTable(&table).ok());
    CHECK(table->num_rows() == static_cast<std::int64_t>(rows.size()));
    CHECK(table->ColumnNames() == ReportWriter::columns());
    std::int64_t absent = 0;
    for (const auto& r : rows) { absent += std::isnan(r.im) ? 1 : 0; }
    CHECK(absent > 0);
    CHECK(table->GetColumnByName("IM")->null_count() == absent);
    CHECK(table->GetColumnByName("Q.Value")->null_count() == 0);
    CHECK(table->GetColumnByName("Decoy")->type()->id() == arrow::Type::INT32);
    const auto md = table->schema()->metadata();
    CHECK(md && md->Contains("odia.identifier"));
    if (md && md->Contains("odia.identifier"))
    {
      const auto j = nlohmann::json::parse(md->Get("odia.identifier").ValueOrDie());
      CHECK(j["search"]["settings"]["seed"] == 42);
    }
  }

  // refine's reader: every gate column present, decoys rejected, minutes kept.
  {
    ODIA::RefineParams p;
    p.write_im = true;
    ODIA::RefineStats st;
    const auto obs = ODIA::LibraryRefiner::readObservations(path, p, st);
    std::cout << "refine reader: " << st.ids_rows << " rows, " << st.ids_passing << " passing, " << st.ids_decoy
              << " decoy, " << st.ids_q_above << " above q, run " << st.run << "\n";
    CHECK(st.ids_rows == rows.size());
    CHECK(st.ids_rows == targets + 300);
    CHECK(st.ids_decoy == 300);
    CHECK(st.ids_passing == pass_all);
    CHECK(st.ids_q_invalid == 0 && st.ids_charge_invalid == 0 && st.ids_sequence_invalid == 0);
    CHECK(st.run == "synthetic_run");
    CHECK(st.gates_bypassed.empty());
    std::size_t checked = 0;
    for (const auto& r : rows)
    {
      if (r.decoy || r.q > 0.01 || r.global_q > 0.01) { continue; }
      const auto it = obs.find(ODIA::LibraryRefiner::key(r.modified_sequence, r.charge));
      CHECK(it != obs.end());
      if (it == obs.end()) { continue; }
      CHECK(std::fabs(it->second.rt - r.rt) < 1e-4);
      CHECK(std::isfinite(r.im) ? std::fabs(it->second.im - r.im) < 1e-5 : std::isnan(it->second.im));
      CHECK(std::fabs(it->second.q - r.q) < 1e-7);
      CHECK(std::fabs(it->second.evidence - r.evidence) < 1e-5);
      ++checked;
    }
    CHECK(checked == pass_all);
  }

#ifdef HAS_TUNE
  // tune's reader, both heads: the same file must pass its cohort checks.
  for (const bool ccs : {false, true})
  {
    ODIA::tune::TuneParams tp;
    tp.report = path;
    tp.head = ccs ? ODIA::tune::HeadKind::CCS : ODIA::tune::HeadKind::RT;
    std::ostringstream log;
    bool ok = true;
    try { ODIA::tune::validateTrainingReport(tp, log); }
    catch (const std::exception& e) { ok = false; std::cerr << "tune reader: " << e.what() << "\n"; }
    CHECK(ok);
    std::cout << "tune reader (" << (ccs ? "ccs" : "rt") << "): " << log.str();
    const std::string expected = "report: " + std::to_string(rows.size()) + " rows, " +
                                 std::to_string(ccs ? pass_tune_ccs : pass_tune_rt) + " observations";
    CHECK(log.str().find(expected) != std::string::npos);
    CHECK(log.str().find("decoy 300") != std::string::npos);
  }
#else
  std::cout << "tune reader: not built (no fine-tuning stage)\n";
#endif

  fs::remove_all(dir);
  if (synth::failures) { std::cerr << synth::failures << " check(s) failed\n"; return 1; }
  std::cout << "identify_report_roundtrip: PASS\n";
  return 0;
}
