// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <odia/search/ReportWriter.h>

#include <odia/AtomicFile.h>

#include <arrow/api.h>
#include <arrow/io/file.h>
#include <parquet/arrow/writer.h>

#include <cmath>
#include <filesystem>
#include <memory>
#include <stdexcept>

namespace ODIA::search
{
  const std::vector<std::string>& ReportWriter::columns()
  {
    static const std::vector<std::string> names = {
      "Run", "Precursor.Id", "Modified.Sequence", "Precursor.Charge", "Precursor.Mz", "Protein.Group", "Decoy",
      "RT", "RT.Start", "RT.Stop", "iRT", "IM", "Q.Value", "Global.Q.Value", "PG.Q.Value", "PEP", "Evidence"};
    return names;
  }

  void ReportWriter::write(const std::string& path, const std::string& run, const std::vector<ReportRow>& rows,
                           const std::vector<std::pair<std::string, std::string>>& metadata)
  {
    if (std::filesystem::exists(path) || std::filesystem::is_symlink(path))
    { throw std::runtime_error("refusing to overwrite existing output: " + path); }

    auto ok = [](const arrow::Status& st) {
      if (!st.ok()) { throw std::runtime_error("identification report: " + st.ToString()); }
    };
    arrow::StringBuilder b_run, b_id, b_seq, b_pg;
    arrow::Int32Builder b_z, b_decoy;
    arrow::DoubleBuilder b_mz, b_rt, b_start, b_stop, b_irt, b_im, b_q, b_gq, b_pgq, b_pep, b_ev;
    // NaN means absent and is written as null, as the library writer does.
    auto number = [&ok](arrow::DoubleBuilder& b, double v) {
      if (std::isfinite(v)) { ok(b.Append(v)); } else { ok(b.AppendNull()); }
    };
    for (const auto& r : rows)
    {
      ok(b_run.Append(run));
      ok(b_id.Append(r.precursor_id));
      ok(b_seq.Append(r.modified_sequence));
      ok(b_z.Append(r.charge));
      number(b_mz, r.precursor_mz);
      ok(b_pg.Append(r.protein_group));
      ok(b_decoy.Append(r.decoy ? 1 : 0));
      number(b_rt, r.rt);
      number(b_start, r.rt_start);
      number(b_stop, r.rt_stop);
      number(b_irt, r.irt);
      number(b_im, r.im);
      number(b_q, r.q);
      number(b_gq, r.global_q);
      number(b_pgq, r.pg_q);
      number(b_pep, r.pep);
      number(b_ev, r.evidence);
    }
    std::vector<std::shared_ptr<arrow::Array>> arrays(columns().size());
    ok(b_run.Finish(&arrays[0]));   ok(b_id.Finish(&arrays[1]));    ok(b_seq.Finish(&arrays[2]));
    ok(b_z.Finish(&arrays[3]));     ok(b_mz.Finish(&arrays[4]));    ok(b_pg.Finish(&arrays[5]));
    ok(b_decoy.Finish(&arrays[6])); ok(b_rt.Finish(&arrays[7]));    ok(b_start.Finish(&arrays[8]));
    ok(b_stop.Finish(&arrays[9]));  ok(b_irt.Finish(&arrays[10]));  ok(b_im.Finish(&arrays[11]));
    ok(b_q.Finish(&arrays[12]));    ok(b_gq.Finish(&arrays[13]));   ok(b_pgq.Finish(&arrays[14]));
    ok(b_pep.Finish(&arrays[15]));  ok(b_ev.Finish(&arrays[16]));

    std::vector<std::shared_ptr<arrow::Field>> fields;
    for (std::size_t c = 0; c < arrays.size(); ++c) { fields.push_back(arrow::field(columns()[c], arrays[c]->type())); }
    std::vector<std::string> keys, values;
    for (const auto& kv : metadata) { keys.push_back(kv.first); values.push_back(kv.second); }
    auto schema = arrow::schema(fields)->WithMetadata(arrow::key_value_metadata(keys, values));
    auto table = arrow::Table::Make(schema, arrays);

    AtomicFile output(path);
    auto outfile = arrow::io::FileOutputStream::Open(output.temporaryPath().string());
    if (!outfile.ok()) { throw std::runtime_error("cannot write identification report: " + path); }
    auto props = parquet::WriterProperties::Builder().compression(parquet::Compression::ZSTD)->build();
    // store_schema() carries the key-value metadata into the file.
    auto arrow_props = parquet::ArrowWriterProperties::Builder().store_schema()->build();
    ok(parquet::arrow::WriteTable(*table, arrow::default_memory_pool(), *outfile, 1 << 20, props, arrow_props));
    ok((*outfile)->Close());
    output.commit();
  }
}
