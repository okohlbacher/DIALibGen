// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <odia/LibraryExport.h>
#include <odia/AtomicFile.h>
#include <odia/DIANNLibraryFile.h>
#include <odia/ModifiedSequence.h>
#include <odia/Trace.h>

#include "Parallel.h"

#include <arrow/api.h>
#include <arrow/io/file.h>
#include <parquet/arrow/writer.h>
#include <parquet/properties.h>

#include <algorithm>
#include <atomic>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string_view>

namespace ODIA
{
  namespace
  {
    struct FormatName { LibraryFormat format; const char* name; const char* extension; };
    constexpr FormatName FORMATS[] = {
      {LibraryFormat::DIALibGenParquet, "dialibgen_parquet", ".parquet"},
      {LibraryFormat::DIANNTSV, "diann_tsv", ".tsv"},
      {LibraryFormat::DIANN2Parquet, "diann2_parquet", ".parquet"},
    };

    std::string lowerExtension(const std::string& filename)
    {
      auto extension = std::filesystem::path(filename).extension().string();
      std::transform(extension.begin(), extension.end(), extension.begin(),
                     [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
      return extension;
    }

    void ok(const arrow::Status& status)
    {
      if (!status.ok()) { throw std::runtime_error("Parquet: " + status.ToString()); }
    }

    /// A member of a protein group, "sp|P02768|ALBU_HUMAN", as accession and
    /// entry name; a member without the db|accession|name form is both.
    std::pair<std::string_view, std::string_view> accessionAndName(std::string_view member)
    {
      const auto bar1 = member.find('|');
      const auto bar2 = bar1 == std::string_view::npos ? bar1 : member.find('|', bar1 + 1);
      if (bar2 == std::string_view::npos) { return {member, member}; }
      return {member.substr(bar1 + 1, bar2 - bar1 - 1), member.substr(bar2 + 1)};
    }

    std::vector<std::string_view> members(std::string_view group)
    {
      std::vector<std::string_view> out;
      std::size_t pos = 0;
      while (pos < group.size())
      {
        std::size_t end = group.find(';', pos);
        if (end == std::string_view::npos) { end = group.size(); }
        if (end > pos) { out.push_back(group.substr(pos, end - pos)); }
        pos = end + 1;
      }
      return out;
    }

    std::string strippedSequence(std::string_view modified)
    {
      std::string out;
      int depth = 0;
      for (const char c : modified)
      {
        if (c == '(' || c == '[') { ++depth; continue; }
        if (c == ')' || c == ']') { --depth; continue; }
        if (depth == 0 && std::isupper(static_cast<unsigned char>(c))) { out.push_back(c); }
      }
      return out;
    }

    // ------------------------------------------------------- DIA-NN 2 Parquet

    std::shared_ptr<arrow::Schema> diann2Schema(const std::vector<std::pair<std::string, std::string>>& metadata)
    {
      const auto s = arrow::utf8();
      const auto i = arrow::int64();
      const auto f = arrow::float32();
      std::vector<std::shared_ptr<arrow::Field>> fields;
      for (const auto& [name, type] : std::vector<std::pair<const char*, std::shared_ptr<arrow::DataType>>>{
             {"Precursor.Id", s}, {"Modified.Sequence", s}, {"Stripped.Sequence", s},
             {"Precursor.Charge", i}, {"Proteotypic", i}, {"Decoy", i}, {"N.Term", i}, {"C.Term", i},
             {"RT", f}, {"IM", f}, {"Q.Value", f}, {"Peptidoform.Q.Value", f}, {"PTM.Site.Confidence", f},
             {"PG.Q.Value", f}, {"Precursor.Mz", f}, {"Product.Mz", f}, {"Relative.Intensity", f},
             {"Fragment.Type", s}, {"Fragment.Charge", i}, {"Fragment.Series.Number", i},
             {"Fragment.Loss.Type", s}, {"Exclude.From.Quant", i},
             {"Protein.Ids", s}, {"Protein.Group", s}, {"Protein.Names", s}, {"Genes", s}, {"Flags", i}})
      { fields.push_back(arrow::field(name, type, /*nullable=*/false)); }
      std::shared_ptr<arrow::KeyValueMetadata> kv;
      if (!metadata.empty())
      {
        kv = std::make_shared<arrow::KeyValueMetadata>();
        for (const auto& [key, value] : metadata) { kv->Append(key, value); }
      }
      return arrow::schema(fields, kv);
    }

    /// One row group: the transitions of precursors [first, last).
    std::shared_ptr<arrow::Table> diann2RowGroup(const Library& library, const std::shared_ptr<arrow::Schema>& schema,
                                                 std::size_t first, std::size_t last, std::size_t& unknown_mods)
    {
      const auto& p = library.precursors();
      const auto& t = library.transitions();
      arrow::StringBuilder id, seq, stripped, ftype, loss, ids, group, names, genes;
      arrow::Int64Builder charge, proteotypic, decoy, n_term, c_term, fcharge, ordinal, exclude, flags;
      arrow::FloatBuilder rt, im, q, pq, ptm, pgq, pmz, qmz, intensity;
      std::size_t rows = 0;
      for (std::size_t i = first; i < last; ++i) { rows += p.transition_count[i]; }
      for (auto* b : {&charge, &proteotypic, &decoy, &n_term, &c_term, &fcharge, &ordinal, &exclude, &flags})
      { ok(b->Reserve(static_cast<std::int64_t>(rows))); }
      for (auto* b : {&rt, &im, &q, &pq, &ptm, &pgq, &pmz, &qmz, &intensity})
      { ok(b->Reserve(static_cast<std::int64_t>(rows))); }

      std::string precursor_id, accessions, entry_names;
      for (std::size_t i = first; i < last; ++i)
      {
        const auto stored = library.strings().get(p.modified_sequence[i]);
        const std::string modified = diannModifiedSequence(stored, &unknown_mods);
        const std::string plain = strippedSequence(modified);
        const int z = p.charge[i];
        // As in the TSV: DIA-NN gives a decoy its target's sequence and charge,
        // so the id alone would merge the pair on reading.
        precursor_id = modified + std::to_string(z) + (p.decoy[i] ? "_decoy" : "");
        accessions.clear(); entry_names.clear();
        const auto list = members(library.strings().get(p.protein_group[i]));
        for (std::size_t k = 0; k < list.size(); ++k)
        {
          const auto [accession, entry] = accessionAndName(list[k]);
          if (k) { accessions += ';'; entry_names += ';'; }
          accessions += accession; entry_names += entry;
        }
        const bool has_rt = !std::isnan(p.irt[i]);
        const bool has_im = i < p.im.size() && !std::isnan(p.im[i]);
        for (std::uint32_t k = 0; k < p.transition_count[i]; ++k)
        {
          const std::uint32_t j = p.transition_begin[i] + k;
          ok(id.Append(precursor_id)); ok(seq.Append(modified)); ok(stripped.Append(plain));
          charge.UnsafeAppend(z); proteotypic.UnsafeAppend(list.size() == 1 ? 1 : 0);
          decoy.UnsafeAppend(p.decoy[i] ? 1 : 0); n_term.UnsafeAppend(0); c_term.UnsafeAppend(0);
          rt.UnsafeAppend(has_rt ? p.irt[i] : 0.0f); im.UnsafeAppend(has_im ? p.im[i] : 0.0f);
          q.UnsafeAppend(0.0f); pq.UnsafeAppend(0.0f); ptm.UnsafeAppend(1.0f); pgq.UnsafeAppend(0.0f);
          pmz.UnsafeAppend(static_cast<float>(fromFixed(p.mz[i])));
          qmz.UnsafeAppend(static_cast<float>(fromFixed(t.product_mz[j])));
          intensity.UnsafeAppend(std::isnan(t.library_intensity[j]) ? 0.0f : t.library_intensity[j]);
          ok(ftype.Append(std::string(toString(t.type[j]))));
          fcharge.UnsafeAppend(t.charge[j]); ordinal.UnsafeAppend(t.ordinal[j]);
          ok(loss.Append(std::string(toString(t.loss[j]))));
          exclude.UnsafeAppend(0);
          ok(ids.Append(accessions)); ok(group.Append(accessions)); ok(names.Append(entry_names));
          ok(genes.Append("")); flags.UnsafeAppend(0);
        }
      }
      std::vector<std::shared_ptr<arrow::Array>> arrays;
      const auto finish = [&](arrow::ArrayBuilder& b) {
        std::shared_ptr<arrow::Array> a;
        ok(b.Finish(&a));
        arrays.push_back(std::move(a));
      };
      for (arrow::ArrayBuilder* b : std::initializer_list<arrow::ArrayBuilder*>{
             &id, &seq, &stripped, &charge, &proteotypic, &decoy, &n_term, &c_term,
             &rt, &im, &q, &pq, &ptm, &pgq, &pmz, &qmz, &intensity,
             &ftype, &fcharge, &ordinal, &loss, &exclude, &ids, &group, &names, &genes, &flags})
      { finish(*b); }
      return arrow::Table::Make(schema, arrays, static_cast<std::int64_t>(rows));
    }

  } // namespace

  std::vector<std::string> LibraryExport::formatNames()
  {
    std::vector<std::string> out{"auto"};
    for (const auto& f : FORMATS) { out.emplace_back(f.name); }
    return out;
  }

  const char* LibraryExport::name(LibraryFormat format)
  {
    for (const auto& f : FORMATS) { if (f.format == format) { return f.name; } }
    return "?";
  }

  LibraryFormat LibraryExport::resolve(const std::string& out, const std::string& requested)
  {
    const std::string extension = lowerExtension(out);
    if (requested.empty() || requested == "auto")
    {
      if (extension == ".parquet") { return LibraryFormat::DIALibGenParquet; }
      if (extension == ".tsv") { return LibraryFormat::DIANNTSV; }
      throw std::invalid_argument("-out must end in .parquet or .tsv (got '" + out + "')");
    }
    for (const auto& f : FORMATS)
    {
      if (requested != f.name) { continue; }
      if (extension != f.extension)
      {
        throw std::invalid_argument("-out_format " + requested + " writes a " + f.extension +
                                    " file; -out is '" + out + "'");
      }
      return f.format;
    }
    std::string known;
    for (const auto& n : formatNames()) { known += (known.empty() ? "" : ", ") + n; }
    throw std::invalid_argument("unknown -out_format '" + requested + "' (known: " + known + ")");
  }

  bool LibraryExport::carriesRecipe(LibraryFormat format)
  {
    return format == LibraryFormat::DIALibGenParquet || format == LibraryFormat::DIANN2Parquet;
  }

  void LibraryExport::storeDIANN2Parquet(const std::string& filename, const Library& library,
                                         const std::vector<std::pair<std::string, std::string>>& metadata,
                                         unsigned threads)
  {
    trace::Scope write_step("write_diann2_parquet", 1);
    write_step.fields().rows = static_cast<std::int64_t>(library.transitionCount());
    const auto schema = diann2Schema(metadata);
    AtomicFile output(filename);
    auto opened = arrow::io::FileOutputStream::Open(output.temporaryPath().string());
    if (!opened.ok()) { throw std::runtime_error("cannot write library: " + filename); }
    const auto outfile = *opened;
    // ZSTD and per-row-group dictionaries for the strings each precursor
    // repeats on every fragment row, as DIALibGen's own layout; DIA-NN 2.0
    // reads the file through Arrow either way.
    const auto properties = parquet::WriterProperties::Builder()
                              .compression(parquet::Compression::ZSTD)
                              ->enable_dictionary()
                              ->build();
    auto writer_result = parquet::arrow::FileWriter::Open(*schema, arrow::default_memory_pool(), outfile, properties,
                                                          parquet::default_arrow_writer_properties());
    if (!writer_result.ok()) { throw std::runtime_error("cannot write Parquet: " + writer_result.status().ToString()); }
    auto writer = std::move(*writer_result);
    // Without store_schema (DIA-NN's own file has no ARROW:schema), the
    // schema's metadata reaches the file only this way.
    if (schema->metadata()) { ok(writer->AddKeyValueMetadata(schema->metadata())); }

    // Row groups of 65,536 precursors (~800k rows), built a wave at a time in
    // parallel and written in order: the file does not depend on threads.
    constexpr std::size_t GROUP = std::size_t(1) << 16;
    const std::size_t n = library.precursorCount();
    const std::size_t groups = (n + GROUP - 1) / GROUP;
    const unsigned workers = std::max(1u, std::min(parallel::workers(threads), 8u));
    std::atomic<std::size_t> unknown_mods{0};
    for (std::size_t wave = 0; wave < groups; wave += workers)
    {
      const std::size_t count = std::min<std::size_t>(workers, groups - wave);
      std::vector<std::shared_ptr<arrow::Table>> tables(count);
      parallel::forEach(count, workers, [&](std::size_t k)
      {
        const std::size_t first = (wave + k) * GROUP;
        std::size_t unknown = 0;
        tables[k] = diann2RowGroup(library, schema, first, std::min(n, first + GROUP), unknown);
        unknown_mods += unknown;
      });
      for (const auto& table : tables)
      {
        if (table->num_rows() == 0) { continue; }
        ok(writer->WriteTable(*table, table->num_rows()));
      }
    }
    if (n == 0)
    {
      // An empty library is still a readable file with the schema.
      std::size_t unknown = 0;
      const auto empty = diann2RowGroup(library, schema, 0, 0, unknown);
      ok(writer->WriteTable(*empty, 1));
    }
    ok(writer->Close());
    ok(outfile->Close());
    output.commit();
    if (unknown_mods.load())
    {
      std::cerr << "warning: " << unknown_mods.load() << " modification tokens have no UniMod accession and "
                   "were written verbatim; DIA-NN may not read those precursors\n";
    }
  }

  void LibraryExport::store(const std::string& filename, const Library& library, LibraryFormat format,
                            const std::vector<std::pair<std::string, std::string>>& metadata, unsigned threads)
  {
    switch (format)
    {
      case LibraryFormat::DIALibGenParquet:
        DIANNLibraryFile::storeParquetCompact(filename, library, DIANNLibraryFile::Fingerprint{}, "", metadata, threads);
        return;
      case LibraryFormat::DIANNTSV: DIANNLibraryFile::storeTSV(filename, library, threads); return;
      case LibraryFormat::DIANN2Parquet: storeDIANN2Parquet(filename, library, metadata, threads); return;
    }
  }

} // namespace ODIA
