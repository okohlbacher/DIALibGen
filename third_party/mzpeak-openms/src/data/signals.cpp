/*

This file is part of the mzpeak.h project.  It is subject to the
license specified in the LICENSE file which can be found in the
top-level directory of this repository.

*/

#include "mzpeak/data/signals.h"

#include <algorithm>
#include <arrow/record_batch.h>
#include <bit>
#include <chrono>
#include <memory>
#include <optional>
#include <parquet/arrow/reader.h>
#include <parquet/statistics.h>
#include <utility>

#include "mzpeak/data/array_index.h"
#include "mzpeak/exception.h"
#include "mzpeak/schema/entity_type.h"
#include "mzpeak/util/executor.h"
#include "mzpeak/util/planner.h"
#include "mzpeak/util/projection.h"

namespace MzPeak::Data {

/******************************************************************************/
struct Signals::Impl {
  Impl(std::unique_ptr<Util::Parquet> parquet)
      : parquet_(std::move(parquet))
      , array_index_(parse_array_index())
  {
  }

  std::shared_ptr<ArrayIndex> parse_array_index() const;

  std::unique_ptr<Util::Parquet> parquet_;
  std::shared_ptr<ArrayIndex> array_index_;
};

/******************************************************************************/
std::shared_ptr<ArrayIndex> Signals::Impl::parse_array_index() const
{
  Util::Parquet::file_metadata_t fmd(parquet_->file_metadata());
  EntityType entity_type = parquet_->index_file().entity_type();

  std::string num_key(entity_type.metadata_count_key());
  std::optional<std::size_t> num_entities(parquet_->kv_size_t(fmd, num_key));

  std::string index_key(entity_type.array_index_name());
  auto index_str(parquet_->kv_string(fmd, index_key));
  if (!index_str.has_value()) throw InvalidFormatError("missing array_index");

  namespace json = boost::json;
  boost::system::error_code ec;
  json::value v(json::parse(index_str.value()));
  if (ec) throw MzPeak::JsonError(ec.message());

  if (v.is_object()) {
    auto ai = std::make_shared<ArrayIndex>(entity_type, v.as_object());
    ai->num_entities(num_entities);
    return ai;
  } else {
    throw JsonError("array_index should be a JSON object");
  }
}

/******************************************************************************/
Signals::Signals(std::unique_ptr<Util::Parquet> parquet)
    : impl_(std::make_unique<Impl>(std::move(parquet)))
{
}

/******************************************************************************/
Signals::~Signals() = default;

/******************************************************************************/
const std::shared_ptr<ArrayIndex>& Signals::array_index() const
{
  return impl_->array_index_;
}

/******************************************************************************/
std::size_t Signals::record_count() const
{
  auto ne(impl_->array_index_->num_entities());
  if (ne.has_value()) return *ne;

  // Fall back to scanning row-group max statistics on the entity-index column.
  // The index is 0-based and contiguous so COUNT = max_index + 1.
  const std::string prefix = impl_->array_index_->prefix();
  const std::string index_col = [&]() -> std::string {
    switch (impl_->array_index_->entity_type().type().value_or(
        Schema::EntityType::Spectrum)) {
    case Schema::EntityType::Chromatogram:
      return "chromatogram_index";
    case Schema::EntityType::WavelengthSpectrum:
      return "wavelength_spectrum_index";
    default:
      return "spectrum_index";
    }
  }();

  auto dest = impl_->parquet_->field(prefix, index_col);
  if (!dest.has_value()) {
    throw ParquetError("record_count: missing " + prefix + "." + index_col +
                       " column");
  }
  const int col = dest->second->absolute_index();
  auto fmd = impl_->parquet_->file_metadata();
  if (fmd->num_row_groups() == 0) return 0;

  int64_t max_signed = -1;
  for (int g = 0; g < fmd->num_row_groups(); ++g) {
    auto chunk = fmd->RowGroup(g)->ColumnChunk(col);
    if (!chunk->is_stats_set()) continue;
    auto stats = chunk->statistics();
    if (!stats || !stats->HasMinMax()) continue;
    auto typed = std::dynamic_pointer_cast<parquet::Int64Statistics>(stats);
    if (!typed) continue;
    max_signed = std::max(max_signed, typed->max());
  }

  if (max_signed < 0) {
    if (fmd->num_rows() == 0) return 0;
    throw ParquetError("record_count: no statistics for " + prefix + "." +
                       index_col);
  }
  // The index column is unsigned 64-bit stored with a signed physical type.
  return static_cast<std::size_t>(std::bit_cast<uint64_t>(max_signed)) + 1;
}

/******************************************************************************/
std::optional<Schema::Column> Signals::column(std::string_view name) const
{
  return impl_->parquet_->field(impl_->array_index_->prefix(), name);
}

/******************************************************************************/
std::optional<Schema::Column> Signals::column(const ArrayIndex::Entry& entry) const
{
  return impl_->array_index_->entry_column(*impl_->parquet_->groups(), entry);
}

/******************************************************************************/
std::optional<Schema::Column> Signals::column(const ArrayIndex::Dimension& dim,
                                              Schema::BufferFormat format) const
{
  return dim.entry_with(format).and_then(
      [this](const auto& entry) { return column(entry); });
}

/******************************************************************************/
const std::shared_ptr<Schema::GroupMap>& Signals::groups() const
{
  return impl_->parquet_->groups();
}

/******************************************************************************/
Util::Query::Builder Signals::index() const
{
  auto entity_type = impl_->array_index_->entity_type();
  auto field_name = entity_type.index_column_name();
  auto index_field = column(field_name);

  if (!index_field.has_value()) {
    throw InvalidFormatError("parquet file is missing the index column: " +
                             field_name);
  }

  return Util::Query::Builder(index_field.value());
}

/******************************************************************************/
std::unique_ptr<Util::Slice>
Signals::select(const std::vector<ArrayIndex::Dimension>& projection,
                const Util::Query& query)
{
  Util::Projection columns;

  for (const auto& dim : projection) {
    for (const auto& entry : dim.entries) {
      if (!entry.needed_for_decoding()) continue;

      auto field =
          impl_->array_index_->entry_column(*impl_->parquet_->groups(), entry);
      if (!field.has_value()) {
        throw InvalidFormatError("array entry not present in schema: " + entry.name);
      } else {
        columns.project(field.value());
      }
    }
  }

#ifdef MZPEAK_READ_COUNTERS
  const auto t_plan0 = std::chrono::steady_clock::now();
#endif
  Util::Planner planner = impl_->parquet_->planner(query);
#ifdef MZPEAK_READ_COUNTERS
  Util::count_ns_plan_ctor(std::chrono::duration_cast<std::chrono::nanoseconds>(
                               std::chrono::steady_clock::now() - t_plan0)
                               .count());
#endif
  auto plan = planner.plan();
#ifdef MZPEAK_READ_COUNTERS
  Util::count_ns_plan(std::chrono::duration_cast<std::chrono::nanoseconds>(
                          std::chrono::steady_clock::now() - t_plan0)
                          .count());
#endif

  Util::Executor executor = impl_->parquet_->executor(columns);
  return executor.execute(plan);
}

/******************************************************************************/
std::optional<std::vector<std::pair<uint64_t, std::unique_ptr<Util::Slice>>>>
Signals::select_range(const std::vector<ArrayIndex::Dimension>& projection,
                      uint64_t first,
                      uint64_t last)
{
  // The same columns select() projects, in the same order.
  std::vector<Schema::Column> fields;
  for (const auto& dim : projection) {
    for (const auto& entry : dim.entries) {
      if (!entry.needed_for_decoding()) continue;
      auto field =
          impl_->array_index_->entry_column(*impl_->parquet_->groups(), entry);
      if (!field.has_value()) {
        throw InvalidFormatError("array entry not present in schema: " + entry.name);
      }
      fields.push_back(field.value());
    }
  }

  // Plus the index column, which tells the rows apart.
  const auto index_name = impl_->array_index_->entity_type().index_column_name();
  const auto index_field = column(index_name);
  if (!index_field.has_value()) {
    throw InvalidFormatError("parquet file is missing the index column: " +
                             index_name);
  }
  const auto key = [](const Schema::Column& c) {
    return c.second->absolute_index();
  };
  const bool index_projected = std::ranges::any_of(
      fields, [&](const auto& f) { return key(f) == key(*index_field); });

  Util::Projection columns;
  for (const auto& f : fields)
    columns.project(f);
  if (!index_projected) columns.project(*index_field);

  const Util::Query query = index().ge(first).and_then(index().le(last));
  Util::Planner planner = impl_->parquet_->planner(query);
  auto plan = planner.plan();

  // Sorted groups held as KeyRuns already say where each spectrum's rows are:
  // split by them, without the index column and without a per-row filter.
  if (!index_projected) {
    Util::Projection decode_only;
    for (const auto& f : fields)
      decode_only.project(f);
    Util::Executor by_runs = impl_->parquet_->executor(decode_only);
    if (auto split = by_runs.execute_runs(plan, first, last)) return split;
  }

  Util::Executor executor = impl_->parquet_->executor(columns);
  std::unique_ptr<Util::Slice> whole = executor.execute(plan);

  std::vector<std::pair<uint64_t, std::unique_ptr<Util::Slice>>> out;
  if (!whole) return out;

  // raw() hands each column over once.  The executor appends one array per
  // field per record batch, so chunk c of every field covers the same rows.
  std::shared_ptr<Util::Slice::Raw> ids = whole->raw(*index_field);
  if (!ids) return out;
  std::vector<std::shared_ptr<Util::Slice::Raw>> raws;
  for (const auto& f : fields) {
    raws.push_back(key(f) == key(*index_field) ? ids : whole->raw(f));
    if (!raws.back() || raws.back()->size() != ids->size()) {
      throw ParquetError(
          "range read: column chunks do not line up with the index column");
    }
  }

  const auto id_at = [](const std::shared_ptr<arrow::Array>& a,
                        int64_t r) -> std::optional<uint64_t> {
    if (a->IsNull(r)) return std::nullopt;
    switch (a->type_id()) {
    case arrow::Type::UINT64:
      return std::static_pointer_cast<arrow::UInt64Array>(a)->Value(r);
    case arrow::Type::INT64:
      return static_cast<uint64_t>(
          std::static_pointer_cast<arrow::Int64Array>(a)->Value(r));
    case arrow::Type::UINT32:
      return std::static_pointer_cast<arrow::UInt32Array>(a)->Value(r);
    case arrow::Type::INT32:
      return static_cast<uint64_t>(
          std::static_pointer_cast<arrow::Int32Array>(a)->Value(r));
    default:
      throw ParquetError("range read: unsupported index column type " +
                         a->type()->ToString());
    }
  };

  // Append rows [begin, begin + length) of chunk c to the slice for spectrum id;
  // a spectrum that continues into the next batch keeps its slice.
  const auto emit = [&](uint64_t id, std::size_t c, int64_t begin, int64_t length) {
    if (out.empty() || out.back().first != id) {
      out.emplace_back(id, std::unique_ptr<Util::Slice>(new Util::Slice(fields)));
    }
    for (std::size_t f = 0; f < fields.size(); ++f) {
      out.back().second->append(fields[f], (*raws[f])[c]->Slice(begin, length));
    }
  };

  for (std::size_t c = 0; c < ids->size(); ++c) {
    const auto& chunk = (*ids)[c];
    int64_t begin = 0;
    std::optional<uint64_t> current;
    for (int64_t r = 0; r < chunk->length(); ++r) {
      const auto id = id_at(chunk, r);
      if (id != current) {
        if (current.has_value()) emit(*current, c, begin, r - begin);
        current = id;
        begin = r;
      }
    }
    if (current.has_value()) emit(*current, c, begin, chunk->length() - begin);
  }

  // Rows sorted by index give one slice per spectrum.  Unsorted rows would
  // split a spectrum into several slices; the format allows them, so report it
  // and let the caller read per spectrum rather than hand out fragments.
  for (std::size_t i = 1; i < out.size(); ++i) {
    if (out[i].first <= out[i - 1].first) return std::nullopt;
  }
  return out;
}

} // namespace MzPeak::Data
