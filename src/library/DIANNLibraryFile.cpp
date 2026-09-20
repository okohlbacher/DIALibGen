// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <odia/DIANNLibraryFile.h>
#include <odia/AtomicFile.h>
#include <odia/TextWriter.h>

#include <arrow/api.h>
#include <arrow/compute/api.h>
#include <arrow/io/file.h>
#include <parquet/arrow/reader.h>
#include <parquet/arrow/writer.h>
#include <parquet/file_reader.h>
#include <iomanip>
#include <filesystem>
#include <arrow/builder.h>
#include <arrow/table.h>

#include <charconv>
#include <algorithm>
#include <cctype>
#include <cstring>
#include <cmath>
#include <fstream>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

namespace ODIA
{
  namespace
  {
    bool isParquet(const std::string& filename)
    {
      auto extension = std::filesystem::path(filename).extension().string();
      std::transform(extension.begin(), extension.end(), extension.begin(),
                     [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
      if (extension == ".parquet") { return true; }
      if (extension == ".tsv") { return false; }
      throw std::runtime_error("unsupported library extension (expected .tsv or .parquet): " + filename);
    }

    double toDouble(std::string_view s, const char* column)
    {
      double v = 0.0;
      if (s.empty()) { return std::nan(""); }
      const auto* first = s.data();
      const auto* last = first + s.size();
      if (*first == '+') { ++first; }
      const auto parsed = std::from_chars(first, last, v, std::chars_format::general);
      if (parsed.ec != std::errc{} || parsed.ptr != last || std::isinf(v))
      { throw std::runtime_error(std::string("invalid numeric value in ") + column + ": " + std::string(s)); }
      return v;
    }

    long toLong(std::string_view s, const char* column)
    {
      if (s.empty()) { return 0; }
      long v = 0;
      const auto* first = s.data();
      const auto* last = s.data() + s.size();
      while (first != last && (*first == ' ' || *first == '+')) { ++first; }
      while (last != first && last[-1] == ' ') { --last; }
      const auto parsed = std::from_chars(first, last, v);
      if (parsed.ec != std::errc{} || parsed.ptr != last)
      { throw std::runtime_error(std::string("invalid integer value in ") + column + ": " + std::string(s)); }
      return v;
    }

    /// Residues in a modified sequence, ignoring bracketed modification names.
    ///
    /// The name must not be scanned for residues: "UniMod" contains an M.
    std::size_t residueCount(std::string_view seq)
    {
      std::size_t n = 0;
      for (std::size_t i = 0; i < seq.size(); ++i)
      {
        if (seq[i] == '(' || seq[i] == '[')
        {
          const char close = seq[i] == '(' ? ')' : ']';
          while (i < seq.size() && seq[i] != close) { ++i; }
          continue;
        }
        if (seq[i] >= 'A' && seq[i] <= 'Z') { ++n; }
      }
      return n;
    }

    /// One parsed row, in the order the builder needs it.
    struct Row
    {
      std::string_view precursor_id;
      std::string_view modified_sequence;
      std::string_view protein_group;
      double precursor_mz = 0, product_mz = 0;
      double rt = 0, im = 0, ccs = 0, intensity = 0;
      long precursor_charge = 0, fragment_charge = 0, ordinal = 0, decoy = 0;
      std::string_view fragment_type, loss_type;
    };

    /// Appends rows to the library, starting a new precursor when the
    /// Precursor.Id changes.
    class Builder
    {
    public:
      explicit Builder(Library& lib) : lib_(lib) {}

      void append(const Row& r)
      {
        if (r.precursor_charge < 1 || r.precursor_charge > 255)
        { throw std::runtime_error("invalid precursor charge: expected an integer between 1 and 255"); }
        if (!have_current_ || r.precursor_id != current_id_ || (r.decoy != 0) != current_decoy_)
        {
          closeCurrent();
          // A precursor's rows must be contiguous, because grouping keys off a
          // change of Precursor.Id. If they are not, the same id reappears and
          // one precursor silently becomes several single-transition ones with
          // a plausible-looking count -- which is what any library that has been
          // sorted by m/z or concatenated will do.
          // Hash rather than store the id: a set<pair<string,bool>> measured
          // 128 B per precursor, ~0.9 GB at the scale D3 is written against.
          std::size_t key = std::hash<std::string_view>{}(r.precursor_id);
          key = key * 31 + (r.decoy != 0 ? 1u : 0u);
          if (!seen_.insert(key).second) { ++reopened_; }
          startPrecursor(r);
        }
        // Reject fragments that cannot exist, rather than storing them.
        //
        // An ordinal past the stored width wraps into a *plausible* ordinal --
        // 300 becomes 44 -- so the row would survive every downstream check as
        // a different, wrong fragment. Ordinal 0, an ordinal at or beyond the
        // peptide length, and a non-positive charge are not fragments at all.
        if (r.ordinal <= 0 || r.ordinal > 255 ||
            r.fragment_charge <= 0 || r.fragment_charge > 127 ||
            static_cast<std::size_t>(r.ordinal) >= residueCount(r.modified_sequence))
        {
          ++out_of_range_;
          return;
        }

        auto& t = lib_.transitions();
        t.product_mz.push_back(toFixed(r.product_mz));
        t.library_intensity.push_back(static_cast<float>(r.intensity));
        t.type.push_back(parseFragmentType(r.fragment_type));
        t.ordinal.push_back(static_cast<std::uint8_t>(r.ordinal));
        t.charge.push_back(static_cast<std::int8_t>(r.fragment_charge));
        t.loss.push_back(parseLossType(r.loss_type));
        ++count_;
      }

      void finish() { closeCurrent(); }

      /// Number of precursors whose rows were not contiguous.
      std::size_t reopened() const { return reopened_; }

      /// Transitions dropped because an ordinal or charge could not be stored.
      std::size_t outOfRange() const { return out_of_range_; }

    private:
      void startPrecursor(const Row& r)
      {
        auto& p = lib_.precursors();
        p.mz.push_back(toFixed(r.precursor_mz));
        p.irt.push_back(static_cast<float>(r.rt));
        p.im.push_back(r.im > 0.0 ? static_cast<float>(r.im) : std::nanf(""));
        // CCS was written by both writers and read by NEITHER reader, so a
        // library round-tripped through ODIA lost it silently: the generator
        // predicts it, storeTSV writes it, and every search that reloaded that
        // file saw an empty column. Measured on a round trip -- CCS 627.119385
        // went in and an empty field came out.
        //
        // It matters because CCS is the ONLY mobility information a generated
        // library carries: IM is left unset (NaN) because 1/K0 needs the drift
        // gas and the instrument calibration to derive. Dropping CCS therefore
        // discards the mobility dimension entirely on reload.
        p.ccs.push_back(r.ccs > 0.0 ? static_cast<float>(r.ccs) : std::nanf(""));
        p.charge.push_back(static_cast<std::uint8_t>(r.precursor_charge));
        p.decoy.push_back(r.decoy != 0 ? std::uint8_t{1} : std::uint8_t{0});
        p.modified_sequence.push_back(lib_.strings().intern(r.modified_sequence));
        p.protein_group.push_back(lib_.strings().intern(r.protein_group));
        p.transition_begin.push_back(static_cast<std::uint32_t>(lib_.transitionCount()));
        p.transition_count.push_back(0);
        // Appending breaks any m/z ordering; leaving the flag set would make
        // lowerBound binary-search an unsorted array and return an index
        // indistinguishable from a real answer.
        lib_.markUnsorted();
        current_id_owned_.assign(r.precursor_id);
        current_id_ = current_id_owned_;
        current_decoy_ = r.decoy != 0;
        have_current_ = true;
        count_ = 0;
      }

      void closeCurrent()
      {
        if (!have_current_) { return; }
        lib_.precursors().transition_count.back() = static_cast<std::uint32_t>(count_);
      }

      Library& lib_;
      std::unordered_set<std::size_t> seen_;
      std::size_t reopened_ = 0;
      std::size_t out_of_range_ = 0;
      std::string current_id_owned_;
      std::string_view current_id_;
      bool have_current_ = false;
      bool current_decoy_ = false;
      std::size_t count_ = 0;
    };
  } // namespace

  void DIANNLibraryFile::completeMobility(Library& library)
  {
    auto& p = library.precursors();
    const std::size_t n = library.precursorCount();
    p.im.resize(n, std::nanf(""));
    p.ccs.resize(n, std::nanf(""));
    for (std::size_t i = 0; i < n; ++i)
    {
      const bool has_im = !std::isnan(p.im[i]);
      const bool has_ccs = !std::isnan(p.ccs[i]);
      if (has_im == has_ccs) { continue; }        // both, or neither: nothing to do
      const double mz = fromFixed(p.mz[i]);
      const int z = p.charge[i];
      if (has_ccs) { p.im[i] = static_cast<float>(mobilityFromCCS(p.ccs[i], mz, z)); }
      else         { p.ccs[i] = static_cast<float>(ccsFromMobility(p.im[i], mz, z)); }
    }
  }

  void DIANNLibraryFile::load(const std::string& filename, Library& library)
  {
    if (isParquet(filename)) { loadParquet(filename, library); }
    else { loadTSV(filename, library); }
    // A library that carries only one of the two mobility quantities gets the
    // other derived, so a consumer never has to know which the producer chose.
    // Measured need: our own generated library reached DIA-NN with 1/K0 absent
    // for every precursor, and DIA-NN's `iIM` was 0 on 100% of a diaPASEF run's
    // identifications (doc/32).
    completeMobility(library);
    // The interning index answered "have I seen this string" while reading and
    // is never consulted again -- handles resolve through the arena's entries.
    // Held for the library's lifetime it is one hash node per distinct string
    // plus a bucket array, for a question nobody asks twice.
    library.strings().releaseLookup();
  }

  void DIANNLibraryFile::loadTSV(const std::string& filename, Library& library)
  {
    std::ifstream in(filename);
    if (!in) { throw std::runtime_error("cannot open library: " + filename); }

    std::string header;
    if (!std::getline(in, header)) { throw std::runtime_error("empty library: " + filename); }
    if (!header.empty() && header.back() == '\r') { header.pop_back(); }
    if (header.starts_with("\xEF\xBB\xBF")) { header.erase(0, 3); }

    std::unordered_map<std::string, int> index;
    {
      int i = 0;
      std::size_t pos = 0;
      while (pos <= header.size())
      {
        const std::size_t tab = header.find('\t', pos);
        const std::size_t end = tab == std::string::npos ? header.size() : tab;
        const auto name = header.substr(pos, end - pos);
        if (!index.emplace(name, i++).second)
        { throw std::runtime_error("duplicate TSV column: " + name); }
        if (tab == std::string::npos) { break; }
        pos = tab + 1;
      }
    }

    auto col = [&](const char* name) {
      const auto it = index.find(name);
      return it == index.end() ? -1 : it->second;
    };

    const int c_id = col(Columns::PRECURSOR_ID);
    const int c_seq = col(Columns::MODIFIED_SEQUENCE);
    const int c_pg = col(Columns::PROTEIN_GROUP);
    const int c_pmz = col(Columns::PRECURSOR_MZ);
    const int c_qmz = col(Columns::PRODUCT_MZ);
    const int c_rt = col(Columns::RT);
    const int c_im = col(Columns::IM);
    const int c_ccs = col(Columns::CCS);
    const int c_int = col(Columns::RELATIVE_INTENSITY);
    const int c_z = col(Columns::PRECURSOR_CHARGE);
    const int c_fz = col(Columns::FRAGMENT_CHARGE);
    const int c_ord = col(Columns::FRAGMENT_SERIES_NUMBER);
    const int c_dec = col(Columns::DECOY);
    const int c_ft = col(Columns::FRAGMENT_TYPE);
    const int c_lt = col(Columns::FRAGMENT_LOSS_TYPE);

    if (c_id < 0 || c_pmz < 0 || c_qmz < 0)
    {
      throw std::runtime_error("not a DIA-NN library (missing Precursor.Id / Precursor.Mz / "
                               "Product.Mz): " + filename);
    }

    Builder builder(library);
    std::string line;
    std::vector<std::string_view> f;
    std::size_t line_number = 1;
    while (std::getline(in, line))
    {
      ++line_number;
      if (!line.empty() && line.back() == '\r') { line.pop_back(); }
      if (line.empty()) { continue; }

      f.clear();
      std::size_t pos = 0;
      while (pos <= line.size())
      {
        const std::size_t tab = line.find('\t', pos);
        const std::size_t end = tab == std::string::npos ? line.size() : tab;
        f.emplace_back(line.data() + pos, end - pos);
        if (tab == std::string::npos) { break; }
        pos = tab + 1;
      }
      if (f.size() != index.size())
      { throw std::runtime_error("TSV row width differs from header at line " + std::to_string(line_number)); }

      auto get = [&](int i) -> std::string_view {
        return (i >= 0 && static_cast<std::size_t>(i) < f.size()) ? f[i] : std::string_view{};
      };

      Row r;
      r.precursor_id = get(c_id);
      r.modified_sequence = get(c_seq);
      r.protein_group = get(c_pg);
      r.precursor_mz = toDouble(get(c_pmz), Columns::PRECURSOR_MZ);
      r.product_mz = toDouble(get(c_qmz), Columns::PRODUCT_MZ);
      r.rt = toDouble(get(c_rt), Columns::RT);
      r.im = toDouble(get(c_im), Columns::IM);
      r.ccs = toDouble(get(c_ccs), Columns::CCS);
      r.intensity = toDouble(get(c_int), Columns::RELATIVE_INTENSITY);
      r.precursor_charge = toLong(get(c_z), "precursor charge");
      r.fragment_charge = toLong(get(c_fz), "fragment charge");
      r.ordinal = toLong(get(c_ord), Columns::FRAGMENT_SERIES_NUMBER);
      r.decoy = toLong(get(c_dec), Columns::DECOY);
      r.fragment_type = get(c_ft);
      r.loss_type = get(c_lt);
      builder.append(r);
    }
    if (in.bad()) { throw std::runtime_error("error reading library: " + filename); }
    builder.finish();
    if (builder.outOfRange())
    {
      std::cerr << "warning: dropped " << builder.outOfRange()
                << " transitions whose fragment ordinal or charge cannot describe "
                   "a fragment of the stated peptide\n";
    }
    if (builder.reopened())
    {
      throw std::runtime_error(
        "library rows are not grouped by precursor: " + std::to_string(builder.reopened())
        + " precursors reappear after another precursor's rows. Sort the file so each "
        "precursor's transitions are contiguous.");
    }
  }

  void DIANNLibraryFile::loadParquet(const std::string& filename, Library& library)
  {
    auto infile = arrow::io::ReadableFile::Open(filename);
    if (!infile.ok()) { throw std::runtime_error("cannot open library: " + filename); }

    auto reader_result = parquet::arrow::OpenFile(*infile, arrow::default_memory_pool());
    if (!reader_result.ok()) { throw std::runtime_error("not a Parquet file: " + filename); }
    std::unique_ptr<parquet::arrow::FileReader> reader = std::move(*reader_result);

    // Bound decoded input memory independently of the full spectral library.
    // The Builder spans batches so a precursor split at a boundary stays intact.
    reader->set_batch_size(65536);
    std::shared_ptr<arrow::RecordBatchReader> batches;
    if (!reader->GetRecordBatchReader(&batches).ok())
    { throw std::runtime_error("cannot read Parquet batches: " + filename); }
    Builder builder(library);
    for (;;)
    {
      std::shared_ptr<arrow::RecordBatch> batch;
      if (!batches->ReadNext(&batch).ok())
      { throw std::runtime_error("cannot read Parquet batch: " + filename); }
      if (!batch) { break; }
      if (batch->num_rows() == 0) { continue; }
      auto result = arrow::Table::FromRecordBatches({batch});
      if (!result.ok()) { throw std::runtime_error("cannot materialize Parquet batch: " + filename); }
      auto table = *result;

      // COMPACT layout? Product.Mz being a list means one row per PRECURSOR with
      // the transitions nested, rather than one row per transition. Detected from
      // the schema rather than from a version field, so a file is readable on its
      // own terms.
      {
        const int qi = table->schema()->GetFieldIndex(Columns::PRODUCT_MZ);
        if (qi >= 0 && table->schema()->field(qi)->type()->id() == arrow::Type::LIST)
        {
          loadParquetCompact(table, library);
          continue;
        }
      }

      // Decode dictionary and large_string columns rather than rejecting them.
      // Parquet writes repeated strings dictionary-encoded by default -- which is
      // the very property that makes this path cheap -- and any library over 2 GB
      // of characters needs large_string. Casting straight to StringArray yields
      // nullptr for both, after which the reader blamed the file for "missing
      // Precursor.Id" when the column was present all along.
      std::vector<std::shared_ptr<arrow::Array>> keep_alive;
      auto str = [&](const char* n) -> std::shared_ptr<arrow::StringArray> {
        const int i = table->schema()->GetFieldIndex(n);
        if (i < 0) { return nullptr; }
        std::shared_ptr<arrow::Array> a = table->column(i)->chunk(0);
        if (a->type_id() != arrow::Type::STRING)
        {
          auto casted = arrow::compute::Cast(arrow::Datum(a), arrow::utf8());
          if (!casted.ok()) { throw std::runtime_error(std::string("cannot decode Parquet string column ") + n + ": " + casted.status().ToString()); }
          a = casted->make_array();
          keep_alive.push_back(a);
        }
        return std::dynamic_pointer_cast<arrow::StringArray>(a);
      };
      auto num = [&](const char* n) -> std::shared_ptr<arrow::Array> {
        const int i = table->schema()->GetFieldIndex(n);
        return i < 0 ? nullptr : table->column(i)->chunk(0);
      };
      auto at = [](const std::shared_ptr<arrow::Array>& a, int64_t i) -> double {
        if (!a || a->IsNull(i)) { return std::numeric_limits<double>::quiet_NaN(); }
        switch (a->type_id())
        {
          case arrow::Type::DOUBLE: return static_cast<const arrow::DoubleArray&>(*a).Value(i);
          case arrow::Type::FLOAT:  return static_cast<const arrow::FloatArray&>(*a).Value(i);
          case arrow::Type::INT64:  return static_cast<double>(static_cast<const arrow::Int64Array&>(*a).Value(i));
          case arrow::Type::INT32:  return static_cast<double>(static_cast<const arrow::Int32Array&>(*a).Value(i));
          case arrow::Type::INT16:  return static_cast<double>(static_cast<const arrow::Int16Array&>(*a).Value(i));
          case arrow::Type::INT8:   return static_cast<double>(static_cast<const arrow::Int8Array&>(*a).Value(i));
          case arrow::Type::UINT64: return static_cast<double>(static_cast<const arrow::UInt64Array&>(*a).Value(i));
          case arrow::Type::UINT32: return static_cast<double>(static_cast<const arrow::UInt32Array&>(*a).Value(i));
          case arrow::Type::UINT16: return static_cast<double>(static_cast<const arrow::UInt16Array&>(*a).Value(i));
          case arrow::Type::UINT8:  return static_cast<double>(static_cast<const arrow::UInt8Array&>(*a).Value(i));
          // Decoy is naturally a bool, and that is what pandas and polars write.
          // Falling through to 0.0 here silently made every decoy a target.
          case arrow::Type::BOOL:   return static_cast<const arrow::BooleanArray&>(*a).Value(i) ? 1.0 : 0.0;
          default: return std::nan("");
        }
      };
      auto sv = [](const std::shared_ptr<arrow::StringArray>& a, int64_t i) -> std::string_view {
        if (!a || a->IsNull(i)) { return {}; }
        return a->GetView(i);
      };

      const auto a_id = str(Columns::PRECURSOR_ID);
      const auto a_seq = str(Columns::MODIFIED_SEQUENCE);
      const auto a_pg = str(Columns::PROTEIN_GROUP);
      const auto a_ft = str(Columns::FRAGMENT_TYPE);
      const auto a_lt = str(Columns::FRAGMENT_LOSS_TYPE);
      const auto a_pmz = num(Columns::PRECURSOR_MZ);
      const auto a_qmz = num(Columns::PRODUCT_MZ);
      const auto a_rt = num(Columns::RT);
      const auto a_im = num(Columns::IM);
      const auto a_ccs = num(Columns::CCS);
      const auto a_int = num(Columns::RELATIVE_INTENSITY);
      const auto a_z = num(Columns::PRECURSOR_CHARGE);
      const auto a_fz = num(Columns::FRAGMENT_CHARGE);
      const auto a_ord = num(Columns::FRAGMENT_SERIES_NUMBER);
      const auto a_dec = num(Columns::DECOY);

      if (!a_id || !a_pmz || !a_qmz)
      {
        throw std::runtime_error("not a DIA-NN library (missing Precursor.Id / Precursor.Mz / "
                                 "Product.Mz): " + filename);
      }

      const int64_t rows = table->num_rows();
      for (int64_t i = 0; i < rows; ++i)
      {
        Row r;
        r.precursor_id = sv(a_id, i);
        r.modified_sequence = sv(a_seq, i);
        r.protein_group = sv(a_pg, i);
        r.fragment_type = sv(a_ft, i);
        r.loss_type = sv(a_lt, i);
        r.precursor_mz = at(a_pmz, i);
        r.product_mz = at(a_qmz, i);
        r.rt = at(a_rt, i);
        r.im = at(a_im, i);
        r.ccs = at(a_ccs, i);
        r.intensity = at(a_int, i);
        // at() yields NaN for a type it cannot decode or a null cell; casting
        // that to long is undefined behaviour, which is exactly what the toFixed
        // guard was added to eliminate two functions away.
        // NaN was guarded; +-inf and any |v| >= 2^63 still reached the cast and
        // UBSan flagged them. On x86-64 they yield INT64_MIN, which narrows to
        // charge 0 -- a wrong value rather than a crash.
        auto as_long = [](double v) -> long {
          if (!std::isfinite(v) || v != std::floor(v)) { return 0L; }
          if (v <= static_cast<double>(std::numeric_limits<long>::min())) { return 0L; }
          if (v >= static_cast<double>(std::numeric_limits<long>::max())) { return 0L; }
          return static_cast<long>(v);
        };
        r.precursor_charge = as_long(at(a_z, i));
        r.fragment_charge = as_long(at(a_fz, i));
        r.ordinal = as_long(at(a_ord, i));
        r.decoy = as_long(at(a_dec, i));
        builder.append(r);
      }
    }
    builder.finish();
    if (builder.outOfRange())
    {
      std::cerr << "warning: dropped " << builder.outOfRange()
                << " transitions whose fragment ordinal or charge cannot describe "
                   "a fragment of the stated peptide\n";
    }
    if (builder.reopened())
    {
      throw std::runtime_error(
        "library rows are not grouped by precursor: " + std::to_string(builder.reopened())
        + " precursors reappear after another precursor's rows. Sort the file so each "
        "precursor's transitions are contiguous.");
    }
  }

  /// 7.28 GB for the human library, and it was the second-largest item in
  /// Phase 1 -- 211.7 s, 26% of the stage, and untimed until now. At 34.4 MB/s
  /// it sat on `ofstream <<`'s measured 25.4 MB/s line and nowhere near the
  /// NVMe's 1,042 MB/s, so the cost is number formatting, not disk. TextWriter
  /// renders the same digits through to_chars into a reusable buffer; see its
  /// header for why the precisions below are passed explicitly, and
  /// `odia_tsv_writers` for the byte-for-byte check against this function's
  /// previous ofstream form.
  void DIANNLibraryFile::store(const std::string& filename, const Library& library)
  {
    if (isParquet(filename)) { storeParquet(filename, library); }
    else { storeTSV(filename, library); }
  }

  std::string DIANNLibraryFile::hashFile(const std::string& path)
  {
    // Empty means the bundled OpenMS model, which is pinned and
    // read-only, so its identity is the install's. Hash it explicitly if a
    // build ever ships more than one.
    if (path.empty()) { return "bundled"; }
    std::ifstream in(path, std::ios::binary);
    if (!in) { throw std::runtime_error("cannot read model: " + path); }
    std::uint64_t h = 14695981039346656037ull;
    std::vector<char> buf(1 << 20);
    while (in)
    {
      in.read(buf.data(), std::streamsize(buf.size()));
      const std::streamsize got = in.gcount();
      for (std::streamsize i = 0; i < got; ++i)
      { h ^= static_cast<unsigned char>(buf[std::size_t(i)]); h *= 1099511628211ull; }
    }
    if (in.bad()) { throw std::runtime_error("error reading model: " + path); }
    std::ostringstream hex;
    hex << std::hex << std::setw(16) << std::setfill('0') << h;
    return hex.str();
  }

  DIANNLibraryFile::Fingerprint
  DIANNLibraryFile::fingerprintFasta(const std::string& fasta)
  {
    std::ifstream in(fasta, std::ios::binary);
    if (!in) { throw std::runtime_error("cannot read FASTA: " + fasta); }

    // FNV-1a over the CONTENT. Not path, size or mtime: the same FASTA staged
    // to /scratch must hit the cache, and an edited one that happens to keep
    // its size and timestamp must miss it. Not cryptographic -- this guards
    // against accident, not against an adversary editing a FASTA to collide.
    std::uint64_t h = 14695981039346656037ull;
    std::uint64_t bytes = 0;
    std::vector<char> buf(1 << 20);
    while (in)
    {
      in.read(buf.data(), std::streamsize(buf.size()));
      const std::streamsize got = in.gcount();
      for (std::streamsize i = 0; i < got; ++i)
      {
        h ^= static_cast<unsigned char>(buf[std::size_t(i)]);
        h *= 1099511628211ull;
      }
      bytes += std::uint64_t(got);
    }
    if (in.bad()) { throw std::runtime_error("error reading FASTA: " + fasta); }
    std::ostringstream hex;
    hex << std::hex << std::setw(16) << std::setfill('0') << h;
    Fingerprint fp;
    fp.fasta_hash = hex.str();
    fp.fasta_bytes = bytes;
    return fp;
  }

  std::string DIANNLibraryFile::readFingerprint(const std::string& filename,
                                                const char* key)
  {
    if (!std::filesystem::exists(filename)) { return {}; }
    try
    {
      auto infile = arrow::io::ReadableFile::Open(filename);
      if (!infile.ok()) { return {}; }
      // Read it back the SAME way it is written -- through the Arrow schema,
      // not the raw Parquet key-value block. store_schema() serialises the
      // whole Arrow schema under one key, so looking for "odia.fingerprint" at
      // the Parquet level finds nothing even though it is present.
      //
      // Schema ONLY, never the table: deciding whether to use a multi-GiB
      // library must not cost reading it.
      auto reader_result = parquet::arrow::OpenFile(*infile, arrow::default_memory_pool());
      if (!reader_result.ok()) { return {}; }
      std::shared_ptr<arrow::Schema> schema;
      if (!(*reader_result)->GetSchema(&schema).ok() || !schema) { return {}; }
      const auto kv = schema->metadata();
      if (!kv) { return {}; }
      const auto got = kv->Get(key);
      if (!got.ok()) { return {}; }
      return *got;
    }
    catch (const std::exception&) { return {}; }
  }

  void DIANNLibraryFile::loadParquetCompact(const std::shared_ptr<arrow::Table>& table_in,
                                            Library& library)
  {
    auto combined = table_in->CombineChunks(arrow::default_memory_pool());
    if (!combined.ok()) { throw std::runtime_error("cannot combine Parquet chunks"); }
    auto table = *combined;

    std::vector<std::shared_ptr<arrow::Array>> keep_alive;
    auto col = [&](const char* n) -> std::shared_ptr<arrow::Array> {
      const int i = table->schema()->GetFieldIndex(n);
      return i < 0 ? nullptr : table->column(i)->chunk(0);
    };
    // Dictionary columns decode to their values; the flat path does the same.
    auto str = [&](const char* n) -> std::shared_ptr<arrow::StringArray> {
      auto a = col(n);
      if (!a) { return nullptr; }
      if (a->type_id() != arrow::Type::STRING)
      {
        auto casted = arrow::compute::Cast(arrow::Datum(a), arrow::utf8());
        if (!casted.ok()) { throw std::runtime_error(std::string("cannot decode Parquet string column ") + n + ": " + casted.status().ToString()); }
        a = casted->make_array();
        keep_alive.push_back(a);
      }
      return std::static_pointer_cast<arrow::StringArray>(a);
    };
    auto num = [&](const std::shared_ptr<arrow::Array>& a, int64_t i) -> double {
      if (!a || a->IsNull(i)) { return std::numeric_limits<double>::quiet_NaN(); }
      switch (a->type_id())
      {
        case arrow::Type::DOUBLE: return static_cast<const arrow::DoubleArray&>(*a).Value(i);
        case arrow::Type::FLOAT:  return static_cast<const arrow::FloatArray&>(*a).Value(i);
        case arrow::Type::UINT8:  return static_cast<const arrow::UInt8Array&>(*a).Value(i);
        case arrow::Type::BOOL:   return static_cast<const arrow::BooleanArray&>(*a).Value(i) ? 1.0 : 0.0;
        case arrow::Type::INT64:  return static_cast<double>(static_cast<const arrow::Int64Array&>(*a).Value(i));
        case arrow::Type::INT32:  return static_cast<const arrow::Int32Array&>(*a).Value(i);
        default: return std::numeric_limits<double>::quiet_NaN();
      }
    };
    auto lst = [&](const char* n) -> std::shared_ptr<arrow::ListArray> {
      auto a = col(n);
      if (!a || a->type_id() != arrow::Type::LIST) { return nullptr; }
      return std::static_pointer_cast<arrow::ListArray>(a);
    };

    const auto a_seq = str(Columns::MODIFIED_SEQUENCE);
    const auto a_pg  = str(Columns::PROTEIN_GROUP);
    const auto a_z   = col(Columns::PRECURSOR_CHARGE);
    const auto a_dec = col(Columns::DECOY);
    const auto a_rt  = col(Columns::RT);
    const auto a_im  = col(Columns::IM);
    const auto a_ccs = col(Columns::CCS);
    const auto a_pmz = col(Columns::PRECURSOR_MZ);
    const auto l_qmz = lst(Columns::PRODUCT_MZ);
    const auto l_int = lst(Columns::RELATIVE_INTENSITY);
    const auto l_ft  = lst(Columns::FRAGMENT_TYPE);
    const auto l_fz  = lst(Columns::FRAGMENT_CHARGE);
    const auto l_ord = lst(Columns::FRAGMENT_SERIES_NUMBER);
    const auto l_loss = lst(Columns::FRAGMENT_LOSS_TYPE);
    if (!a_seq || !a_z || !a_pmz || !l_qmz || !l_int)
    { throw std::runtime_error("compact Parquet library is missing required columns"); }

    // Validate nested input before typed access: each precursor owns a separate
    // list in every column, and malformed lengths must never index another row.
    auto requireList = [&](const char* name, const std::shared_ptr<arrow::ListArray>& a,
                           arrow::Type::type type) {
      if (!a)
      {
        if (col(name)) { throw std::runtime_error(std::string("compact column is not a list: ") + name); }
        return;
      }
      if (a->values()->type_id() != type)
      { throw std::runtime_error(std::string("unsupported compact list element type: ") + name); }
      if (a->null_count() || a->values()->null_count())
      { throw std::runtime_error(std::string("null compact list or element: ") + name); }
      for (int64_t i = 0; i < table->num_rows(); ++i)
      {
        if (a->value_length(i) != l_qmz->value_length(i))
        { throw std::runtime_error("compact transition list lengths disagree"); }
      }
    };
    requireList(Columns::PRODUCT_MZ, l_qmz, arrow::Type::DOUBLE);
    requireList(Columns::RELATIVE_INTENSITY, l_int, arrow::Type::FLOAT);
    requireList(Columns::FRAGMENT_TYPE, l_ft, arrow::Type::UINT8);
    requireList(Columns::FRAGMENT_CHARGE, l_fz, arrow::Type::UINT8);
    requireList(Columns::FRAGMENT_SERIES_NUMBER, l_ord, arrow::Type::UINT8);
    requireList(Columns::FRAGMENT_LOSS_TYPE, l_loss, arrow::Type::UINT8);
    if (a_seq->null_count()) { throw std::runtime_error("null compact modified sequence"); }
    for (int64_t i = 0; i < table->num_rows(); ++i)
    {
      const double charge = num(a_z, i);
      if (!std::isfinite(charge) || charge != std::floor(charge) || charge < 1 || charge > 255)
      { throw std::runtime_error("invalid compact precursor charge"); }
    }

    const auto v_qmz = std::static_pointer_cast<arrow::DoubleArray>(l_qmz->values());
    const auto v_int = std::static_pointer_cast<arrow::FloatArray>(l_int->values());
    const auto v_ft  = l_ft ? std::static_pointer_cast<arrow::UInt8Array>(l_ft->values()) : nullptr;
    const auto v_fz  = l_fz ? std::static_pointer_cast<arrow::UInt8Array>(l_fz->values()) : nullptr;
    const auto v_ord = l_ord ? std::static_pointer_cast<arrow::UInt8Array>(l_ord->values()) : nullptr;
    const auto v_loss = l_loss ? std::static_pointer_cast<arrow::UInt8Array>(l_loss->values()) : nullptr;

    auto& p = library.precursors();
    auto& t = library.transitions();
    const int64_t n = table->num_rows();

    for (int64_t i = 0; i < n; ++i)
    {
      const auto seq = a_seq->GetView(i);
      p.mz.push_back(toFixed(num(a_pmz, i)));
      p.irt.push_back(static_cast<float>(num(a_rt, i)));
      const double im = num(a_im, i);
      p.im.push_back(im > 0.0 ? static_cast<float>(im) : std::nanf(""));
      const double ccs = num(a_ccs, i);
      p.ccs.push_back(ccs > 0.0 ? static_cast<float>(ccs) : std::nanf(""));
      p.charge.push_back(static_cast<std::uint8_t>(num(a_z, i)));
      p.decoy.push_back(a_dec && !a_dec->IsNull(i) && num(a_dec, i) != 0.0 ? std::uint8_t{1} : std::uint8_t{0});
      p.modified_sequence.push_back(library.strings().intern(
        std::string_view(seq.data(), seq.size())));
      const auto pg = a_pg && !a_pg->IsNull(i) ? a_pg->GetView(i) : std::string_view{};
      p.protein_group.push_back(library.strings().intern(
        std::string_view(pg.data(), pg.size())));

      p.transition_begin.push_back(static_cast<std::uint32_t>(t.product_mz.size()));
      const int64_t b = l_qmz->value_offset(i), e = l_qmz->value_offset(i + 1);
      for (int64_t k = b; k < e; ++k)
      {
        t.product_mz.push_back(toFixed(v_qmz->Value(k)));
        const auto offset = k - b;
        t.library_intensity.push_back(v_int->Value(l_int->value_offset(i) + offset));
        const auto type = v_ft ? v_ft->Value(l_ft->value_offset(i) + offset) : std::uint8_t{0};
        const auto charge = v_fz ? v_fz->Value(l_fz->value_offset(i) + offset) : std::uint8_t{1};
        const auto ordinal = v_ord ? v_ord->Value(l_ord->value_offset(i) + offset) : std::uint8_t{0};
        if (type > static_cast<std::uint8_t>(FragmentType::Precursor))
        { throw std::runtime_error("invalid compact fragment type"); }
        if (charge < 1 || charge > 127) { throw std::runtime_error("invalid compact fragment charge"); }
        if (v_ord && (ordinal < 1 || ordinal >= residueCount(seq)))
        { throw std::runtime_error("invalid compact fragment ordinal"); }
        t.type.push_back(static_cast<FragmentType>(type));
        t.charge.push_back(static_cast<std::int8_t>(charge));
        t.ordinal.push_back(ordinal);
        const auto loss = v_loss ? v_loss->Value(l_loss->value_offset(i) + offset) : std::uint8_t{0};
        if (loss > static_cast<std::uint8_t>(LossType::Other)) { throw std::runtime_error("invalid compact fragment loss"); }
        // Older compact-v1 files omitted the all-noloss column.
        t.loss.push_back(static_cast<LossType>(loss));
      }
      p.transition_count.push_back(static_cast<std::uint32_t>(e - b));
    }
    library.markUnsorted();
  }

  void DIANNLibraryFile::storeParquetCompact(const std::string& filename,
                                             const Library& library,
                                             const Fingerprint& fp,
                                             const std::string& config_json)
  {
    const auto& p = library.precursors();
    const auto& t = library.transitions();
    const std::size_t n = library.precursorCount();

    auto ok = [](const arrow::Status& st) {
      if (!st.ok()) { throw std::runtime_error("parquet build: " + st.ToString()); }
    };

    // Precursor-level: ONE value each, not one per transition.
    arrow::StringDictionaryBuilder b_id, b_seq, b_pg;
    arrow::UInt8Builder b_z, b_dec;
    arrow::FloatBuilder b_rt, b_im, b_ccs;
    arrow::DoubleBuilder b_pmz;

    // Transition-level: a LIST per precursor. The list offsets replace the
    // repeated precursor key entirely -- there is no foreign key to store.
    auto pool = arrow::default_memory_pool();
    auto qmz_v = std::make_shared<arrow::DoubleBuilder>(pool);
    auto int_v = std::make_shared<arrow::FloatBuilder>(pool);
    auto ft_v  = std::make_shared<arrow::UInt8Builder>(pool);
    auto fz_v  = std::make_shared<arrow::UInt8Builder>(pool);
    auto ord_v = std::make_shared<arrow::UInt8Builder>(pool);
    auto loss_v = std::make_shared<arrow::UInt8Builder>(pool);
    arrow::ListBuilder b_qmz(pool, qmz_v), b_int(pool, int_v),
                       b_ft(pool, ft_v), b_fz(pool, fz_v), b_ord(pool, ord_v), b_loss(pool, loss_v);

    for (std::size_t i = 0; i < n; ++i)
    {
      const auto seq = library.strings().get(p.modified_sequence[i]);
      const auto pg = library.strings().get(p.protein_group[i]);
      const std::string seq_s(seq.data(), seq.size());
      const std::string pg_s(pg.data(), pg.size());
      const int z = static_cast<int>(p.charge[i]);
      const std::string id = seq_s + std::to_string(z) + (p.decoy[i] ? "_decoy" : "");

      ok(b_id.Append(id)); ok(b_seq.Append(seq_s)); ok(b_pg.Append(pg_s));
      ok(b_z.Append(static_cast<std::uint8_t>(z)));
      ok(b_dec.Append(static_cast<std::uint8_t>(p.decoy[i])));
      // NaN means ABSENT and is written as null. A zero ion mobility is a
      // VALUE that gates every peak against 0 and rejects them all.
      if (!std::isnan(p.irt[i])) { ok(b_rt.Append(p.irt[i])); } else { ok(b_rt.AppendNull()); }
      if (i < p.im.size() && !std::isnan(p.im[i])) { ok(b_im.Append(p.im[i])); }
      else { ok(b_im.AppendNull()); }
      if (i < p.ccs.size() && !std::isnan(p.ccs[i])) { ok(b_ccs.Append(p.ccs[i])); }
      else { ok(b_ccs.AppendNull()); }
      ok(b_pmz.Append(fromFixed(p.mz[i])));

      ok(b_qmz.Append()); ok(b_int.Append()); ok(b_ft.Append());
      ok(b_fz.Append()); ok(b_ord.Append()); ok(b_loss.Append());
      const std::uint32_t begin = p.transition_begin[i];
      for (std::uint32_t k = 0; k < p.transition_count[i]; ++k)
      {
        const std::uint32_t j = begin + k;
        ok(qmz_v->Append(fromFixed(t.product_mz[j])));
        ok(int_v->Append(t.library_intensity[j]));
        // The enum's own value, not its spelling: "noloss" and "y" cost bytes
        // and a parse on every read.
        ok(ft_v->Append(static_cast<std::uint8_t>(t.type[j])));
        ok(fz_v->Append(static_cast<std::uint8_t>(t.charge[j])));
        ok(ord_v->Append(static_cast<std::uint8_t>(t.ordinal[j])));
        ok(loss_v->Append(static_cast<std::uint8_t>(t.loss[j])));
      }
    }

    std::shared_ptr<arrow::Array> a_id, a_seq, a_pg, a_z, a_dec, a_rt, a_im,
                                  a_ccs, a_pmz, a_qmz, a_int, a_ft, a_fz, a_ord, a_loss;
    ok(b_id.Finish(&a_id));   ok(b_seq.Finish(&a_seq)); ok(b_pg.Finish(&a_pg));
    ok(b_z.Finish(&a_z));     ok(b_dec.Finish(&a_dec)); ok(b_rt.Finish(&a_rt));
    ok(b_im.Finish(&a_im));   ok(b_ccs.Finish(&a_ccs)); ok(b_pmz.Finish(&a_pmz));
    ok(b_qmz.Finish(&a_qmz)); ok(b_int.Finish(&a_int)); ok(b_ft.Finish(&a_ft));
    ok(b_fz.Finish(&a_fz));   ok(b_ord.Finish(&a_ord)); ok(b_loss.Finish(&a_loss));

    auto schema = arrow::schema({
      arrow::field(Columns::PRECURSOR_ID, a_id->type()),
      arrow::field(Columns::MODIFIED_SEQUENCE, a_seq->type()),
      arrow::field(Columns::PROTEIN_GROUP, a_pg->type()),
      arrow::field(Columns::PRECURSOR_CHARGE, arrow::uint8()),
      arrow::field(Columns::DECOY, arrow::uint8()),
      arrow::field(Columns::RT, arrow::float32()),
      arrow::field(Columns::IM, arrow::float32()),
      arrow::field(Columns::CCS, arrow::float32()),
      arrow::field(Columns::PRECURSOR_MZ, arrow::float64()),
      arrow::field(Columns::PRODUCT_MZ, arrow::list(arrow::float64())),
      arrow::field(Columns::RELATIVE_INTENSITY, arrow::list(arrow::float32())),
      arrow::field(Columns::FRAGMENT_TYPE, arrow::list(arrow::uint8())),
      arrow::field(Columns::FRAGMENT_CHARGE, arrow::list(arrow::uint8())),
      arrow::field(Columns::FRAGMENT_SERIES_NUMBER, arrow::list(arrow::uint8())),
      arrow::field(Columns::FRAGMENT_LOSS_TYPE, arrow::list(arrow::uint8())),
    });
    std::vector<std::string> md_keys{
      "odia.fingerprint", "odia.target_fingerprint", "odia.layout",
      "odia.decoy_semantics", "odia.decoy_method",
      "odia.fasta_fnv1a64", "odia.fasta_bytes", "odia.params"};
    schema = schema->WithMetadata(arrow::key_value_metadata(
      md_keys,
      // decoy_semantics is stated because a consumer CANNOT infer it and one
      // already got it wrong: OpenSwathAssayGenerator recomputed fragments from
      // Modified.Sequence, found the decoys identical to their targets, and
      // deduplicated them into a PQP with zero decoys. A decoy row here carries
      // its TARGET's sequence by design -- DIA-NN's convention -- with only the
      // fragment m/z shifted, so the sequence does NOT generate the fragments.
      {fp.key(), fp.targetKey(), "compact-v1",
       "target-sequence-with-shifted-fragments", fp.decoy_method,
       fp.fasta_hash, std::to_string(fp.fasta_bytes), fp.params}));
    // The recipe travels with the library. Appended rather than folded into the
    // list above so the existing keys keep their positions for any reader that
    // indexes them.
    if (!config_json.empty())
    {
      auto md = schema->metadata()->Copy();
      md->Append("odia.config_json", config_json);
      schema = schema->WithMetadata(md);
    }

    auto table = arrow::Table::Make(schema,
      {a_id, a_seq, a_pg, a_z, a_dec, a_rt, a_im, a_ccs, a_pmz,
       a_qmz, a_int, a_ft, a_fz, a_ord, a_loss});

    AtomicFile output(filename);
    auto outfile = arrow::io::FileOutputStream::Open(output.temporaryPath().string());
    if (!outfile.ok()) { throw std::runtime_error("cannot write library: " + filename); }
    auto props = parquet::WriterProperties::Builder()
                   .compression(parquet::Compression::ZSTD)
                   ->enable_dictionary()
                   ->build();
    auto arrow_props = parquet::ArrowWriterProperties::Builder().store_schema()->build();
    const auto st = parquet::arrow::WriteTable(*table, pool, *outfile, 1 << 20,
                                               props, arrow_props);
    if (!st.ok()) { throw std::runtime_error("cannot write Parquet: " + st.ToString()); }
    ok((*outfile)->Close());
    output.commit();
  }

  void DIANNLibraryFile::storeParquet(const std::string& filename, const Library& library)
  { storeParquet(filename, library, Fingerprint{}); }

  void DIANNLibraryFile::storeParquet(const std::string& filename, const Library& library,
                                      const Fingerprint& fp)
  {
    const auto& p = library.precursors();
    const auto& t = library.transitions();

    // One row per TRANSITION, matching the TSV exactly, because loadParquet and
    // loadTSV populate the same structures from the same column names. A
    // narrower per-precursor schema would be smaller still but would not round
    // trip through the existing reader.
    // DICTIONARY builders: store an INDEX per row and each distinct string
    // once, which is what the in-memory library already does with its
    // StringArena and what this writer was throwing away.
    //
    // Measured on the parity library: 1,588,688 distinct strings totalling
    // 32,920,618 bytes (31.4 MiB), materialised across 119,088,506 transition
    // rows at 69.0 characters each -- Precursor.Id 20.3, Modified.Sequence
    // 19.3, Protein.Group 22.3 -- for 7.65 GiB of string bytes. A 250x
    // amplification of 31 MiB of actual content, because a precursor's
    // sequence is rewritten for all ~12 of its transitions.
    //
    // That is what overflowed the 32-bit offset limit (2147483648 bytes of
    // Precursor.Id alone) and cost 17.9 GB of RAM during the write. Parquet
    // would have dictionary-compressed the FILE anyway; the blow-up was
    // entirely in the Arrow builder before the file was ever touched.
    //
    // loadParquet already decodes dictionary columns -- "Parquet writes
    // repeated strings dictionary-encoded by default, which is the very
    // property that makes this path cheap" -- so nothing on the read side
    // changes.
    arrow::StringDictionaryBuilder b_id, b_seq, b_pg, b_ft, b_lt;
    // Widths that MATCH THE IN-MEMORY TYPES, which is both lossless and
    // minimal. charge/decoy/fragment charge/ordinal are uint8 in the library;
    // irt, im, ccs and library_intensity are float. Writing them as int32 and
    // float64 was pure upcasting: 28 bytes per row of padding across
    // 119,088,506 rows, about 3.3 GB of buffer for no information.
    //
    // m/z stays float64 DELIBERATELY. The library holds it as uint32
    // fixed-point at 1e-5 Th precisely because that beats float32 -- 0.0025 ppm
    // at m/z 2000 against float32's flat 0.03-0.05 ppm (Library.h:20-23) -- so
    // narrowing to float32 would store LESS precision than memory holds and
    // break the exact round trip. Storing the raw fixed-point integer instead
    // would be 4 bytes and exact, but changes what the column MEANS: the
    // reader would return 97199181 where it expects 971.99181. That is a
    // format-version change, not a width change, and is left as one.
    arrow::UInt8Builder b_z, b_dec, b_fz, b_ord;
    arrow::FloatBuilder b_rt, b_im, b_ccs, b_int;
    arrow::DoubleBuilder b_pmz, b_qmz;

    auto ok = [](const arrow::Status& st) {
      if (!st.ok()) { throw std::runtime_error("parquet build: " + st.ToString()); }
    };

    for (std::size_t i = 0; i < library.precursorCount(); ++i)
    {
      const auto seq = library.strings().get(p.modified_sequence[i]);
      const auto pg = library.strings().get(p.protein_group[i]);
      const std::string seq_s(seq.data(), seq.size());
      const std::string pg_s(pg.data(), pg.size());
      const int z = static_cast<int>(p.charge[i]);
      // The suffix is what makes a decoy's id distinct from its target's: they
      // share sequence, charge and precursor m/z, so <sequence><charge> alone
      // collides and a reload silently merges the pair. See storeTSV.
      const std::string id = seq_s + std::to_string(z) + (p.decoy[i] ? "_decoy" : "");
      // NaN means ABSENT and is written as null, not as 0. A zero ion mobility
      // is a VALUE that gates every peak against 0 and rejects them all; the
      // TSV writer emits 0 here and only survives because the reader maps 0
      // back to NaN. Parquet has real nulls, so this path does not need that
      // round-trip coincidence.
      const bool has_im = i < p.im.size() && !std::isnan(p.im[i]);
      const bool has_ccs = i < p.ccs.size() && !std::isnan(p.ccs[i]);
      const bool has_rt = !std::isnan(p.irt[i]);

      const std::uint32_t begin = p.transition_begin[i];
      for (std::uint32_t k = 0; k < p.transition_count[i]; ++k)
      {
        const std::uint32_t j = begin + k;
        ok(b_id.Append(id));
        ok(b_seq.Append(seq_s));
        ok(b_pg.Append(pg_s));
        ok(b_z.Append(static_cast<std::uint8_t>(z)));
        ok(b_dec.Append(static_cast<std::uint8_t>(p.decoy[i])));
        if (has_rt) { ok(b_rt.Append(p.irt[i])); } else { ok(b_rt.AppendNull()); }
        if (has_im) { ok(b_im.Append(p.im[i])); } else { ok(b_im.AppendNull()); }
        if (has_ccs) { ok(b_ccs.Append(p.ccs[i])); } else { ok(b_ccs.AppendNull()); }
        ok(b_pmz.Append(fromFixed(p.mz[i])));
        ok(b_qmz.Append(fromFixed(t.product_mz[j])));
        ok(b_int.Append(t.library_intensity[j]));
        ok(b_ft.Append(std::string(toString(t.type[j]))));
        ok(b_fz.Append(static_cast<std::uint8_t>(t.charge[j])));
        ok(b_ord.Append(static_cast<std::uint8_t>(t.ordinal[j])));
        ok(b_lt.Append(std::string(toString(t.loss[j]))));
      }
    }

    std::shared_ptr<arrow::Array> a_id, a_seq, a_pg, a_ft, a_lt, a_z, a_dec,
                                  a_fz, a_ord, a_rt, a_im, a_ccs, a_pmz, a_qmz, a_int;
    ok(b_id.Finish(&a_id));   ok(b_seq.Finish(&a_seq)); ok(b_pg.Finish(&a_pg));
    ok(b_ft.Finish(&a_ft));   ok(b_lt.Finish(&a_lt));   ok(b_z.Finish(&a_z));
    ok(b_dec.Finish(&a_dec)); ok(b_fz.Finish(&a_fz));   ok(b_ord.Finish(&a_ord));
    ok(b_rt.Finish(&a_rt));   ok(b_im.Finish(&a_im));   ok(b_ccs.Finish(&a_ccs));
    ok(b_pmz.Finish(&a_pmz)); ok(b_qmz.Finish(&a_qmz)); ok(b_int.Finish(&a_int));

    auto schema = arrow::schema({
      arrow::field(Columns::PRECURSOR_ID, a_id->type()),
      arrow::field(Columns::MODIFIED_SEQUENCE, a_seq->type()),
      arrow::field(Columns::PRECURSOR_CHARGE, arrow::uint8()),
      arrow::field(Columns::DECOY, arrow::uint8()),
      arrow::field(Columns::RT, arrow::float32()),
      arrow::field(Columns::IM, arrow::float32()),
      arrow::field(Columns::CCS, arrow::float32()),
      arrow::field(Columns::PRECURSOR_MZ, arrow::float64()),
      arrow::field(Columns::PRODUCT_MZ, arrow::float64()),
      arrow::field(Columns::RELATIVE_INTENSITY, arrow::float32()),
      arrow::field(Columns::FRAGMENT_TYPE, a_ft->type()),
      arrow::field(Columns::FRAGMENT_CHARGE, arrow::uint8()),
      arrow::field(Columns::FRAGMENT_SERIES_NUMBER, arrow::uint8()),
      arrow::field(Columns::FRAGMENT_LOSS_TYPE, a_lt->type()),
      arrow::field(Columns::PROTEIN_GROUP, a_pg->type()),
    });
    if (!fp.params.empty() || !fp.fasta_hash.empty())
    {
      // Human-readable siblings alongside the compared key, so a library found
      // on disk months later can be explained without running anything.
      schema = schema->WithMetadata(arrow::key_value_metadata(
        {"odia.fingerprint", "odia.fasta_fnv1a64", "odia.fasta_bytes", "odia.params"},
        {fp.key(), fp.fasta_hash, std::to_string(fp.fasta_bytes), fp.params}));
    }
    auto table = arrow::Table::Make(schema,
      {a_id, a_seq, a_z, a_dec, a_rt, a_im, a_ccs, a_pmz, a_qmz, a_int,
       a_ft, a_fz, a_ord, a_lt, a_pg});

    AtomicFile output(filename);
    auto outfile = arrow::io::FileOutputStream::Open(output.temporaryPath().string());
    if (!outfile.ok()) { throw std::runtime_error("cannot write library: " + filename); }
    // Dictionary encoding is the whole point: the sequence and protein group
    // repeat once per transition, twelve times per precursor.
    auto props = parquet::WriterProperties::Builder()
                   .compression(parquet::Compression::ZSTD)
                   ->enable_dictionary()
                   ->build();
    // store_schema(), or the Arrow schema metadata -- which is where the
    // fingerprint lives -- is silently dropped and every cache lookup misses.
    // Verified: without it the written file reports no key-value metadata at
    // all, from either the file level or the Arrow schema.
    auto arrow_props = parquet::ArrowWriterProperties::Builder().store_schema()->build();
    const auto st = parquet::arrow::WriteTable(*table, arrow::default_memory_pool(),
                                               *outfile, 1 << 20, props, arrow_props);
    if (!st.ok()) { throw std::runtime_error("cannot write Parquet: " + st.ToString()); }
    ok((*outfile)->Close());
    output.commit();
  }

  void DIANNLibraryFile::storeTSV(const std::string& filename, const Library& library)
  {
    AtomicFile output(filename);
    TextWriter out(output.temporaryPath().string());

    for (const char* column : {Columns::PRECURSOR_ID, Columns::MODIFIED_SEQUENCE,
                               Columns::PRECURSOR_CHARGE, Columns::DECOY,
                               Columns::RT, Columns::IM,
                               Columns::PRECURSOR_MZ, Columns::PRODUCT_MZ,
                               Columns::RELATIVE_INTENSITY, Columns::FRAGMENT_TYPE,
                               Columns::FRAGMENT_CHARGE, Columns::FRAGMENT_SERIES_NUMBER,
                               Columns::FRAGMENT_LOSS_TYPE, Columns::PROTEIN_GROUP})
    {
      out.put(column);
      out.put('\t');
    }
    out.put(Columns::CCS);
    out.put('\n');

    const auto& p = library.precursors();
    const auto& t = library.transitions();
    for (std::size_t i = 0; i < library.precursorCount(); ++i)
    {
      const auto seq = library.strings().get(p.modified_sequence[i]);
      const auto pg = library.strings().get(p.protein_group[i]);
      const int z = p.charge[i];
      if (seq.find_first_of("\t\r\n\0", 0, 4) != std::string_view::npos ||
          pg.find_first_of("\t\r\n\0", 0, 4) != std::string_view::npos)
      { throw std::runtime_error("TSV sequence or protein group contains a tab, newline or NUL"); }
      const bool has_im = i < p.im.size() && !std::isnan(p.im[i]);
      const bool has_ccs = i < p.ccs.size() && !std::isnan(p.ccs[i]);
      const std::uint32_t begin = p.transition_begin[i];
      for (std::uint32_t k = 0; k < p.transition_count[i]; ++k)
      {
        const std::uint32_t j = begin + k;
        // DIA-NN gives a decoy the same Modified.Sequence, charge and
        // Precursor.Mz as its target, so <sequence><charge> alone is not unique.
        // Without the suffix, reloading our own output merged each such pair and
        // re-labelled the decoy's transitions as target evidence -- 33,390
        // precursors became 33,386, silently, with no transitions lost.
        out.put(seq); out.integer(z); out.put(p.decoy[i] ? "_decoy" : ""); out.put('\t');
        out.put(seq); out.put('\t'); out.integer(z); out.put('\t');
        out.integer(p.decoy[i]); out.put('\t');
        // An unpredicted retention time is written as an empty field, not as
        // "nan": no TSV consumer accepts the latter, and it propagates into
        // anything that reads the library back.
        if (!std::isnan(p.irt[i])) { out.number(p.irt[i], 9); }
        out.put('\t');
        if (has_im) { out.number(p.im[i], 9); }
        out.put('\t');
        out.number(fromFixed(p.mz[i]), 10); out.put('\t');
        out.number(fromFixed(t.product_mz[j]), 10); out.put('\t');
        if (!std::isnan(t.library_intensity[j])) { out.number(t.library_intensity[j], 9); }
        out.put('\t');
        out.put(toString(t.type[j])); out.put('\t');
        out.integer(t.charge[j]); out.put('\t');
        out.integer(t.ordinal[j]); out.put('\t');
        out.put(toString(t.loss[j])); out.put('\t');
        out.put(pg); out.put('\t');
        // Empty rather than "nan" when absent, as for RT: no TSV consumer
        // accepts the latter, and a 0 here would read as a real cross-section.
        if (has_ccs) { out.number(p.ccs[i], 9); }
        out.put('\n');
      }
    }

    // A file checked only at open reports success on a full disk, an exceeded
    // quota or a broken mount, leaving a truncated library behind.
    out.close();
    output.commit();
  }

} // namespace ODIA
