// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <odia/DIANNLibraryFile.h>
#include <odia/ModifiedSequence.h>
#include <odia/Trace.h>
#include <odia/AtomicFile.h>
#include <odia/TextWriter.h>

#include "Parallel.h"
#include "ParseDouble.h"

#include <arrow/api.h>
#include <arrow/compute/api.h>
#include <arrow/io/file.h>
#include <arrow/io/interfaces.h>
#include <arrow/util/thread_pool.h>
#include <parquet/arrow/reader.h>
#include <parquet/arrow/writer.h>
#include <parquet/file_reader.h>
#include <iomanip>
#include <filesystem>
#include <arrow/builder.h>
#include <arrow/table.h>

#include <array>
#include <deque>
#include <atomic>
#include <charconv>
#include <algorithm>
#include <cctype>
#include <cstring>
#include <cmath>
#include <fstream>
#include <initializer_list>
#include <cstdio>
#include <iostream>
#include <locale>
#include <sstream>
#include <stdexcept>
#include <type_traits>
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
      if (!onlyDecimalCharacters(s))
      {
        std::string token(s);
        std::transform(token.begin(), token.end(), token.begin(),
                       [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
        if (token == "nan" || token == "+nan" || token == "-nan") { return std::nan(""); }
        throw std::runtime_error(std::string("invalid numeric value in ") + column + ": " + std::string(s));
      }
      // Correctly rounded on every platform; see ParseDouble.h.
      if (!parseDouble(s, v) || std::isinf(v))
      { throw std::runtime_error(std::string("invalid numeric value in ") + column + ": " + std::string(s)); }
      return v;
    }

    long toLong(std::string_view s, const char* column)
    {
      if (s.empty()) { return 0; }
      long v = 0;
      const auto* first = s.data();
      const auto* last = s.data() + s.size();
      if (*first == '+' && s.size() > 1 && s[1] != '-') { ++first; }
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
      static constexpr std::size_t kUnknown = static_cast<std::size_t>(-1);

      std::string_view precursor_id;
      std::string_view modified_sequence;
      std::string_view protein_group;
      double precursor_mz = 0, product_mz = 0;
      double rt = 0, im = 0, ccs = 0, intensity = 0;
      long precursor_charge = 0, fragment_charge = 0, ordinal = 0, decoy = 0;
      FragmentType fragment_type = FragmentType::Unknown;
      LossType loss_type = LossType::None;
      /// residueCount(modified_sequence), when the reader already has it.
      std::size_t residues = kUnknown;
    };

    /// The precursor-charge rule, checked where a row is parsed: the TSV
    /// reader parses on several threads and must report the first bad line
    /// a serial read would have met, whichever check that line fails.
    void checkPrecursorCharge(long charge)
    {
      if (charge < 1 || charge > 255)
      { throw std::runtime_error("invalid precursor charge: expected an integer between 1 and 255"); }
    }

    /// Appends rows to the library, starting a new precursor when the
    /// Precursor.Id changes.
    class Builder
    {
    public:
      explicit Builder(Library& lib) : lib_(lib) {}

      void append(const Row& r)
      {
        checkPrecursorCharge(r.precursor_charge);
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
            static_cast<std::size_t>(r.ordinal) >= residues(r))
        {
          ++out_of_range_;
          return;
        }

        auto& t = lib_.transitions();
        Library::checkTransitionCapacity(t.product_mz.size(), 1);
        t.product_mz.push_back(toFixed(r.product_mz));
        t.library_intensity.push_back(static_cast<float>(r.intensity));
        t.type.push_back(r.fragment_type);
        t.ordinal.push_back(static_cast<std::uint8_t>(r.ordinal));
        t.charge.push_back(static_cast<std::int8_t>(r.fragment_charge));
        t.loss.push_back(r.loss_type);
        ++count_;
      }

      void finish() { closeCurrent(); }

      /// Number of precursors whose rows were not contiguous.
      std::size_t reopened() const { return reopened_; }

      /// Transitions dropped because an ordinal or charge could not be stored.
      std::size_t outOfRange() const { return out_of_range_; }

    private:
      /// Residues of the row's sequence, counted once per distinct sequence
      /// run rather than once per transition.
      std::size_t residues(const Row& r)
      {
        if (r.residues != Row::kUnknown) { return r.residues; }
        if (!have_residues_ || r.modified_sequence != residue_sequence_)
        {
          residue_sequence_.assign(r.modified_sequence);
          residue_count_ = residueCount(r.modified_sequence);
          have_residues_ = true;
        }
        return residue_count_;
      }

      void startPrecursor(const Row& r)
      {
        Library::checkTransitionCapacity(lib_.transitionCount());
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
      std::string residue_sequence_;
      std::size_t residue_count_ = 0;
      bool have_residues_ = false;
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

  void DIANNLibraryFile::load(const std::string& filename, Library& library, unsigned threads)
  {
    if (isParquet(filename)) { loadParquet(filename, library, threads); }
    else { loadTSV(filename, library, threads); }
    // A library that carries only one of the two mobility quantities gets the
    // other derived, so a consumer never has to know which the producer chose.
    // Measured need: our own generated library reached DIA-NN with 1/K0 absent
    // for every precursor, and DIA-NN's `iIM` was 0 on 100% of a diaPASEF run's
    // identifications.
    completeMobility(library);
    // The interning index answered "have I seen this string" while reading and
    // is never consulted again -- handles resolve through the arena's entries.
    // Held for the library's lifetime it is one hash node per distinct string
    // plus a bucket array, for a question nobody asks twice.
    library.strings().releaseLookup();
  }

  namespace
  {
    /// Column positions of a TSV library, -1 for an absent column.
    struct TsvLayout
    {
      std::size_t width = 0;
      int id = -1, seq = -1, pg = -1, pmz = -1, qmz = -1, rt = -1, im = -1, ccs = -1, intensity = -1;
      int z = -1, fz = -1, ord = -1, dec = -1, ft = -1, lt = -1;
    };

    /// Bytes of text one parse task covers: 4 MiB, or ODIA_TSV_RANGE_BYTES
    /// (tests: many ranges and reads from a small file).
    std::size_t tsvRangeBytes()
    {
      static const std::size_t bytes = [] {
        if (const char* v = std::getenv("ODIA_TSV_RANGE_BYTES"); v != nullptr && *v != '\0')
        {
          char* end = nullptr;
          const long long x = std::strtoll(v, &end, 10);
          if (end != v && x > 0) { return static_cast<std::size_t>(x); }
        }
        return std::size_t(4) << 20;
      }();
      return bytes;
    }

    /// Whole lines of a TSV library, parsed on their own: the rows, the
    /// synthesized ids they point into, and the first line that failed.
    struct TsvRange
    {
      const char* begin = nullptr;
      const char* end = nullptr;
      std::vector<Row> rows;
      std::string ids;
      /// Lines in the range, empty ones included: what getline would count.
      std::size_t lines = 0;
      /// The first failure: its exception, or a row width that differs from
      /// the header (whose message needs the file's line number), at
      /// 0-based line `error_line` of the range. The rows before it are kept.
      std::exception_ptr error;
      bool width_error = false;
      std::size_t error_line = 0;
    };

    /// Parses one range into Rows, exactly as the serial reader parsed each
    /// line: every field of every row, precursor-level ones included, so a
    /// malformed value on any row is still reported. Never throws.
    void parseTsvRange(TsvRange& range, const TsvLayout& c)
    {
      range.rows.clear();
      range.ids.clear();
      range.lines = 0;
      range.error = nullptr;
      range.width_error = false;
      std::vector<std::string_view> f;
      f.reserve(c.width);
      // Synthesized ids live in range.ids, which may reallocate while it
      // grows; rows record offsets and get their views at the end.
      std::vector<std::pair<std::size_t, std::size_t>> id_spans;
      std::string_view residue_sequence;
      std::size_t residues = 0;
      bool have_residues = false;

      try
      {
        const char* p = range.begin;
        while (p < range.end)
        {
          const auto* nl = static_cast<const char*>(std::memchr(p, '\n', static_cast<std::size_t>(range.end - p)));
          const char* stop = nl != nullptr ? nl : range.end;
          std::string_view line(p, static_cast<std::size_t>(stop - p));
          p = nl != nullptr ? nl + 1 : range.end;
          ++range.lines;
          if (!line.empty() && line.back() == '\r') { line.remove_suffix(1); }
          if (line.empty()) { continue; }

          f.clear();
          std::size_t pos = 0;
          while (pos <= line.size())
          {
            const std::size_t tab = line.find('\t', pos);
            const std::size_t end = tab == std::string_view::npos ? line.size() : tab;
            f.push_back(line.substr(pos, end - pos));
            if (tab == std::string_view::npos) { break; }
            pos = tab + 1;
          }
          if (f.size() != c.width)
          {
            range.width_error = true;
            range.error_line = range.lines - 1;
            break;
          }

          auto get = [&](int i) -> std::string_view {
            return (i >= 0 && static_cast<std::size_t>(i) < f.size()) ? f[static_cast<std::size_t>(i)] : std::string_view{};
          };

          Row r;
          r.modified_sequence = get(c.seq);
          r.protein_group = get(c.pg);
          r.precursor_mz = toDouble(get(c.pmz), DIANNLibraryFile::Columns::PRECURSOR_MZ);
          r.product_mz = toDouble(get(c.qmz), DIANNLibraryFile::Columns::PRODUCT_MZ);
          r.rt = toDouble(get(c.rt), DIANNLibraryFile::Columns::RT);
          r.im = toDouble(get(c.im), DIANNLibraryFile::Columns::IM);
          r.ccs = toDouble(get(c.ccs), DIANNLibraryFile::Columns::CCS);
          r.intensity = toDouble(get(c.intensity), DIANNLibraryFile::Columns::RELATIVE_INTENSITY);
          r.precursor_charge = toLong(get(c.z), "precursor charge");
          r.fragment_charge = toLong(get(c.fz), "fragment charge");
          r.ordinal = toLong(get(c.ord), DIANNLibraryFile::Columns::FRAGMENT_SERIES_NUMBER);
          {
            // DIA-NN also writes and accepts True/False here.
            const auto d = get(c.dec);
            if (!d.empty() && (d[0] == 't' || d[0] == 'T')) { r.decoy = 1; }
            else if (!d.empty() && (d[0] == 'f' || d[0] == 'F')) { r.decoy = 0; }
            else { r.decoy = toLong(d, DIANNLibraryFile::Columns::DECOY); }
          }
          r.fragment_type = parseFragmentType(get(c.ft));
          r.loss_type = parseLossType(get(c.lt));
          // Once per run of rows with the same sequence -- a precursor's
          // transitions -- instead of once per transition.
          if (!have_residues || r.modified_sequence != residue_sequence)
          {
            residue_sequence = r.modified_sequence;
            residues = residueCount(r.modified_sequence);
            have_residues = true;
          }
          r.residues = residues;
          if (c.id >= 0) { r.precursor_id = get(c.id); }
          else
          {
            // The id storeTSV would have written, so a library without the
            // column groups its rows exactly as one with it. A precursor's
            // rows share one copy.
            const bool same = !range.rows.empty() && !id_spans.empty() &&
                              range.rows.back().modified_sequence == r.modified_sequence &&
                              range.rows.back().precursor_charge == r.precursor_charge &&
                              (range.rows.back().decoy != 0) == (r.decoy != 0);
            if (same) { id_spans.push_back(id_spans.back()); }
            else
            {
              const std::size_t at = range.ids.size();
              range.ids.append(r.modified_sequence);
              range.ids += std::to_string(r.precursor_charge);
              if (r.decoy != 0) { range.ids += "_decoy"; }
              id_spans.emplace_back(at, range.ids.size() - at);
            }
          }
          checkPrecursorCharge(r.precursor_charge);
          range.rows.push_back(r);
        }
      }
      catch (...)
      {
        range.error = std::current_exception();
        range.error_line = range.lines - 1;
      }
      if (c.id < 0)
      {
        for (std::size_t i = 0; i < range.rows.size(); ++i)
        { range.rows[i].precursor_id = std::string_view(range.ids).substr(id_spans[i].first, id_spans[i].second); }
      }
    }

    /// The last '\n' in [begin, end), or nullptr.
    const char* lastNewline(const char* begin, const char* end)
    {
      for (const char* p = end; p != begin; --p)
      {
        if (p[-1] == '\n') { return p - 1; }
      }
      return nullptr;
    }

    /// A block of a TSV library's text, cut into ranges and parsed: whole
    /// lines in [0, complete), the start of the next line after that.
    struct TsvBlock
    {
      std::vector<char> text;
      std::size_t have = 0;
      std::size_t complete = 0;
      std::size_t used = 0;           // ranges in use
      std::vector<TsvRange> ranges;
    };

    /// Reads the next block -- the partial line @p previous ended with, then
    /// as much of @p in as fits -- and parses its ranges on @p workers
    /// threads. Leaves `have` 0 at the end of the file.
    void readTsvBlock(std::istream& in, const std::string& filename, TsvBlock& block,
                      const TsvBlock* previous, const TsvLayout& layout, unsigned workers,
                      std::size_t range_bytes)
    {
      block.have = 0;
      block.complete = 0;
      block.used = 0;
      if (previous != nullptr && previous->have > previous->complete)
      {
        const std::size_t tail = previous->have - previous->complete;
        if (block.text.size() < tail * 2) { block.text.resize(std::max(block.text.size(), tail) * 2); }
        std::memcpy(block.text.data(), previous->text.data() + previous->complete, tail);
        block.have = tail;
      }
      for (;;)
      {
        bool eof = !in;
        if (!eof && block.have < block.text.size())
        {
          in.read(block.text.data() + block.have, static_cast<std::streamsize>(block.text.size() - block.have));
          block.have += static_cast<std::size_t>(in.gcount());
          if (!in)
          {
            if (in.bad()) { throw std::runtime_error("error reading library: " + filename); }
            eof = true;
          }
        }
        if (block.have == 0) { return; }
        if (eof) { block.complete = block.have; break; }
        const char* nl = lastNewline(block.text.data(), block.text.data() + block.have);
        if (nl != nullptr) { block.complete = static_cast<std::size_t>(nl - block.text.data()) + 1; break; }
        // A line longer than the whole buffer: read more of it.
        block.text.resize(block.text.size() * 2);
      }

      const char* const end = block.text.data() + block.complete;
      for (const char* p = block.text.data(); p < end; ++block.used)
      {
        const char* cut = end;
        if (static_cast<std::size_t>(end - p) > range_bytes)
        {
          const char* from = p + range_bytes - 1;
          const auto* nl = static_cast<const char*>(std::memchr(from, '\n', static_cast<std::size_t>(end - from)));
          cut = nl != nullptr ? nl + 1 : end;
        }
        if (block.used == block.ranges.size()) { block.ranges.emplace_back(); }
        block.ranges[block.used].begin = p;
        block.ranges[block.used].end = cut;
        p = cut;
      }
      parallel::forEach(block.used, workers, [&](std::size_t i) { parseTsvRange(block.ranges[i], layout); });
    }
  } // namespace

  void DIANNLibraryFile::loadTSV(const std::string& filename, Library& library, unsigned threads)
  {
    std::ifstream in(filename, std::ios::binary);
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

    // First present name wins. The TSVColumns spelling comes first for what
    // storeTSV writes now, the dotted Columns spelling for what it wrote up to
    // 0.11.1, then DIA-NN's other aliases for libraries from elsewhere.
    auto col = [&](std::initializer_list<const char*> names) {
      for (const char* name : names)
      {
        const auto it = index.find(name);
        if (it != index.end()) { return it->second; }
      }
      return -1;
    };
    using T = TSVColumns;
    const int c_id = col({T::PRECURSOR_ID});
    const int c_seq = col({T::MODIFIED_PEPTIDE, Columns::MODIFIED_SEQUENCE,
                           "FullUniModPeptideName", "ModifiedPeptideSequence"});
    // Protein.Group holds DIALibGen's verbatim group; a library from elsewhere
    // names its proteins only in DIA-NN's accession column.
    const int c_pg = col({T::PROTEIN_GROUP, "Protein.Ids", T::PROTEIN_IDS, "UniprotId",
                          "UniProtIds", "ProteinId", "ProteinID"});
    const int c_pmz = col({T::PRECURSOR_MZ, Columns::PRECURSOR_MZ, "Q1"});
    const int c_qmz = col({T::PRODUCT_MZ, Columns::PRODUCT_MZ, "FragmentMz", "Q3"});
    const int c_rt = col({T::RT, Columns::RT, "Tr_recalibrated", "NormalizedRetentionTime",
                          "RetentionTime"});
    const int c_im = col({T::IM, Columns::IM, "PrecursorIonMobility", "Ion.Mobility"});
    const int c_ccs = col({T::CCS});
    const int c_int = col({T::RELATIVE_INTENSITY, Columns::RELATIVE_INTENSITY,
                           "RelativeIntensity", "RelativeFragmentIntensity"});
    const int c_z = col({T::PRECURSOR_CHARGE, Columns::PRECURSOR_CHARGE});
    const int c_fz = col({T::FRAGMENT_CHARGE, Columns::FRAGMENT_CHARGE});
    const int c_ord = col({T::FRAGMENT_SERIES_NUMBER, Columns::FRAGMENT_SERIES_NUMBER,
                           "FragmentNumber"});
    const int c_dec = col({T::DECOY, "decoy"});
    const int c_ft = col({T::FRAGMENT_TYPE, Columns::FRAGMENT_TYPE});
    const int c_lt = col({T::FRAGMENT_LOSS_TYPE, Columns::FRAGMENT_LOSS_TYPE});

    if (c_pmz < 0 || c_qmz < 0)
    {
      throw std::runtime_error("not a DIA-NN library (missing PrecursorMz / ProductMz): " + filename);
    }
    if (c_seq < 0) { throw std::runtime_error("missing Modified.Sequence (or ModifiedPeptide) in library: " + filename); }

    TsvLayout layout;
    layout.width = index.size();
    layout.id = c_id; layout.seq = c_seq; layout.pg = c_pg; layout.pmz = c_pmz; layout.qmz = c_qmz;
    layout.rt = c_rt; layout.im = c_im; layout.ccs = c_ccs; layout.intensity = c_int; layout.z = c_z;
    layout.fz = c_fz; layout.ord = c_ord; layout.dec = c_dec; layout.ft = c_ft; layout.lt = c_lt;

    // Parsing is the cost -- splitting, six numbers and four integers per
    // row -- and it is row-local, so the text is read in blocks, cut into
    // ranges at newlines, and the ranges are parsed on all workers. The rows
    // then enter the library in file order on one thread, which is all
    // Builder::append needs to group, intern and number them exactly as a
    // serial read does. A range keeps the rows before its first bad line, so
    // the error and everything appended before it are the serial reader's.
    //
    // With more than one thread, appending block k (one thread) overlaps
    // reading and parsing block k+1 (the others); appending is about half of
    // what is left once parsing is parallel.
    const unsigned workers = parallel::workers(threads);
    const bool overlap = workers > 1;
    const unsigned parsers = overlap ? workers - 1 : 1;
    const std::size_t range_bytes = tsvRangeBytes();
    const std::size_t ranges_per_block = std::clamp<std::size_t>(std::size_t(workers) * 2, 2, 64);
    TsvBlock blocks[2];
    for (auto& block : blocks) { block.text.resize(range_bytes * ranges_per_block); }

    Builder builder(library);
    std::size_t lines_before = 1;  // the header
    const auto append = [&](const TsvBlock& block) {
      for (std::size_t i = 0; i < block.used; ++i)
      {
        const auto& range = block.ranges[i];
        for (const auto& row : range.rows) { builder.append(row); }
        if (range.width_error)
        {
          throw std::runtime_error("TSV row width differs from header at line " +
                                   std::to_string(lines_before + range.error_line + 1));
        }
        if (range.error) { std::rethrow_exception(range.error); }
        lines_before += range.lines;
      }
    };

    readTsvBlock(in, filename, blocks[0], nullptr, layout, parsers, range_bytes);
    for (std::size_t k = 0; blocks[k % 2].have != 0; ++k)
    {
      TsvBlock& current = blocks[k % 2];
      TsvBlock& next = blocks[(k + 1) % 2];
      if (!overlap)
      {
        append(current);
        readTsvBlock(in, filename, next, &current, layout, parsers, range_bytes);
        continue;
      }
      std::exception_ptr append_error;
      std::thread appender;
      try { appender = std::thread([&] { try { append(current); } catch (...) { append_error = std::current_exception(); } }); }
      catch (const std::system_error&) { append(current); }  // no thread to spare: in turn
      // The block before is the one a serial read fails in first.
      std::exception_ptr read_error;
      try { readTsvBlock(in, filename, next, &current, layout, parsers, range_bytes); }
      catch (...) { read_error = std::current_exception(); }
      if (appender.joinable()) { appender.join(); }
      if (append_error) { std::rethrow_exception(append_error); }
      if (read_error) { std::rethrow_exception(read_error); }
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

  void DIANNLibraryFile::loadParquet(const std::string& filename, Library& library, unsigned threads)
  {
    auto infile = arrow::io::ReadableFile::Open(filename);
    if (!infile.ok()) { throw std::runtime_error("cannot open library: " + filename); }

    // The default reader decodes a batch's columns one after another on the
    // calling thread. use_threads spreads them over Arrow's CPU pool, which
    // useThreads() sizes to -threads; pre_buffer fetches each row group's
    // column chunks in a few coalesced reads instead of one per page.
    parquet::ArrowReaderProperties properties = parquet::default_arrow_reader_properties();
    properties.set_use_threads(parallel::workers(threads) > 1);
    properties.set_pre_buffer(true);
    std::unique_ptr<parquet::arrow::FileReader> reader;
    {
      parquet::arrow::FileReaderBuilder open;
      if (!open.Open(*infile).ok() ||
          !open.memory_pool(arrow::default_memory_pool())->properties(properties)->Build(&reader).ok() || !reader)
      { throw std::runtime_error("not a Parquet file: " + filename); }
    }

    // Size the library once from the footer instead of doubling every array
    // its way up: rows are transitions in the flat layout; in the compact one
    // they are precursors, and the values of the Product.Mz list are the
    // transitions. Capacity only -- the content does not depend on it.
    {
      const auto metadata = reader->parquet_reader()->metadata();
      std::shared_ptr<arrow::Schema> schema;
      if (metadata && reader->GetSchema(&schema).ok() && schema)
      {
        const auto rows = static_cast<std::size_t>(std::max<std::int64_t>(0, metadata->num_rows()));
        const int qi = schema->GetFieldIndex(Columns::PRODUCT_MZ);
        const bool compact = qi >= 0 && schema->field(qi)->type()->id() == arrow::Type::LIST;
        std::size_t precursors = 0, transitions = rows;
        if (compact)
        {
          precursors = rows;
          transitions = 0;
          const auto* columns = metadata->schema();
          for (int j = 0; j < columns->num_columns(); ++j)
          {
            const auto path = columns->Column(j)->path()->ToDotVector();
            if (path.empty() || path.front() != Columns::PRODUCT_MZ) { continue; }
            for (int g = 0; g < metadata->num_row_groups(); ++g)
            { transitions += static_cast<std::size_t>(std::max<std::int64_t>(0, metadata->RowGroup(g)->ColumnChunk(j)->num_values())); }
            break;
          }
        }
        const std::size_t limit = std::numeric_limits<std::uint32_t>::max();
        const std::size_t have = library.transitionCount();
        transitions = have >= limit ? have : std::min(have + transitions, limit);
        library.reserve(library.precursorCount() + precursors, transitions);
      }
    }

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
      if (!a_seq) { throw std::runtime_error("missing Modified.Sequence in library: " + filename); }

      const int64_t rows = table->num_rows();
      for (int64_t i = 0; i < rows; ++i)
      {
        Row r;
        r.precursor_id = sv(a_id, i);
        r.modified_sequence = sv(a_seq, i);
        r.protein_group = sv(a_pg, i);
        r.fragment_type = parseFragmentType(sv(a_ft, i));
        r.loss_type = parseLossType(sv(a_lt, i));
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
    trace::Scope hashing("hash_model", 1);
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
    trace::Scope hashing("hash_fasta", 1);
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
    // loadParquet hands over one record batch, i.e. one chunk per column, so
    // nothing needs combining; anything else is combined first.
    std::shared_ptr<arrow::Table> table = table_in;
    for (const auto& column : table_in->columns())
    {
      if (column->num_chunks() == 1) { continue; }
      auto combined = table_in->CombineChunks(arrow::default_memory_pool());
      if (!combined.ok()) { throw std::runtime_error("cannot combine Parquet chunks"); }
      table = *combined;
      break;
    }

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

      const int64_t b = l_qmz->value_offset(i), e = l_qmz->value_offset(i + 1);
      const std::size_t residues = v_ord ? residueCount(seq) : 0;
      Library::checkTransitionCapacity(t.product_mz.size(), static_cast<std::size_t>(e - b));
      p.transition_begin.push_back(static_cast<std::uint32_t>(t.product_mz.size()));
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
        if (v_ord && (ordinal < 1 || ordinal >= residues))
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

  namespace
  {
    auto okStatus = [](const arrow::Status& st) {
      if (!st.ok()) { throw std::runtime_error("parquet build: " + st.ToString()); }
    };

    /// One precursor row to write: its fields in a library's precursor arrays,
    /// its transitions in @p transitions. The whole library (begin and count
    /// from the CSR) and a block of MS2 results held apart from it (a
    /// streaming generate) are both written through this.
    struct WriteRow
    {
      std::uint32_t precursor;
      const Library::TransitionArrays* transitions;
      std::uint32_t begin;
      std::uint32_t count;
    };

    std::vector<WriteRow> libraryRows(const Library& library, std::size_t first, std::size_t last)
    {
      const auto& p = library.precursors();
      std::vector<WriteRow> rows;
      rows.reserve(last - first);
      for (std::size_t i = first; i < last; ++i)
      {
        rows.push_back({static_cast<std::uint32_t>(i), &library.transitions(), p.transition_begin[i],
                        p.transition_count[i]});
      }
      return rows;
    }

    // Precursor-level: ONE value each, not one per transition.
    //
    // Plain strings, not dictionary builders. The flat writer below needs its
    // dictionary builders (one row per transition repeats every string ~12x),
    // but here each Precursor.Id occurs once. A dictionary builder produced ONE
    // global dictionary of every string in the library, and Parquet stores an
    // Arrow dictionary as the dictionary page of EVERY row group, so every
    // reader decoded every string of the library once per row group. With
    // plain strings, Parquet dictionary-encodes each row group only until the
    // dictionary page reaches dictionary_pagesize_limit (1 MiB by default) and
    // writes the rest of that row group as PLAIN.
    //
    // Arrays are finished per row group (compactChunkRows()). A plain
    // StringArray holds at most 2 GiB of characters, so a chunk fails loudly
    // (CapacityError, and no partial file) if a string column averages 2 KiB or
    // more per row. Real libraries hold at most ~30 MB per string column and
    // chunk, ~1.4 % of it.
    /// Rows per row group: 2^20, or ODIA_ROW_GROUP_ROWS (tests: several row
    /// groups from a small library).
    std::size_t compactChunkRows()
    {
      static const std::size_t rows = [] {
        if (const char* v = std::getenv("ODIA_ROW_GROUP_ROWS"); v != nullptr && *v != '\0')
        {
          char* end = nullptr;
          const long long x = std::strtoll(v, &end, 10);
          if (end != v && x > 0) { return static_cast<std::size_t>(x); }
        }
        return std::size_t(1) << 20;
      }();
      return rows;
    }
    constexpr std::size_t kCompactColumns = 15;
    using CompactChunk = std::array<std::shared_ptr<arrow::Array>, kCompactColumns>;

    /// The arrays of @p count rows, one per column in the order of
    /// compactSchema()'s fields; the static_assert there checks the count,
    /// ChunkedArray/RecordBatch the types.
    CompactChunk buildCompactChunk(const Library& library, const WriteRow* rows, std::size_t count)
    {
      auto pool = arrow::default_memory_pool();
      const auto& p = library.precursors();
      arrow::StringBuilder b_id(pool), b_seq(pool), b_pg(pool);
      arrow::UInt8Builder b_z(pool), b_dec(pool);
      arrow::FloatBuilder b_rt(pool), b_im(pool), b_ccs(pool);
      arrow::DoubleBuilder b_pmz(pool);

      // Transition-level: a LIST per precursor. The list offsets replace the
      // repeated precursor key entirely -- there is no foreign key to store.
      auto qmz_v = std::make_shared<arrow::DoubleBuilder>(pool);
      auto int_v = std::make_shared<arrow::FloatBuilder>(pool);
      auto ft_v  = std::make_shared<arrow::UInt8Builder>(pool);
      auto fz_v  = std::make_shared<arrow::UInt8Builder>(pool);
      auto ord_v = std::make_shared<arrow::UInt8Builder>(pool);
      auto loss_v = std::make_shared<arrow::UInt8Builder>(pool);
      arrow::ListBuilder b_qmz(pool, qmz_v), b_int(pool, int_v),
                         b_ft(pool, ft_v), b_fz(pool, fz_v), b_ord(pool, ord_v), b_loss(pool, loss_v);
      arrow::ArrayBuilder* const builders[] = {&b_id, &b_seq, &b_pg, &b_z, &b_dec, &b_rt, &b_im, &b_ccs, &b_pmz,
                                               &b_qmz, &b_int, &b_ft, &b_fz, &b_ord, &b_loss};
      static_assert(std::extent_v<decltype(builders)> == kCompactColumns, "one builder per column");
      const auto ok = okStatus;

      for (std::size_t r = 0; r < count; ++r)
      {
        const auto& row = rows[r];
        const std::size_t i = row.precursor;
        const auto& t = *row.transitions;
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
        for (std::uint32_t k = 0; k < row.count; ++k)
        {
          const std::uint32_t j = row.begin + k;
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

      CompactChunk arrays;
      for (std::size_t c = 0; c < kCompactColumns; ++c)
      {
        if (builders[c]->length() != static_cast<std::int64_t>(count))
        { throw std::logic_error("compact Parquet columns differ in length"); }
        ok(builders[c]->Finish(&arrays[c]));
      }
      return arrays;
    }

    std::shared_ptr<arrow::Schema> compactSchema(const DIANNLibraryFile::Fingerprint& fp,
                                                 const std::string& config_json,
                                                 const std::vector<std::pair<std::string, std::string>>& extra_metadata)
    {
      using Columns = DIANNLibraryFile::Columns;
      const std::shared_ptr<arrow::Field> fields[] = {
        arrow::field(Columns::PRECURSOR_ID, arrow::utf8()),
        arrow::field(Columns::MODIFIED_SEQUENCE, arrow::utf8()),
        arrow::field(Columns::PROTEIN_GROUP, arrow::utf8()),
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
      };
      static_assert(std::extent_v<decltype(fields)> == kCompactColumns, "one builder per schema field");
      auto schema = arrow::schema(arrow::FieldVector(std::begin(fields), std::end(fields)));
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
      if (!extra_metadata.empty())
      {
        auto md = schema->metadata()->Copy();
        for (const auto& [key, value] : extra_metadata) { md->Append(key, value); }
        schema = schema->WithMetadata(md);
      }
      return schema;
    }

    std::shared_ptr<parquet::WriterProperties> compactWriterProperties()
    {
      return parquet::WriterProperties::Builder()
               .compression(parquet::Compression::ZSTD)
               ->enable_dictionary()
               ->build();
    }

    /// The row groups through Arrow's buffered writer, which encodes and
    /// compresses a row group's columns in parallel (Arrow's CPU pool, sized
    /// to @p writers); the same bytes as WriteTable with compactChunkRows().
    std::unique_ptr<parquet::arrow::FileWriter>
    openCompactWriter(const std::shared_ptr<arrow::Schema>& schema,
                      const std::shared_ptr<arrow::io::FileOutputStream>& outfile, unsigned writers)
    {
      if (writers > 1 && ::arrow::GetCpuThreadPoolCapacity() != static_cast<int>(writers))
      { okStatus(::arrow::SetCpuThreadPoolCapacity(static_cast<int>(writers))); }
      auto arrow_props = parquet::ArrowWriterProperties::Builder().store_schema()->set_use_threads(writers > 1)->build();
      auto opened = parquet::arrow::FileWriter::Open(*schema, arrow::default_memory_pool(), outfile,
                                                     compactWriterProperties(), arrow_props);
      if (!opened.ok()) { throw std::runtime_error("cannot write Parquet: " + opened.status().ToString()); }
      return std::move(*opened);
    }

    /// One row group. Buffered: Arrow encodes the columns in parallel and
    /// holds the whole encoded row group until the next one starts. Direct:
    /// column after column, each streamed to the file as it is encoded, so
    /// only one column's pages are held -- what the streaming writer uses,
    /// since it works while MS2 holds its peak memory. The bytes are the same.
    void writeCompactRowGroup(parquet::arrow::FileWriter& writer, const std::shared_ptr<arrow::Schema>& schema,
                              const CompactChunk& chunk, bool buffered = true)
    {
      const std::vector<std::shared_ptr<arrow::Array>> arrays(chunk.begin(), chunk.end());
      const auto rows = arrays.front()->length();
      if (buffered)
      {
        const auto batch = arrow::RecordBatch::Make(schema, rows, arrays);
        okStatus(writer.NewBufferedRowGroup());
        okStatus(writer.WriteRecordBatch(*batch));
        return;
      }
      okStatus(writer.WriteTable(*arrow::Table::Make(schema, arrays, rows), rows));
    }

    /// Builds the row groups of library rows [first, last) (a multiple of
    /// compactChunkRows() from first) in parallel, then writes them in order.
    void writeCompactLibraryRows(parquet::arrow::FileWriter& writer, const std::shared_ptr<arrow::Schema>& schema,
                                 const Library& library, std::size_t first, std::size_t last, unsigned threads,
                                 bool buffered = true)
    {
      const std::size_t chunk_count = (last - first + compactChunkRows() - 1) / compactChunkRows();
      std::vector<CompactChunk> built(chunk_count);
      parallel::forEach(chunk_count, parallel::workers(threads), [&](std::size_t c)
      {
        const std::size_t a = first + c * compactChunkRows(), b = std::min(last, a + compactChunkRows());
        const auto rows = libraryRows(library, a, b);
        built[c] = buildCompactChunk(library, rows.data(), rows.size());
      });
      for (auto& chunk : built) { writeCompactRowGroup(writer, schema, chunk, buffered); chunk = CompactChunk{}; }
    }
  } // namespace

  void DIANNLibraryFile::useThreads(unsigned threads)
  {
    const auto workers = static_cast<int>(parallel::workers(threads));
    if (::arrow::GetCpuThreadPoolCapacity() != workers) { okStatus(::arrow::SetCpuThreadPoolCapacity(workers)); }
    const int io = std::min(workers, 8);
    if (::arrow::io::GetIOThreadPoolCapacity() != io) { okStatus(::arrow::io::SetIOThreadPoolCapacity(io)); }
  }

  void DIANNLibraryFile::storeParquetCompact(const std::string& filename,
                                             const Library& library,
                                             const Fingerprint& fp,
                                             const std::string& config_json,
                                             const std::vector<std::pair<std::string, std::string>>& extra_metadata,
                                             unsigned threads)
  {
    const std::size_t n = library.precursorCount();
    auto pool = arrow::default_memory_pool();
    const auto schema = compactSchema(fp, config_json, extra_metadata);

    trace::Scope write_step("write_table", 1);
    write_step.fields().rows = static_cast<std::int64_t>(n);
    AtomicFile output(filename);
    auto outfile = arrow::io::FileOutputStream::Open(output.temporaryPath().string());
    if (!outfile.ok()) { throw std::runtime_error("cannot write library: " + filename); }
    const unsigned writers = parallel::workers(threads);
    if (writers <= 1 || n == 0)
    {
      // Chunks are independent: the table is the same as from one pass.
      std::vector<arrow::ArrayVector> chunks(kCompactColumns);
      for (std::size_t first = 0; first < n; first += compactChunkRows())
      {
        const auto rows = libraryRows(library, first, std::min(n, first + compactChunkRows()));
        const auto arrays = buildCompactChunk(library, rows.data(), rows.size());
        for (std::size_t c = 0; c < kCompactColumns; ++c) { chunks[c].push_back(arrays[c]); }
      }
      std::vector<std::shared_ptr<arrow::ChunkedArray>> columns;
      for (std::size_t c = 0; c < kCompactColumns; ++c)
      {
        auto column = arrow::ChunkedArray::Make(chunks[c], schema->field(static_cast<int>(c))->type());
        okStatus(column.status());
        columns.push_back(*column);
      }
      auto table = arrow::Table::Make(schema, columns, static_cast<std::int64_t>(n));
      auto arrow_props = parquet::ArrowWriterProperties::Builder().store_schema()->build();
      const auto st = parquet::arrow::WriteTable(*table, pool, *outfile, static_cast<std::int64_t>(compactChunkRows()),
                                                 compactWriterProperties(), arrow_props);
      if (!st.ok()) { throw std::runtime_error("cannot write Parquet: " + st.ToString()); }
    }
    else
    {
      auto writer = openCompactWriter(schema, *outfile, writers);
      writeCompactLibraryRows(*writer, schema, library, 0, n, threads);
      okStatus(writer->Close());
    }
    okStatus((*outfile)->Close());
    write_step.end();
    trace::Scope commit_step("commit", 1);
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

  namespace
  {
    void writeTsvHeader(TextWriter& out)
    {
      using T = DIANNLibraryFile::TSVColumns;
      bool first_column = true;
      for (const char* column : {T::PRECURSOR_ID, T::MODIFIED_PEPTIDE, T::PRECURSOR_CHARGE, T::DECOY,
                                 T::RT, T::IM, T::PRECURSOR_MZ, T::PRODUCT_MZ,
                                 T::RELATIVE_INTENSITY, T::FRAGMENT_TYPE, T::FRAGMENT_CHARGE,
                                 T::FRAGMENT_SERIES_NUMBER, T::FRAGMENT_LOSS_TYPE,
                                 T::PROTEIN_IDS, T::PROTEIN_NAMES, T::PROTEOTYPIC,
                                 T::PROTEIN_GROUP, T::CCS})
      {
        if (!first_column) { out.put('\t'); }
        out.put(column);
        first_column = false;
      }
      out.put('\n');
    }

    /// DIA-NN's protein columns from DIALibGen's group "sp|P1|N1_HUMAN;sp|P2|N2_HUMAN":
    /// accessions "P1;P2", entry names "N1_HUMAN;N2_HUMAN". A member without the
    /// db|accession|name form stands for itself in both. Stateless, so
    /// precursors can be rendered on several threads.
    std::size_t tsvProteins(std::string_view group, std::string& ids, std::string& names)
    {
      ids.clear(); names.clear();
      std::size_t members = 0;
      std::size_t pos = 0;
      while (pos < group.size())
      {
        std::size_t end = group.find(';', pos);
        if (end == std::string_view::npos) { end = group.size(); }
        const auto member = group.substr(pos, end - pos);
        pos = end + 1;
        if (member.empty()) { continue; }
        std::string_view accession = member, name = member;
        const auto bar1 = member.find('|');
        const auto bar2 = bar1 == std::string_view::npos ? bar1 : member.find('|', bar1 + 1);
        if (bar2 != std::string_view::npos)
        {
          accession = member.substr(bar1 + 1, bar2 - bar1 - 1);
          name = member.substr(bar2 + 1);
        }
        if (members++) { ids += ';'; names += ';'; }
        ids += accession; names += name;
      }
      return members;
    }

    /// One precursor's rows, into a TextWriter or a TextBuffer alike.
    ///
    /// The fields every row of a precursor repeats -- the leading ones up to
    /// its m/z and the trailing ones from the protein columns on -- are
    /// formatted once into @p head and @p tail and copied per row: number
    /// formatting is the writer's cost, and these are four of its seven
    /// numbers. The same digits, so the same bytes.
    template <class Sink>
    void writeTsvPrecursor(Sink& sink, const Library& library, const WriteRow& row,
                           std::atomic<std::size_t>& unknown_mods, TextBuffer& head, TextBuffer& tail)
    {
      const auto& p = library.precursors();
      const auto& t = *row.transitions;
      const std::size_t i = row.precursor;
      const auto stored = library.strings().get(p.modified_sequence[i]);
      const auto pg = library.strings().get(p.protein_group[i]);
      if (stored.find_first_of("\t\r\n\0", 0, 4) != std::string_view::npos ||
          pg.find_first_of("\t\r\n\0", 0, 4) != std::string_view::npos)
      { throw std::runtime_error("TSV sequence or protein group contains a tab, newline or NUL"); }
      // DIA-NN parses UniMod accessions; OpenMS names such as Carbamidomethyl
      // are not what either version expects.
      std::size_t unknown = 0;
      const std::string seq = diannModifiedSequence(stored, &unknown);
      if (unknown) { unknown_mods += unknown; }
      std::string ids, names;
      const std::size_t members = tsvProteins(pg, ids, names);
      const int z = p.charge[i];
      const bool has_im = i < p.im.size() && !std::isnan(p.im[i]);
      const bool has_ccs = i < p.ccs.size() && !std::isnan(p.ccs[i]);

      head.clear();
      // DIA-NN gives a decoy the same Modified.Sequence, charge and
      // Precursor.Mz as its target, so <sequence><charge> alone is not unique.
      // Without the suffix, reloading our own output merged each such pair and
      // re-labelled the decoy's transitions as target evidence -- 33,390
      // precursors became 33,386, silently, with no transitions lost.
      head.put(seq); head.integer(z); head.put(p.decoy[i] ? "_decoy" : ""); head.put('\t');
      head.put(seq); head.put('\t'); head.integer(z); head.put('\t');
      head.integer(p.decoy[i]); head.put('\t');
      // An unpredicted retention time is written as an empty field, not as
      // "nan": no TSV consumer accepts the latter, and it propagates into
      // anything that reads the library back.
      if (!std::isnan(p.irt[i])) { head.number(p.irt[i], 9); }
      head.put('\t');
      if (has_im) { head.number(p.im[i], 9); }
      head.put('\t');
      head.number(fromFixed(p.mz[i]), 10); head.put('\t');

      tail.clear();
      tail.put(ids); tail.put('\t');
      tail.put(names); tail.put('\t');
      tail.put(members == 1 ? "1" : "0"); tail.put('\t');
      tail.put(pg); tail.put('\t');
      // Empty rather than "nan" when absent, as for RT: no TSV consumer
      // accepts the latter, and a 0 here would read as a real cross-section.
      if (has_ccs) { tail.number(p.ccs[i], 9); }
      tail.put('\n');

      for (std::uint32_t k = 0; k < row.count; ++k)
      {
        const std::uint32_t j = row.begin + k;
        sink.put(head.text());
        sink.number(fromFixed(t.product_mz[j]), 10); sink.put('\t');
        if (!std::isnan(t.library_intensity[j])) { sink.number(t.library_intensity[j], 9); }
        sink.put('\t');
        sink.put(toString(t.type[j])); sink.put('\t');
        sink.integer(t.charge[j]); sink.put('\t');
        sink.integer(t.ordinal[j]); sink.put('\t');
        sink.put(toString(t.loss[j])); sink.put('\t');
        sink.put(tail.text());
      }
    }

    /// Rendering numbers is the cost (see store()), so blocks of precursors are
    /// rendered in parallel and written in order: the same bytes as one pass.
    /// The rendered text held at once is bounded by a fixed number of
    /// precursors (WAVE), not by the thread count: four tasks per thread grew
    /// it to ~4 GB at 96 threads.
    void writeTsvRows(TextWriter& out, const Library& library, const std::vector<WriteRow>& rows,
                      unsigned threads, std::atomic<std::size_t>& unknown_mods)
    {
      const unsigned workers = parallel::workers(threads);
      const std::size_t n = rows.size();
      if (workers <= 1)
      {
        TextBuffer head, tail;
        for (const auto& row : rows) { writeTsvPrecursor(out, library, row, unknown_mods, head, tail); }
        return;
      }
      constexpr std::size_t TASK = 2048;
      constexpr std::size_t WAVE = 128 * TASK;   // ~262k precursors, ~0.7 GB of text
      std::vector<TextBuffer> rendered;
      for (std::size_t first = 0; first < n; first += WAVE)
      {
        const std::size_t count = std::min(WAVE, n - first);
        rendered.assign((count + TASK - 1) / TASK, TextBuffer{});
        parallel::forEach(rendered.size(), workers, [&](std::size_t task)
        {
          TextBuffer head, tail;
          const std::size_t end = first + std::min(count, (task + 1) * TASK);
          for (std::size_t i = first + task * TASK; i < end; ++i)
          { writeTsvPrecursor(rendered[task], library, rows[i], unknown_mods, head, tail); }
        });
        for (auto& block : rendered) { out.put(block.text()); block = TextBuffer{}; }
      }
    }

    void warnUnknownMods(std::size_t unknown_mods)
    {
      if (unknown_mods)
      {
        std::cerr << "warning: " << unknown_mods << " modification tokens have no UniMod accession and "
                     "were written verbatim; DIA-NN may not read those precursors\n";
      }
    }
  } // namespace

  void DIANNLibraryFile::storeTSV(const std::string& filename, const Library& library,
                                   unsigned threads)
  {
    // Formatting (in parallel blocks) and the ordered, buffered writes.
    trace::Scope format_step("format_rows", 1);
    format_step.fields().rows = static_cast<std::int64_t>(library.transitionCount());
    AtomicFile output(filename);
    TextWriter out(output.temporaryPath().string());
    writeTsvHeader(out);
    std::atomic<std::size_t> unknown_mods{0};
    writeTsvRows(out, library, libraryRows(library, 0, library.precursorCount()), threads, unknown_mods);
    format_step.end();
    // A file checked only at open reports success on a full disk, an exceeded
    // quota or a broken mount, leaving a truncated library behind.
    trace::Scope close_step("close", 1);
    out.close();
    close_step.end();
    trace::Scope commit_step("commit", 1);
    output.commit();
    warnUnknownMods(unknown_mods.load());
  }

  // ------------------------------------------------------------ StreamWriter

  struct DIANNLibraryFile::StreamWriter::Impl
  {
    const Library& library;
    unsigned threads;
    bool parquet;
    AtomicFile output;
    // Parquet
    std::shared_ptr<arrow::Schema> schema;
    std::shared_ptr<arrow::io::FileOutputStream> outfile;
    std::unique_ptr<parquet::arrow::FileWriter> writer;
    std::deque<Library::PrecursorRows> pending;
    std::size_t pending_offset = 0;   ///< rows of pending.front() already written
    std::size_t pending_rows = 0;     ///< rows in pending, minus pending_offset
    // TSV
    std::unique_ptr<TextWriter> text;
    std::atomic<std::size_t> unknown_mods{0};

    std::size_t written = 0;
    bool finished = false;

    Impl(const std::string& filename, const Library& lib, unsigned t)
      : library(lib), threads(t), parquet(filename.ends_with(".parquet")), output(filename) {}
  };

  DIANNLibraryFile::StreamWriter::StreamWriter(const std::string& filename, const Library& library,
                                               const Fingerprint& fp, const std::string& config_json,
                                               const std::vector<std::pair<std::string, std::string>>& extra_metadata,
                                               unsigned threads)
    : impl_(std::make_unique<Impl>(filename, library, threads))
  {
    auto& m = *impl_;
    if (m.parquet)
    {
      m.schema = compactSchema(fp, config_json, extra_metadata);
      auto opened = arrow::io::FileOutputStream::Open(m.output.temporaryPath().string());
      if (!opened.ok()) { throw std::runtime_error("cannot write library: " + filename); }
      m.outfile = *opened;
      // Direct row groups, one column at a time (writeCompactRowGroup).
      m.writer = openCompactWriter(m.schema, m.outfile, 1);
    }
    else
    {
      m.text = std::make_unique<TextWriter>(m.output.temporaryPath().string());
      writeTsvHeader(*m.text);
    }
  }

  DIANNLibraryFile::StreamWriter::~StreamWriter() = default;

  std::size_t DIANNLibraryFile::StreamWriter::written() const { return impl_->written; }

  void DIANNLibraryFile::StreamWriter::add(Library::PrecursorRows rows)
  {
    auto& m = *impl_;
    if (m.finished) { throw std::logic_error("StreamWriter::add after finish"); }
    if (rows.precursor.empty()) { return; }
    trace::Scope streamed("stream_write", 1);
    if (!m.parquet)
    {
      std::vector<WriteRow> list;
      list.reserve(rows.precursor.size());
      for (std::size_t k = 0; k < rows.precursor.size(); ++k)
      { list.push_back({rows.precursor[k], &rows.transitions, rows.begin[k], rows.count[k]}); }
      writeTsvRows(*m.text, m.library, list, m.threads, m.unknown_mods);
      m.written += list.size();
      streamed.fields().rows = static_cast<std::int64_t>(list.size());
      return;
    }
    // Parquet: whole row groups only, so the file is cut where the one-pass
    // writer cuts it.
    m.pending_rows += rows.precursor.size();
    m.pending.push_back(std::move(rows));
    std::int64_t done = 0;
    while (m.pending_rows >= compactChunkRows())
    {
      // A row group can span several blocks: it is built before any of them
      // is dropped.
      std::vector<WriteRow> list;
      list.reserve(compactChunkRows());
      std::size_t block = 0, offset = m.pending_offset;
      while (list.size() < compactChunkRows())
      {
        auto& b = m.pending[block];
        const std::size_t take = std::min(compactChunkRows() - list.size(), b.precursor.size() - offset);
        for (std::size_t k = offset; k < offset + take; ++k)
        { list.push_back({b.precursor[k], &b.transitions, b.begin[k], b.count[k]}); }
        offset += take;
        if (offset == b.precursor.size()) { ++block; offset = 0; }
      }
      const auto chunk = buildCompactChunk(m.library, list.data(), list.size());
      writeCompactRowGroup(*m.writer, m.schema, chunk, false);
      for (std::size_t k = 0; k < block; ++k) { m.pending.pop_front(); }
      m.pending_offset = offset;
      m.pending_rows -= compactChunkRows();
      m.written += compactChunkRows();
      done += static_cast<std::int64_t>(compactChunkRows());
    }
    streamed.fields().rows = done;
  }

  void DIANNLibraryFile::StreamWriter::finish(const Library& final_library)
  {
    auto& m = *impl_;
    if (m.finished) { throw std::logic_error("StreamWriter::finish twice"); }
    m.finished = true;
    const std::size_t n = final_library.precursorCount();
    if (n < m.written) { throw std::logic_error("the final library is shorter than what was streamed"); }
    // Rows held for an incomplete row group are taken from the final library,
    // where they are the same rows: everything streamed is its first rows.
    m.pending.clear();
    trace::Scope write_step("write_tail", 1);
    write_step.fields().rows = static_cast<std::int64_t>(n - m.written);
    if (m.parquet)
    {
      writeCompactLibraryRows(*m.writer, m.schema, final_library, m.written, n, m.threads, false);
      okStatus(m.writer->Close());
      okStatus(m.outfile->Close());
    }
    else
    {
      writeTsvRows(*m.text, final_library, libraryRows(final_library, m.written, n), m.threads, m.unknown_mods);
      m.text->close();
    }
    write_step.end();
    trace::Scope commit_step("commit", 1);
    m.output.commit();
    if (!m.parquet) { warnUnknownMods(m.unknown_mods.load()); }
  }

} // namespace ODIA
