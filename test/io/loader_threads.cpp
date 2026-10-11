// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

// The library loaders parse on several threads; the library they build must
// be the one a single thread builds, bit for bit, and an error must name the
// line a single thread stops at. Also: the number parser that replaced
// istringstream returns the stream's doubles.

#include <odia/DIANNLibraryFile.h>
#include "../../src/library/ParseDouble.h"

#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <locale>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace ODIA;
namespace fs = std::filesystem;

namespace
{
  void check(bool condition, const std::string& message) { if (!condition) { throw std::runtime_error(message); } }

  /// FNV-1a over everything a Library holds: every array, element by element
  /// in its bit pattern, and the text of every interned string.
  class Hash
  {
  public:
    void bytes(const void* data, std::size_t n)
    {
      const auto* p = static_cast<const unsigned char*>(data);
      for (std::size_t i = 0; i < n; ++i) { h_ ^= p[i]; h_ *= 1099511628211ull; }
    }
    template <class T> void array(const std::vector<T>& v)
    {
      const std::uint64_t n = v.size();
      bytes(&n, sizeof n);
      if (!v.empty()) { bytes(v.data(), v.size() * sizeof(T)); }
    }
    void text(std::string_view s)
    {
      const std::uint64_t n = s.size();
      bytes(&n, sizeof n);
      bytes(s.data(), s.size());
    }
    std::uint64_t value() const { return h_; }
  private:
    std::uint64_t h_ = 14695981039346656037ull;
  };

  std::uint64_t contentHash(const Library& lib)
  {
    Hash h;
    const auto& p = lib.precursors();
    h.array(p.mz); h.array(p.irt); h.array(p.im); h.array(p.ccs); h.array(p.charge); h.array(p.decoy);
    h.array(p.transition_begin); h.array(p.transition_count);
    h.array(p.modified_sequence); h.array(p.protein_group);  // the handles, i.e. the interning order
    for (std::size_t i = 0; i < lib.precursorCount(); ++i)
    {
      h.text(lib.strings().get(p.modified_sequence[i]));
      h.text(lib.strings().get(p.protein_group[i]));
    }
    const auto& t = lib.transitions();
    h.array(t.product_mz); h.array(t.library_intensity); h.array(t.type); h.array(t.ordinal);
    h.array(t.charge); h.array(t.loss);
    return h.value();
  }

  /// A library with the variety a real one has: modified sequences, shared
  /// protein groups, decoys, absent RT/IM/CCS, 1-20 transitions of every
  /// type, charge and loss, and m/z and intensities with full precision.
  Library synthetic(std::size_t precursors)
  {
    Library lib;
    auto& p = lib.precursors(); auto& t = lib.transitions();
    std::mt19937_64 rng(20261011);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    const std::string residues = "ACDEFGHIKLMNPQRSTVWY";
    const char* mods[] = {"", "M(UniMod:35)", "C(UniMod:4)", "S(UniMod:21)"};
    const LossType losses[] = {LossType::None, LossType::Water, LossType::Ammonia, LossType::Phospho};
    for (std::size_t i = 0; i < precursors; ++i)
    {
      std::string seq;
      const std::size_t len = 7 + rng() % 18;
      for (std::size_t k = 0; k < len; ++k) { seq += residues[rng() % residues.size()]; }
      seq += mods[rng() % 4];
      seq += 'K';
      p.modified_sequence.push_back(lib.strings().intern(seq));
      // One line far longer than a read buffer at a small ODIA_TSV_RANGE_BYTES.
      p.protein_group.push_back(lib.strings().intern(
        i == 1234 ? std::string(20000, 'G') : "P" + std::to_string(rng() % 3000) + ";Q" + std::to_string(i % 17)));
      p.mz.push_back(toFixed(350.0 + 1000.0 * unit(rng)));
      p.irt.push_back(i % 13 == 0 ? std::nanf("") : static_cast<float>(-20.0 + 160.0 * unit(rng)));
      const bool mobile = i % 9 != 0;
      p.im.push_back(mobile ? static_cast<float>(0.6 + unit(rng)) : std::nanf(""));
      p.ccs.push_back(mobile ? static_cast<float>(280.0 + 500.0 * unit(rng)) : std::nanf(""));
      p.charge.push_back(static_cast<std::uint8_t>(1 + rng() % 4));
      p.decoy.push_back(static_cast<std::uint8_t>(i % 2));
      p.transition_begin.push_back(static_cast<std::uint32_t>(t.product_mz.size()));
      const std::size_t n = 1 + rng() % 20;
      for (std::size_t k = 0; k < n; ++k)
      {
        t.product_mz.push_back(toFixed(100.0 + 1500.0 * unit(rng)));
        t.library_intensity.push_back(static_cast<float>(unit(rng)));
        t.type.push_back(rng() % 2 ? FragmentType::Y : FragmentType::B);
        t.charge.push_back(static_cast<std::int8_t>(1 + rng() % 2));
        t.ordinal.push_back(static_cast<std::uint8_t>(1 + rng() % (len - 1)));
        t.loss.push_back(losses[rng() % 4]);
      }
      p.transition_count.push_back(static_cast<std::uint32_t>(n));
    }
    return lib;
  }

  std::uint64_t loadHash(const fs::path& file, unsigned threads)
  {
    Library lib;
    DIANNLibraryFile::load(file.string(), lib, threads);
    check(lib.precursorCount() > 0, file.string() + ": nothing loaded");
    return contentHash(lib);
  }

  std::string loadError(const fs::path& file, unsigned threads)
  {
    try { Library lib; DIANNLibraryFile::load(file.string(), lib, threads); }
    catch (const std::exception& e) { return e.what(); }
    return "(no error)";
  }

  std::vector<std::string> lines(const fs::path& file)
  {
    std::ifstream in(file, std::ios::binary);
    std::vector<std::string> out;
    for (std::string line; std::getline(in, line);) { out.push_back(line); }
    return out;
  }

  void writeLines(const fs::path& file, const std::vector<std::string>& rows, const char* eol)
  {
    std::ofstream out(file, std::ios::binary);
    for (const auto& row : rows) { out << row << eol; }
  }

  std::string field(const std::string& line, std::size_t index)
  {
    std::size_t start = 0;
    for (std::size_t i = 0; i < index; ++i) { start = line.find('\t', start) + 1; }
    return line.substr(start, line.find('\t', start) - start);
  }

  std::string setField(const std::string& line, std::size_t index, const std::string& value)
  {
    std::size_t start = 0;
    for (std::size_t i = 0; i < index; ++i) { start = line.find('\t', start) + 1; }
    const std::size_t end = line.find('\t', start);
    return line.substr(0, start) + value + (end == std::string::npos ? std::string() : line.substr(end));
  }

  std::size_t column(const std::string& header, const std::string& name)
  {
    std::size_t index = 0, start = 0;
    for (;;)
    {
      const std::size_t tab = header.find('\t', start);
      if (header.substr(start, tab - start) == name) { return index; }
      check(tab != std::string::npos, "no column " + name);
      start = tab + 1;
      ++index;
    }
  }

  /// The same library through every loader at 1, 4 and all threads; a
  /// one-thread load is the serial reference.
  int loadersIdentical(const fs::path& dir)
  {
    fs::create_directories(dir);
    const auto library = synthetic(30000);
    const fs::path tsv = dir / "library.tsv", compact = dir / "compact.parquet", flat = dir / "flat.parquet";
    DIANNLibraryFile::storeTSV(tsv.string(), library);
    DIANNLibraryFile::storeParquetCompact(compact.string(), library, {});
    DIANNLibraryFile::storeParquet(flat.string(), library, {});

    // A TSV as other tools write it: CRLF, blank lines, no Precursor.Id (the
    // id is synthesized), True/False decoys, a '+' sign, "nan", and rows whose
    // ordinal cannot exist (dropped, and counted).
    auto text = lines(tsv);
    const std::string header = text.front();
    const std::size_t c_id = column(header, "Precursor.Id"), c_dec = column(header, "Decoy");
    const std::size_t c_rt = column(header, "iRT"), c_ord = column(header, "FragmentSeriesNumber");
    const std::size_t c_pmz = column(header, "PrecursorMz");
    std::vector<std::string> foreign;
    for (std::size_t i = 0; i < text.size(); ++i)
    {
      std::string row = text[i];
      if (i > 0)
      {
        row = setField(row, c_dec, field(row, c_dec) == "1" ? "True" : "False");
        if (i % 7 == 0) { row = setField(row, c_pmz, "+" + field(row, c_pmz)); }
        if (i % 11 == 0) { row = setField(row, c_rt, "nan"); }
        if (i % 101 == 0) { row = setField(row, c_ord, "250"); }
      }
      // Drop the Precursor.Id column.
      std::size_t start = 0;
      for (std::size_t k = 0; k < c_id; ++k) { start = row.find('\t', start) + 1; }
      const std::size_t end = row.find('\t', start);
      row.erase(start, end == std::string::npos ? std::string::npos : end - start + 1);
      foreign.push_back(row);
      if (i % 997 == 0) { foreign.emplace_back(); }
    }
    const fs::path other = dir / "foreign.tsv";
    writeLines(other, foreign, "\r\n");

    for (const auto& file : {tsv, compact, flat, other})
    {
      const auto reference = loadHash(file, 1);
      for (const unsigned threads : {4u, 0u})
      {
        check(loadHash(file, threads) == reference,
              file.filename().string() + ": " + std::to_string(threads) + " threads differ from one");
      }
      std::cout << file.filename().string() << ' ' << std::hex << std::setw(16) << std::setfill('0')
                << reference << std::dec << '\n';
    }
    // The two Parquet layouts carry the same values in the same types.
    check(loadHash(compact, 1) == loadHash(flat, 1), "compact and flat Parquet load differently");

    // Errors: the first bad line of the file wins, whichever check it fails
    // and whichever range parsed it.
    const std::size_t c_z = column(header, "PrecursorCharge");
    auto broken = text;
    const std::size_t charge_line = text.size() / 3, number_line = text.size() / 2, width_line = text.size() - 5;
    broken[charge_line] = setField(broken[charge_line], c_z, "300");
    broken[number_line] = setField(broken[number_line], c_pmz, "12.3.4");
    broken[width_line] += "\textra";
    const fs::path bad = dir / "bad.tsv";
    const std::vector<std::pair<std::string, std::size_t>> expected = {
      {"invalid precursor charge: expected an integer between 1 and 255", charge_line},
      {"invalid numeric value in Precursor.Mz: 12.3.4", number_line},
      {"TSV row width differs from header at line " + std::to_string(width_line + 1), width_line},
    };
    for (const auto& [message, line] : expected)
    {
      writeLines(bad, broken, "\n");
      for (const unsigned threads : {1u, 4u, 0u})
      {
        const auto got = loadError(bad, threads);
        check(got == message, "threads " + std::to_string(threads) + ": expected \"" + message + "\", got \"" + got + "\"");
      }
      broken[line] = text[line];  // fixed: the next one is now the first
    }
    // The fixed file loads as the original.
    writeLines(bad, broken, "\n");
    check(loadHash(bad, 4) == loadHash(tsv, 1), "a repaired file loads differently");
    return 0;
  }

  /// ODIA::parseDouble against a classic-locale istringstream, the
  /// parser it replaced: the same accept/reject and the same bits.
  int parseDoubleMatchesStream()
  {
    auto stream = [](const std::string& s, double& v) {
      if (s.empty() || s.find_first_not_of("0123456789+-.eE") != std::string::npos) { return false; }
      std::istringstream in{s};
      in.imbue(std::locale::classic());
      in >> std::noskipws >> v;
      return in && in.peek() == std::char_traits<char>::eof() && std::isfinite(v);
    };
    std::vector<std::string> cases = {
      "0", "-0", "+0", "1", "+1", "-1", "1.", ".5", "-.5", "+.5", "1e5", "1E5", "1e+05", "1e-05", "00.100",
      "123456789012345678901234567890", "0.1", "0.30000000000000004", "1.7976931348623157e308",
      "2.2250738585072014e-308", "1e309", "-1e309", "1e", "e5", ".", "+", "-", "", "+-1", "-+1", "--1", "++1",
      "1-2", "1e5e3", "1..2", "1.2.3", "12.3.4", " 1", "1 ", "nan", "inf", "0x10", "1,5"};
    std::mt19937_64 rng(7);
    for (int i = 0; i < 200000; ++i)
    {
      // Shortest round-trip, fixed and scientific forms, normal range only:
      // how a stream treats an underflow differs between C++ libraries.
      const double x = std::ldexp(std::uniform_real_distribution<double>(0.5, 1.0)(rng),
                                  static_cast<int>(rng() % 1600) - 800) * (rng() % 2 ? 1 : -1);
      std::ostringstream a, b, c;
      a.imbue(std::locale::classic()); b.imbue(std::locale::classic()); c.imbue(std::locale::classic());
      a << std::setprecision(17) << x;
      b << std::fixed << std::setprecision(static_cast<int>(rng() % 12)) << std::fmod(x, 1e6);
      c << std::scientific << std::setprecision(static_cast<int>(rng() % 20)) << x;
      cases.push_back(a.str()); cases.push_back(b.str()); cases.push_back(c.str());
      // Arbitrary digit strings: more digits than a double holds.
      std::string d;
      const int digits = 1 + static_cast<int>(rng() % 30);
      for (int k = 0; k < digits; ++k) { d += static_cast<char>('0' + rng() % 10); }
      d.insert(rng() % (d.size() + 1), ".");
      if (rng() % 3 == 0) { d += "e" + std::to_string(static_cast<int>(rng() % 80) - 40); }
      cases.push_back(d);
    }
    for (const auto& s : cases)
    {
      double expected = 0.0, got = 0.0;
      const bool ok_stream = stream(s, expected);
      const bool ok_parse = parseDouble(s, got) && std::isfinite(got);
      check(ok_stream == ok_parse, "\"" + s + "\": stream " + (ok_stream ? "accepts" : "rejects") +
                                   ", parseDouble " + (ok_parse ? "accepts" : "rejects"));
      if (ok_stream)
      {
        check(std::memcmp(&expected, &got, sizeof(double)) == 0, "\"" + s + "\": different double");
      }
    }
    std::cout << cases.size() << " strings, parser " << (ODIA_FLOAT_FROM_CHARS ? "from_chars" : "strtod_l") << '\n';
    return 0;
  }
} // namespace

int main(int argc, char** argv)
{
  try
  {
    check(argc >= 2, "usage: loader_threads loaders-identical DIR | parse-double");
    const std::string mode(argv[1]);
    if (mode == "parse-double") { return parseDoubleMatchesStream(); }
    check(mode == "loaders-identical" && argc >= 3, "unknown mode");
    return loadersIdentical(argv[2]);
  }
  catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
