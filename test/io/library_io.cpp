#include <odia/DIANNLibraryFile.h>
#include "../../src/library/Parallel.h"
#include <parquet/file_reader.h>
#include <parquet/metadata.h>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <numeric>
#include <stdexcept>

using namespace ODIA;
namespace fs = std::filesystem;
void check(bool condition, const std::string& message) { if (!condition) { throw std::runtime_error(message); } }
Library fixture()
{
  Library lib;
  auto& p = lib.precursors(); auto& t = lib.transitions();
  const float absent = std::numeric_limits<float>::quiet_NaN();
  p.mz = {toFixed(600.12345), toFixed(600.12345), toFixed(450.98765)};
  p.irt = {15.123456f, 15.123456f, absent};
  p.im = {1.2345678f, 1.2345678f, absent}; p.ccs = {423.12345f, 423.12345f, absent};
  p.charge = {2, 2, 3}; p.decoy = {0, 1, 0};
  for (const auto* s : {"AC(UniMod:4)DEFGK", "AC(UniMod:4)DEFGK", "PEPTIDER"}) { p.modified_sequence.push_back(lib.strings().intern(s)); }
  for (const auto* s : {"P001;P002", "P001;P002", "P003"}) { p.protein_group.push_back(lib.strings().intern(s)); }
  p.transition_begin = {0, 2, 4}; p.transition_count = {2, 2, 1};
  t.product_mz = {toFixed(300.12345), toFixed(400.56789), toFixed(310.12345), toFixed(410.56789), toFixed(220.98765)};
  t.library_intensity = {1.0f, 0.31415927f, 1.0f, 0.31415927f, 0.25f};
  t.type = {FragmentType::B, FragmentType::Y, FragmentType::B, FragmentType::Y, FragmentType::Y};
  t.charge = {1, 2, 1, 2, 1}; t.ordinal = {2, 4, 2, 4, 3};
  t.loss = {LossType::Water, LossType::None, LossType::Ammonia, LossType::None, LossType::None};
  return lib;
}
// One precursor per row and a 2^20-row writer chunk: n = 2^20 + 1 puts the
// last precursor alone in a second row group. Unique ids, as in a real library.
Library chunked(std::size_t n)
{
  Library lib;
  auto& p = lib.precursors(); auto& t = lib.transitions();
  const std::string residues = "ACDEFGHIKLMNPQRSTVWY";
  for (std::size_t i = 0; i < n; ++i)
  {
    std::string seq = "PEP";
    for (std::size_t k = i; k > 0; k /= residues.size()) { seq += residues[k % residues.size()]; }
    seq += 'K';
    p.modified_sequence.push_back(lib.strings().intern(seq));
    p.protein_group.push_back(lib.strings().intern("P" + std::to_string(i % 5000)));
    p.mz.push_back(toFixed(400.0 + 0.001 * static_cast<double>(i % 100000)));
    p.irt.push_back(i % 7 == 0 ? std::nanf("") : static_cast<float>(i % 1000) / 10.0f);
    // IM and CCS absent together: load() derives either one from the other.
    p.im.push_back(i % 11 == 0 ? std::nanf("") : 0.8f + static_cast<float>(i % 500) / 1000.0f);
    p.ccs.push_back(i % 11 == 0 ? std::nanf("") : 300.0f + static_cast<float>(i % 400));
    p.charge.push_back(static_cast<std::uint8_t>(2 + i % 3)); p.decoy.push_back(static_cast<std::uint8_t>(i % 2));
    p.transition_begin.push_back(static_cast<std::uint32_t>(i)); p.transition_count.push_back(1);
    t.product_mz.push_back(toFixed(200.0 + 0.25 * static_cast<double>(i % 997)));
    t.library_intensity.push_back(1.0f); t.type.push_back(FragmentType::Y);
    t.charge.push_back(1); t.ordinal.push_back(1); t.loss.push_back(LossType::None);
  }
  return lib;
}
std::string bytes(const fs::path& path)
{
  std::ifstream in(path, std::ios::binary);
  return {std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>()};
}
std::vector<std::int64_t> rowGroups(const fs::path& path)
{
  const auto metadata = parquet::ParquetFileReader::OpenFile(path.string())->metadata();
  std::vector<std::int64_t> rows;
  for (int g = 0; g < metadata->num_row_groups(); ++g) { rows.push_back(metadata->RowGroup(g)->num_rows()); }
  return rows;
}
void compare(const Library& actual, const Library& expected)
{
  const auto& a = actual.precursors(); const auto& b = expected.precursors();
  check(a.mz == b.mz && a.charge == b.charge && a.decoy == b.decoy, "precursor identity/mass changed");
  check(a.transition_begin == b.transition_begin && a.transition_count == b.transition_count, "precursor transition ranges changed");
  auto floats = [](const auto& x, const auto& y) {
    if (x.size() != y.size()) { return false; }
    for (std::size_t i = 0; i < x.size(); ++i)
    { if (!(x[i] == y[i] || (std::isnan(x[i]) && std::isnan(y[i])))) { return false; } }
    return true;
  };
  check(floats(a.irt, b.irt) && floats(a.im, b.im) && floats(a.ccs, b.ccs), "RT/IM/CCS or missing values changed");
  for (std::size_t i = 0; i < actual.precursorCount(); ++i)
  {
    check(actual.strings().get(a.modified_sequence[i]) == expected.strings().get(b.modified_sequence[i]), "modified sequence changed");
    check(actual.strings().get(a.protein_group[i]) == expected.strings().get(b.protein_group[i]), "protein group changed");
  }
  const auto& x = actual.transitions(); const auto& y = expected.transitions();
  check(x.product_mz == y.product_mz && x.library_intensity == y.library_intensity, "transition mass/intensity changed");
  check(x.type == y.type && x.charge == y.charge && x.ordinal == y.ordinal && x.loss == y.loss, "fragment annotation changed");
}
int main(int argc, char** argv)
{
  try
  {
    check(argc >= 3, "usage: library_io fixtures|load|roundtrip PATH");
    const std::string mode(argv[1]); const fs::path path(argv[2]);
    if (mode == "load")
    {
      Library lib; DIANNLibraryFile::load(path.string(), lib);
      std::cout << lib.precursorCount() << ' ' << lib.transitionCount() << ' ' << lib.decoyCount() << ' '
                << lib.invalidMzCount() << ' ' << lib.invalidMzTransitionCount() << '\n';
      return 0;
    }
    // The compact writer itself, not pyarrow, across its chunk boundary and empty.
    if (mode == "compact-chunks" || mode == "compact-empty")
    {
      const auto expected = chunked(mode == "compact-chunks" ? (std::size_t(1) << 20) + 1 : 0);
      DIANNLibraryFile::storeParquetCompact(path.string(), expected, {});
      const auto rows = rowGroups(path);
      if (mode == "compact-chunks")
      { check(rows == std::vector<std::int64_t>{1 << 20, 1}, "compact row groups are not [1048576, 1]"); }
      else { check(std::accumulate(rows.begin(), rows.end(), std::int64_t{0}) == 0, "empty compact library has rows"); }
      Library loaded; DIANNLibraryFile::load(path.string(), loaded);
      check(loaded.precursorCount() == expected.precursorCount(), "compact precursor count changed");
      compare(loaded, expected);
      return 0;
    }
    // The writers render rows (TSV) and build row groups (Parquet) on several
    // threads; the files must be the bytes a single thread writes. Three row
    // groups, the last a partial one, and more TSV blocks than threads.
    if (mode == "threads-identical")
    {
      fs::create_directories(path);
      const auto library = chunked((std::size_t(2) << 20) + 12345);
      for (const unsigned threads : {1u, 3u, 8u})
      {
        const std::string suffix = std::to_string(threads);
        DIANNLibraryFile::storeParquetCompact((path / ("t" + suffix + ".parquet")).string(), library, {}, "{}", {}, threads);
        DIANNLibraryFile::storeTSV((path / ("t" + suffix + ".tsv")).string(), library, threads);
      }
      for (const char* ext : {".parquet", ".tsv"})
      {
        const auto reference = bytes(path / (std::string("t1") + ext));
        check(!reference.empty(), std::string("empty ") + ext);
        for (const char* other : {"t3", "t8"})
        { check(bytes(path / (other + std::string(ext))) == reference, std::string(other) + ext + " differs from one thread"); }
      }
      check(rowGroups(path / "t8.parquet").size() == 3, "expected three row groups");
      // A failure reports the row a serial writer would have hit first.
      auto bad = chunked((std::size_t(2) << 20) + 12345);
      for (const std::size_t row : {std::size_t(1500000), std::size_t(70000)})
      { bad.precursors().protein_group[row] = bad.strings().intern("bad\tgroup " + std::to_string(row)); }
      try { DIANNLibraryFile::storeTSV((path / "bad.tsv").string(), bad, 8); check(false, "a tab was written"); }
      catch (const std::runtime_error& e) { check(std::string(e.what()).find("contains a tab") != std::string::npos, e.what()); }
      check(!fs::exists(path / "bad.tsv"), "a failed TSV write left a file");
      return 0;
    }
    // parallel::forEach: every index once, and the lowest failure wins.
    if (mode == "parallel-for-each")
    {
      for (const unsigned threads : {1u, 2u, 7u, 64u})
      {
        std::vector<int> hits(10007, 0);
        parallel::forEach(hits.size(), threads, [&](std::size_t i) { ++hits[i]; });
        check(std::all_of(hits.begin(), hits.end(), [](int h) { return h == 1; }), "an index ran other than once");
        try
        {
          parallel::forEach(hits.size(), threads, [&](std::size_t i)
          { if (i == 9000 || i == 4242 || i == 5000) { throw std::runtime_error(std::to_string(i)); } });
          check(false, "no exception");
        }
        catch (const std::runtime_error& e) { check(std::string(e.what()) == "4242", std::string("rethrew ") + e.what()); }
      }
      return 0;
    }
    auto original = fixture();
    if (mode == "write")
    {
      if (path.filename().string().starts_with("compact"))
      { DIANNLibraryFile::storeParquetCompact(path.string(), original, {}); }
      else { DIANNLibraryFile::store(path.string(), original); }
      return 0;
    }
    if (mode == "missing-values")
    {
      original.precursors().im.clear(); original.precursors().ccs.clear();
      original.transitions().library_intensity[0] = std::nanf("");
      DIANNLibraryFile::storeTSV(path.string(), original); return 0;
    }
    // OpenMS modification names and UniProt-style protein groups, as a
    // generated library holds them; the TSV must carry DIA-NN's spelling.
    if (mode == "openms-names")
    {
      auto& p = original.precursors();
      p.modified_sequence = {original.strings().intern("AC(Carbamidomethyl)DEFGK"),
                             original.strings().intern("AC(Carbamidomethyl)DEFGK"),
                             original.strings().intern(".(Acetyl)PEPM(Oxidation)TIDER")};
      const auto group = original.strings().intern("sp|P00001|ONE_HUMAN;tr|P00002|TWO_HUMAN");
      p.protein_group[0] = group; p.protein_group[1] = group;
      DIANNLibraryFile::storeTSV(path.string(), original); return 0;
    }
    if (mode == "invalid-text")
    {
      original.precursors().protein_group[1] = original.strings().intern("bad\tgroup\n");
      DIANNLibraryFile::storeTSV(path.string(), original); return 0;
    }
    if (mode == "roundtrip") { Library loaded; DIANNLibraryFile::load(path.string(), loaded); compare(loaded, original); return 0; }
    check(mode == "fixtures", "unknown operation"); fs::create_directories(path);
    DIANNLibraryFile::Fingerprint fp{"abcdef", 42, "v4:recipe", "v4:targets", "mutate"};
    DIANNLibraryFile::storeTSV((path / "library.tsv").string(), original);
    DIANNLibraryFile::storeParquet((path / "flat.parquet").string(), original, fp);
    DIANNLibraryFile::storeParquetCompact((path / "compact.parquet").string(), original, fp, "{\"example\":true}");
    for (const auto* name : {"output.tsv", "output.parquet", "compact.parquet"})
    {
      bool refused = false;
      try
      {
        const auto missing = (path / "missing-parent" / name).string();
        if (std::string_view(name) == "compact.parquet") { DIANNLibraryFile::storeParquetCompact(missing, original, fp); }
        else { DIANNLibraryFile::store(missing, original); }
      }
      catch (const std::runtime_error&) { refused = true; }
      check(refused, "writer accepted an inaccessible output path");
    }
    for (const auto* name : {"flat.parquet", "compact.parquet"})
    { check(DIANNLibraryFile::readFingerprint((path / name).string()) == fp.key(), "fingerprint lost"); }
    check(DIANNLibraryFile::readFingerprint((path / "compact.parquet").string(), "odia.target_fingerprint") == fp.targetKey(), "target fingerprint lost");
    check(DIANNLibraryFile::readFingerprint((path / "missing").string()).empty(), "missing file reported fingerprint");
    check(DIANNLibraryFile::readFingerprint((path / "library.tsv").string()).empty(), "TSV reported fingerprint");
    const auto fasta = path / "fingerprint.fasta";
    for (const auto& [content, expected] : {std::pair{"", "cbf29ce484222325"}, std::pair{"a", "af63dc4c8601ec8c"}})
    {
      { std::ofstream out(fasta); out << content; }
      check(DIANNLibraryFile::hashFile(fasta.string()) == expected, "FNV-1a hash differs from known vector");
      check(DIANNLibraryFile::fingerprintFasta(fasta.string()).fasta_hash == expected, "FASTA FNV-1a differs from known vector");
    }
    { std::ofstream out(fasta); out << ">P\nPEPTIDER\n"; }
    const auto old = DIANNLibraryFile::fingerprintFasta(fasta.string());
    fs::copy_file(fasta, path / "moved.fasta", fs::copy_options::overwrite_existing);
    check(old.key() == DIANNLibraryFile::fingerprintFasta((path / "moved.fasta").string()).key(), "fingerprint depends on file path");
    { std::ofstream out(fasta); out << ">P\nPEPTIDEK\n"; }
    check(old.key() != DIANNLibraryFile::fingerprintFasta(fasta.string()).key(), "same-size content change reused fingerprint");
    auto lib = fixture();
    check(lib.dropDecoys() == 1 && lib.precursorCount() == 2 && lib.transitionCount() == 3, "decoy removal corrupts arrays");
    bool refused = false;
    try { (void)lib.lowerBound(450); } catch (const std::exception&) { refused = true; }
    check(refused, "unsorted library accepted binary search");
    lib.sortByPrecursorMz();
    check(lib.lowerBound(450.99) == 1 && lib.transitions().product_mz.front() == toFixed(220.98765), "sort did not preserve transitions");
    auto optional = fixture(); optional.precursors().im.clear(); optional.precursors().ccs.clear();
    auto subset = optional.subsetByIndex({2, 0});
    check(subset.precursors().im.size() == 2 && std::isnan(subset.precursors().im[1]), "subset lost missing IM");
    optional.sortByPrecursorMz();
    check(std::isnan(optional.precursors().im[0]) && std::isnan(optional.precursors().ccs[0]), "sort lost missing mobility");
    optional.precursors().im.clear(); optional.precursors().ccs.resize(1);
    optional.dropDecoys();
    check(optional.precursors().im.size() == 2 && std::isnan(optional.precursors().im[1]) && std::isnan(optional.precursors().ccs[1]), "decoy removal lost missing mobility");
    return 0;
  }
  catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
