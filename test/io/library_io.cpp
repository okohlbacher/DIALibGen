#include <odia/DIANNLibraryFile.h>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
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
