// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

// An mzPeak run and an mzML of the same data load into the same OpenSWATH maps.
//
//   identify_mzpeak_run <run.mzpeak> <run.mzML> [threads] [normal|cache] [--report]
//
// Both go through Identifier::loadRun -- the mzML through stock SwathFile, the
// archive through the vendored mzpeak-openms reader and the same SwathFile
// consumers -- and every map is compared: its window (m/z bounds and centre,
// 1/K0 band, MS1 or not), and every spectrum's RT, MS level and arrays (m/z,
// intensity, 1/K0 with its name), value for value. Exit 0 when the arrays are
// identical and the scalars (RT, window bounds) agree within 1e-12 relative
// (the mzML's text, parsed by OpenMS), 1 otherwise; each kind of difference
// is counted and its largest size printed. With --report it exits 0 either
// way: the differences are the output (the full-run equivalence checks).

#include <odia/search/Identifier.h>
#include <odia/search/MzPeakRun.h>

#include <OpenMS/OPENSWATHALGO/DATAACCESS/ISpectrumAccess.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <vector>

using namespace ODIA::search;

namespace
{
  class Loader : public Identifier
  {
  public:
    explicit Loader(SearchParams p) : Identifier(std::move(p), [](const std::string& m) { std::cout << "  info: " << m << "\n"; },
                                                 [](const std::string& m) { std::cout << "  warn: " << m << "\n"; }) {}
    using Identifier::loadRun;
  };

  struct Tally
  {
    std::uint64_t compared = 0, differing = 0;
    double largest = 0.0, largest_relative = 0.0;
    void see(double a, double b)
    {
      ++compared;
      if (a == b || (std::isnan(a) && std::isnan(b))) { return; }
      ++differing;
      largest = std::max(largest, std::abs(a - b));
      largest_relative = std::max(largest_relative, std::abs(a - b) / std::max(std::abs(a), std::abs(b)));
    }
  };

  /// Scalars may differ in the last bits: an mzML carries them as text, and
  /// OpenMS's number parser does not always round the 17 digits an mzML
  /// writer prints to the nearest double. Peak arrays are binary in both and
  /// must be equal.
  constexpr double scalar_relative_tolerance = 1e-12;
  bool exact(const std::string& name) { return name == "m/z" || name == "intensity" || name == "1/K0"; }

  std::map<std::string, Tally> tallies;
  std::uint64_t structural = 0;

  void structure(const std::string& what)
  {
    if (structural++ < 20) { std::cout << "  differs: " << what << "\n"; }
  }

  void compareArrays(const std::string& name, const OpenSwath::BinaryDataArrayPtr& a, const OpenSwath::BinaryDataArrayPtr& b,
                     const std::string& where)
  {
    if (!a && !b) { return; }
    if (!a || !b) { structure(where + ": " + name + " array only in " + (a ? "mzPeak" : "mzML")); return; }
    if (a->data.size() != b->data.size())
    {
      structure(where + ": " + name + " " + std::to_string(a->data.size()) + " vs " + std::to_string(b->data.size()) + " values");
      return;
    }
    Tally& t = tallies[name];
    for (std::size_t k = 0; k < a->data.size(); ++k) { t.see(a->data[k], b->data[k]); }
  }
}

int main(int argc, char** argv)
{
  if (argc < 3)
  {
    std::cerr << "usage: identify_mzpeak_run <run.mzpeak> <run.mzML> [threads] [normal|cache] [--report]\n";
    return 2;
  }
  if (!mzpeak::available())
  {
    std::cout << "identify_mzpeak_run: built without mzPeak support; skipped\n";
    return 77;
  }
  const std::string mzpeak_path = argv[1], mzml_path = argv[2];
  const int threads = argc > 3 ? std::stoi(argv[3]) : 2;
  const std::string mode = argc > 4 ? argv[4] : "normal";
  const bool report = argc > 5 && std::string(argv[5]) == "--report";
  try
  {
    SearchParams p;
    p.threads = threads;
    p.readoptions = mode == "cache" ? ReadMode::Cache : ReadMode::Normal;
    const std::filesystem::path scratch = std::filesystem::temp_directory_path();
    p.cache_dir = scratch.string();

    auto timed = [&](const std::string& path) {
      Loader loader(p);
      const auto t0 = std::chrono::steady_clock::now();
      RunData r = loader.loadRun(path);
      const double s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
      std::cout << "loaded " << path << " in " << std::fixed << std::setprecision(1) << s << " s: " << r.maps.size() << " maps, "
                << r.read_mode << ", ion mobility " << (r.ion_mobility ? "yes" : "no") << "\n";
      return r;
    };
    RunData a = timed(mzpeak_path);
    RunData b = timed(mzml_path);

    if (a.ion_mobility != b.ion_mobility || a.ms1_ion_mobility != b.ms1_ion_mobility) { structure("ion-mobility flags"); }
    if (a.maps.size() != b.maps.size())
    { structure("maps: " + std::to_string(a.maps.size()) + " vs " + std::to_string(b.maps.size())); }
    std::uint64_t spectra = 0;
    for (std::size_t m = 0; m < std::min(a.maps.size(), b.maps.size()); ++m)
    {
      const auto& x = a.maps[m];
      const auto& y = b.maps[m];
      const std::string where = "map " + std::to_string(m);
      if (x.ms1 != y.ms1) { structure(where + ": MS1 flag"); continue; }
      tallies["window lower"].see(x.lower, y.lower);
      tallies["window upper"].see(x.upper, y.upper);
      tallies["window center"].see(x.center, y.center);
      tallies["window 1/K0 lower"].see(x.imLower, y.imLower);
      tallies["window 1/K0 upper"].see(x.imUpper, y.imUpper);
      const std::size_t n = x.sptr ? x.sptr->getNrSpectra() : 0, n2 = y.sptr ? y.sptr->getNrSpectra() : 0;
      if (n != n2) { structure(where + ": " + std::to_string(n) + " vs " + std::to_string(n2) + " spectra"); continue; }
      for (std::size_t k = 0; k < n; ++k)
      {
        const auto mx = x.sptr->getSpectrumMetaById(static_cast<int>(k));
        const auto my = y.sptr->getSpectrumMetaById(static_cast<int>(k));
        tallies["spectrum RT"].see(mx.RT, my.RT);
        if (mx.ms_level != my.ms_level) { structure(where + " spectrum " + std::to_string(k) + ": MS level"); }
        const OpenSwath::SpectrumPtr sx = x.sptr->getSpectrumById(static_cast<int>(k));
        const OpenSwath::SpectrumPtr sy = y.sptr->getSpectrumById(static_cast<int>(k));
        const std::string at = where + " spectrum " + std::to_string(k);
        compareArrays("m/z", sx->getMZArray(), sy->getMZArray(), at);
        static int shown = 0;
        if (shown < 3 && sx->getDriftTimeArray() && sy->getDriftTimeArray() &&
            sx->getDriftTimeArray()->data.size() == sy->getDriftTimeArray()->data.size() &&
            sx->getDriftTimeArray()->data != sy->getDriftTimeArray()->data)
        {
          // The first points where the 1/K0 order differs: (m/z, intensity, 1/K0) on both sides.
          ++shown;
          const auto& xi = sx->getDriftTimeArray()->data;
          std::size_t first = 0;
          while (first < xi.size() && xi[first] == sy->getDriftTimeArray()->data[first]) { ++first; }
          std::cout << "  " << at << ": first 1/K0 difference at point " << first << "\n";
          for (std::size_t q = first >= 2 ? first - 2 : 0; q < std::min(xi.size(), first + 4); ++q)
          {
            std::cout << std::setprecision(10) << "    " << sx->getMZArray()->data[q] << " " << sx->getIntensityArray()->data[q] << " "
                      << xi[q] << "  |  " << sy->getMZArray()->data[q] << " " << sy->getIntensityArray()->data[q] << " "
                      << sy->getDriftTimeArray()->data[q] << "\n";
          }
        }
        compareArrays("intensity", sx->getIntensityArray(), sy->getIntensityArray(), at);
        compareArrays("1/K0", sx->getDriftTimeArray(), sy->getDriftTimeArray(), at);
        if (sx->getDataArrays().size() != sy->getDataArrays().size()) { structure(at + ": number of arrays"); }
        else
        {
          for (std::size_t d = 0; d < sx->getDataArrays().size(); ++d)
          {
            if (sx->getDataArrays()[d]->description != sy->getDataArrays()[d]->description)
            { structure(at + ": array name '" + sx->getDataArrays()[d]->description + "' vs '" + sy->getDataArrays()[d]->description + "'"); }
          }
        }
        ++spectra;
      }
    }
    std::cout << "compared " << a.maps.size() << " maps, " << spectra << " spectra\n";
    std::uint64_t differing = 0, beyond = 0;
    for (const auto& [name, t] : tallies)
    {
      std::cout << "  " << std::left << std::setw(20) << name << std::right << " " << t.compared << " values, " << t.differing
                << " differ" << (t.differing ? ", largest by " : "");
      if (t.differing)
      {
        std::cout << std::scientific << std::setprecision(3) << t.largest << " (relative " << t.largest_relative << ")" << std::fixed;
      }
      std::cout << "\n";
      differing += t.differing;
      if (t.differing && (exact(name) || t.largest_relative > scalar_relative_tolerance)) { ++beyond; }
    }
    std::cout << "  structural differences: " << structural << "\n";
    for (const RunData* r : {&a, &b})
    {
      std::error_code ec;
      if (!r->scratch.empty()) { std::filesystem::remove_all(r->scratch, ec); }
    }
    const bool same = differing == 0 && structural == 0;
    const bool equivalent = beyond == 0 && structural == 0;
    std::cout << (same ? "IDENTICAL" : equivalent ? "EQUIVALENT (scalars within 1e-12 relative)" : "DIFFERENT") << "\n";
    return equivalent || report ? 0 : 1;
  }
  catch (const std::exception& e)
  {
    std::cerr << "identify_mzpeak_run: " << e.what() << "\n";
    return 1;
  }
}
