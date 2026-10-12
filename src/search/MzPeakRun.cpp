// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

// mzPeak runs for -run, through the vendored mzpeak-openms library. This file
// includes the library's headers and no OpenMS header (see MzPeakRun.h).

#include <odia/search/MzPeakRun.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <iomanip>
#include <locale>
#include <sstream>
#include <filesystem>
#include <optional>
#include <stdexcept>
#include <system_error>

#ifdef DIALIBGEN_WITH_MZPEAK
#include <mzpeak/open.h>
#include <mzpeak/spectra.h>
#include <mzpeak/spectrum.h>
#include <mzpeak/spectrum_metadata.h>
#endif

namespace fs = std::filesystem;

namespace ODIA::search::mzpeak
{
  bool available()
  {
#ifdef DIALIBGEN_WITH_MZPEAK
    return true;
#else
    return false;
#endif
  }

  bool isMzPeak(const std::string& path)
  {
    std::error_code ec;
    if (fs::is_directory(path, ec)) { return fs::is_regular_file(fs::path(path) / "mzpeak_index.json", ec); }
    std::string ext = fs::path(path).extension().string();
    std::transform(ext.begin(), ext.end(), ext.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return ext == ".mzpeak";
  }

#ifdef DIALIBGEN_WITH_MZPEAK
  namespace
  {
    /// A float32 field as the decimal it was written from: the shortest text
    /// that reads back as the same float, read as a double. This is the value
    /// an mzML of the archive carries (its writers print float32 fields that
    /// way, and OpenMS parses the text as a double), and the value the
    /// instrument reported (a float32 isolation target 1276.05078125 is the
    /// 1276.0508 the method states).
    /// Text to double in the "C" locale; false unless all of @p s is a number.
    bool parseDouble(const std::string& s, double& out)
    {
      std::istringstream in(s);
      in.imbue(std::locale::classic());
      in >> out;
      return !in.fail() && in.peek() == std::char_traits<char>::eof();
    }

    double decimalOf(float v)
    {
      if (!std::isfinite(v)) { return static_cast<double>(v); }
      for (int digits = 1; digits <= 9; ++digits)
      {
        std::ostringstream text;
        text.imbue(std::locale::classic());
        text << std::setprecision(digits) << v;
        double back = 0.0;
        if (parseDouble(text.str(), back) && static_cast<float>(back) == v) { return back; }
      }
      return static_cast<double>(v);
    }

    /// The window's 1/K0 band of a selected ion: the typed columns, or the
    /// parameters mzpeak-convert writes (MZP:1000006 / MZP:1000007, "isolation
    /// window inverse reduced ion mobility lower/upper limit").
    void imBand(const MzPeak::SelectedIonInfo& ion, Window& w)
    {
      std::optional<double> lo = ion.ion_mobility_lower_limit, hi = ion.ion_mobility_upper_limit;
      for (const auto& p : ion.parameters)
      {
        if (!p.value) { continue; }
        const std::string acc = p.accession.value_or(""), name = p.name.value_or("");
        double v = 0.0;
        if (!parseDouble(*p.value, v)) { continue; }
        if (!lo && (acc == "MZP:1000006" || name == "isolation window inverse reduced ion mobility lower limit")) { lo = v; }
        if (!hi && (acc == "MZP:1000007" || name == "isolation window inverse reduced ion mobility upper limit")) { hi = v; }
      }
      if (lo && hi)
      {
        w.has_im_limits = true;
        w.im_lower = *lo;
        w.im_upper = *hi;
      }
      if (ion.ion_mobility_value)
      {
        w.has_im = true;
        w.im = *ion.ion_mobility_value;
      }
    }

    SpectrumInfo infoOf(const MzPeak::Spectrum& sp)
    {
      const MzPeak::SpectrumMetadata& md = sp.metadata();
      SpectrumInfo info;
      info.index = md.index;
      info.id = md.id;
      info.ms_level = md.ms_level.value_or(0);
      info.has_rt = md.retention_time.has_value();
      info.rt_s = md.retention_time.value_or(0.0);
      info.points = std::max(md.number_of_peaks.value_or(0), md.number_of_data_points.value_or(0));
      // The selected ion of each precursor. An archive whose precursor rows
      // carry no precursor_index (mzpeak-convert 0.18's ims-compact frames)
      // gets all of a frame's selected ions joined to its first precursor by
      // the library; there, the k-th ion belongs to the k-th window.
      std::vector<const MzPeak::SelectedIonInfo*> ions(md.precursors.size(), nullptr);
      std::size_t holders = 0;
      for (const auto& p : md.precursors) { holders += p.selected_ions.empty() ? 0 : 1; }
      if (md.precursors.size() > 1 && holders == 1 && md.precursors.front().selected_ions.size() == md.precursors.size())
      {
        for (std::size_t k = 0; k < ions.size(); ++k) { ions[k] = &md.precursors.front().selected_ions[k]; }
      }
      else
      {
        for (std::size_t k = 0; k < ions.size(); ++k)
        {
          if (!md.precursors[k].selected_ions.empty()) { ions[k] = &md.precursors[k].selected_ions.front(); }
        }
      }
      for (std::size_t k = 0; k < md.precursors.size(); ++k)
      {
        const auto& p = md.precursors[k];
        Window w;
        const auto& iw = p.isolation_window;
        if (iw.target_mz) { w.target = decimalOf(*iw.target_mz); }
        w.mz = w.target;
        if (iw.lower_offset && iw.upper_offset)
        {
          w.has_offsets = true;
          w.lower_offset = decimalOf(*iw.lower_offset);
          w.upper_offset = decimalOf(*iw.upper_offset);
        }
        if (ions[k] != nullptr)
        {
          const auto& ion = *ions[k];
          // OpenMS takes the selected ion's m/z as the precursor's when it has
          // one (PeakFileOptions::getPrecursorMZSelectedIon, the default).
          if (ion.selected_ion_mz) { w.mz = *ion.selected_ion_mz; }
          imBand(ion, w);
        }
        info.windows.push_back(w);
      }
      // A per-window spectrum written from an mzML keeps the window's band as
      // the spectrum's own "ion mobility lower/upper limit" (userParams there,
      // which OpenMS reads as meta values of those names).
      if (info.windows.size() == 1 && !info.windows.front().has_im_limits)
      {
        std::optional<double> lo, hi;
        for (const auto& p : md.parameters)
        {
          double v = 0.0;
          if (!p.value || !p.name || !parseDouble(*p.value, v)) { continue; }
          if (*p.name == "ion mobility lower limit") { lo = v; }
          if (*p.name == "ion mobility upper limit") { hi = v; }
        }
        if (lo && hi)
        {
          info.windows.front().has_im_limits = true;
          info.windows.front().im_lower = *lo;
          info.windows.front().im_upper = *hi;
        }
      }
      return info;
    }

    class SpectraCursor final : public Reader::Cursor
    {
    public:
      explicit SpectraCursor(const MzPeak::Index& index) : spectra_(index.spectra(MzPeak::MetadataDetail::Full)) {}
      void read(std::uint64_t index, Peaks& out) override
      {
        const MzPeak::Spectrum sp = spectra_[index];   // bound first: mz() and intensity() live as long as sp
        out.mz = sp.mz();
        out.intensity = sp.intensity();
        out.mobility = sp.ion_mobility_array();
      }

    private:
      MzPeak::Spectra spectra_;
    };

    std::string what(const std::exception& e) { return e.what(); }
  }

  struct Reader::Impl
  {
    std::unique_ptr<MzPeak::Index> index;
    std::vector<SpectrumInfo> spectra;
  };

  Reader::Reader(const std::string& path) : impl_(std::make_unique<Impl>())
  {
    try
    {
      impl_->index = std::make_unique<MzPeak::Index>(MzPeak::open(fs::path(path)));
      // Full: the 1/K0 bands of diaPASEF windows are parameters of the
      // selected ion, which the leaner levels leave out. The map is cached
      // by the archive and shared with every cursor.
      MzPeak::Spectra all = impl_->index->spectra(MzPeak::MetadataDetail::Full);
      const std::size_t n = all.size();
      impl_->spectra.reserve(n);
      for (std::size_t i = 0; i < n; ++i)
      {
        const MzPeak::Spectrum sp = all[i];
        impl_->spectra.push_back(infoOf(sp));
        if (impl_->spectra.back().index != i)
        {
          throw std::runtime_error("spectrum " + std::to_string(i) + " carries index " +
                                   std::to_string(impl_->spectra.back().index) + "; the spectra are not numbered 0..n-1");
        }
      }
    }
    catch (const std::exception& e)
    {
      throw std::runtime_error("cannot read the mzPeak archive " + path + ": " + what(e));
    }
  }

  Reader::~Reader() = default;

  const std::vector<SpectrumInfo>& Reader::spectra() const { return impl_->spectra; }

  std::string Reader::version() const { return impl_->index->version(); }

  std::unique_ptr<Reader::Cursor> Reader::cursor() const { return std::make_unique<SpectraCursor>(*impl_->index); }

#else

  struct Reader::Impl {};

  Reader::Reader(const std::string& path)
  {
    throw std::runtime_error("cannot read " + path + ": this build of DIALibGen has no mzPeak support (DIALIBGEN_WITH_MZPEAK=OFF); "
                             "convert the run to mzML (mzpeak-convert --to mzml)");
  }
  Reader::~Reader() = default;
  const std::vector<SpectrumInfo>& Reader::spectra() const { static const std::vector<SpectrumInfo> none; return none; }
  std::string Reader::version() const { return {}; }
  std::unique_ptr<Reader::Cursor> Reader::cursor() const { return nullptr; }

#endif
}
