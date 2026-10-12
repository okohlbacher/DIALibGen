// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

/// Reading a DIA run from an mzPeak archive for -run (experimental).
///
/// The archive is read with the vendored mzpeak-openms library
/// (third_party/mzpeak-openms), which this header does not expose: the
/// library's public headers need a newer Boost than OpenMS's headers on some
/// platforms, so no translation unit includes both. RunStages.cpp turns the
/// spectra handed out here into OpenMS spectra and feeds them into the same
/// SwathFile consumers an mzML run goes through.
///
/// What is read: per spectrum its MS level, retention time (seconds), native
/// id and isolation windows (precursor m/z, offsets, and the window's 1/K0
/// band where the archive has one); per peak m/z, intensity and, on an
/// ion-mobility run, 1/K0.

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace ODIA::search::mzpeak
{
  /// Whether this build reads mzPeak (DIALIBGEN_WITH_MZPEAK).
  bool available();

  /// Whether @p path names an mzPeak run: a file ending in .mzpeak (any
  /// case), or a directory holding mzpeak_index.json (an unpacked archive).
  bool isMzPeak(const std::string& path);

  /// One isolation window of a spectrum.
  struct Window
  {
    /// The precursor m/z OpenMS reads from an mzML of the same archive: the
    /// selected ion's m/z where the archive has one, else the isolation
    /// window's target.
    double mz = 0.0;
    double target = 0.0;   ///< the isolation window target m/z (float32 in the archive)
    double lower_offset = 0.0;
    double upper_offset = 0.0;
    bool has_offsets = false;
    /// The window's 1/K0 band (diaPASEF), when the archive records it.
    bool has_im_limits = false;
    double im_lower = -1.0;
    double im_upper = -1.0;
    /// The selected ion's 1/K0 (the band's midpoint on diaPASEF), when present.
    bool has_im = false;
    double im = 0.0;
  };

  /// What the run's spectra are, without their peaks.
  struct SpectrumInfo
  {
    std::uint64_t index = 0;
    std::string id;
    int ms_level = 0;
    double rt_s = 0.0;
    bool has_rt = false;
    std::vector<Window> windows;   ///< empty for MS1
    /// Points the metadata announces (peaks or data points, whichever the
    /// spectrum is stored as); 0 when it says neither.
    std::uint64_t points = 0;
  };

  /// A spectrum's peaks; mobility is empty unless the archive stores 1/K0 per peak.
  struct Peaks
  {
    std::vector<double> mz;
    std::vector<float> intensity;
    std::vector<double> mobility;
  };

  class Reader
  {
  public:
    /// Opens the archive and reads every spectrum's metadata; throws
    /// std::runtime_error with the library's reason when it cannot.
    explicit Reader(const std::string& path);
    ~Reader();
    Reader(const Reader&) = delete;
    Reader& operator=(const Reader&) = delete;

    const std::vector<SpectrumInfo>& spectra() const;
    /// The archive's format version (mzpeak_index.json).
    std::string version() const;

    /// A reader of peaks for ONE thread: the library's spectra are not
    /// thread-safe, so every worker reads through its own cursor.
    class Cursor
    {
    public:
      virtual ~Cursor() = default;
      /// The peaks of spectrum @p index (an index of spectra()), into @p out.
      virtual void read(std::uint64_t index, Peaks& out) = 0;
    };
    std::unique_ptr<Cursor> cursor() const;

  private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
  };
}
