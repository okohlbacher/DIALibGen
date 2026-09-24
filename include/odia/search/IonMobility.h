// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

/// Ion mobility (diaPASEF) in the built-in search: the pieces that are not a
/// stock OpenSWATH call. See docs/design/built-in-identification.md, M2.
///
///   * mobilityApex: where a precursor's fragments co-locate in 1/K0, measured
///     in the spectra of its isolation window around its elution apex over the
///     WHOLE 1/K0 range of the window. An estimate that does not start from
///     the library's 1/K0 cannot be pulled towards it;
///   * mobilityAt: mobilityApex at one run time in every window holding a
///     precursor's m/z. The ion-mobility calibration measures its seeds with
///     it, and the report the 1/K0 of every reported peak group;
///   * libraryMobility: a precursor's library 1/K0 (the IM column, else from CCS);
///   * automaticImWindow: search:im_window 0, the RT window's rule on the
///     calibration's 1/K0 residuals, clamped;
///   * reportedMobility: OpenSWATH's own 1/K0 of a peak group (its MS2 value
///     cross-checked by its MS1 value), read INSIDE the extraction window and
///     so pulled towards its centre: a diagnostic, not what the report says;
///   * windowOf: OpenSWATH's own rule for which diaPASEF window extracts a
///     precursor (strictly inside in m/z and 1/K0, then the most central in
///     1/K0), so the calibration can say how many seeds the library's 1/K0
///     would have put in the wrong window;
///   * assignWindow: windowOf, or for a precursor whose calibrated 1/K0 lies
///     just outside every window at its m/z, the nearest window its 1/K0
///     extraction range reaches into.
///
/// Everything here reads a peak or a precursor, never a label: a target and
/// its decoy share precursor m/z and 1/K0, so they get the same window and
/// the same 1/K0 extraction range.
#pragma once

#include <odia/Library.h>

#include <OpenMS/OPENSWATHALGO/DATAACCESS/DataStructures.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/SwathMap.h>

#include <cstddef>
#include <limits>
#include <vector>

namespace ODIA::search
{
  struct MobilityApex
  {
    double im = std::numeric_limits<double>::quiet_NaN();   ///< 1/K0 of the apex; NaN = none
    std::size_t fragments = 0;   ///< fragments with signal within mobility_support of the apex
    double intensity = 0.0;      ///< summed intensity within mobility_support of the apex
    std::size_t windows = 0;     ///< mobilityAt: isolation windows read (0 = none holds the m/z)
  };

  /// Bin width, smoothing and support of mobilityApex, 1/K0 units.
  constexpr double mobility_bin = 0.002;
  constexpr double mobility_smoothing = 0.004;   ///< Gaussian SD of the per-fragment mobilograms
  constexpr double mobility_support = 0.01;      ///< half-width around the apex that counts as "at the apex"
  /// Spectra per window mobilityAt reads: the one closest to the run time and one on each side.
  constexpr int mobility_probe_spectra = 3;

  /// Every peak of @p spectra within +-@p ppm_half of a fragment m/z in
  /// @p fragment_mz, with 1/K0 in [@p im_low, @p im_high], adds its intensity
  /// to that fragment's mobilogram (bins of mobility_bin, smoothed with a
  /// Gaussian of mobility_smoothing). Each fragment's mobilogram is scaled to
  /// a maximum of 1 -- every fragment gets one vote, so one intense
  /// interference cannot outvote the others -- and the votes are summed; the
  /// apex is the highest bin (the lowest 1/K0 on a tie), refined to the
  /// intensity-weighted mean 1/K0 of the matched peaks within
  /// mobility_support of it. Spectra without a 1/K0 array are skipped.
  MobilityApex mobilityApex(const std::vector<OpenSwath::SpectrumPtr>& spectra, const std::vector<double>& fragment_mz,
                            double ppm_half, double im_low, double im_high);

  /// mobilityApex over the mobility_probe_spectra spectra closest to @p rt_s
  /// of EVERY MS2 map with lower < @p mz < upper, over the union of those
  /// maps' 1/K0 ranges: where the precursor's fragments co-locate at that
  /// time, whatever its library or calibrated 1/K0 says and whatever the
  /// extraction window. windows = the maps read; NaN when none holds @p mz.
  /// Reads the maps through their own accessors: give each thread its own
  /// (lightClone) when calling it concurrently.
  MobilityApex mobilityAt(const std::vector<OpenSwath::SwathMap>& maps, double mz, double rt_s,
                          const std::vector<double>& fragment_mz, double ppm_half);

  /// Precursor @p i's library 1/K0: the IM column when finite and positive,
  /// else converted from a finite CCS; NaN when the library has neither.
  double libraryMobility(const Library& library, std::size_t i);

  /// search:im_window 0 from the absolute residuals of the ion-mobility
  /// calibration: 2 x im_window_padding x max(their 0.99 quantile,
  /// im_window_normal_quantile x 1.4826 x their median), clamped to
  /// [im_window_min, im_window_max]. Also returns the quantile and the robust SD.
  double automaticImWindow(const std::vector<double>& abs_residuals, double& q99, double& robust_sd);

  /// OpenSWATH's 1/K0 of a peak group: @p ms2 (im_drift) when it is a
  /// mobility, unless @p ms1 (im_ms1_drift) is one too and differs by more
  /// than SearchParams::im_ms1_agreement. NaN otherwise. Both are read inside
  /// the 1/K0 extraction window around the CALIBRATED library value, so this
  /// is pulled towards it (slope about 0.75 against independent measurements
  /// at a 0.11 window, 0.55 at 0.05); the report re-measures instead
  /// (mobilityAt) and keeps this as a diagnostic.
  double reportedMobility(double ms2, double ms1);

  /// OpenSWATH's diaPASEF window for a precursor at (@p mz, @p im): the MS2
  /// map with lower < mz < upper and imLower < im < imUpper whose 1/K0 centre
  /// is closest to @p im (the first on a tie, as OpenSwathWorkflow); -1 when
  /// none holds it.
  int windowOf(const std::vector<OpenSwath::SwathMap>& maps, double mz, double im);

  /// The window a precursor at (@p mz, @p im) is extracted from with a 1/K0
  /// extraction range of +-@p half_width, and the 1/K0 its TRANSITIONS must
  /// carry for OpenSWATH to assign it there (@p assign_im; the compound keeps
  /// @p im, which centres the extraction range). windowOf when a window holds
  /// @p im (assign_im = im). Otherwise the MS2 map with lower < mz < upper
  /// whose 1/K0 range lies closest to @p im (the first on a tie), when that
  /// distance is below @p half_width -- the range reaches into it -- with
  /// assign_im just inside its limits. -1 when neither: nothing of the range
  /// was acquired at that m/z. A precursor only: a target and its decoy share
  /// m/z and 1/K0, and so the answer.
  int assignWindow(const std::vector<OpenSwath::SwathMap>& maps, double mz, double im, double half_width, double& assign_im);
}
