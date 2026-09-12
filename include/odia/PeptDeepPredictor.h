// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <odia/PeptDeepEncoder.h>

#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace OpenMS { class AASequence; }

namespace ODIA
{

  /// Runs AlphaPeptDeep's PeptDeep models through ONNX Runtime.
  ///
  /// ODIA drives the session itself rather than using OpenMS's PeptDeep
  /// binding, for two reasons that are not preference: the binding rejects
  /// modified peptides outright, and its headers are compiled into libOpenMS.so
  /// but never installed, so no external project can include them. The model
  /// weights still come from OpenMS -- WITH_ONNX=ON downloads them against
  /// pinned checksums -- so nothing is vendored here.
  class PeptDeepPredictor
  {
  public:
    /// Which execution provider the session ended up using.
    enum class Provider { CPU, CUDA };

    /// @param model_path a PeptDeep .onnx file.
    /// @param prefer_gpu attempt CUDA first, falling back to CPU (D4).
    /// @param intra_op_threads 0 leaves it to ONNX Runtime.
    /// @param sessions how many independent sessions to run batches across.
    ///
    /// @note **Prefer sessions over intra_op_threads on CPU.** These are
    /// recurrent networks whose per-op tensors are too small to spread across
    /// many threads. Measured on 16 pinned cores with the MS2 model, one
    /// session: intra-op 1 gives 169 peptides/s, 8 gives 375 (the peak, only
    /// 2.2x), 16 gives 318 and 32 gives 231 -- past 8 threads it goes
    /// backwards. Sixteen single-threaded sessions on the same cores give
    /// 1760/s, 4.7x the best intra-op configuration, and it keeps scaling to
    /// 4585/s at 64. Each session costs a copy of the weights and an arena,
    /// so the count is capped by the caller rather than defaulted to the core
    /// count.
    /// @param gpu_device which CUDA device. -1 (default) PROBES: it tries each
    ///        visible device in turn and keeps the first that actually
    ///        constructs a session. Probing rather than trusting device 0 is
    ///        not defensive programming -- GPU 0 is dead on both of this
    ///        project's GPU nodes, so a hard-coded 0 falls back to CPU and the
    ///        build silently takes 25 minutes instead of 4.
    explicit PeptDeepPredictor(const std::string& model_path,
                               bool prefer_gpu = true,
                               int intra_op_threads = 0,
                               int sessions = 1,
                               int gpu_device = -1);
    ~PeptDeepPredictor();

    /// Which CUDA device the session ended up on, or -1 on CPU.
    int device() const;

    /// How many sessions batches are actually spread across (at least 1).
    std::size_t sessionCount() const;

    PeptDeepPredictor(const PeptDeepPredictor&) = delete;
    PeptDeepPredictor& operator=(const PeptDeepPredictor&) = delete;

    /// The provider actually in use. Requesting CUDA on a machine without it is
    /// not an error -- the same binary has to work on both -- so callers that
    /// care must ask rather than assume.
    Provider provider() const;

    /// Largest number of peptides submitted to the model in one call.
    ///
    /// Not a tuning knob. Without a cap, a length group is submitted whole:
    /// 2 M peptides with a tryptic length distribution put 234,593 in the
    /// largest group and peaked at 26.7 GiB, and under a memory limit it did
    /// not degrade but aborted inside the BLAS allocator with no indication of
    /// which stage failed. Capping also fixes the batch size, which the result
    /// depends on at the last bit.
    static constexpr std::size_t MAX_BATCH_ROWS = 2048;

    /// Peptides that could not be encoded, by index into the input, with the
    /// reason. One unencodable peptide must not discard a whole run's work.
    struct Failure
    {
      std::size_t index;
      std::string reason;
    };

    /// Predict normalised iRT, one value per peptide, in the input order.
    ///
    /// Peptides are grouped by length internally: padding is not inert, because
    /// index 0 is one-hot encoded and no model applies a padding mask, so a
    /// mixed-length batch would give every peptide a different answer from the
    /// one it gets alone.
    /// @param failures if given, receives the peptides that could not be
    ///        encoded; their entries in the result are left NaN. If not given,
    ///        an unencodable peptide throws.
    std::vector<float> predictRT(const std::vector<OpenMS::AASequence>& peptides,
                                 std::vector<Failure>* failures = nullptr);

    /// One predicted fragment spectrum: (nAA - 1) positions x 8 channels.
    ///
    /// Channel order comes from alphabase's sort_charged_frag_types, which is
    /// sorted(no_loss) + sorted(loss):
    ///   0 b_z1  1 b_z2  2 y_z1  3 y_z2
    ///   4 b_modloss_z1  5 b_modloss_z2  6 y_modloss_z1  7 y_modloss_z2
    struct Spectrum
    {
      std::size_t positions = 0;
      static constexpr std::size_t CHANNELS = 8;
      std::vector<float> intensities;   ///< [positions, CHANNELS], row-major

      float at(std::size_t position, std::size_t channel) const
      {
        return intensities[position * CHANNELS + channel];
      }
    };

    /// Predict fragment intensities. @p charges is one per peptide.
    ///
    /// @param failures if given, receives the peptides that could not be
    ///        encoded -- an unsupported modification, or a charge below 1 --
    ///        and their entries in the result are left with positions == 0 and
    ///        no intensities. If not given, any of those throws.
    ///
    /// A failed entry is an empty spectrum, so Spectrum::at() must not be
    /// called on one; check positions first. Unlike predictRT there is no NaN
    /// to carry the distinction, because there is no value to carry it in.
    std::vector<Spectrum> predictMS2(const std::vector<OpenMS::AASequence>& peptides,
                                     const std::vector<int>& charges,
                                     float nce = 30.0f,
                                     const std::string& instrument = "QE",
                                     std::vector<Failure>* failures = nullptr);

    /// Predict collision cross-section in square angstroms, one per peptide.
    ///
    /// This is CCS, not ion mobility. Converting it to the 1/K0 a timsTOF
    /// reports needs the Mason-Schamp relation with the drift gas mass and a
    /// calibration constant, and those live in alphabase, which is not on this
    /// machine -- inventing them would produce plausible numbers in the wrong
    /// units. OpenMS's own PeptDeepCCSInference likewise stops at CCS.
    ///
    /// @param failures as for predictMS2; unfilled entries are left NaN.
    std::vector<float> predictCCS(const std::vector<OpenMS::AASequence>& peptides,
                                  const std::vector<int>& charges,
                                  std::vector<Failure>* failures = nullptr);

  private:
    /// @param session which session to submit on. Sessions are independent, so
    /// concurrent calls are safe as long as no two use the same index and no
    /// two write overlapping @p group indices -- both of which hold because
    /// the chunks partition the input.
    void runBatch_(std::size_t session, const PeptDeepEncoder::Batch& batch,
                   const std::vector<std::size_t>& group, std::vector<float>& out);
    void runMS2Batch_(std::size_t session, const PeptDeepEncoder::Batch& batch,
                      const std::vector<std::size_t>& group, std::vector<Spectrum>& out);
    void runCCSBatch_(std::size_t session, const PeptDeepEncoder::Batch& batch,
                      const std::vector<std::size_t>& group, std::vector<float>& out);

    /// Splits @p count peptides into length-homogeneous chunks of at most
    /// MAX_BATCH_ROWS and runs @p body over them across sessionCount()
    /// threads. The chunking is identical to the serial order, so results do
    /// not depend on how many sessions ran them.
    void forEachChunk_(const std::vector<OpenMS::AASequence>& peptides,
                       const std::function<void(std::size_t session,
                                                const std::vector<std::size_t>& group)>& body);

    struct Impl;
    std::unique_ptr<Impl> impl_;
  };

} // namespace ODIA
