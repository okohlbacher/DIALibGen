// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <odia/PeptDeepEncoder.h>

#include <cstdint>
#include <filesystem>
#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace OpenMS { class AASequence; }
namespace Ort { struct Session; }

namespace ODIA
{

  /// How the models use a GPU. Only consulted when a GPU is requested.
  struct GpuOptions
  {
    /// Apple: run the models through ONNX Runtime's CoreML execution
    /// provider (macOS; -machine:device coreml or mps) instead of CUDA.
    /// Opt-in: CoreML has no LSTM operator, so RT and CCS stay on the CPU
    /// with a warning, and the MS2 library agrees with the CPU's to float
    /// rounding rather than exactly. See usage.md#apple-gpu.
    bool coreml = false;
    /// CoreML compute units: "CPUAndGPU" (default: measured fastest),
    /// "ALL" (Core ML's choice, may use the Neural Engine), "CPUAndNeuralEngine"
    /// or "CPUOnly".
    std::string coreml_units = "CPUAndGPU";
    /// CoreML: residue counts whose batch shapes are compiled (or loaded
    /// from the cache) in the background as soon as the predictor exists,
    /// several at a time, so the builds overlap the phases before the first
    /// prediction instead of running one after the other inside it (a cold
    /// cache cost a whole-proteome MS2 phase 66 s of serialised compiles).
    /// generate passes the peptide-length window; 0: nothing ahead of use.
    int prebuild_min_residues = 0;
    int prebuild_max_residues = 0;
    /// Throw instead of running on the CPU when CUDA was requested but no
    /// CUDA session could be created (missing cuDNN, no visible GPU, a
    /// CPU-only ONNX Runtime).
    bool require = false;
    /// TensorFloat-32 math in the CUDA provider. Off by default: TF32 rounds
    /// matmul and convolution inputs to 10 mantissa bits, so a TF32 library
    /// differs from the CPU library (about 1 % of fragment sets, CCS by up
    /// to ~1.9 A^2), while with it off the GPU agrees with the CPU to float
    /// rounding. On is faster (one H100, whole proteome: ~25 s less).
    bool tf32 = false;
    /// CUDA sessions per device. Each has its own stream, so the host-side
    /// encoding of one batch overlaps the GPU work of another; measured on an
    /// H100, 2-4 give RT and CCS +30-45 % and 8 is slower again.
    int sessions_per_device = 4;
    /// CUDA devices to spread the sessions over, round robin. Empty: the
    /// first visible device that constructs a session.
    std::vector<int> devices;
    /// Where callers collect "requested CUDA, got the CPU" warnings for the
    /// output's provenance; may be null. The predictor itself only reports
    /// them through fallbackReason().
    std::vector<std::string>* warnings = nullptr;
  };

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
    enum class Provider { CPU, CUDA, CoreML };

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
    /// 4585/s at 64. CPU sessions share one copy of the weights, and with
    /// CPU_CHUNK_ROWS batches each arena is small, so generation runs one
    /// session per thread, up to the physical cores.
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
                               int gpu_device = -1,
                               const GpuOptions& gpu = {});
    ~PeptDeepPredictor();

    /// Why a GPU was requested but this predictor runs on the CPU (or skips a
    /// listed GPU), or empty.
    ///
    /// Set only when a GPU exists and could not be used (for example the
    /// CUDA provider found no cuDNN, or CoreML keeps the model's LSTM on the
    /// CPU): a machine without a GPU, or a CPU-only ONNX Runtime, is the
    /// expected CPU case and leaves this empty. With GpuOptions::require the
    /// constructor throws instead.
    const std::string& fallbackReason() const;
    /// CoreML: what ONNX Runtime kept on the CPU although the model runs on
    /// CoreML (cheap ops such as OneHot), for the provider line; otherwise
    /// empty.
    const std::string& providerNote() const;
    /// The CUDA devices sessions run on, one entry per session; empty on CPU.
    const std::vector<int>& sessionDevices() const;


    /// Which CUDA device the session ended up on, or -1 on CPU and CoreML.
    int device() const;

    /// How many sessions batches are actually spread across (at least 1).
    /// On CoreML: how many batches run concurrently; the sessions there are
    /// one per batch shape, shared by the workers (coremlSession_).
    std::size_t sessionCount() const;

    /// CoreML: waits for the background shape builds (GpuOptions::
    /// prebuild_*) and returns how many shapes this predictor compiled or
    /// loaded and the seconds they took (summed over the builders). 0 on
    /// the other providers. Prediction calls never need this: a worker that
    /// needs a shape still being built waits for that one alone.
    std::pair<std::size_t, double> waitForShapes();

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

    /// On the CPU, rows per batch within a MAX_BATCH_ROWS chunk (see
    /// forEachChunk_). Small batches let every session stay busy to the end of
    /// a call and keep each session's ONNX arena small (about 53 MB for MS2 at
    /// 64 rows against 1.9 GB at 2048), which is what makes one session per
    /// core affordable. The result does not change (forEachChunk_).
    static constexpr std::size_t CPU_CHUNK_ROWS = 64;

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
    /// MAX_BATCH_ROWS (on the CPU cut further into CPU_CHUNK_ROWS pieces) and
    /// runs @p body over them across sessionCount() threads, costliest first.
    /// The chunking does not depend on the session count, so neither do the
    /// results.
    void forEachChunk_(const std::vector<OpenMS::AASequence>& peptides,
                       const std::function<void(std::size_t session,
                                                const std::vector<std::size_t>& group)>& body);

    /// Hands the float weights of the optimised graph at @p optimized to
    /// every later CPU session as one shared, read-only copy, with one
    /// pre-packed-weights container. False if they could not be read, in
    /// which case nothing is shared.
    bool shareWeights_(const std::filesystem::path& optimized);

    /// The session a batch runs on: session @p worker on the CPU and CUDA;
    /// on CoreML the static-shape session for @p rows x @p length, built on
    /// first use (any worker needing it meanwhile waits for it).
    Ort::Session& sessionFor_(std::size_t worker, std::int64_t rows, std::int64_t length);
    /// CoreML: rows padded up to the batch bucket @p rows belongs to, so that
    /// every batch is one of a few compiled shapes: MAX_BATCH_ROWS alone up
    /// to 21 residues, where a block has many full chunks per length and the
    /// one partial chunk pads 5-11 %; 512, 1024 and MAX_BATCH_ROWS from 22
    /// residues on, where the rarer lengths padded 15-38 % (measured on a
    /// 5.8 M-precursor library: 11.5 % padded rows over the whole phase,
    /// 8.8 % with the finer buckets, at the cost of two more compiled shapes
    /// per long length). Returns @p batch itself when nothing is padded.
    const PeptDeepEncoder::Batch& padded_(const PeptDeepEncoder::Batch& batch, PeptDeepEncoder::Batch& storage) const;

    struct Impl;
    std::unique_ptr<Impl> impl_;
  };

} // namespace ODIA
