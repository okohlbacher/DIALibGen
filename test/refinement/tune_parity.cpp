// Copyright (c) 2026, Oliver Kohlbacher and the DIALibRefine authors.
// SPDX-License-Identifier: BSD-3-Clause

/// Two things the whole trainer rests on, tested against the stock ONNX:
///  1. PARITY: the libtorch transcription loaded from the ONNX predicts what
///     ONNX Runtime predicts (via DIALibGen's PeptDeepPredictor), for both
///     heads, on real peptides of several lengths with a modification.
///  2. ROUND TRIP: load -> store writes back a byte-identical file, and a
///     perturbed model written back and re-loaded returns the perturbation
///     (so the gate permutation is its own inverse and every tensor is
///     reached).
///  3. METADATA: provenance survives write/read, unrelated and unknown fields
///     remain byte-identical, and adding metadata does not move weight offsets.
///
///   tune_parity <rt.onnx> <ccs.onnx>

#include <odia/tune/OnnxWeights.h>
#include <odia/AtomicFile.h>
#include <cstring>
#include <limits>
#include <odia/PeptDeepEncoder.h>
#include <odia/PeptDeepPredictor.h>

#include <OpenMS/CHEMISTRY/AASequence.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <random>

using namespace ODIA;
using namespace ODIA::tune;

namespace
{
  int failures = 0;
  void check(bool ok, const std::string& what)
  {
    std::cout << (ok ? "  ok   " : "  FAIL ") << what << "\n";
    if (!ok) { ++failures; }
  }

  std::vector<OpenMS::AASequence> peptides()
  {
    std::vector<OpenMS::AASequence> v;
    for (const char* s : {"PEPTIDEK", "LGEHNIDVLEGNEQFINAAK", "AC(UniMod:4)DEFGHIK", "M(UniMod:35)VLSDGK",
                          "TTPSYVAFTDTER", "QWERTYIPASDFGHK", "GLVLIAFSQYLQQC(UniMod:4)PFDEHVK", "SVAAAR"})
    { v.push_back(OpenMS::AASequence::fromString(s)); }
    return v;
  }

  /// Predict with the libtorch model, one length group at a time (as training does).
  std::vector<float> torchPredict(Head& model, const std::vector<OpenMS::AASequence>& peps, const std::vector<int>& charges)
  {
    torch::NoGradGuard ng;
    model->eval();
    std::vector<float> out(peps.size(), NAN);
    for (const auto& group : PeptDeepEncoder::groupByLength(peps))
    {
      std::vector<OpenMS::AASequence> gp; std::vector<int> gz;
      for (auto i : group) { gp.push_back(peps[i]); gz.push_back(charges[i]); }
      auto b = PeptDeepEncoder::encode(gp, gz, 30.0f, "Lumos");
      const auto rows = static_cast<std::int64_t>(b.rows), L = static_cast<std::int64_t>(b.sequence_length);
      auto aa = torch::from_blob(b.aa_indices.data(), {rows, L}, torch::kInt64).clone();
      auto mx = torch::from_blob(b.mod_x.data(), {rows, L, MOD_FEATURES}, torch::kFloat32).clone();
      auto ch = torch::from_blob(b.charges.data(), {rows, 1}, torch::kFloat32).clone();
      auto y = model->forward(aa, mx, ch).contiguous();
      for (std::size_t k = 0; k < group.size(); ++k) { out[group[k]] = y[static_cast<std::int64_t>(k)].item<float>(); }
    }
    return out;
  }

  double maxAbsDiff(const std::vector<float>& a, const std::vector<float>& b)
  {
    double m = 0;
    for (std::size_t i = 0; i < a.size(); ++i)
    {
      const double d = std::fabs(static_cast<double>(a[i]) - static_cast<double>(b[i]));
      if (!std::isfinite(d)) { return INFINITY; }   // a NaN anywhere is a failure, not a zero
      m = std::max(m, d);
    }
    return m;
  }

  void metadataRoundtrip(const OnnxFile& stock, Head& model)
  {
    struct TempFile
    {
      std::string path;
      std::filesystem::path dir;
      TempFile()
      {
        dir = std::filesystem::temp_directory_path() / ("dialibgen-metadata-" + std::to_string(std::random_device{}()));
        if (!std::filesystem::create_directory(dir)) { throw std::runtime_error("cannot create metadata test directory"); }
        path = (dir / "model.onnx").string();
      }
      ~TempFile() { std::error_code ec; std::filesystem::remove(path, ec); std::filesystem::remove(dir, ec); }
    } temp;

    // Hand-encoded metadata with unknown nested fields, followed by unknown
    // top-level fields of each supported wire type. Put them BEFORE the graph
    // so rewriting metadata would invalidate cached initializer offsets.
    const std::string vendor = "\x72\x16\x0a\x06" "vendor" "\x12\x03" "old" "\x48\x81\x01\x1a\x04" "keep";
    const std::string unknown = "\xc0\x3e\x81\x01\xc9\x3e" "01234567" "\xd5\x3e" "abcd" "\xda\x3e\x04" "keep";
    const std::string prefix = vendor + unknown;
    OnnxFile fixture = stock;
    fixture.bytes.insert(fixture.bytes.begin(), prefix.begin(), prefix.end());
    fixture.write(temp.path);
    auto f = OnnxFile::read(temp.path);
    check(f.metadata().at("vendor") == "old", "existing metadata is read");
    f.write(temp.path);
    check(OnnxFile::read(temp.path).bytes == fixture.bytes, "unmodified metadata and unknown fields round-trip byte-identically");
    f.setMetadata("vendor", "old");
    f.write(temp.path);
    check(OnnxFile::read(temp.path).bytes == fixture.bytes, "setting the existing metadata value is byte-identical");

    const std::string key = "org.openms.dialibgen.training";
    const std::string provenance = "{\"tool\":\"DIALibGen\",\"note\":\"" + std::string(300, 'x') + "\"}";
    auto expected_metadata = f.metadata();
    expected_metadata[key] = provenance;
    f.setMetadata(key, provenance);
    storeWeights(f, model);
    check(f.bytes == fixture.bytes, "pending metadata preserves initializer offsets and untouched weights");
    f.write(temp.path);
    auto added = OnnxFile::read(temp.path);
    check(added.metadata() == expected_metadata, "provenance and existing metadata survive write/read");
    check(added.bytes.size() > fixture.bytes.size() && std::equal(fixture.bytes.begin(), fixture.bytes.end(), added.bytes.begin()),
          "adding provenance preserves every existing byte");

    added.setMetadata("vendor", "new value");
    added.write(temp.path);
    auto replaced = OnnxFile::read(temp.path);
    const std::string changed_vendor = "\x72\x1c\x0a\x06" "vendor" "\x48\x81\x01\x1a\x04" "keep" "\x12\x09" "new value";
    std::vector<std::uint8_t> expected_bytes(changed_vendor.begin(), changed_vendor.end());
    expected_bytes.insert(expected_bytes.end(), added.bytes.begin() + vendor.size(), added.bytes.end());
    check(replaced.metadata().at("vendor") == "new value" && replaced.bytes == expected_bytes,
          "replacing metadata preserves its unknown nested fields and all unrelated bytes");

    replaced.setMetadata(key, "{\"tool\":\"DIALibGen\",\"best_epoch\":2}");
    replaced.write(temp.path);
    auto updated = OnnxFile::read(temp.path);
    const std::string wire(updated.bytes.begin(), updated.bytes.end());
    const auto first = wire.find(key);
    check(updated.metadata().at(key) == "{\"tool\":\"DIALibGen\",\"best_epoch\":2}" &&
          first != std::string::npos && wire.find(key, first + key.size()) == std::string::npos,
          "a later tuning run replaces provenance without duplicating its key");
    const auto before_store = updated.bytes;
    storeWeights(updated, model);
    updated.write(temp.path);
    check(updated.bytes == before_store && OnnxFile::read(temp.path).bytes == before_store,
          "metadata-bearing model still supports byte-identical weight round-trip");
  }

  void rejectsNonfinite(const std::string& path, const std::string& kind)
  {
    OnnxFile broken = OnnxFile::read(path);
    const float nan = std::numeric_limits<float>::quiet_NaN();
    for (const auto& tensor : broken.initializers)
    {
      if (tensor.data_type != 1) { continue; }
      for (std::size_t i = 0; i + sizeof(nan) <= tensor.raw_length; i += sizeof(nan))
      { std::memcpy(broken.bytes.data() + tensor.raw_offset + i, &nan, sizeof(nan)); }
    }
    AtomicFile temporary(std::filesystem::temp_directory_path() / "nonfinite-test.onnx");
    broken.write(temporary.temporaryPath().string());
    const auto peps = peptides(); // Multiple lengths ensure independent workers fail concurrently.
    for (const int sessions : {1, 4})
    {
      bool refused = false;
      try
      {
        PeptDeepPredictor predictor(temporary.temporaryPath().string(), false, 1, sessions);
        if (kind == "rt") { (void)predictor.predictRT(peps); }
        else if (kind == "ccs") { (void)predictor.predictCCS(peps, std::vector<int>(peps.size(), 2)); }
        else { (void)predictor.predictMS2(peps, std::vector<int>(peps.size(), 2), 30, "QE"); }
      }
      catch (const std::exception& e) { refused = std::string(e.what()).find("non-finite") != std::string::npos; }
      check(refused, kind + ": non-finite model output rejected with " + std::to_string(sessions) + " sessions");
    }
  }

  void parity(const std::string& path, bool ccs)
  {
    // Stage markers, flushed: where a crash happens is the whole diagnosis on
    // a platform that cannot be debugged interactively.
    auto stage = [&](const char* what) { std::cout << "  ..   " << (ccs ? "ccs" : "rt") << ": " << what << std::endl; };
    const auto peps = peptides();
    std::vector<int> charges; for (std::size_t i = 0; i < peps.size(); ++i) { charges.push_back(2 + static_cast<int>(i % 3)); }
    stage("reading the ONNX");
    OnnxFile f = OnnxFile::read(path);
    stage("building the libtorch model");
    Head model(std::make_shared<HeadImpl>(ccs));
    stage("loading the weights");
    const std::size_t n = loadWeights(f, model);
    check(n == 21, std::string(ccs ? "ccs" : "rt") + ": 21 initializers loaded (" + std::to_string(n) + ")");
    stage("libtorch forward");
    auto ours = torchPredict(model, peps, charges);
    stage("ONNX Runtime forward");
    PeptDeepPredictor ort(path, /*prefer_gpu=*/false, /*intra_op_threads=*/1);
    auto ref = ccs ? ort.predictCCS(peps, charges) : ort.predictRT(peps);
    const double d = maxAbsDiff(ours, ref);
    // RT is rt_norm (~0..1), CCS is A^2 (~300-600): tolerances scaled accordingly
    const double tol = ccs ? 5e-3 : 1e-5;
    char buf[256]; std::snprintf(buf, sizeof buf, "%s parity vs ONNX Runtime: max |diff| = %.3g (tol %.0e)", ccs ? "ccs" : "rt", d, tol);
    check(d <= tol, buf);
    for (std::size_t i = 0; i < 2; ++i) { std::printf("       %-28s z%d torch %.6f ort %.6f\n", peps[i].toString().c_str(), charges[i], ours[i], ref[i]); }

    // Round trip 1: unchanged model writes back byte-identical
    OnnxFile g = f;
    storeWeights(g, model);
    check(g.bytes == f.bytes, std::string(ccs ? "ccs" : "rt") + ": store(load(x)) is byte-identical");
    metadataRoundtrip(f, model);

    // Round trip 2: perturb every trainable parameter, write back, reload into a fresh model, compare
    {
      torch::NoGradGuard ng;
      for (auto& p : model->parameters()) { p.add_(0.001 * torch::arange(p.numel(), torch::kFloat32).reshape(p.sizes()) / static_cast<double>(std::max<std::int64_t>(1, p.numel()))); }
    }
    OnnxFile h = f;
    storeWeights(h, model);
    Head again(std::make_shared<HeadImpl>(ccs));
    loadWeights(h, again);
    double worst = 0;
    auto a = model->named_parameters(); auto b = again->named_parameters();
    for (const auto& kv : a)
    {
      const double w = (kv.value() - *b.find(kv.key())).abs().max().item<double>();
      worst = std::isfinite(w) ? std::max(worst, w) : INFINITY;
    }
    check(worst == 0.0, std::string(ccs ? "ccs" : "rt") + ": perturbed model survives store->load exactly (worst " + std::to_string(worst) + ")");
    check(h.bytes != f.bytes, std::string(ccs ? "ccs" : "rt") + ": perturbed bytes differ from the original");
  }
}

int main(int argc, char** argv)
{
  if (argc != 3) { std::cerr << "usage: tune_parity <rt.onnx> <ccs.onnx>\n"; return 2; }
  torch::set_num_threads(1);
  // Diagnostics: DLR_NO_MKLDNN=1 routes conv/LSTM away from oneDNN.
  if (const char* e = std::getenv("DLR_NO_MKLDNN"); e && *e && *e != '0') { at::globalContext().setUserEnabledMkldnn(false); std::cout << "  (oneDNN disabled)\n"; }
  try
  {
    parity(argv[1], false);
    parity(argv[2], true);
    rejectsNonfinite(argv[1], "rt");
    rejectsNonfinite(argv[2], "ccs");
    rejectsNonfinite((std::filesystem::path(argv[1]).parent_path() / "peptdeep_ms2_dynamic.onnx").string(), "ms2");
  }
  catch (const std::exception& e) { std::cerr << "exception: " << e.what() << "\n"; return 1; }
  std::cout << (failures ? "FAILED" : "PASSED") << " (" << failures << " failures)\n";
  return failures ? 1 : 0;
}
