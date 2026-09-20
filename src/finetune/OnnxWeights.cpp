// Copyright (c) 2026, Oliver Kohlbacher and the DIALibRefine authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <odia/tune/OnnxWeights.h>

#include <cstring>
#include <fstream>
#include <functional>
#include <stdexcept>

namespace ODIA::tune
{
  namespace
  {
    struct Walker
    {
      const std::vector<std::uint8_t>& b;
      std::uint64_t varint(std::size_t& pos) const
      {
        std::uint64_t r = 0; int sh = 0;
        for (;;)
        {
          if (pos >= b.size()) { throw std::runtime_error("truncated ONNX (varint)"); }
          if (sh > 63) { throw std::runtime_error("malformed ONNX (varint longer than 10 bytes)"); }
          const std::uint8_t x = b[pos++];
          r |= static_cast<std::uint64_t>(x & 0x7f) << sh; sh += 7;
          if (!(x & 0x80)) { return r; }
        }
      }
      /// Also expose each complete wire field's range so unknown fields can be copied verbatim.
      void walkRaw(std::size_t pos, std::size_t end,
                   const std::function<void(int, int, std::uint64_t, std::size_t, std::size_t, std::size_t, std::size_t)>& cb) const
      {
        while (pos < end)
        {
          const std::size_t start = pos;
          const std::uint64_t key = varint(pos);
          if (pos > end) { throw std::runtime_error("malformed ONNX (field key runs past its message)"); }
          const int field = static_cast<int>(key >> 3), wt = static_cast<int>(key & 7);
          std::uint64_t value = 0;
          std::size_t off = 0, len = 0;
          if (wt == 0) { value = varint(pos); }
          else if (wt == 1) { if (end - pos < 8) { throw std::runtime_error("truncated ONNX (fixed64)"); } off = pos; len = 8; pos += len; }
          else if (wt == 5) { if (end - pos < 4) { throw std::runtime_error("truncated ONNX (fixed32)"); } off = pos; len = 4; pos += len; }
          else if (wt == 2)
          {
            const std::uint64_t ln = varint(pos);
            if (pos > end || ln > end - pos) { throw std::runtime_error("truncated ONNX (length-delimited field runs past its message)"); }
            off = pos; len = static_cast<std::size_t>(ln); pos += len;
          }
          else { throw std::runtime_error("ONNX: unsupported wire type"); }
          if (pos > end) { throw std::runtime_error("truncated ONNX (field runs past its message)"); }
          cb(field, wt, value, off, len, start, pos);
        }
      }
      void walk(std::size_t pos, std::size_t end,
                const std::function<void(int, int, std::uint64_t, std::size_t, std::size_t)>& cb) const
      {
        walkRaw(pos, end, [&](int field, int wt, std::uint64_t value, std::size_t off, std::size_t len, std::size_t, std::size_t)
        { cb(field, wt, value, off, len); });
      }
      std::string str(std::size_t off, std::size_t len) const { return std::string(reinterpret_cast<const char*>(b.data() + off), len); }
    };

    std::pair<std::string, std::string> metadataEntry(const Walker& w, std::size_t off, std::size_t len)
    {
      std::pair<std::string, std::string> entry;
      w.walk(off, off + len, [&](int field, int wt, std::uint64_t, std::size_t soff, std::size_t slen)
      {
        if (field == 1 && wt == 2) { entry.first = w.str(soff, slen); }
        else if (field == 2 && wt == 2) { entry.second = w.str(soff, slen); }
      });
      return entry;
    }

    void appendVarint(std::string& bytes, std::uint64_t value)
    {
      while (value >= 0x80) { bytes.push_back(static_cast<char>((value & 0x7f) | 0x80)); value >>= 7; }
      bytes.push_back(static_cast<char>(value));
    }

    void appendString(std::string& bytes, int field, const std::string& value)
    {
      appendVarint(bytes, (static_cast<std::uint64_t>(field) << 3) | 2);
      appendVarint(bytes, value.size());
      bytes += value;
    }
  }

  OnnxFile OnnxFile::read(const std::string& path)
  {
    OnnxFile f;
    std::ifstream in(path, std::ios::binary);
    if (!in) { throw std::runtime_error("cannot open ONNX file: " + path); }
    f.bytes.assign(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
    Walker w{f.bytes};
    w.walk(0, f.bytes.size(), [&](int field, int wt, std::uint64_t, std::size_t off, std::size_t len)
    {
      if (field != 7 || wt != 2) { return; }                   // ModelProto.graph
      w.walk(off, off + len, [&](int gf, int gwt, std::uint64_t, std::size_t goff, std::size_t glen)
      {
        if (gwt != 2) { return; }
        if (gf == 1)                                             // GraphProto.node
        {
          Node n;
          w.walk(goff, goff + glen, [&](int nf, int nwt, std::uint64_t, std::size_t noff, std::size_t nlen)
          {
            if (nwt != 2) { return; }
            if (nf == 1) { n.inputs.push_back(w.str(noff, nlen)); }
            else if (nf == 4) { n.op_type = w.str(noff, nlen); }
          });
          f.nodes.push_back(std::move(n));
        }
        else if (gf == 5)                                        // GraphProto.initializer
        {
          Tensor t;
          w.walk(goff, goff + glen, [&](int tf, int twt, std::uint64_t v, std::size_t toff, std::size_t tlen)
          {
            if (tf == 1 && twt == 0) { t.dims.push_back(static_cast<std::int64_t>(v)); }
            else if (tf == 1 && twt == 2)                        // packed dims
            { std::size_t p = toff; while (p < toff + tlen) { t.dims.push_back(static_cast<std::int64_t>(w.varint(p))); } }
            else if (tf == 2 && twt == 0) { t.data_type = static_cast<std::int32_t>(v); }
            else if (tf == 8 && twt == 2) { t.name = w.str(toff, tlen); }
            else if (tf == 9 && twt == 2) { t.raw_offset = toff; t.raw_length = tlen; }
            else if (tf == 4 && twt == 2) { throw std::runtime_error("ONNX initializer uses float_data, not raw_data: " + t.name); }
          });
          f.initializers.push_back(std::move(t));
        }
        else if (gf == 11)                                       // GraphProto.input (ValueInfoProto)
        {
          w.walk(goff, goff + glen, [&](int vf, int vwt, std::uint64_t, std::size_t voff, std::size_t vlen)
          { if (vf == 1 && vwt == 2) { f.graph_inputs.push_back(w.str(voff, vlen)); } });
        }
      });
    });
    if (f.initializers.empty()) { throw std::runtime_error("no initializers in " + path); }
    return f;
  }

  void OnnxFile::write(const std::string& path) const
  {
    std::ofstream out(path, std::ios::binary);
    auto raw = [&](std::size_t start, std::size_t end)
    { out.write(reinterpret_cast<const char*>(bytes.data() + start), static_cast<std::streamsize>(end - start)); };
    if (metadata_updates_.empty()) { raw(0, bytes.size()); }
    else
    {
      Walker w{bytes};
      auto remaining = metadata_updates_;
      auto entry = [&](const std::string& value)
      {
        std::string field;
        appendString(field, 14, value);  // ModelProto.metadata_props
        out.write(field.data(), static_cast<std::streamsize>(field.size()));
      };
      w.walkRaw(0, bytes.size(), [&](int field, int wt, std::uint64_t, std::size_t off, std::size_t len, std::size_t start, std::size_t end)
      {
        if (field != 14 || wt != 2) { raw(start, end); return; }
        const auto [key, value] = metadataEntry(w, off, len);
        const auto update = metadata_updates_.find(key);
        if (update == metadata_updates_.end()) { raw(start, end); return; }
        remaining.erase(key);
        if (update->second == value) { raw(start, end); return; }
        std::string replaced;
        w.walkRaw(off, off + len, [&](int ef, int ewt, std::uint64_t, std::size_t, std::size_t, std::size_t first, std::size_t last)
        { if (ef != 2 || ewt != 2) { replaced += w.str(first, last - first); } });
        appendString(replaced, 2, update->second);
        entry(replaced);
      });
      for (const auto& [key, value] : remaining)
      {
        std::string added;
        appendString(added, 1, key); appendString(added, 2, value);
        entry(added);
      }
    }
    out.close();
    if (!out) { throw std::runtime_error("cannot write " + path); }
  }

  std::map<std::string, std::string> OnnxFile::metadata() const
  {
    std::map<std::string, std::string> result;
    Walker w{bytes};
    w.walk(0, bytes.size(), [&](int field, int wt, std::uint64_t, std::size_t off, std::size_t len)
    {
      if (field == 14 && wt == 2)
      {
        const auto [key, value] = metadataEntry(w, off, len);
        result[key] = value;
      }
    });
    for (const auto& [key, value] : metadata_updates_) { result[key] = value; }
    return result;
  }

  void OnnxFile::setMetadata(const std::string& key, const std::string& value)
  {
    if (key.empty()) { throw std::invalid_argument("ONNX metadata key must not be empty"); }
    metadata_updates_[key] = value;
  }

  bool OnnxFile::hasInitializer(const std::string& name) const
  { for (const auto& t : initializers) { if (t.name == name) { return true; } } return false; }

  const OnnxFile::Tensor& OnnxFile::initializer(const std::string& name) const
  {
    for (const auto& t : initializers) { if (t.name == name) { return t; } }
    throw std::runtime_error("ONNX has no initializer named " + name);
  }

  namespace
  {
    torch::Tensor asTensor(const OnnxFile& f, const OnnxFile::Tensor& t)
    {
      if (t.data_type != 1) { throw std::runtime_error("initializer is not float32: " + t.name); }
      std::uint64_t n = 1;
      for (auto d : t.dims)
      {
        if (d < 0 || (d > 0 && n > (std::uint64_t{1} << 40) / static_cast<std::uint64_t>(d))) { throw std::runtime_error("implausible initializer shape: " + t.name); }
        n *= static_cast<std::uint64_t>(d);
      }
      if (n * 4 != t.raw_length || t.raw_offset + t.raw_length > f.bytes.size()) { throw std::runtime_error("raw_data length mismatch: " + t.name); }
      auto out = torch::empty(t.dims, torch::kFloat32);
      std::memcpy(out.data_ptr<float>(), f.bytes.data() + t.raw_offset, t.raw_length);
      return out;
    }

    void putTensor(OnnxFile& f, const OnnxFile::Tensor& t, const torch::Tensor& src)
    {
      auto c = src.detach().to(torch::kCPU, torch::kFloat32).contiguous();
      if (static_cast<std::size_t>(c.numel()) * 4 != t.raw_length)
      { throw std::runtime_error("write-back size mismatch for " + t.name); }
      std::memcpy(f.bytes.data() + t.raw_offset, c.data_ptr<float>(), t.raw_length);
    }

    /// torch gate rows [i,f,g,o] -> ONNX [i,o,f,c]
    torch::Tensor gatesToOnnx(const torch::Tensor& w)
    { auto b = w.chunk(4, 0); return torch::cat({b[0], b[3], b[1], b[2]}, 0); }
    /// ONNX [i,o,f,c] -> torch [i,f,g,o]
    torch::Tensor gatesToTorch(const torch::Tensor& w)
    { auto b = w.chunk(4, 0); return torch::cat({b[0], b[2], b[3], b[1]}, 0); }

    struct Names
    {
      std::vector<std::string> lstm_w, lstm_r, lstm_b;   // per layer
      std::vector<std::string> matmul;                   // in graph order
    };
    Names positional(const OnnxFile& f)
    {
      Names n;
      for (const auto& node : f.nodes)
      {
        if (node.op_type == "LSTM" && node.inputs.size() >= 4)
        { n.lstm_w.push_back(node.inputs[1]); n.lstm_r.push_back(node.inputs[2]); n.lstm_b.push_back(node.inputs[3]); }
        else if (node.op_type == "MatMul" && node.inputs.size() == 2 && f.hasInitializer(node.inputs[1]))
        { n.matmul.push_back(node.inputs[1]); }
      }
      if (n.lstm_w.size() != 2) { throw std::runtime_error("expected 2 LSTM nodes, found " + std::to_string(n.lstm_w.size())); }
      if (n.matmul.size() != 2) { throw std::runtime_error("expected 2 MatMul constants, found " + std::to_string(n.matmul.size())); }
      return n;
    }

    /// The name-matched tensors: {torch parameter name, ONNX initializer name}.
    std::vector<std::pair<std::string, std::string>> named(const Head& model)
    {
      const std::string enc = model->ccs ? "ccs_encoder" : "rt_encoder";
      const std::string dec = model->ccs ? "ccs_decoder" : "rt_decoder";
      return {
        {"encoder.cnn_short.weight", enc + ".input_cnn.cnn_short.weight"}, {"encoder.cnn_short.bias", enc + ".input_cnn.cnn_short.bias"},
        {"encoder.cnn_medium.weight", enc + ".input_cnn.cnn_medium.weight"}, {"encoder.cnn_medium.bias", enc + ".input_cnn.cnn_medium.bias"},
        {"encoder.cnn_long.weight", enc + ".input_cnn.cnn_long.weight"}, {"encoder.cnn_long.bias", enc + ".input_cnn.cnn_long.bias"},
        {"encoder.rnn_h0", enc + ".hidden_nn.rnn_h0"}, {"encoder.rnn_c0", enc + ".hidden_nn.rnn_c0"},
        {"dec0.weight", dec + ".nn.0.weight"}, {"dec0.bias", dec + ".nn.0.bias"},
        {"prelu.weight", dec + ".nn.1.weight"},
        {"dec2.weight", dec + ".nn.2.weight"}, {"dec2.bias", dec + ".nn.2.bias"},
      };
    }

    /// A handle to the named parameter (shares storage; copy_ on it edits the model).
    torch::Tensor param(Head& model, const std::string& name)
    {
      auto params = model->named_parameters();
      auto* p = params.find(name);
      if (!p) { throw std::runtime_error("model has no parameter " + name); }
      return *p;
    }
  }

  std::size_t loadWeights(const OnnxFile& f, Head& model)
  {
    torch::NoGradGuard ng;
    std::size_t count = 0;
    for (const auto& [tname, oname] : named(model))
    {
      auto src = asTensor(f, f.initializer(oname));
      auto dst = param(model, tname);
      if (src.sizes() != dst.sizes()) { throw std::runtime_error("shape mismatch " + oname + " vs " + tname); }
      dst.copy_(src); ++count;
    }
    const Names n = positional(f);
    // mod_nn (2,103) is stored transposed (103,2); attn (1,256) as (256,1)
    auto copyExact = [&](const std::string& tname, const torch::Tensor& src)
    {
      auto dst = param(model, tname);
      if (src.sizes() != dst.sizes()) { throw std::runtime_error("shape mismatch loading " + tname + ": ONNX " + std::to_string(src.numel()) + " values vs model " + std::to_string(dst.numel())); }
      dst.copy_(src);
    };
    copyExact("encoder.mod_nn.weight", asTensor(f, f.initializer(n.matmul[0])).t()); ++count;
    copyExact("encoder.attn.weight", asTensor(f, f.initializer(n.matmul[1])).t()); ++count;
    for (int l = 0; l < 2; ++l)
    {
      auto W = asTensor(f, f.initializer(n.lstm_w[l]));   // (2, 512, in)
      auto R = asTensor(f, f.initializer(n.lstm_r[l]));   // (2, 512, 128)
      auto B = asTensor(f, f.initializer(n.lstm_b[l]));   // (2, 1024)
      if (W.dim() != 3 || W.size(0) != 2 || W.size(1) != 4 * HIDDEN || R.sizes() != torch::IntArrayRef({2, 4 * HIDDEN, HIDDEN}) || B.sizes() != torch::IntArrayRef({2, 8 * HIDDEN}))
      { throw std::runtime_error("LSTM layer " + std::to_string(l) + " has an unexpected W/R/B shape"); }
      for (int d = 0; d < 2; ++d)
      {
        const std::string sfx = "_l" + std::to_string(l) + (d ? "_reverse" : "");
        copyExact("encoder.rnn.weight_ih" + sfx, gatesToTorch(W[d]));
        copyExact("encoder.rnn.weight_hh" + sfx, gatesToTorch(R[d]));
        copyExact("encoder.rnn.bias_ih" + sfx, gatesToTorch(B[d].slice(0, 0, 4 * HIDDEN)));
        copyExact("encoder.rnn.bias_hh" + sfx, gatesToTorch(B[d].slice(0, 4 * HIDDEN, 8 * HIDDEN)));
      }
      count += 3;
    }
    return count;
  }

  std::size_t storeWeights(OnnxFile& f, Head& model)
  {
    torch::NoGradGuard ng;
    std::size_t count = 0;
    for (const auto& [tname, oname] : named(model)) { putTensor(f, f.initializer(oname), param(model, tname)); ++count; }
    const Names n = positional(f);
    putTensor(f, f.initializer(n.matmul[0]), param(model, "encoder.mod_nn.weight").t()); ++count;
    putTensor(f, f.initializer(n.matmul[1]), param(model, "encoder.attn.weight").t()); ++count;
    for (int l = 0; l < 2; ++l)
    {
      std::vector<torch::Tensor> W, R, B;
      for (int d = 0; d < 2; ++d)
      {
        const std::string sfx = "_l" + std::to_string(l) + (d ? "_reverse" : "");
        W.push_back(gatesToOnnx(param(model, "encoder.rnn.weight_ih" + sfx)));
        R.push_back(gatesToOnnx(param(model, "encoder.rnn.weight_hh" + sfx)));
        B.push_back(torch::cat({gatesToOnnx(param(model, "encoder.rnn.bias_ih" + sfx)),
                                gatesToOnnx(param(model, "encoder.rnn.bias_hh" + sfx))}, 0));
      }
      putTensor(f, f.initializer(n.lstm_w[l]), torch::stack(W));
      putTensor(f, f.initializer(n.lstm_r[l]), torch::stack(R));
      putTensor(f, f.initializer(n.lstm_b[l]), torch::stack(B));
      count += 3;
    }
    return count;
  }
}
