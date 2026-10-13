// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

#include "CoreMLSupport.h"
#include "OnnxWire.h"

#include <algorithm>
#include <array>
#include <cstring>
#include <map>
#include <sstream>

namespace ODIA::coreml
{
  using onnxwire::Reader;

  namespace
  {
    // ONNX field numbers (onnx.proto3), the only schema knowledge here.
    constexpr std::uint32_t MODEL_GRAPH = 7;
    constexpr std::uint32_t GRAPH_NODE = 1;
    constexpr std::uint32_t GRAPH_INITIALIZER = 5;
    constexpr std::uint32_t GRAPH_INPUT = 11;
    constexpr std::uint32_t NODE_INPUT = 1, NODE_OUTPUT = 2, NODE_NAME = 3, NODE_OP_TYPE = 4, NODE_ATTRIBUTE = 5, NODE_DOMAIN = 7;
    constexpr std::uint32_t ATTR_NAME = 1, ATTR_F = 2, ATTR_TYPE = 20;
    constexpr std::uint64_t ATTR_TYPE_FLOAT = 1;
    constexpr std::uint32_t TENSOR_DIMS = 1, TENSOR_DATA_TYPE = 2, TENSOR_NAME = 8, TENSOR_RAW = 9, TENSOR_FLOAT_DATA = 4;
    constexpr std::uint64_t TENSOR_FLOAT = 1;
    constexpr std::uint32_t VALUE_NAME = 1, VALUE_TYPE = 2;
    constexpr std::uint32_t TYPE_TENSOR = 1;
    constexpr std::uint32_t TENSOR_TYPE_SHAPE = 2;
    constexpr std::uint32_t SHAPE_DIM = 1;
    constexpr std::uint32_t DIM_PARAM = 2;

    /// The scalar float initializers of a graph: name -> value.
    std::map<std::string, float> scalarInitializers(std::string_view graph)
    {
      std::map<std::string, float> out;
      Reader g{graph.data(), graph.data() + graph.size()};
      std::uint32_t field = 0, type = 0;
      std::string_view payload;
      std::uint64_t value = 0;
      while (g.p < g.end)
      {
        if (!g.next(field, type, payload, value)) { return {}; }
        if (field != GRAPH_INITIALIZER || type != 2) { continue; }
        std::string name;
        std::uint64_t elements = 1, data_type = 0;
        bool have = false;
        float f = 0.0f;
        Reader t{payload.data(), payload.data() + payload.size()};
        std::string_view inner;
        while (t.p < t.end)
        {
          if (!t.next(field, type, inner, value)) { return {}; }
          if (field == TENSOR_DIMS && type == 0) { elements *= value; }
          else if (field == TENSOR_DIMS && type == 2)
          {
            Reader d{inner.data(), inner.data() + inner.size()};
            while (d.p < d.end) { if (!d.varint(value)) { return {}; } elements *= value; }
          }
          else if (field == TENSOR_DATA_TYPE && type == 0) { data_type = value; }
          else if (field == TENSOR_NAME && type == 2) { name = std::string(inner); }
          else if (field == TENSOR_RAW && type == 2 && inner.size() == 4) { std::memcpy(&f, inner.data(), 4); have = true; }
          else if (field == TENSOR_FLOAT_DATA && type == 5) { std::memcpy(&f, inner.data(), 4); have = true; }
          else if (field == TENSOR_FLOAT_DATA && type == 2 && inner.size() == 4) { std::memcpy(&f, inner.data(), 4); have = true; }
        }
        if (have && data_type == TENSOR_FLOAT && elements == 1 && !name.empty()) { out[name] = f; }
      }
      return out;
    }

    /// The node re-emitted as LeakyRelu(alpha) if it is a PRelu with a
    /// scalar slope; otherwise the original bytes. True if rewritten.
    bool rewriteNode(std::string_view node, const std::map<std::string, float>& scalars, std::string& out)
    {
      std::string op, domain;
      std::vector<std::string_view> inputs, outputs;
      std::string_view name;
      Reader n{node.data(), node.data() + node.size()};
      std::uint32_t field = 0, type = 0;
      std::string_view payload;
      std::uint64_t value = 0;
      while (n.p < n.end)
      {
        if (!n.next(field, type, payload, value)) { return false; }
        if (type != 2) { continue; }
        if (field == NODE_INPUT) { inputs.push_back(payload); }
        else if (field == NODE_OUTPUT) { outputs.push_back(payload); }
        else if (field == NODE_NAME) { name = payload; }
        else if (field == NODE_OP_TYPE) { op = std::string(payload); }
        else if (field == NODE_DOMAIN) { domain = std::string(payload); }
      }
      if (op != "PRelu" || !(domain.empty() || domain == "ai.onnx") || inputs.size() != 2 || outputs.size() != 1)
      { return false; }
      const auto slope = scalars.find(std::string(inputs[1]));
      if (slope == scalars.end()) { return false; }
      std::string attribute;
      onnxwire::putBytesField(attribute, ATTR_NAME, "alpha");
      onnxwire::putFloatField(attribute, ATTR_F, slope->second);
      onnxwire::putVarintField(attribute, ATTR_TYPE, ATTR_TYPE_FLOAT);
      onnxwire::putBytesField(out, NODE_INPUT, inputs[0]);
      onnxwire::putBytesField(out, NODE_OUTPUT, outputs[0]);
      if (!name.empty()) { onnxwire::putBytesField(out, NODE_NAME, name); }
      onnxwire::putBytesField(out, NODE_OP_TYPE, "LeakyRelu");
      onnxwire::putBytesField(out, NODE_ATTRIBUTE, attribute);
      return true;
    }
  }

  std::size_t rewritePRelu(std::string& model)
  {
    std::uint32_t field = 0, type = 0;
    std::string_view payload;
    std::uint64_t value = 0;
    std::string_view graph;
    {
      Reader m{model.data(), model.data() + model.size()};
      while (m.p < m.end)
      {
        if (!m.next(field, type, payload, value)) { return 0; }
        if (field == MODEL_GRAPH && type == 2) { graph = payload; }
      }
    }
    if (graph.empty()) { return 0; }
    const auto scalars = scalarInitializers(graph);
    if (scalars.empty()) { return 0; }

    // Pass 1: the nodes, noting which slopes the rewrite leaves unused (a
    // slope another node still reads stays; ONNX Runtime would otherwise
    // warn about every unused initializer, once per session).
    std::string nodes;
    nodes.reserve(graph.size());
    std::size_t rewritten = 0;
    std::map<std::string, std::size_t> uses;   // initializer name -> reads left
    {
      Reader g{graph.data(), graph.data() + graph.size()};
      while (g.p < g.end)
      {
        if (!g.next(field, type, payload, value)) { return 0; }
        if (field != GRAPH_NODE || type != 2) { continue; }
        std::string node;
        if (rewriteNode(payload, scalars, node))
        {
          onnxwire::putBytesField(nodes, GRAPH_NODE, node);
          ++rewritten;
          Reader n{payload.data(), payload.data() + payload.size()};
          std::string_view inner;
          while (n.p < n.end) { if (!n.next(field, type, inner, value)) { return 0; } }
          continue;
        }
        onnxwire::putRaw(nodes, field, type, payload, value);
        Reader n{payload.data(), payload.data() + payload.size()};
        std::string_view inner;
        while (n.p < n.end)
        {
          if (!n.next(field, type, inner, value)) { return 0; }
          if (field == NODE_INPUT && type == 2) { ++uses[std::string(inner)]; }
        }
      }
    }
    if (rewritten == 0) { return 0; }
    // Pass 2: everything else, in order, minus the now unused scalar slopes.
    std::string new_graph;
    new_graph.reserve(graph.size());
    bool nodes_written = false;
    {
      Reader g{graph.data(), graph.data() + graph.size()};
      while (g.p < g.end)
      {
        if (!g.next(field, type, payload, value)) { return 0; }
        if (field == GRAPH_NODE && type == 2)
        {
          if (!nodes_written) { new_graph += nodes; nodes_written = true; }
          continue;
        }
        if (field == GRAPH_INITIALIZER && type == 2)
        {
          std::string name;
          Reader t{payload.data(), payload.data() + payload.size()};
          std::string_view inner;
          while (t.p < t.end)
          {
            if (!t.next(field, type, inner, value)) { return 0; }
            if (field == TENSOR_NAME && type == 2) { name = std::string(inner); }
          }
          if (scalars.count(name) != 0 && uses[name] == 0) { continue; }
          onnxwire::putRaw(new_graph, GRAPH_INITIALIZER, 2, payload, 0);
          continue;
        }
        onnxwire::putRaw(new_graph, field, type, payload, value);
      }
    }

    std::string new_model;
    new_model.reserve(model.size() + 64);
    Reader m{model.data(), model.data() + model.size()};
    while (m.p < m.end)
    {
      if (!m.next(field, type, payload, value)) { return 0; }
      if (field == MODEL_GRAPH && type == 2) { onnxwire::putBytesField(new_model, MODEL_GRAPH, new_graph); }
      else { onnxwire::putRaw(new_model, field, type, payload, value); }
    }
    model.swap(new_model);
    return rewritten;
  }

  bool freeDimensionNames(std::string_view model, std::string& rows, std::string& length)
  {
    std::uint32_t field = 0, type = 0;
    std::string_view payload, graph, input;
    std::uint64_t value = 0;
    Reader m{model.data(), model.data() + model.size()};
    while (m.p < m.end)
    {
      if (!m.next(field, type, payload, value)) { return false; }
      if (field == MODEL_GRAPH && type == 2) { graph = payload; }
    }
    Reader g{graph.data(), graph.data() + graph.size()};
    while (g.p < g.end && input.empty())
    {
      if (!g.next(field, type, payload, value)) { return false; }
      if (field == GRAPH_INPUT && type == 2) { input = payload; }   // the first input
    }
    if (input.empty()) { return false; }
    // ValueInfoProto.type -> TypeProto.tensor_type -> Tensor.shape -> dim[]
    std::string_view level = input;
    for (const std::uint32_t want : {VALUE_TYPE, TYPE_TENSOR, TENSOR_TYPE_SHAPE})
    {
      std::string_view found;
      Reader r{level.data(), level.data() + level.size()};
      while (r.p < r.end)
      {
        if (!r.next(field, type, payload, value)) { return false; }
        if (field == want && type == 2) { found = payload; }
      }
      if (found.empty()) { return false; }
      level = found;
    }
    std::vector<std::string> names;
    Reader s{level.data(), level.data() + level.size()};
    while (s.p < s.end)
    {
      if (!s.next(field, type, payload, value)) { return false; }
      if (field != SHAPE_DIM || type != 2) { continue; }
      std::string name;
      Reader d{payload.data(), payload.data() + payload.size()};
      std::string_view inner;
      while (d.p < d.end)
      {
        if (!d.next(field, type, inner, value)) { return false; }
        if (field == DIM_PARAM && type == 2) { name = std::string(inner); }
      }
      names.push_back(name);
    }
    if (names.size() != 2 || names[0].empty() || names[1].empty() || names[0] == names[1]) { return false; }
    rows = names[0];
    length = names[1];
    return true;
  }

  Placement parsePlacement(const std::string& log)
  {
    // The session log, one message per line:
    //   Node(s) placed on [CoreMLExecutionProvider]. Number of nodes: 2
    //    <fused name> (<fused name>)
    //   Node(s) placed on [CPUExecutionProvider]. Number of nodes: 3
    //    OneHot (/meta_nn/OneHot)
    // or "All nodes placed on [X]. Number of nodes: N" for a single provider.
    Placement p;
    std::istringstream in(log);
    std::string line, current;
    while (std::getline(in, line))
    {
      const auto placed = line.find("placed on [");
      if (placed != std::string::npos)
      {
        const auto close = line.find(']', placed);
        current = close == std::string::npos ? "" : line.substr(placed + 11, close - placed - 11);
        p.found = true;
        continue;
      }
      // "  OpType (name)": the op type is the first token.
      const auto first = line.find_first_not_of(' ');
      if (first == std::string::npos || current.empty()) { continue; }
      const auto space = line.find(' ', first);
      const std::string op = line.substr(first, space == std::string::npos ? std::string::npos : space - first);
      if (op.empty() || line.find('(', first) == std::string::npos) { continue; }
      if (current == "CoreMLExecutionProvider") { ++p.coreml_nodes; }
      else if (current == "CPUExecutionProvider") { p.cpu_ops.push_back(op); }
    }
    return p;
  }

  bool heavyOp(std::string_view op)
  {
    static constexpr std::array<std::string_view, 19> heavy = {
      "LSTM", "GRU", "RNN", "MatMul", "FusedMatMul", "Gemm", "Conv", "ConvTranspose", "Einsum",
      "Attention", "MultiHeadAttention", "LayerNormalization", "SkipLayerNormalization",
      "EmbedLayerNormalization", "Softmax", "Gelu", "FastGelu", "BiasGelu", "QuickGelu"};
    return std::find(heavy.begin(), heavy.end(), op) != heavy.end();
  }

  std::uint64_t fnv64(std::string_view bytes)
  {
    std::uint64_t h = 1469598103934665603ULL;
    for (const unsigned char c : bytes)
    {
      h ^= c;
      h *= 1099511628211ULL;
    }
    return h;
  }
} // namespace ODIA::coreml
