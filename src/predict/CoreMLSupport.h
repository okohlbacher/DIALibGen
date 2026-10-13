// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

/// What running the PeptDeep models through ONNX Runtime's CoreML execution
/// provider needs besides the provider itself: a graph it can compile, static
/// shapes (the session options do that), and a reading of where ONNX Runtime
/// placed the nodes, since the provider silently keeps on the CPU whatever
/// it cannot express.
namespace ODIA::coreml
{
  /// Rewrites every PRelu whose slope is a scalar initializer into the
  /// equivalent LeakyRelu, in place; returns how many. The PeptDeep decoders
  /// use nn.PReLU(1), and the CoreML provider's MLProgram builder emits an
  /// alpha of the wrong rank for it once the shapes are static ("Invalid
  /// 'alpha' shape"), which fails the whole compile; LeakyRelu with the same
  /// slope is the same function and a native CoreML operator. False from the
  /// return value's point of view (0) when the bytes do not parse.
  std::size_t rewritePRelu(std::string& model);

  /// The names of the two free dimensions of the first graph input, which
  /// the PeptDeep models declare as [batch_size, seq_length]: ONNX Runtime
  /// fixes them per session by name. False if the input has no such names.
  bool freeDimensionNames(std::string_view model, std::string& rows, std::string& length);

  /// ONNX Runtime's node placement, read from the session's verbose log.
  struct Placement
  {
    std::size_t coreml_nodes = 0;         ///< fused CoreML partitions
    std::vector<std::string> cpu_ops;     ///< op types of the nodes left on the CPU, in graph order
    bool found = false;                   ///< the log had a placement at all
  };
  Placement parsePlacement(const std::string& verbose_log);

  /// The ops that carry the models' arithmetic. One of them on the CPU means
  /// the model does not run on CoreML in any sense worth the name, whatever
  /// the partition count says; the others (OneHot, shape ops) are cheap and
  /// routinely left on the CPU by design.
  bool heavyOp(std::string_view op_type);

  /// FNV-1a over the bytes: names the compiled-model cache per model content.
  std::uint64_t fnv64(std::string_view bytes);
} // namespace ODIA::coreml
