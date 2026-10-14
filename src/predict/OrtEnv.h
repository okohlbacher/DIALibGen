// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <onnxruntime_cxx_api.h>

#include <string>

namespace ODIA::ort
{
  /// The process's ONNX Runtime environment. ONNX Runtime keeps ONE
  /// environment per process and reference-counts it; the logger is the first
  /// creator's, so every session goes through this one, whose logger prints
  /// warnings and errors as the default does and can hand a session's verbose
  /// log to the code that asked for it (LogCapture).
  Ort::Env& env();

  /// Releases env() at exit before the static destructors of a provider
  /// library loaded so far run (OrtEnv.cpp). Called once the CUDA provider
  /// has been appended to a session's options; calls after the first do
  /// nothing.
  void releaseBeforeProviders();

  /// Collects the verbose messages of the sessions whose log id is @p id
  /// (SessionOptions::SetLogId) for as long as it lives. Only the node
  /// placement lines are kept: that is what a caller wants to know, and a
  /// verbose session log is otherwise thousands of lines.
  class LogCapture
  {
  public:
    explicit LogCapture(std::string id);
    ~LogCapture();
    LogCapture(const LogCapture&) = delete;
    LogCapture& operator=(const LogCapture&) = delete;
    const std::string& id() const { return id_; }
    std::string text() const;

  private:
    std::string id_;
  };
} // namespace ODIA::ort
