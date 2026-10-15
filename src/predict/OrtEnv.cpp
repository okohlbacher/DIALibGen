// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

#include "OrtEnv.h"

#include <odia/DeviceRoute.h>

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <map>
#include <mutex>

namespace ODIA::ort
{
  namespace
  {
    struct Captures
    {
      std::mutex mutex;
      std::map<std::string, std::string> buffers;   // log id -> placement lines
    };
    Captures& captures()
    {
      static Captures c;
      return c;
    }

    void ORT_API_CALL logMessage(void* /*param*/, OrtLoggingLevel severity, const char* /*category*/,
                                 const char* logid, const char* location, const char* message) noexcept
    {
      const bool placement = location != nullptr && std::strstr(location, "VerifyEachNodeIsAssignedToAnEp") != nullptr;
      if (severity == ORT_LOGGING_LEVEL_VERBOSE)
      {
        if (!placement || logid == nullptr) { return; }
        auto& c = captures();
        const std::lock_guard<std::mutex> lock(c.mutex);
        if (const auto it = c.buffers.find(logid); it != c.buffers.end())
        {
          it->second.append(message);
          it->second.push_back('\n');
        }
        return;
      }
      if (severity < ORT_LOGGING_LEVEL_WARNING) { return; }
      // ONNX Runtime's warnings are for diagnosis (-debug 1). Its errors
      // show: some (a CUDA provider library that cannot load one of its own
      // dependencies) are logged and then reported only as a generic
      // failure.
      if (severity == ORT_LOGGING_LEVEL_WARNING && ODIA::diagnosticsLevel() == 0) { return; }
      // Two notices every CoreML session logs at warning level, which the
      // predictor reports in its own words (the provider line, and a WARNING
      // when it matters): how many nodes the provider took, and that some
      // nodes are on the CPU.
      if (placement && severity == ORT_LOGGING_LEVEL_WARNING) { return; }
      if (location != nullptr && std::strstr(location, "coreml_execution_provider.cc") != nullptr &&
          severity == ORT_LOGGING_LEVEL_WARNING)
      { return; }
      const char level = severity == ORT_LOGGING_LEVEL_WARNING ? 'W' : severity == ORT_LOGGING_LEVEL_ERROR ? 'E' : 'F';
      std::cerr << "[" << level << ":onnxruntime:" << (logid != nullptr ? logid : "") << ", "
                << (location != nullptr ? location : "") << "] " << (message != nullptr ? message : "") << "\n";
    }
  }

  Ort::Env& env()
  {
    static Ort::Env instance(ORT_LOGGING_LEVEL_WARNING, "odia", logMessage, nullptr);
    return instance;
  }

  void releaseBeforeProviders()
  {
    // Released by the static's destructor at exit, the environment unloads
    // ONNX Runtime's provider libraries. The CUDA provider library is loaded
    // after the environment exists, so its static destructors run first: its
    // static kernel registry was already freed when the unload released it a
    // second time, a write into freed heap memory, which glibc caught as
    // "corrupted double-linked list" after most GPU runs of 0.11.4. An exit
    // handler registered once the library is loaded runs before them and
    // releases the environment while the library is intact; the destructor
    // then has nothing left to release. (Never releasing it is no way out:
    // ONNX Runtime 1.22, conda-forge's on macOS x64, deletes an environment
    // still held at exit during its own static destruction, and aborts.)
    // env() first: the handler must run before the static's destructor.
    static std::once_flag once;
    std::call_once(once, [] {
      env();
      std::atexit([] { Ort::GetApi().ReleaseEnv(env().release()); });
    });
  }

  LogCapture::LogCapture(std::string id) : id_(std::move(id))
  {
    auto& c = captures();
    const std::lock_guard<std::mutex> lock(c.mutex);
    c.buffers[id_].clear();
  }

  LogCapture::~LogCapture()
  {
    auto& c = captures();
    const std::lock_guard<std::mutex> lock(c.mutex);
    c.buffers.erase(id_);
  }

  std::string LogCapture::text() const
  {
    auto& c = captures();
    const std::lock_guard<std::mutex> lock(c.mutex);
    const auto it = c.buffers.find(id_);
    return it == c.buffers.end() ? std::string() : it->second;
  }
} // namespace ODIA::ort
