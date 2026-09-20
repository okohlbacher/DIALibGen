// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <filesystem>
#include <random>
#include <stdexcept>
#include <string>
#include <system_error>
#include <utility>
#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#else
#include <cerrno>
#include <sys/stat.h>
#endif

namespace ODIA
{
  /// Write temporaryPath(), close it successfully, then commit(). A failed or
  /// interrupted write never changes the destination; exception unwinding also
  /// removes the temporary file. An exclusive directory reserves its name.
  class AtomicFile
  {
  public:
    explicit AtomicFile(std::filesystem::path destination) : destination_(std::move(destination))
    {
      auto parent = destination_.parent_path();
      if (parent.empty()) { parent = "."; }
      std::random_device random;
      for (unsigned attempt = 0; attempt < 16; ++attempt)
      {
        directory_ = parent / (".dialibgen-tmp-" + std::to_string(random()) + "-" + std::to_string(random()));
        temporary_ = directory_ / destination_.filename();
        std::error_code ec;
#ifdef _WIN32
        if (std::filesystem::create_directory(directory_, ec)) { return; }
#else
        if (::mkdir(directory_.c_str(), 0700) == 0) { return; }
        if (errno != EEXIST) { ec = std::error_code(errno, std::generic_category()); }
#endif
        if (ec) { throw std::runtime_error("cannot stage output " + destination_.string() + ": " + ec.message()); }
      }
      throw std::runtime_error("cannot reserve a temporary output for " + destination_.string());
    }

    AtomicFile(const AtomicFile&) = delete;
    AtomicFile& operator=(const AtomicFile&) = delete;
    ~AtomicFile()
    {
      std::error_code ignored;
      std::filesystem::remove(temporary_, ignored);
      std::filesystem::remove(directory_, ignored);
    }

    const std::filesystem::path& temporaryPath() const { return temporary_; }

    void commit()
    {
      std::error_code ec;
#ifdef _WIN32
      if (!MoveFileExW(temporary_.c_str(), destination_.c_str(), MOVEFILE_REPLACE_EXISTING | MOVEFILE_WRITE_THROUGH))
      { ec = std::error_code(static_cast<int>(GetLastError()), std::system_category()); }
#else
      std::filesystem::rename(temporary_, destination_, ec);
#endif
      if (ec) { throw std::runtime_error("cannot replace output " + destination_.string() + ": " + ec.message()); }
    }

  private:
    std::filesystem::path destination_, directory_, temporary_;
  };
}
