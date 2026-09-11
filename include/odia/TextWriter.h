// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <charconv>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <vector>

namespace ODIA
{

  /// Buffered text output whose numbers are rendered with std::to_chars.
  ///
  /// Both TSV writers were formatting-bound, not disk-bound. Measured on this
  /// cluster's node-local NVMe:
  ///
  ///     std::ofstream operator<<        25.4 MB/s
  ///     fwrite of pre-rendered bytes  1,042 MB/s      41x
  ///
  /// and both writers sat on the ofstream line, nowhere near the disk line:
  /// 34.4 MB/s for the 7.28 GB library TSV (211.7 s, 26% of Phase 1) and
  /// ~24 MB/s for the 7.16 GB chromatogram TSV (~301 s, 33% of Phase 2). The
  /// cost is `num_put` dispatch and locale, once per field.
  ///
  /// **The numbers must not move.** `defaultfloat` at precision P is specified
  /// as printf's `%.Pg`, and `to_chars` with `chars_format::general` and the
  /// same precision is specified the same way -- so `number()` takes the
  /// precision EXPLICITLY at every call site. to_chars' default is
  /// shortest-round-trip, which prints up to 17 digits where the stream printed
  /// 6, so an omitted argument would not be a rounding difference but a
  /// different file. Checked two ways: `odia_tsv_writers` diffs both real
  /// output files byte for byte against an ofstream reference, and its
  /// `formatting` case compares the two renderings directly over ~4.8 M values
  /// -- random bit patterns, denormals, infinities, NaNs and realistic m/z,
  /// intensity and retention-time magnitudes -- at precisions 6, 9 and 10.
  /// Zero disagree.
  ///
  /// Failure is reported at `close()`, not at open and not in the destructor. A
  /// stream checked only at open reports success on a full disk, an exceeded
  /// quota or a broken mount, and leaves a truncated file behind.
  class TextWriter
  {
  public:
    explicit TextWriter(const std::string& path) : path_(path)
    {
      file_ = std::fopen(path.c_str(), "wb");
      if (file_ == nullptr) { throw std::runtime_error("cannot write: " + path); }
      buffer_.resize(CAPACITY);
    }

    TextWriter(const TextWriter&) = delete;
    TextWriter& operator=(const TextWriter&) = delete;

    /// Closes without reporting. Only reached when close() was not called,
    /// which is the exception path; the error there is the one being thrown.
    ~TextWriter()
    {
      if (file_ != nullptr) { std::fclose(file_); }
    }

    void put(char c)
    {
      need(1);
      buffer_[used_++] = c;
    }

    void put(std::string_view s)
    {
      // A protein group has no bound, so a field larger than the buffer goes
      // straight out rather than growing it.
      if (s.size() >= CAPACITY)
      {
        flush();
        write(s.data(), s.size());
        return;
      }
      need(s.size());
      std::memcpy(buffer_.data() + used_, s.data(), s.size());
      used_ += s.size();
    }

    void integer(long long v)
    {
      need(SCRATCH);
      const auto r = std::to_chars(buffer_.data() + used_, buffer_.data() + CAPACITY, v);
      check(r.ec);
      used_ = static_cast<std::size_t>(r.ptr - buffer_.data());
    }

    /// @p v as `std::ostream` would write it under `defaultfloat` with
    /// `precision(digits)` -- i.e. printf's `%.<digits>g`.
    ///
    /// Takes a double rather than a float because the stream did: `num_put` has
    /// no float overload, so `out << someFloat` formats the promoted double,
    /// and rendering the float directly is a different set of digits.
    void number(double v, int digits)
    {
      need(SCRATCH);
      const auto r = std::to_chars(buffer_.data() + used_, buffer_.data() + CAPACITY,
                                   v, std::chars_format::general, digits);
      check(r.ec);
      used_ = static_cast<std::size_t>(r.ptr - buffer_.data());
    }

    /// Flush and close, reporting a write that failed. Idempotent.
    void close()
    {
      if (file_ == nullptr) { return; }
      flush();
      const bool bad = std::ferror(file_) != 0;
      const bool closed = std::fclose(file_) == 0;
      file_ = nullptr;
      if (bad || !closed)
      {
        throw std::runtime_error("failed while writing: " + path_);
      }
    }

  private:
    /// 1 MiB: large enough that one fwrite covers thousands of rows, small
    /// enough to stay in cache. The gain is in not calling num_put, not in the
    /// size of the write, so this is not a tuned number.
    static constexpr std::size_t CAPACITY = 1u << 20;
    /// Headroom one number can need. `%.10g` of a double is at most ~24 bytes.
    static constexpr std::size_t SCRATCH = 64;

    void need(std::size_t n)
    {
      if (used_ + n > CAPACITY) { flush(); }
    }

    void flush()
    {
      if (used_ == 0) { return; }
      write(buffer_.data(), used_);
      used_ = 0;
    }

    void write(const char* data, std::size_t n)
    {
      // A short write is not retried: on a regular file it means the write
      // failed, and ferror() is checked at close().
      std::fwrite(data, 1, n, file_);
    }

    static void check(std::errc ec)
    {
      if (ec != std::errc{})
      {
        // Unreachable while SCRATCH exceeds what any format below can produce.
        // Stated rather than assumed: to_chars reports a short buffer by
        // writing nothing, so an unchecked failure is a silently missing field.
        throw std::runtime_error("number did not fit its buffer");
      }
    }

    std::string path_;
    std::FILE* file_ = nullptr;
    std::vector<char> buffer_;
    std::size_t used_ = 0;
  };

} // namespace ODIA
