// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <cstdint>
#include <cstring>
#include <string>
#include <string_view>

namespace ODIA::onnxwire
{
  /// The protobuf wire format, just enough of it to walk an ONNX file: the
  /// predictor reads the weights it shares between CPU sessions this way,
  /// and rewrites one operator for CoreML. Nothing here knows ONNX's schema
  /// beyond the field numbers its callers name.
  struct Reader
  {
    const char* p;
    const char* end;

    bool varint(std::uint64_t& v)
    {
      v = 0;
      for (int shift = 0; shift < 64 && p < end; shift += 7)
      {
        const auto byte = static_cast<std::uint8_t>(*p++);
        v |= static_cast<std::uint64_t>(byte & 0x7f) << shift;
        if (!(byte & 0x80)) { return true; }
      }
      return false;
    }
    /// Next field: number, wire type and, for length-delimited fields, the
    /// payload; for varints the value. Fixed-width fields are skipped and
    /// their bytes exposed through @p payload too, so a caller can copy them.
    bool next(std::uint32_t& field, std::uint32_t& type, std::string_view& payload, std::uint64_t& value)
    {
      std::uint64_t key = 0;
      if (p >= end || !varint(key)) { return false; }
      field = static_cast<std::uint32_t>(key >> 3);
      type = static_cast<std::uint32_t>(key & 7);
      switch (type)
      {
        case 0: return varint(value);
        case 1: if (end - p < 8) { return false; } payload = std::string_view(p, 8); p += 8; return true;
        case 5: if (end - p < 4) { return false; } payload = std::string_view(p, 4); p += 4; return true;
        case 2:
        {
          std::uint64_t n = 0;
          if (!varint(n) || n > static_cast<std::uint64_t>(end - p)) { return false; }
          payload = std::string_view(p, static_cast<std::size_t>(n));
          p += n;
          return true;
        }
        default: return false;   // groups are not used by ONNX
      }
    }
  };

  inline void putVarint(std::string& out, std::uint64_t v)
  {
    while (v >= 0x80)
    {
      out.push_back(static_cast<char>((v & 0x7f) | 0x80));
      v >>= 7;
    }
    out.push_back(static_cast<char>(v));
  }
  inline void putVarintField(std::string& out, std::uint32_t field, std::uint64_t v)
  {
    putVarint(out, (static_cast<std::uint64_t>(field) << 3) | 0);
    putVarint(out, v);
  }
  inline void putBytesField(std::string& out, std::uint32_t field, std::string_view bytes)
  {
    putVarint(out, (static_cast<std::uint64_t>(field) << 3) | 2);
    putVarint(out, bytes.size());
    out.append(bytes.data(), bytes.size());
  }
  inline void putFloatField(std::string& out, std::uint32_t field, float f)
  {
    putVarint(out, (static_cast<std::uint64_t>(field) << 3) | 5);
    char raw[4];
    std::memcpy(raw, &f, 4);
    out.append(raw, 4);
  }
  /// Re-emits a field exactly as read (any wire type).
  inline void putRaw(std::string& out, std::uint32_t field, std::uint32_t type, std::string_view payload, std::uint64_t value)
  {
    if (type == 0) { putVarintField(out, field, value); return; }
    if (type == 2) { putBytesField(out, field, payload); return; }
    putVarint(out, (static_cast<std::uint64_t>(field) << 3) | type);
    out.append(payload.data(), payload.size());
  }
} // namespace ODIA::onnxwire
