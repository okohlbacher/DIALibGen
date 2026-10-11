// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <charconv>
#include <cstddef>
#include <string_view>
#include <system_error>

// ODIA_FLOAT_FROM_CHARS comes from CMake, which compiles a std::from_chars(double)
// call for the actual target -- on macOS for CMAKE_OSX_DEPLOYMENT_TARGET, where
// Apple's libc++ either lacks the overload or marks it unavailable below the
// macOS release that ships it. Without CMake's answer, the library's own
// feature-test macro decides.
#ifndef ODIA_FLOAT_FROM_CHARS
#  if defined(__cpp_lib_to_chars) && __cpp_lib_to_chars >= 201611L
#    define ODIA_FLOAT_FROM_CHARS 1
#  else
#    define ODIA_FLOAT_FROM_CHARS 0
#  endif
#endif

#if !ODIA_FLOAT_FROM_CHARS
#  include <cerrno>
#  include <cstdlib>
#  include <locale.h>
#  include <string>
#  if defined(__APPLE__)
#    include <xlocale.h>
#  endif
#endif

namespace ODIA
{
  /// Parses all of @p s as a decimal floating-point number, correctly rounded.
  ///
  /// Accepts only the characters `0-9 + - . e E` (no whitespace, no `inf`,
  /// `nan` or hex), an optional leading sign included, and returns false for
  /// anything else, for trailing characters, and for a value outside the
  /// range of double (overflow, or underflow the library reports as such).
  ///
  /// Where the standard library has floating-point `std::from_chars`
  /// (libstdc++ 11 and later, MSVC) that is what parses. Elsewhere -- Apple's
  /// libc++ on the deployment target -- it is `strtod_l` in the "C" locale,
  /// which is what `std::istringstream >> double` called underneath on that
  /// library, with the same ERANGE rejection, so the doubles are the ones the
  /// stream produced; both are correctly rounded, so they are also the ones
  /// `from_chars` produces elsewhere. It replaces a stream per number, about
  /// 1 µs each, which a whole-proteome TSV pays some 400 million times.
  inline bool parseDouble(std::string_view s, double& v)
  {
    if (s.empty() || s.find_first_not_of("0123456789+-.eE") != std::string_view::npos) { return false; }
#if ODIA_FLOAT_FROM_CHARS
    const char* first = s.data();
    const char* const last = first + s.size();
    // from_chars takes no '+'; a stream and strtod do. "+-1" stays invalid.
    if (*first == '+' && s.size() > 1 && s[1] != '-') { ++first; }
    double parsed = 0.0;
    const auto r = std::from_chars(first, last, parsed, std::chars_format::general);
    if (r.ec != std::errc{} || r.ptr != last) { return false; }
    v = parsed;
    return true;
#else
    // One locale per thread: Apple's strtod_l takes the locale's lock (in
    // localeconv_l) on every call, and with one shared locale the parsing
    // threads of a TSV load spent most of their time waiting on it.
    struct CLocale
    {
      locale_t handle = newlocale(LC_ALL_MASK, "C", static_cast<locale_t>(0));
      ~CLocale() { if (handle != static_cast<locale_t>(0)) { freelocale(handle); } }
    };
    static thread_local const CLocale c_locale;
    // strtod needs a terminator. Library fields are short; the heap only for
    // something pathological.
    char small[64];
    std::string large;
    const char* text = small;
    if (s.size() < sizeof(small))
    {
      s.copy(small, s.size());
      small[s.size()] = '\0';
    }
    else
    {
      large.assign(s);
      text = large.c_str();
    }
    char* end = nullptr;
    const int saved = errno;
    errno = 0;
    const double parsed = strtod_l(text, &end, c_locale.handle);
    const bool range_error = errno == ERANGE;
    errno = saved;
    if (end != text + s.size() || range_error) { return false; }
    v = parsed;
    return true;
#endif
  }
} // namespace ODIA
