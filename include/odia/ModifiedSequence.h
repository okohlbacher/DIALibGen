// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <cstddef>
#include <string>
#include <string_view>

namespace ODIA
{
  /// Canonical form of a modified sequence, so a DIA-NN report joins an ODIA
  /// library.
  ///
  /// NOT cosmetic. DIA-NN writes `C(UniMod:4)`; a library generated through
  /// OpenMS writes `C(Carbamidomethyl)`. On an alkylated sample that is every
  /// cysteine-containing precursor -- measured: 3,444 of 37,193 in one run, exactly
  /// the cysteine population -- and a verbatim join drops all of them silently.
  ///
  /// Brackets are parsed BALANCED, so `K(Label:13C(6)15N(2))` is one token.
  /// A token that is neither a known name nor a `UniMod:` accession is passed
  /// through verbatim and, when @p unknown is given, counted there.
  ///
  /// The TSV writer uses the same function, because DIA-NN reads modifications
  /// in UniMod notation: OpenMS names such as `(Carbamidomethyl)` are written
  /// as `(UniMod:4)`.
  std::string canonicalModifiedSequence(std::string_view seq, std::size_t* unknown = nullptr);

} // namespace ODIA
