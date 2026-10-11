// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <OpenMS/CHEMISTRY/AASequence.h>

#include <string_view>

namespace ODIA
{
  /// AASequence::fromString(@p sequence), without OpenMS's locks in the
  /// common case.
  ///
  /// fromString looks every modified residue up through
  /// ResidueDB::getModifiedResidue, which holds the process-wide named
  /// critical section `ResidueDB`: with 96 threads parsing, every one of them
  /// was found spinning on it (gdb, rt/unique_parse). Here a sequence made
  /// only of one-letter residues, each optionally followed by one
  /// modification in round brackets ("PEPC(Carbamidomethyl)TM(Oxidation)K"),
  /// is assembled from one-residue sequences that fromString built once per
  /// distinct token and that are then only read: the same canonical Residue
  /// pointers, in the same order, as fromString's. Anything else -- terminal
  /// modifications, dot notation, square brackets, nested brackets, a token
  /// fromString would not resolve to one internal residue -- goes to
  /// fromString itself, errors included.
  ///
  /// Thread-safe. ODIA_PARSE_CHECK=1 parses every sequence both ways and
  /// throws std::logic_error if the residues or termini differ.
  OpenMS::AASequence parsePeptide(std::string_view sequence);
} // namespace ODIA
