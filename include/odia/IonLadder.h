// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/CONCEPT/Constants.h>

#include <vector>

namespace ODIA
{
  /// Fragment ion m/z of one peptide, bit for bit what OpenMS returns for
  /// peptide.getPrefix(i).getMZ(z, type) (a, b and c ions) and
  /// peptide.getSuffix(i).getMZ(z, type) (x, y and z ions).
  ///
  /// Those build a prefix or suffix sequence per ordinal (a heap allocation
  /// each), and every call re-derives the terminal offset from its empirical
  /// formula. Here the offsets are computed once per process and shared
  /// read-only, and the residue masses once per peptide. The sums are formed
  /// with the same terms in the same order as AASequence::getMonoWeight --
  /// proton mass times charge, the terminal modification of the fragment's
  /// own terminus (N for a/b/c, C for x/y/z), residues first to last, the ion
  /// offset, then the division by the charge -- so the result does not change
  /// by a bit. In particular a suffix ion is summed from its own first
  /// residue, not as a running suffix sum, whose different association would
  /// change the last bit. test/search/identify_ion_ladder.cpp compares every
  /// ion type against OpenMS.
  ///
  /// OpenMS's own path is kept where it would throw or log instead of
  /// returning a mass: a peptide with an unknown residue X (no mass), an
  /// ordinal outside 1 .. size() - 1 and charge 0. The ladder reads the
  /// peptide it was built from; pass that peptide to every call.
  class IonLadder
  {
  public:
    explicit IonLadder(const OpenMS::AASequence& peptide)
    {
      static const OpenMS::Residue* const unknown = OpenMS::ResidueDB::getInstance()->getResidue("X");
      mass_.reserve(peptide.size());
      for (OpenMS::Size k = 0; k < peptide.size(); ++k)
      {
        if (&peptide[k] == unknown) { exact_ = false; }
        mass_.push_back(peptide[k].getMonoWeight(OpenMS::Residue::Internal));
      }
    }

    double b(const OpenMS::AASequence& peptide, OpenMS::Size i, int z) const
    {
      return prefix(peptide, i, z, OpenMS::Residue::BIon, offsets().b);
    }

    double y(const OpenMS::AASequence& peptide, OpenMS::Size i, int z) const
    {
      return suffix(peptide, i, z, OpenMS::Residue::YIon, offsets().y);
    }

    /// Any of the a, b, c, x, y and z ions; another type goes OpenMS's way
    /// (a prefix for NTerminal and the N-terminal ion types, else a suffix).
    double ion(const OpenMS::AASequence& peptide, OpenMS::Size i, int z, OpenMS::Residue::ResidueType type) const
    {
      using OpenMS::Residue;
      const Offsets& o = offsets();
      switch (type)
      {
        case Residue::AIon: return prefix(peptide, i, z, type, o.a);
        case Residue::BIon: return prefix(peptide, i, z, type, o.b);
        case Residue::CIon: return prefix(peptide, i, z, type, o.c);
        case Residue::XIon: return suffix(peptide, i, z, type, o.x);
        case Residue::YIon: return suffix(peptide, i, z, type, o.y);
        case Residue::ZIon: return suffix(peptide, i, z, type, o.z);
        case Residue::NTerminal: return peptide.getPrefix(i).getMZ(z, type);
        default: return peptide.getSuffix(i).getMZ(z, type);
      }
    }

  private:
    struct Offsets
    {
      double a, b, c, x, y, z;
    };

    /// Built once, thread-safely, then only read.
    static const Offsets& offsets()
    {
      using OpenMS::Residue;
      static const Offsets o{Residue::getInternalToAIon().getMonoWeight(), Residue::getInternalToBIon().getMonoWeight(),
                             Residue::getInternalToCIon().getMonoWeight(), Residue::getInternalToXIon().getMonoWeight(),
                             Residue::getInternalToYIon().getMonoWeight(), Residue::getInternalToZIon().getMonoWeight()};
      return o;
    }

    bool exact(OpenMS::Size i, int z) const { return exact_ && i > 0 && i < mass_.size() && z != 0; }

    double prefix(const OpenMS::AASequence& peptide, OpenMS::Size i, int z, OpenMS::Residue::ResidueType type, double offset) const
    {
      if (!exact(i, z)) { return peptide.getPrefix(i).getMZ(z, type); }
      double m = OpenMS::Constants::PROTON_MASS_U * z;
      if (const auto* mod = peptide.getNTerminalModification()) { m += mod->getDiffMonoMass(); }
      for (OpenMS::Size k = 0; k < i; ++k) { m += mass_[k]; }
      return (m + offset) / z;
    }

    double suffix(const OpenMS::AASequence& peptide, OpenMS::Size i, int z, OpenMS::Residue::ResidueType type, double offset) const
    {
      if (!exact(i, z)) { return peptide.getSuffix(i).getMZ(z, type); }
      double m = OpenMS::Constants::PROTON_MASS_U * z;
      if (const auto* mod = peptide.getCTerminalModification()) { m += mod->getDiffMonoMass(); }
      for (OpenMS::Size k = mass_.size() - i; k < mass_.size(); ++k) { m += mass_[k]; }
      return (m + offset) / z;
    }

    std::vector<double> mass_;
    bool exact_ = true;
  };
}
