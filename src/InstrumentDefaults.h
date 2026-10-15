// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <array>
#include <filesystem>
#include <map>
#include <optional>
#include <string>
#include <vector>

/// The generation defaults that depend on the instrument, read from YAML
/// files (docs/usage.md#instrument-defaults): default.yaml for every
/// instrument, and one file per instrument (astral.yaml, timstof.yaml, ...)
/// for what differs. An explicit option (command line, INI, -config) always
/// wins over both.
namespace ODIA::instruments
{
  using Window = std::array<double, 2>;

  /// What one MS2 model slot implies when no instrument file says otherwise.
  struct ModelDefaults
  {
    std::optional<double> nce;
    std::optional<bool> ion_mobility;
  };

  /// One instrument's file.
  struct InstrumentFile
  {
    std::filesystem::path path;
    std::string hash;                 ///< FNV-1a64 of the file's bytes
    std::string name;                 ///< as the file names it, e.g. "Astral"
    std::vector<std::string> aliases;
    std::string model_instrument;     ///< the MS2 model's slot: QE, Lumos, timsTOF, SciexTOF, ThermoTOF
    std::optional<bool> ion_mobility;
    std::optional<double> nce;
    std::optional<Window> precursor_mz, fragment_mz;
    std::optional<std::vector<int>> precursor_charges;

    /// Whether @p instrument is this file's name or one of its aliases (case,
    /// ' ', '-' and '_' ignored; a leading "Orbitrap" and a trailing model
    /// number too, as for the built-in names).
    bool names(const std::string& instrument) const;
  };

  /// default.yaml: the instrument when none is named, and every value an
  /// instrument file leaves out.
  struct DefaultsFile
  {
    std::filesystem::path path;
    std::string hash;
    std::string instrument;
    Window precursor_mz{}, fragment_mz{};
    std::vector<int> precursor_charges;
    std::map<std::string, ModelDefaults> models;   ///< by MS2 model slot
  };

  /// Strict readers: unknown keys, wrong types and values out of range are
  /// errors (std::runtime_error, "<file>: <key>: <problem>").
  InstrumentFile readInstrumentFile(const std::filesystem::path& file);
  DefaultsFile readDefaultsFile(const std::filesystem::path& file);

  /// The directories searched, first wins: $DIALIBGEN_INSTRUMENT_DIR (a
  /// list, ':' or ';' on Windows), the user's (~/.config/DIALibGen/instruments,
  /// ~/Library/Application Support/DIALibGen/instruments on macOS,
  /// %APPDATA%\DIALibGen\instruments on Windows), then the first of
  /// @p data_dirs holding instruments/default.yaml (the shipped files).
  std::vector<std::filesystem::path> searchPath(const std::vector<std::filesystem::path>& data_dirs);

  /// Every *.yaml file of the search path, read at once. A file in an earlier
  /// directory hides one of the same name in a later one (a user's
  /// astral.yaml replaces the shipped one).
  class Catalog
  {
  public:
    explicit Catalog(const std::vector<std::filesystem::path>& dirs);
    const DefaultsFile& defaults() const { return defaults_; }
    /// The file that names @p instrument, the earliest directory first; two
    /// files of one directory naming it are an error. nullptr if none does.
    const InstrumentFile* find(const std::string& instrument) const;
    const std::vector<InstrumentFile>& files() const { return files_; }
    const std::vector<std::filesystem::path>& dirs() const { return dirs_; }

  private:
    std::vector<std::filesystem::path> dirs_;
    DefaultsFile defaults_;
    std::vector<InstrumentFile> files_;
    std::vector<std::size_t> rank_;   ///< the directory of each file
  };
}
