// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

#include "InstrumentDefaults.h"
#include <odia/DIANNLibraryFile.h>
#include <odia/PeptDeepEncoder.h>
#include <yaml-cpp/yaml.h>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <set>
#include <stdexcept>

namespace fs = std::filesystem;

namespace ODIA::instruments
{
  namespace
  {
    std::string fold(const std::string& name)
    {
      std::string f;
      for (const char c : name)
      {
        if (c == '-' || c == '_' || c == ' ') { continue; }
        f.push_back(static_cast<char>(std::toupper(static_cast<unsigned char>(c))));
      }
      return f;
    }

    /// The folded name without a leading "Orbitrap" and a trailing model number.
    std::string stripped(std::string f)
    {
      if (f.rfind("ORBITRAP", 0) == 0) { f.erase(0, 8); }
      while (!f.empty() && std::isdigit(static_cast<unsigned char>(f.back()))) { f.pop_back(); }
      return f;
    }

    [[noreturn]] void fail(const fs::path& file, const std::string& key, const std::string& problem)
    {
      throw std::runtime_error(file.filename().string() + ": " + (key.empty() ? "" : key + ": ") + problem);
    }

    YAML::Node load(const fs::path& file)
    {
      YAML::Node root;
      try { root = YAML::LoadFile(file.string()); }
      catch (const YAML::Exception& e) { fail(file, "", std::string("not valid YAML (") + e.what() + ")"); }
      if (!root.IsMap()) { fail(file, "", "must be a mapping of keys to values"); }
      return root;
    }

    void onlyKeys(const fs::path& file, const YAML::Node& node, const std::set<std::string>& allowed,
                  const std::string& where = "")
    {
      for (const auto& item : node)
      {
        const std::string key = item.first.as<std::string>();
        if (!allowed.count(key))
        {
          std::string list;
          for (const auto& k : allowed) { list += (list.empty() ? "" : ", ") + k; }
          fail(file, where + key, "unknown key (known: " + list + ")");
        }
      }
    }

    std::string text(const fs::path& file, const YAML::Node& n, const std::string& key)
    {
      if (!n.IsScalar() || n.Scalar().empty()) { fail(file, key, "must be a non-empty string"); }
      return n.Scalar();
    }

    bool boolean(const fs::path& file, const YAML::Node& n, const std::string& key)
    {
      bool v = false;
      if (!n.IsScalar() || !YAML::convert<bool>::decode(n, v)) { fail(file, key, "must be true or false"); }
      return v;
    }

    double number(const fs::path& file, const YAML::Node& n, const std::string& key, double lo, double hi)
    {
      double v = 0;
      if (!n.IsScalar() || !YAML::convert<double>::decode(n, v) || !std::isfinite(v))
      { fail(file, key, "must be a number"); }
      if (!(v > lo && v <= hi)) { fail(file, key, "must be greater than " + std::to_string(static_cast<long long>(lo)) +
                                                  " and at most " + std::to_string(static_cast<long long>(hi))); }
      return v;
    }

    Window window(const fs::path& file, const YAML::Node& n, const std::string& key)
    {
      if (!n.IsSequence() || n.size() != 2) { fail(file, key, "must be [min, max]"); }
      Window w{};
      for (std::size_t i = 0; i < 2; ++i)
      {
        double v = 0;
        if (!n[i].IsScalar() || !YAML::convert<double>::decode(n[i], v) || !std::isfinite(v) || v < 0 || v > 100000)
        { fail(file, key, "values must be numbers between 0 and 100000"); }
        w[i] = v;
      }
      if (!(w[0] < w[1])) { fail(file, key, "needs min < max"); }
      return w;
    }

    std::vector<int> charges(const fs::path& file, const YAML::Node& n, const std::string& key)
    {
      if (!n.IsSequence() || n.size() == 0) { fail(file, key, "must be a non-empty list of charges"); }
      std::vector<int> out;
      for (const auto& item : n)
      {
        int z = 0;
        double d = 0;
        if (!item.IsScalar() || !YAML::convert<double>::decode(item, d) || d != std::floor(d) ||
            !YAML::convert<int>::decode(item, z) || z < 1 || z > 8)
        { fail(file, key, "charges must be whole numbers from 1 to 8"); }
        if (std::find(out.begin(), out.end(), z) != out.end()) { fail(file, key, "charges must not repeat"); }
        out.push_back(z);
      }
      return out;
    }

    std::string modelSlot(const fs::path& file, const YAML::Node& n, const std::string& key)
    {
      const std::string named = text(file, n, key);
      const std::string slot = PeptDeepEncoder::canonicalInstrument(named);
      if (slot.empty()) { fail(file, key, "'" + named + "' is not an MS2 model instrument (QE, Lumos, timsTOF, SciexTOF, ThermoTOF)"); }
      return slot;
    }
  }

  bool InstrumentFile::names(const std::string& instrument) const
  {
    const std::string query = fold(instrument), plain = stripped(query);
    for (const std::string* n : {&name})
    { if (fold(*n) == query || fold(*n) == plain) { return true; } }
    for (const auto& alias : aliases)
    { if (fold(alias) == query || fold(alias) == plain) { return true; } }
    return false;
  }

  InstrumentFile readInstrumentFile(const fs::path& file)
  {
    const YAML::Node root = load(file);
    onlyKeys(file, root, {"name", "aliases", "model_instrument", "ion_mobility", "nce", "precursor_mz", "fragment_mz",
                          "precursor_charges"});
    InstrumentFile f;
    f.path = file;
    f.hash = DIANNLibraryFile::hashFile(file.string());
    if (!root["name"]) { fail(file, "name", "missing (the instrument's name)"); }
    f.name = text(file, root["name"], "name");
    if (!root["model_instrument"]) { fail(file, "model_instrument", "missing (the MS2 model's instrument: QE, Lumos, timsTOF, SciexTOF or ThermoTOF)"); }
    f.model_instrument = modelSlot(file, root["model_instrument"], "model_instrument");
    if (const auto n = root["aliases"])
    {
      if (!n.IsSequence()) { fail(file, "aliases", "must be a list of names"); }
      for (const auto& a : n) { f.aliases.push_back(text(file, a, "aliases")); }
    }
    if (const auto n = root["ion_mobility"]) { f.ion_mobility = boolean(file, n, "ion_mobility"); }
    if (const auto n = root["nce"]) { f.nce = number(file, n, "nce", 0, 100); }
    if (const auto n = root["precursor_mz"]) { f.precursor_mz = window(file, n, "precursor_mz"); }
    if (const auto n = root["fragment_mz"]) { f.fragment_mz = window(file, n, "fragment_mz"); }
    if (const auto n = root["precursor_charges"]) { f.precursor_charges = charges(file, n, "precursor_charges"); }
    return f;
  }

  DefaultsFile readDefaultsFile(const fs::path& file)
  {
    const YAML::Node root = load(file);
    onlyKeys(file, root, {"instrument", "precursor_mz", "fragment_mz", "precursor_charges", "models"});
    DefaultsFile d;
    d.path = file;
    d.hash = DIANNLibraryFile::hashFile(file.string());
    for (const char* key : {"instrument", "precursor_mz", "fragment_mz", "precursor_charges"})
    { if (!root[key]) { fail(file, key, "missing"); } }
    d.instrument = text(file, root["instrument"], "instrument");
    d.precursor_mz = window(file, root["precursor_mz"], "precursor_mz");
    d.fragment_mz = window(file, root["fragment_mz"], "fragment_mz");
    d.precursor_charges = charges(file, root["precursor_charges"], "precursor_charges");
    if (const auto models = root["models"])
    {
      if (!models.IsMap()) { fail(file, "models", "must map MS2 model instruments to {nce, ion_mobility}"); }
      for (const auto& item : models)
      {
        const std::string slot = modelSlot(file, item.first, "models");
        const std::string where = "models." + slot;
        if (!item.second.IsMap()) { fail(file, where, "must be a mapping with nce and ion_mobility"); }
        onlyKeys(file, item.second, {"nce", "ion_mobility"}, where + ".");
        ModelDefaults m;
        if (const auto n = item.second["nce"]) { m.nce = number(file, n, where + ".nce", 0, 100); }
        if (const auto n = item.second["ion_mobility"]) { m.ion_mobility = boolean(file, n, where + ".ion_mobility"); }
        if (!d.models.emplace(slot, m).second) { fail(file, where, "listed twice"); }
      }
    }
    return d;
  }

  std::vector<fs::path> searchPath(const std::vector<fs::path>& data_dirs)
  {
    std::vector<fs::path> out;
#ifdef _WIN32
    const char separator = ';';
#else
    const char separator = ':';
#endif
    if (const char* env = std::getenv("DIALIBGEN_INSTRUMENT_DIR"); env != nullptr && *env != '\0')
    {
      std::string list = env;
      std::size_t start = 0;
      while (start <= list.size())
      {
        const std::size_t end = std::min(list.find(separator, start), list.size());
        if (end > start) { out.emplace_back(list.substr(start, end - start)); }
        start = end + 1;
      }
    }
#if defined(_WIN32)
    if (const char* appdata = std::getenv("APPDATA"); appdata != nullptr && *appdata != '\0')
    { out.push_back(fs::path(appdata) / "DIALibGen" / "instruments"); }
#elif defined(__APPLE__)
    if (const char* home = std::getenv("HOME"); home != nullptr && *home != '\0')
    { out.push_back(fs::path(home) / "Library" / "Application Support" / "DIALibGen" / "instruments"); }
#else
    if (const char* xdg = std::getenv("XDG_CONFIG_HOME"); xdg != nullptr && *xdg != '\0')
    { out.push_back(fs::path(xdg) / "DIALibGen" / "instruments"); }
    else if (const char* home = std::getenv("HOME"); home != nullptr && *home != '\0')
    { out.push_back(fs::path(home) / ".config" / "DIALibGen" / "instruments"); }
#endif
    std::error_code ec;
    for (const auto& dir : data_dirs)
    {
      if (fs::exists(dir / "instruments" / "default.yaml", ec)) { out.push_back(dir / "instruments"); break; }
    }
    return out;
  }

  Catalog::Catalog(const std::vector<fs::path>& dirs) : dirs_(dirs)
  {
    std::set<std::string> seen;   // file names: an earlier directory hides a later one
    bool have_defaults = false;
    std::error_code ec;
    for (std::size_t rank = 0; rank < dirs.size(); ++rank)
    {
      if (!fs::is_directory(dirs[rank], ec)) { continue; }
      std::vector<fs::path> entries;
      for (fs::directory_iterator it(dirs[rank], ec), end; !ec && it != end; it.increment(ec))
      {
        const std::string ext = it->path().extension().string();
        if (it->is_regular_file(ec) && (ext == ".yaml" || ext == ".yml")) { entries.push_back(it->path()); }
      }
      std::sort(entries.begin(), entries.end());
      for (const auto& file : entries)
      {
        const std::string name = file.filename().string();
        if (!seen.insert(name).second) { continue; }
        if (file.stem() == "default")
        {
          if (have_defaults) { continue; }   // default.yaml and default.yml: the first
          defaults_ = readDefaultsFile(file);
          have_defaults = true;
          continue;
        }
        files_.push_back(readInstrumentFile(file));
        rank_.push_back(rank);
      }
    }
    if (!have_defaults)
    {
      std::string list;
      for (const auto& d : dirs) { list += (list.empty() ? "" : ", ") + d.string(); }
      throw std::runtime_error("no default.yaml among the instrument defaults (searched: " +
                               (list.empty() ? std::string("nothing") : list) + "); set DIALIBGEN_DATA_DIR to the "
                               "directory holding instruments/default.yaml");
    }
  }

  const InstrumentFile* Catalog::find(const std::string& instrument) const
  {
    const InstrumentFile* found = nullptr;
    std::size_t found_rank = 0;
    for (std::size_t i = 0; i < files_.size(); ++i)
    {
      if (!files_[i].names(instrument)) { continue; }
      if (found == nullptr) { found = &files_[i]; found_rank = rank_[i]; continue; }
      if (rank_[i] == found_rank)
      {
        throw std::runtime_error("instrument '" + instrument + "' is named by two files of one directory: " +
                                 found->path.string() + " and " + files_[i].path.string());
      }
    }
    return found;
  }
}
