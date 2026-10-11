// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <odia/PeptideParser.h>

#include <OpenMS/CHEMISTRY/Residue.h>

#include <array>
#include <atomic>
#include <cstdlib>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

using namespace OpenMS;

namespace ODIA
{
  namespace
  {
    struct TokenHash
    {
      using is_transparent = void;
      std::size_t operator()(std::string_view s) const noexcept { return std::hash<std::string_view>{}(s); }
    };

    /// What each token resolved to: a one-residue sequence, or null when the
    /// token has to go through fromString in context. Immutable once
    /// published; a new token publishes a copy with it added.
    struct Table
    {
      std::unordered_map<std::string, const AASequence*, TokenHash, std::equal_to<>> token;
      std::array<const AASequence*, 26> letter{};   ///< unmodified 'A'..'Z', a copy of `token`'s
      std::array<bool, 26> known{};
    };

    std::mutex g_learn;
    std::vector<std::unique_ptr<AASequence>> g_residues;   ///< owned for the process
    std::vector<std::unique_ptr<Table>> g_tables;          ///< every version, kept: readers may hold any
    std::atomic<const Table*> g_table{nullptr};

    /// fromString on the token alone. One residue, no terminal modification:
    /// the residue fromString puts there in any position, because an internal
    /// modification is looked up by residue and name only.
    const AASequence* resolve(const std::string& token)
    {
      try
      {
        auto one = std::make_unique<AASequence>(AASequence::fromString(token));
        if (one->size() != 1 || one->hasNTerminalModification() || one->hasCTerminalModification()) { return nullptr; }
        g_residues.push_back(std::move(one));
        return g_residues.back().get();
      }
      catch (...) { return nullptr; }
    }

    const AASequence* learn(std::string_view token)
    {
      const std::lock_guard<std::mutex> lock(g_learn);
      const Table* now = g_table.load(std::memory_order_acquire);
      if (now != nullptr)
      {
        if (auto it = now->token.find(token); it != now->token.end()) { return it->second; }
      }
      auto next = now != nullptr ? std::make_unique<Table>(*now) : std::make_unique<Table>();
      const std::string key(token);
      const AASequence* value = resolve(key);
      next->token.emplace(key, value);
      if (key.size() == 1) { next->letter[key[0] - 'A'] = value; next->known[key[0] - 'A'] = true; }
      g_tables.push_back(std::move(next));
      g_table.store(g_tables.back().get(), std::memory_order_release);
      return value;
    }

    const AASequence* lookup(std::string_view token)
    {
      const Table* t = g_table.load(std::memory_order_acquire);
      if (t != nullptr)
      {
        if (token.size() == 1 && t->known[token[0] - 'A']) { return t->letter[token[0] - 'A']; }
        if (auto it = t->token.find(token); it != t->token.end()) { return it->second; }
      }
      return learn(token);
    }

    /// The lock-free form, or false when fromString has to do it.
    bool fastParse(std::string_view s, AASequence& out)
    {
      std::size_t i = 0;
      while (i < s.size())
      {
        const char c = s[i];
        if (c < 'A' || c > 'Z') { return false; }
        std::size_t end = i + 1;
        if (end < s.size() && s[end] == '(')
        {
          const std::size_t close = s.find_first_of("()", end + 1);
          if (close == std::string_view::npos || s[close] != ')') { return false; }   // nested or unclosed
          end = close + 1;
        }
        const AASequence* one = lookup(s.substr(i, end - i));
        if (one == nullptr) { return false; }
        out += *one;
        i = end;
      }
      return !s.empty();
    }

    bool checking()
    {
      static const bool on = [] {
        const char* v = std::getenv("ODIA_PARSE_CHECK");
        return v != nullptr && std::string_view(v) == "1";
      }();
      return on;
    }
  } // namespace

  AASequence parsePeptide(std::string_view sequence)
  {
    AASequence fast;
    if (!fastParse(sequence, fast)) { return AASequence::fromString(std::string(sequence)); }
    if (checking())
    {
      const AASequence ref = AASequence::fromString(std::string(sequence));
      bool same = ref.size() == fast.size() && ref.getNTerminalModification() == fast.getNTerminalModification() &&
                  ref.getCTerminalModification() == fast.getCTerminalModification();
      for (Size k = 0; same && k < ref.size(); ++k) { same = &ref[k] == &fast[k]; }
      if (!same) { throw std::logic_error("parsePeptide differs from AASequence::fromString for " + std::string(sequence)); }
    }
    return fast;
  }
} // namespace ODIA
