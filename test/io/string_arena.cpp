// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

// StringArena's index: handles in first-seen
// order, deduplication, the precomputed-hash overload, growth, reserve and
// releaseLookup's semantics, checked against a reference std::unordered_map.
#include <odia/Library.h>

#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

using ODIA::StringArena;

namespace
{
  void check(bool condition, const std::string& message)
  { if (!condition) { throw std::runtime_error(message); } }

  /// Interns @p strings (with repeats) into @p arena and into a reference
  /// map assigning handles in first-seen order; every handle must agree.
  void againstReference(StringArena& arena, const std::vector<std::string>& strings, bool precomputed,
                        std::unordered_map<std::string, std::uint32_t>& reference, std::uint32_t& next)
  {
    for (const auto& s : strings)
    {
      const auto h = precomputed ? arena.intern(s, StringArena::hashOf(s)) : arena.intern(s);
      const auto [it, fresh] = reference.emplace(s, next);
      if (fresh) { ++next; }
      check(h == it->second, "handle of '" + s + "' is " + std::to_string(h) + ", expected " + std::to_string(it->second));
      check(arena.get(h) == s, "get() of '" + s + "' returned '" + std::string(arena.get(h)) + "'");
    }
    check(arena.size() == next, "size() " + std::to_string(arena.size()) + ", expected " + std::to_string(next));
  }
}

int main()
{
  try
  {
    std::mt19937_64 rng(42);
    std::vector<std::string> strings;
    // Peptide-like strings with many repeats, the empty string, and strings
    // that share prefixes and lengths (tag collisions are then likely enough
    // to exercise the character comparison).
    for (int i = 0; i < 200000; ++i)
    {
      const auto n = rng() % 50000;
      strings.push_back(n == 0 ? std::string() : "PEPTIDE" + std::to_string(n) + "K");
    }
    for (const bool precomputed : {false, true})
    {
      StringArena arena;
      std::unordered_map<std::string, std::uint32_t> reference;
      std::uint32_t next = 0;
      againstReference(arena, strings, precomputed, reference, next);   // grows many times
      againstReference(arena, strings, !precomputed, reference, next);  // all repeats now
    }
    {
      StringArena reserved;
      reserved.reserve(100000, 1u << 22);
      std::unordered_map<std::string, std::uint32_t> reference;
      std::uint32_t next = 0;
      againstReference(reserved, strings, true, reference, next);
    }
    {
      // releaseLookup: strings interned afterwards are stored again, as by
      // the unordered_map index it replaced, and deduplicated among themselves.
      StringArena arena;
      const auto a = arena.intern("AAA");
      arena.releaseLookup();
      arena.releaseLookup();   // twice is safe
      const auto b = arena.intern("AAA");
      const auto c = arena.intern("AAA");
      check(a == 0 && b == 1 && c == 1, "releaseLookup: handles " + std::to_string(a) + "," + std::to_string(b) + "," + std::to_string(c));
      check(arena.get(a) == "AAA" && arena.get(b) == "AAA", "releaseLookup: strings lost");
      check(arena.get(7) == std::string_view(), "get() of an unknown handle is not empty");
      check(arena.footprintBytes() >= arena.bytes(), "footprintBytes() below bytes()");
    }
    {
      // A string seen before is found without growing the index.
      StringArena arena;
      for (int i = 0; i < 8; ++i) { arena.intern("S" + std::to_string(i)); }
      const auto before = arena.footprintBytes();
      check(arena.intern("S3") == 3, "repeat lookup returned a new handle");
      check(arena.footprintBytes() == before, "a repeat lookup grew the index");
    }
    {
      // Moves leave the source empty and usable; the target keeps its index.
      StringArena a;
      a.intern("old");
      a.releaseLookup();
      StringArena b(std::move(a));
      const auto x = a.intern("x");
      a.reserve(32, 0);
      const auto y = a.intern("x");
      check(x == 0 && y == 0, "moved-from arena: handles " + std::to_string(x) + "," + std::to_string(y));
      check(b.get(0) == "old" && b.intern("new") == 1 && b.intern("new") == 1, "moved-to arena lost its strings");
      StringArena c;
      c.intern("c0");
      c = std::move(b);
      check(c.size() == 2 && c.get(1) == "new" && b.size() == 0, "move assignment");
    }
    std::cout << "string_arena: OK\n";
    return 0;
  }
  catch (const std::exception& e)
  {
    std::cerr << "string_arena: FAIL: " << e.what() << "\n";
    return 1;
  }
}
