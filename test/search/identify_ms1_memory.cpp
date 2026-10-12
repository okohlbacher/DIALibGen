// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

// search:ms1_in_memory: the decision rule (auto holds the MS1 map when the
// free memory is at least twice its estimated size; on and off are obeyed),
// that the free memory can be read on the platforms the tool ships for, and,
// with glibc, that memory freed but kept by the allocator is given back.
// That on and off write the same report is identify_e2e's.

#include <odia/search/MemoryCheck.h>

#if defined(__linux__)
#  include <malloc.h>
#endif

#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <vector>
#include <stdexcept>
#include <string>

namespace
{
  int failures = 0;

  void check(bool ok, const std::string& what)
  {
    if (!ok) { ++failures; std::cerr << "FAIL: " << what << '\n'; }
    else { std::cout << "ok   " << what << '\n'; }
  }
}

int main()
{
  using namespace ODIA::search;
  constexpr std::uint64_t GB = 1000000000ull;

  check(decideMs1InMemory(Ms1InMemory::Auto, 10 * GB, 20 * GB).hold, "auto holds at exactly twice the estimate");
  check(decideMs1InMemory(Ms1InMemory::Auto, 10 * GB, 64 * GB).hold, "auto holds with more than twice the estimate free");
  check(!decideMs1InMemory(Ms1InMemory::Auto, 10 * GB, 20 * GB - 1).hold, "auto streams just below twice the estimate");
  check(!decideMs1InMemory(Ms1InMemory::Auto, 10 * GB, 0).hold, "auto streams when the free memory is unknown");
  check(decideMs1InMemory(Ms1InMemory::On, 10 * GB, 1 * GB).hold, "on holds whatever the free memory");
  check(decideMs1InMemory(Ms1InMemory::On, 10 * GB, 0).hold, "on holds when the free memory is unknown");
  check(!decideMs1InMemory(Ms1InMemory::Off, 1 * GB, 1000 * GB).hold, "off never holds");
  for (const Ms1InMemory m : {Ms1InMemory::Auto, Ms1InMemory::On, Ms1InMemory::Off})
  { check(!decideMs1InMemory(m, 0, 100 * GB).hold, std::string("nothing to hold without MS1 spectra (") + toString(m) + ")"); }

  const Ms1MemoryDecision d = decideMs1InMemory(Ms1InMemory::Auto, 12300000000ull, 20 * GB);
  check(d.reason.find("12.3 GB") != std::string::npos && d.reason.find("20.0 GB free") != std::string::npos &&
          d.reason.find("24.6 GB") != std::string::npos,
        "the reason names the estimate, the free memory and what auto needed: " + d.reason);

  check(parseMs1InMemory("auto") == Ms1InMemory::Auto && parseMs1InMemory("on") == Ms1InMemory::On &&
          parseMs1InMemory("off") == Ms1InMemory::Off,
        "auto, on and off parse");
  bool refused = false;
  try { parseMs1InMemory("yes"); }
  catch (const std::invalid_argument&) { refused = true; }
  check(refused, "anything else is refused");

  const std::uint64_t free_bytes = availableMemoryBytes();
  std::cout << "     free memory: " << static_cast<double>(free_bytes) / 1e9 << " GB\n";
#if defined(__linux__) || defined(__APPLE__) || defined(_WIN32)
  check(free_bytes > 0, "the free memory is known on this platform");
#endif

#if defined(__GLIBC__)
  {
    // What extraction leaves behind: many large blocks from the heap (glibc's
    // dynamic mmap threshold puts them there once a block that size was
    // freed; fixed here), freed below a block still in use, so free() cannot
    // shrink the heap and the pages stay resident.
    mallopt(M_MMAP_THRESHOLD, 64 << 20);
    mallopt(M_TRIM_THRESHOLD, 1 << 30);
    constexpr std::size_t blocks = 256, block = 1 << 20;
    std::vector<char*> held;
    for (std::size_t k = 0; k < blocks; ++k)
    {
      held.push_back(static_cast<char*>(std::malloc(block)));
      std::memset(held.back(), 1, block);
    }
    char* pin = static_cast<char*>(std::malloc(4096));
    std::memset(pin, 1, 4096);
    for (char* p : held) { std::free(p); }
    const std::uint64_t released = releaseFreedMemory();
    std::free(pin);
    check(released >= blocks * block / 2, "freed heap memory kept by glibc is given back: " + memorySize(released) + " of " +
                                            memorySize(blocks * block));
  }
#else
  check(releaseFreedMemory() == 0, "nothing to give back without glibc");
#endif

  if (failures) { std::cerr << failures << " check(s) failed\n"; return 1; }
  return 0;
}
