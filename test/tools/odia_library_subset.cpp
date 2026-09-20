// Library::subsetByIndex must OUTLIVE the library it came from.
//
// `StringArena::Entry` holds a raw `const char*` into the arena's own blocks.
// The arena used to be copied wholesale into the subset, which allocated fresh
// blocks and then copied the entries verbatim -- so every handle in the subset
// resolved into the PARENT's storage. That is correct for exactly as long as
// the parent is alive, and the only caller at the time held it alive, so the
// suite never saw it.
//
// The `-min_library_fragments` filter is written `library = library.subsetByIndex(keep)`.
// The move-assign destroys the parent, the subset's handles dangle, and the
// next `strings().get()` segfaults inside memcpy -- which is what it did, on a
// 9.6M-precursor run, 52 s of library load in, after the filter reported the
// correct 366,084 targets dropped.
//
// So the contract under test is a lifetime, not a value: build a library,
// subset it, DESTROY the parent, and only then read every string back.

#include <odia/Library.h>

#include <cstdlib>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

int main()
{
  using namespace ODIA;

  const std::size_t limit = std::numeric_limits<std::uint32_t>::max();
  Library::checkTransitionCapacity(limit - 1, 1); // The boundary itself remains representable.
  bool refused = false;
  Library oversized;
  try { oversized.reserve(1, limit + std::size_t{1}); }
  catch (const std::length_error& e)
  { refused = std::string(e.what()).find("32-bit transition limit") != std::string::npos; }
  if (!refused || oversized.precursors().mz.capacity() || oversized.transitions().product_mz.capacity())
  { std::cerr << "oversized transition reserve was not refused before allocation\n"; return EXIT_FAILURE; }
  refused = false;
  try { Library::checkTransitionCapacity(limit, std::numeric_limits<std::size_t>::max()); }
  catch (const std::length_error&) { refused = true; }
  if (!refused) { std::cerr << "overflowing append capacity was accepted\n"; return EXIT_FAILURE; }

  // Enough distinct strings to span more than one arena block, so the test
  // exercises the block-walking path rather than a single-block special case.
  constexpr std::size_t N = 4096;

  std::vector<std::string> expect_seq(N), expect_pg(N);
  std::vector<std::size_t> keep;

  auto parent = std::make_unique<Library>();
  {
    auto& p = parent->precursors();
    auto& t = parent->transitions();
    for (std::size_t i = 0; i < N; ++i)
    {
      // Long enough that 4096 of them cannot share one 1 MiB block.
      expect_seq[i] = "PEPTIDESEQUENCE" + std::to_string(i) + std::string(200, 'A');
      expect_pg[i] = "sp|PROTEIN" + std::to_string(i) + "|" + std::string(200, 'B');
      p.mz.push_back(400.0 + static_cast<double>(i));
      p.irt.push_back(static_cast<float>(i));
      p.im.push_back(1.0f);
      p.ccs.push_back(300.0f);
      p.charge.push_back(2);
      p.decoy.push_back(i % 2 ? 1 : 0);
      p.modified_sequence.push_back(parent->strings().intern(expect_seq[i]));
      p.protein_group.push_back(parent->strings().intern(expect_pg[i]));
      p.transition_begin.push_back(static_cast<std::uint32_t>(t.product_mz.size()));
      p.transition_count.push_back(3);
      for (int k = 0; k < 3; ++k)
      {
        t.product_mz.push_back(200.0 + k);
        t.library_intensity.push_back(1.0f);
        t.type.push_back(FragmentType::Y);
        t.ordinal.push_back(static_cast<std::uint8_t>(k + 1));
        t.charge.push_back(1);
        t.loss.push_back(LossType::None);
      }
      if (i % 3 == 0) { keep.push_back(i); }
    }
  }

  Library subset = parent->subsetByIndex(keep);

  // THE POINT OF THE TEST. Everything below reads through the subset alone.
  parent.reset();

  if (subset.precursorCount() != keep.size())
  {
    std::cerr << "subset has " << subset.precursorCount() << " precursors, expected "
              << keep.size() << "\n";
    return EXIT_FAILURE;
  }

  const auto& p = subset.precursors();
  for (std::size_t j = 0; j < keep.size(); ++j)
  {
    const std::size_t i = keep[j];
    const std::string seq(subset.strings().get(p.modified_sequence[j]));
    const std::string pg(subset.strings().get(p.protein_group[j]));
    if (seq != expect_seq[i] || pg != expect_pg[i])
    {
      std::cerr << "precursor " << j << " (source " << i << ") resolved to\n"
                << "  sequence " << seq.substr(0, 40) << "...\n"
                << "  expected " << expect_seq[i].substr(0, 40) << "...\n";
      return EXIT_FAILURE;
    }
    if (p.transition_count[j] != 3)
    {
      std::cerr << "precursor " << j << " kept " << p.transition_count[j]
                << " transitions, expected 3\n";
      return EXIT_FAILURE;
    }
  }

  // The subset must carry only what it uses: re-interning is what makes the
  // filter a memory win rather than a memory copy.
  if (subset.strings().size() != 2 * keep.size())
  {
    std::cerr << "subset arena holds " << subset.strings().size()
              << " strings, expected " << 2 * keep.size() << "\n";
    return EXIT_FAILURE;
  }

  std::cout << "subset of " << keep.size() << " precursors outlived its parent; "
            << subset.strings().size() << " strings re-interned\n";
  return EXIT_SUCCESS;
}
