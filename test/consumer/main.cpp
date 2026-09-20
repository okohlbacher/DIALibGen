#include <odia/Library.h>
#ifdef HAS_REFINE
#include <odia/LibraryRefiner.h>
#endif
#ifdef HAS_TUNE
#include <odia/tune/Trainer.h>
#endif
#include <iostream>
#include <string>

int main()
{
  ODIA::Library library;
  if (library.decoyCount() != 0) return 1;
#ifdef HAS_REFINE
  if (ODIA::canonicalModifiedSequence("AC(UniMod:4)DEFGK") != "AC(UniMod:4)DEFGK") return 2;
#endif
#ifdef HAS_TUNE
  if (std::string(ODIA::tune::headName(ODIA::tune::HeadKind::RT)) != "rt") return 3;
#endif
  std::cout << "Installed DIALibGen libraries linked and ran\n";
}
