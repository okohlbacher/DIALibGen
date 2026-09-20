#include <odia/AtomicFile.h>
#include <fstream>
#include <iostream>
#include <iterator>

namespace fs = std::filesystem;
void check(bool ok) { if (!ok) { throw std::runtime_error("atomic file contract failed"); } }
std::string read(const fs::path& p)
{ std::ifstream in(p); return {std::istreambuf_iterator<char>(in), {}}; }
int main(int argc, char** argv)
{
  try
  {
    check(argc == 2);
    const fs::path work(argv[1]); fs::create_directories(work);
    const auto target = work / "result.txt";
    { std::ofstream out(target); out << "original"; }
    fs::path temporary;
    try
    {
      ODIA::AtomicFile first(target), second(target);
      temporary = first.temporaryPath();
      check(temporary != second.temporaryPath());
      { std::ofstream out(temporary); out << "partial"; }
      throw std::runtime_error("cancelled");
    }
    catch (const std::runtime_error& e) { check(std::string(e.what()) == "cancelled"); }
    check(read(target) == "original" && !fs::exists(temporary.parent_path()));
    {
      ODIA::AtomicFile output(target);
      { std::ofstream out(output.temporaryPath()); out << "complete"; }
      output.commit();
    }
    check(read(target) == "complete");
#ifndef _WIN32
    const auto previous_umask = ::umask(0);
    {
      ODIA::AtomicFile output(target);
      ::umask(previous_umask);
      check(fs::status(output.temporaryPath().parent_path()).permissions() == fs::perms::owner_all);
    }
    const auto link = work / "result-link.txt";
    fs::remove(link); fs::create_symlink(target.filename(), link);
    {
      ODIA::AtomicFile output(link);
      { std::ofstream out(output.temporaryPath()); out << "new link entry"; }
      output.commit();
    }
    check(!fs::is_symlink(link) && read(link) == "new link entry" && read(target) == "complete");
#endif
    const auto directory = work / "existing-directory"; fs::create_directories(directory);
    { std::ofstream out(directory / "keep"); out << "untouched"; }
    bool refused = false;
    try
    {
      ODIA::AtomicFile output(directory); temporary = output.temporaryPath();
      { std::ofstream out(temporary); out << "replacement"; }
      output.commit();
    }
    catch (const std::runtime_error&) { refused = true; }
    check(refused && read(directory / "keep") == "untouched" && !fs::exists(temporary.parent_path()));
    for (const auto& entry : fs::directory_iterator(work))
    { check(!entry.path().filename().string().starts_with(".dialibgen-tmp-")); }
    std::cout << "atomic output preserves existing files on cancellation and commit failure\n";
  }
  catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
