#include <odia/TextWriter.h>
#include <bit>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <locale>
#include <random>
#include <sstream>
#include <stdexcept>

namespace fs = std::filesystem;
void check(bool condition, const char* what) { if (!condition) { throw std::runtime_error(what); } }
std::string read(const fs::path& p) { std::ifstream in(p, std::ios::binary); return {std::istreambuf_iterator<char>(in), {}}; }

int main(int argc, char** argv)
{
  try
  {
    check(argc == 2, "work directory required");
    const fs::path dir(argv[1]); fs::create_directories(dir);
    const auto path = dir / "writer.txt";
    std::ostringstream expected; expected.imbue(std::locale::classic());
    ODIA::TextWriter actual(path.string());
    auto number = [&](double v, int precision) {
      actual.number(v, precision); actual.put('\n');
      expected << std::setprecision(precision) << std::defaultfloat << v << '\n';
    };
    std::mt19937_64 random(317);
    for (int precision : {6, 9, 10})
    {
      for (double v : {0.0, -0.0, 1.0, -1.0, 1e-30, 42949.67295,
                      std::numeric_limits<double>::min(), std::numeric_limits<double>::max(),
                      std::numeric_limits<double>::denorm_min(), std::numeric_limits<double>::infinity(),
                      -std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()})
      { number(v, precision); }
      for (int i = 0; i < 2000; ++i)
      {
        const double v = std::bit_cast<double>(random());
        if (std::isfinite(v)) { number(v, precision); }
      }
    }
    for (long long n : {0LL, -1LL, std::numeric_limits<long long>::min(), std::numeric_limits<long long>::max()})
    { actual.integer(n); actual.put('\n'); expected << n << '\n'; }
    for (std::size_t size : {std::size_t(0), std::size_t(1048575), std::size_t(1048576), std::size_t(1048593)})
    {
      const std::string value(size, static_cast<char>('a' + size % 26));
      actual.put(value); actual.put('|'); expected << value << '|';
    }
    actual.close(); actual.close();
    check(read(path) == expected.str(), "buffered text differs from independent ostream formatting");
    bool refused = false;
    try { ODIA::TextWriter bad((dir / "missing" / "output").string()); }
    catch (const std::runtime_error&) { refused = true; }
    check(refused, "failed open was not reported");
#ifndef _WIN32
    if (fs::exists("/dev/full"))
    {
      refused = false;
      try { ODIA::TextWriter full("/dev/full"); full.put("buffered failure"); full.close(); }
      catch (const std::runtime_error&) { refused = true; }
      check(refused, "buffered close failure was not reported");
      refused = false;
      try { ODIA::TextWriter full("/dev/full"); full.put(std::string(1048577, 'x')); full.close(); }
      catch (const std::runtime_error&) { refused = true; }
      check(refused, "direct write failure was not reported");
    }
#endif
    fs::remove(path);
    return 0;
  }
  catch (const std::exception& e) { std::fprintf(stderr, "%s\n", e.what()); return 1; }
}
