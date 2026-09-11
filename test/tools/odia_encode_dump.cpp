// Copyright (c) 2026, Oliver Kohlbacher and the DIALibraryGenerator authors.
// SPDX-License-Identifier: BSD-3-Clause
//
// Dumps ODIA's PeptDeep encoding as JSON so it can be diffed against
// test/peptdeep_reference.py, which implements the same spec independently.

#include <odia/PeptDeepEncoder.h>

#include <OpenMS/CHEMISTRY/AASequence.h>

#include <iostream>
#include <vector>

int main(int argc, char** argv)
{
  if (argc < 2)
  {
    std::cerr << "usage: odia_encode_dump <modified-sequence>...\n"
                 "Several sequences of equal length are encoded as one batch, so\n"
                 "every row can be checked -- not just row 0.\n";
    return 2;
  }
  try
  {
    std::vector<OpenMS::AASequence> peptides;
    for (int i = 1; i < argc; ++i)
    {
      peptides.push_back(OpenMS::AASequence::fromString(argv[i]));
    }
    const auto batch = ODIA::PeptDeepEncoder::encode(peptides);
    const auto width = ODIA::PEPTDEEP_MOD_ELEMENTS.size();
    const auto length = batch.sequence_length;

    std::cout << "{\"rows\": " << batch.rows
              << ", \"sequence_length\": " << length << ", \"peptides\": [";
    for (std::size_t row = 0; row < batch.rows; ++row)
    {
      std::cout << (row ? ", " : "") << "{\"aa_indices\": [";
      for (std::size_t i = 0; i < length; ++i)
      {
        std::cout << (i ? ", " : "") << batch.aa_indices[row * length + i];
      }
      std::cout << "], \"mod_x\": {";
      bool first_row = true;
      for (std::size_t r = 0; r < length; ++r)
      {
        std::string entries;
        for (std::size_t c = 0; c < width; ++c)
        {
          const float v = batch.mod_x[(row * length + r) * width + c];
          if (v != 0.0f)
          {
            if (!entries.empty()) { entries += ", "; }
            entries += "\"" + std::to_string(c) + "\": " + std::to_string(v);
          }
        }
        if (!entries.empty())
        {
          std::cout << (first_row ? "" : ", ") << "\"" << r << "\": {" << entries << "}";
          first_row = false;
        }
      }
      std::cout << "}}";
    }
    std::cout << "]}\n";
  }
  catch (const std::exception& e)
  {
    std::cerr << "error: " << e.what() << "\n";
    return 1;
  }
  return 0;
}
