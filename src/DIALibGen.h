// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <nlohmann/json.hpp>
#include <set>
#include <string>

class DIALibGen final : public OpenMS::TOPPBase
{
public:
  DIALibGen();
  static std::string bundledModelDir();
protected:
  void registerOptionsAndFlags_() override;
  void registerGenerationOptions_();
  void registerRefinementOptions_();
  ExitCodes main_(int argc, const char** argv) override;
  ExitCodes generate_();
  ExitCodes refine_(bool tune_only);
  nlohmann::json generationValue_(const std::string& key, const nlohmann::json& value);
  std::set<std::string> supplied_;
};
