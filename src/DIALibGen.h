// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <odia/PeptDeepPredictor.h>
#include <nlohmann/json.hpp>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace ODIA::search { struct SearchParams; }

class DIALibGen final : public OpenMS::TOPPBase
{
public:
  DIALibGen();
  static std::string bundledModelDir();
  /// --help with -mode X on the command line: only X's options are shown
  /// (the others are registered as advanced for this invocation; --helphelp
  /// shows everything). Set by main() before the tool registers its options.
  static void setHelpMode(const std::string& mode) { help_mode_ = mode; }
  /// The options on @p argv that the tool does not know (TOPP's rule: a
  /// token is an option when it starts with '-' and a letter, or '--' and a
  /// letter), each with the closest known name or "" when nothing is close.
  static std::vector<std::pair<std::string, std::string>> unknownOptions(int argc, const char** argv);
  /// Option names removed in 0.13 (deprecated in 0.12) and what replaces
  /// each: a command line or an INI that still gives one is refused with a
  /// line naming the replacement (TOPP alone would reject an INI's unknown
  /// name without saying what it became).
  struct Removed { const char* old_name; const char* use; };
  static const std::vector<Removed>& removedOptions();
protected:
  void registerOptionsAndFlags_() override;
  void registerGenerationOptions_();
  void registerRefinementOptions_();
  ExitCodes main_(int argc, const char** argv) override;
  ExitCodes generate_();
  ExitCodes refine_(bool tune_only);
  /// The search: options, validated (throws std::invalid_argument).
  ODIA::search::SearchParams searchParams_();
  ExitCodes convert_();
  /// -mode prepare: compiles the CoreML shapes of the bundled MS2 model
  /// into the cache, for the peptide-length window.
  ExitCodes prepare_();
  /// GPU settings for prediction from -machine:*; warnings collects CPU fallbacks.
  ODIA::GpuOptions gpuOptions_(std::vector<std::string>* warnings);
  nlohmann::json generationValue_(const std::string& key, const nlohmann::json& value);
  /// A warning in the tool's log, as "WARNING: <text>" (the prefix that the
  /// library's own warnings carry, and that scripts look for).
  void warn_(const std::string& text) const { writeLogWarn_("WARNING: " + text); }
  /// Whether @p name is registered as advanced: everything outside the
  /// short list of the mode's main options (DIALibGen.cpp, kVisible).
  static bool advanced_(const std::string& name);
  /// Given on the command line, or in the INI with a value other than the
  /// registered default (a -write_ini file lists every option).
  bool given_(const std::string& name) const;
  /// The directory of the stock models: -models, else $DIALIBGEN_MODEL_DIR,
  /// else the bundled models.
  std::string modelDirectory_() const;
  std::set<std::string> supplied_;        ///< on the command line or in the INI
  std::set<std::string> argv_supplied_;   ///< on the command line
  std::set<std::string> refinement_options_;
  std::string mode_;                      ///< the mode in effect (-mode, or inferred)
  bool old_ini_ = false;                  ///< the INI predates 0.12 (its defaults differ)
  /// A training setting (tune_heads, tune_keep_free_cysteine_offset, filter:,
  /// cohort:, train:, stop:, machine:train_threads, machine:seed,
  /// machine:no_cudnn): the command line or INI, else the recipe's "tuning"
  /// object (-config), else the default.
  OpenMS::ParamValue trainingSetting_(const std::string& name) const;
  /// The training settings as a recipe's "tuning" object.
  nlohmann::json trainingRecipe_() const;
  nlohmann::json tuning_recipe_ = nlohmann::json::object();
  static std::string help_mode_;
};
