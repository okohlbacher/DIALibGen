#include <odia/tune/Trainer.h>
#include <functional>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

int main()
{
  using ODIA::tune::TuneParams;
  using Case = std::pair<std::string, std::function<void(TuneParams&)>>;
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double inf = std::numeric_limits<double>::infinity();
  const std::vector<Case> cases = {
    {"filter:q_value", [](auto& p) { p.q_value = -0.1; }},
    {"filter:q_value", [](auto& p) { p.q_value = 1.1; }},
    {"filter:q_value", [=](auto& p) { p.q_value = nan; }},
    {"cohort:train_frac", [=](auto& p) { p.train_frac = inf; }},
    {"cohort:train_frac", [](auto& p) { p.train_frac = -0.1; }},
    {"filter:rt_spread_max", [](auto& p) { p.rt_spread_max = -1; }},
    {"filter:rt_max_minutes", [=](auto& p) { p.rt_max_minutes = nan; }},
    {"train:lr", [](auto& p) { p.lr = 0; }},
    {"train:lr", [=](auto& p) { p.lr = nan; }},
    {"train:lr", [=](auto& p) { p.lr = inf; }},
    {"train:epochs", [](auto& p) { p.epochs = 0; }},
    {"train:warmup", [](auto& p) { p.warmup = p.epochs + 1; }},
    {"batch size", [](auto& p) { p.batch_size = 0; }},
    {"evaluation interval", [](auto& p) { p.eval_every = 0; }},
    {"threads", [](auto& p) { p.threads = -1; }},
    {"stopping epochs", [](auto& p) { p.patience = -1; }},
    {"stop:rel_tol", [=](auto& p) { p.rel_tol = nan; }},
    {"stop:rel_tol", [](auto& p) { p.rel_tol = 1.1; }},
    {"stop:abs_tol", [](auto& p) { p.abs_tol = -1; }},
    {"stop:max_seconds", [=](auto& p) { p.max_seconds = inf; }},
    {"filter:min_charge", [](auto& p) { p.min_charge = 0; }},
    {"filter:min_charge", [](auto& p) { p.min_charge = 9; }},
    {"not both", [](auto& p) { p.train_size = 10; p.train_frac = 0.5; }},
  };
  for (const auto& [expected, change] : cases)
  {
    TuneParams p;
    p.model_in = "nonexistent-stock.onnx";
    p.model_out = "nonexistent-tuned.onnx";
    change(p);
    std::ostringstream log;
    try { (void)ODIA::tune::finetune(p, log); }
    catch (const std::exception& e)
    {
      if (std::string(e.what()).find(expected) != std::string::npos) { continue; }
      std::cerr << expected << ": wrong rejection: " << e.what() << '\n';
      return 1;
    }
    std::cerr << expected << ": invalid training controls were accepted\n";
    return 1;
  }
  std::cout << "PASS: invalid training controls refused before loading data or models\n";
}
