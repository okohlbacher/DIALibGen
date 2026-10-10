// ortbench: throughput / memory harness for the PeptDeep ONNX models, with the
// same session settings as DIALibGen (intra-op 1, ORT_ENABLE_ALL, CPU arena on,
// one session per worker thread). perf/memory-cache, bid10.
//
//   ortbench --model M --kind ms2|ccs|rt --threads N [--rows R] [--len nAA|all]
//            [--seconds S] [--share 0|1] [--cpus a,b,c,...] [--label L]
//
// Each worker owns one session (or all share one PrepackedWeightsContainer with
// --share 1), builds one synthetic input of R rows x nAA residues, warms up once,
// waits on a barrier, then calls Run until the deadline. --cpus pins worker i to
// cpus[i]. --len all cycles through nAA 7..30 (arena high-water mode).
// Output: one TSV line per worker plus one "total" line (header with --header).
#include <onnxruntime_cxx_api.h>

#include <pthread.h>
#include <sched.h>
#include <time.h>
#include <unistd.h>

#include <atomic>
#include <barrier>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <memory>
#include <random>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

static double now() { timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + 1e-9 * t.tv_nsec; }
static double tcpu() { timespec t; clock_gettime(CLOCK_THREAD_CPUTIME_ID, &t); return t.tv_sec + 1e-9 * t.tv_nsec; }
static long rss_kb()
{
  std::ifstream f("/proc/self/statm"); long a = 0, b = 0; f >> a >> b; return b * (sysconf(_SC_PAGESIZE) / 1024);
}

struct Input
{
  std::int64_t rows, L;  // L = nAA + 2
  std::vector<std::int64_t> aa, inst;
  std::vector<float> mod, charge, nce;
};

static Input makeInput(std::int64_t rows, std::int64_t naa, unsigned seed)
{
  Input in; in.rows = rows; in.L = naa + 2;
  std::mt19937 rng(seed);
  in.aa.assign(rows * in.L, 0);
  in.mod.assign(rows * in.L * 109, 0.0f);
  for (std::int64_t r = 0; r < rows; ++r)
    for (std::int64_t k = 0; k < naa; ++k) in.aa[r * in.L + k + 1] = 1 + rng() % 20;
  in.charge.assign(rows, 0.2f + 0.1f * (rng() % 3));
  in.nce.assign(rows, 0.30f);
  in.inst.assign(rows, 0);
  return in;
}

int main(int argc, char** argv)
{
  std::string model, kind = "ms2", label = "-", lenarg = "16", cpus;
  int threads = 4, rows = 64, share = 0, header = 0, shared_session = 0; double seconds = 20, start_at = 0;
  for (int i = 1; i < argc; ++i)
  {
    std::string a = argv[i]; auto v = [&] { return std::string(argv[++i]); };
    if (a == "--model") model = v(); else if (a == "--kind") kind = v();
    else if (a == "--threads") threads = std::stoi(v()); else if (a == "--rows") rows = std::stoi(v());
    else if (a == "--len") lenarg = v(); else if (a == "--seconds") seconds = std::stod(v());
    else if (a == "--share") share = std::stoi(v()); else if (a == "--cpus") cpus = v();
    else if (a == "--label") label = v(); else if (a == "--header") header = 1;
    else if (a == "--shared-session") shared_session = std::stoi(v()); else if (a == "--start-at") start_at = std::stod(v());
  }
  std::vector<int> cpu;
  { std::stringstream ss(cpus); std::string t; while (std::getline(ss, t, ',')) if (!t.empty()) cpu.push_back(std::stoi(t)); }
  const bool all = lenarg == "all";
  const int naa0 = all ? 7 : std::stoi(lenarg);

  Ort::Env env{ORT_LOGGING_LEVEL_WARNING, "ortbench"};
  Ort::SessionOptions opt;
  opt.SetIntraOpNumThreads(1);
  opt.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_ALL);
  std::unique_ptr<Ort::PrepackedWeightsContainer> pp;
  if (share) pp = std::make_unique<Ort::PrepackedWeightsContainer>();
  Ort::MemoryInfo mem = Ort::MemoryInfo::CreateCpu(OrtArenaAllocator, OrtMemTypeDefault);

  const long rss0 = rss_kb();
  const double tc0 = now();
  std::vector<std::unique_ptr<Ort::Session>> sess;
  for (int t = 0; t < (shared_session ? 1 : threads); ++t)
    sess.push_back(pp ? std::make_unique<Ort::Session>(env, model.c_str(), opt, *pp)
                      : std::make_unique<Ort::Session>(env, model.c_str(), opt));
  const double create_s = now() - tc0;
  const long rss_sessions = rss_kb();

  std::vector<std::string> in_names_s, out_names_s;
  Ort::AllocatorWithDefaultOptions alloc;
  for (size_t i = 0; i < sess[0]->GetInputCount(); ++i) in_names_s.emplace_back(sess[0]->GetInputNameAllocated(i, alloc).get());
  out_names_s.emplace_back(sess[0]->GetOutputNameAllocated(0, alloc).get());
  std::vector<const char*> in_names, out_names{out_names_s[0].c_str()};
  for (auto& n : in_names_s) in_names.push_back(n.c_str());
  const size_t nin = in_names.size();

  auto run = [&](Ort::Session& s, Input& in) {
    std::vector<Ort::Value> v;
    std::int64_t aa_shape[2] = {in.rows, in.L}, mod_shape[3] = {in.rows, in.L, 109}, meta[2] = {in.rows, 1}, ins[1] = {in.rows};
    v.push_back(Ort::Value::CreateTensor<std::int64_t>(mem, in.aa.data(), in.aa.size(), aa_shape, 2));
    v.push_back(Ort::Value::CreateTensor<float>(mem, in.mod.data(), in.mod.size(), mod_shape, 3));
    if (nin >= 3) v.push_back(Ort::Value::CreateTensor<float>(mem, in.charge.data(), in.charge.size(), meta, 2));
    if (nin >= 5)
    {
      v.push_back(Ort::Value::CreateTensor<float>(mem, in.nce.data(), in.nce.size(), meta, 2));
      v.push_back(Ort::Value::CreateTensor<std::int64_t>(mem, in.inst.data(), in.inst.size(), ins, 1));
    }
    auto r = s.Run(Ort::RunOptions{nullptr}, in_names.data(), v.data(), v.size(), out_names.data(), 1);
    return r.front().GetTensorTypeAndShapeInfo().GetElementCount();
  };

  std::vector<double> rate(threads), cpu_s(threads), wall_s(threads); std::vector<long> done(threads), calls(threads);
  std::vector<int> ran_on(threads, -1);
  std::barrier sync(threads + 1);
  std::atomic<double> deadline{0};
  long rss_warm = 0;
  std::vector<std::thread> pool;
  for (int t = 0; t < threads; ++t)
    pool.emplace_back([&, t] {
      if (!cpu.empty())
      { cpu_set_t s; CPU_ZERO(&s); CPU_SET(cpu[t % cpu.size()], &s); pthread_setaffinity_np(pthread_self(), sizeof(s), &s); }
      std::vector<Input> ins;
      if (all) { for (int n = 7; n <= 30; ++n) ins.push_back(makeInput(rows, n, 1234 + 97 * t + n)); }
      else ins.push_back(makeInput(rows, naa0, 1234 + 97 * t));
      Ort::Session& my = *sess[shared_session ? 0 : t];
      for (auto& in : ins) run(my, in);  // warm-up (all lengths in --len all)
      sync.arrive_and_wait();                  // main reads RSS after warm-up
      sync.arrive_and_wait();                  // start
      const double end = deadline.load();
      const double w0 = now(), c0 = tcpu();
      long n = 0, k = 0; size_t li = 0;
      while (now() < end) { run(my, ins[li]); n += rows; ++k; li = (li + 1) % ins.size(); }
      wall_s[t] = now() - w0; cpu_s[t] = tcpu() - c0; done[t] = n; calls[t] = k; rate[t] = n / wall_s[t];
      ran_on[t] = sched_getcpu();
    });
  sync.arrive_and_wait();
  rss_warm = rss_kb();
  if (start_at > 0)
  {  // common wall-clock start for concurrent instances (CLOCK_REALTIME epoch seconds)
    timespec rt; clock_gettime(CLOCK_REALTIME, &rt);
    const double wait = start_at - (rt.tv_sec + 1e-9 * rt.tv_nsec);
    if (wait > 0) std::this_thread::sleep_for(std::chrono::duration<double>(wait));
  }
  deadline = now() + seconds;
  sync.arrive_and_wait();
  for (auto& p : pool) p.join();
  const long rss_end = rss_kb();

  if (header)
    std::printf("label\tkind\trows\tlen\tthreads\tshare\tshared_session\tworker\tcpu\trows_per_s\tcalls\twall_s\tcpu_s\tcreate_s\trss0_kb\trss_sessions_kb\trss_warm_kb\trss_end_kb\n");
  double tot = 0, cs = 0;
  for (int t = 0; t < threads; ++t)
  {
    tot += rate[t]; cs += cpu_s[t];
    std::printf("%s\t%s\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%.2f\t%ld\t%.3f\t%.3f\t%.3f\t%ld\t%ld\t%ld\t%ld\n", label.c_str(), kind.c_str(), rows,
                lenarg.c_str(), threads, share, shared_session, t, ran_on[t], rate[t], calls[t], wall_s[t], cpu_s[t], create_s, rss0, rss_sessions, rss_warm, rss_end);
  }
  std::printf("%s\t%s\t%d\t%s\t%d\t%d\t%d\ttotal\t-1\t%.2f\t-1\t%.3f\t%.3f\t%.3f\t%ld\t%ld\t%ld\t%ld\n", label.c_str(), kind.c_str(), rows,
              lenarg.c_str(), threads, share, shared_session, tot, seconds, cs, create_s, rss0, rss_sessions, rss_warm, rss_end);
  return 0;
}
