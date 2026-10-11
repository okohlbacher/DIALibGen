// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

#include "ThreadPool.h"

#include <odia/Trace.h>

#include <algorithm>
#include <condition_variable>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <memory>
#include <mutex>
#include <string>
#include <string_view>
#include <system_error>
#include <thread>
#include <vector>

#if defined(__linux__)
#include <pthread.h>
#include <sched.h>
#endif

namespace ODIA::parallel
{
  namespace
  {
#if defined(__linux__)
    /// CPUs of every NUMA node that this process may run on (its affinity
    /// mask, e.g. under taskset), nodes without such a CPU left out.
    std::vector<cpu_set_t> nodeCpuSets()
    {
      std::vector<cpu_set_t> sets;
      cpu_set_t allowed;
      CPU_ZERO(&allowed);
      if (sched_getaffinity(0, sizeof allowed, &allowed) != 0) { return sets; }
      for (int node = 0; node < 1024; ++node)
      {
        std::ifstream in("/sys/devices/system/node/node" + std::to_string(node) + "/cpulist");
        if (!in) { if (node > 64) { break; } continue; }
        std::string list;
        std::getline(in, list);
        cpu_set_t set;
        CPU_ZERO(&set);
        std::size_t pos = 0;
        while (pos < list.size())
        {
          const std::size_t comma = std::min(list.find(',', pos), list.size());
          const std::string range = list.substr(pos, comma - pos);
          const std::size_t dash = range.find('-');
          try
          {
            const int a = std::stoi(range.substr(0, dash));
            const int b = dash == std::string::npos ? a : std::stoi(range.substr(dash + 1));
            for (int c = a; c <= b && c < CPU_SETSIZE; ++c) { if (CPU_ISSET(c, &allowed)) { CPU_SET(c, &set); } }
          }
          catch (const std::exception&) {}
          pos = comma + 1;
        }
        if (CPU_COUNT(&set) > 0) { sets.push_back(set); }
      }
      return sets;
    }
#endif

    struct Worker
    {
      std::size_t id = 0;
      std::mutex mutex;
      std::condition_variable wake;
      std::function<void()> task;   ///< set by the pool, cleared by the worker
      bool busy = false;            ///< claimed; guarded by Pool::mutex
      std::thread thread;
    };

    class Pool
    {
    public:
      static Pool& instance()
      {
        // Never destroyed: workers sleep until the process ends.
        static Pool* pool = new Pool;
        return *pool;
      }

      /// Up to @p want free workers, lowest or highest ids first, starting
      /// new ones as needed; fewer if no thread can be started.
      std::vector<Worker*> claim(std::size_t want, bool low)
      {
        std::vector<Worker*> got;
        if (want == 0) { return got; }
        const std::lock_guard<std::mutex> lock(mutex_);
        const std::size_t size = workers_.size();
        for (std::size_t k = 0; k < size && got.size() < want; ++k)
        {
          Worker* w = workers_[low ? k : size - 1 - k].get();
          if (!w->busy) { w->busy = true; got.push_back(w); }
        }
        while (got.size() < want && !exhausted_)
        {
          auto w = std::make_unique<Worker>();
          w->id = workers_.size();
          w->busy = true;
          Worker* raw = w.get();
          try { w->thread = std::thread([this, raw] { loop(raw); }); }
          catch (const std::system_error&) { exhausted_ = true; break; }
          workers_.push_back(std::move(w));
          got.push_back(raw);
        }
        return got;
      }

      void assign(Worker* w, std::function<void()> task)
      {
        {
          const std::lock_guard<std::mutex> lock(w->mutex);
          w->task = std::move(task);
        }
        w->wake.notify_one();
      }

      /// Back to the free workers: before the caller is told the work is
      /// done, so that its next team finds this worker free.
      void release(Worker* w)
      {
        const std::lock_guard<std::mutex> lock(mutex_);
        w->busy = false;
      }

      std::size_t size()
      {
        const std::lock_guard<std::mutex> lock(mutex_);
        return workers_.size();
      }

    private:
      Pool()
      {
#if defined(__linux__)
        const char* pin = std::getenv("ODIA_POOL_PIN");
        if (pin == nullptr || std::string_view(pin) != "0") { nodes_ = nodeCpuSets(); }
        if (nodes_.size() < 2) { nodes_.clear(); }   // one node: nothing to place
#endif
      }

      void loop(Worker* w)
      {
#if defined(__linux__)
        if (!nodes_.empty())
        { pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &nodes_[w->id % nodes_.size()]); }
#endif
        std::unique_lock<std::mutex> lock(w->mutex);
        for (;;)
        {
          w->wake.wait(lock, [&] { return static_cast<bool>(w->task); });
          std::function<void()> task = std::move(w->task);
          w->task = nullptr;
          lock.unlock();
          task();   // never throws: team() and async() catch inside
          trace::flushThread();
          lock.lock();
        }
      }

      std::mutex mutex_;
      std::vector<std::unique_ptr<Worker>> workers_;
      bool exhausted_ = false;   ///< a thread could not be started; no more tries
#if defined(__linux__)
      std::vector<cpu_set_t> nodes_;
#endif
    };

  } // namespace

  std::size_t team(std::size_t n, const std::function<void(std::size_t)>& fn, bool low)
  {
    if (n <= 1) { if (n == 1) { fn(0); } return 1; }
    Pool& pool = Pool::instance();
    const std::vector<Worker*> workers = pool.claim(n - 1, low);
    std::mutex mutex;
    std::condition_variable done;
    std::size_t remaining = workers.size();
    std::exception_ptr error;
    for (std::size_t k = 0; k < workers.size(); ++k)
    {
      Worker* w = workers[k];
      pool.assign(w, [&, w, k] {
        std::exception_ptr e;
        try { fn(k + 1); }
        catch (...) { e = std::current_exception(); }
        pool.release(w);
        const std::lock_guard<std::mutex> lock(mutex);
        if (e && !error) { error = e; }
        if (--remaining == 0) { done.notify_one(); }
      });
    }
    try { fn(0); }
    catch (...)
    {
      const std::lock_guard<std::mutex> lock(mutex);
      if (!error) { error = std::current_exception(); }
    }
    std::unique_lock<std::mutex> lock(mutex);
    done.wait(lock, [&] { return remaining == 0; });
    if (error) { std::rethrow_exception(error); }
    return workers.size() + 1;
  }

  std::future<void> async(std::function<void()> task)
  {
    auto promise = std::make_shared<std::promise<void>>();
    std::future<void> result = promise->get_future();
    Pool& pool = Pool::instance();
    const std::vector<Worker*> workers = pool.claim(1, false);
    if (workers.empty())
    {
      try { task(); promise->set_value(); }
      catch (...) { promise->set_exception(std::current_exception()); }
      return result;
    }
    Worker* w = workers.front();
    pool.assign(w, [&pool, w, promise, task = std::move(task)] {
      std::exception_ptr e;
      try { task(); }
      catch (...) { e = std::current_exception(); }
      pool.release(w);
      if (e) { promise->set_exception(e); } else { promise->set_value(); }
    });
    return result;
  }

  std::size_t poolSize() { return Pool::instance().size(); }
} // namespace ODIA::parallel
