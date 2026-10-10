// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <odia/Trace.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <filesystem>
#include <functional>
#include <limits>
#include <map>
#include <unordered_map>
#include <mutex>
#include <string>
#include <string_view>
#include <thread>
#include <utility>
#include <vector>

#if defined(__linux__) || defined(__APPLE__)
#include <pthread.h>
#include <sys/resource.h>
#include <time.h>
#include <unistd.h>
#endif
#if defined(__linux__)
#include <dirent.h>
#include <dlfcn.h>
#include <fcntl.h>
#include <sched.h>
#include <sys/syscall.h>
#include <sys/utsname.h>
#endif
#if defined(__GLIBC__)
#include <gnu/libc-version.h>
#include <malloc.h>
#endif

namespace ODIA::trace
{
  namespace
  {
    using Clock = std::chrono::steady_clock;

    struct Event
    {
      double t0, t1, cpu;
      std::int64_t tid;
      int level;
      const char* phase;
      const char* name;
      Fields fields;
    };

    /// Process-wide state. Never destroyed: worker threads hand over their
    /// events from thread_local destructors, whose order relative to static
    /// destruction is not something to depend on.
    struct State
    {
      Clock::time_point epoch = Clock::now();
      std::string dir;
      std::atomic<const char*> phase{"-"};
      std::mutex mutex;                     // guards events and notes
      std::vector<Event> events;
      std::vector<std::pair<std::string, std::string>> notes;
      bool session_active = false;
    };

    State& state()
    {
      static State* s = new State;
      return *s;
    }

    /// The calling thread's events, handed over when the thread exits.
    struct Buffer
    {
      std::vector<Event> events;
      Buffer() { events.reserve(1024); }
      ~Buffer() { flush(); }
      void flush()
      {
        if (events.empty()) { return; }
        auto& s = state();
        const std::lock_guard<std::mutex> lock(s.mutex);
        s.events.insert(s.events.end(), events.begin(), events.end());
        events.clear();
      }
    };

    Buffer& buffer()
    {
      thread_local Buffer b;
      return b;
    }

    std::int64_t currentTid() noexcept
    {
      thread_local std::int64_t tid = [] {
#if defined(__linux__)
        return static_cast<std::int64_t>(::syscall(SYS_gettid));
#elif defined(__APPLE__)
        std::uint64_t id = 0;
        pthread_threadid_np(nullptr, &id);
        return static_cast<std::int64_t>(id);
#else
        return static_cast<std::int64_t>(std::hash<std::thread::id>{}(std::this_thread::get_id()) & 0x7fffffff);
#endif
      }();
      return tid;
    }

    void push(const Event& e) noexcept
    {
      try { buffer().events.push_back(e); }
      catch (...) {}   // tracing must never take the run down
    }

    [[maybe_unused]] long envLong(const char* name, long fallback)
    {
      const char* v = std::getenv(name);
      if (v == nullptr || *v == '\0') { return fallback; }
      char* end = nullptr;
      const long x = std::strtol(v, &end, 10);
      return (end != v && x >= 0) ? x : fallback;
    }

    std::string sanitize(std::string s)
    {
      for (char& c : s) { if (c == '\t' || c == '\n' || c == '\r') { c = ' '; } }
      return s;
    }

    FILE* openOut(const std::string& name)
    {
      const auto path = (std::filesystem::path(state().dir) / name).string();
      FILE* f = std::fopen(path.c_str(), "w");
      if (f != nullptr) { std::setvbuf(f, nullptr, _IOFBF, 1 << 20); }
      else { std::fprintf(stderr, "[trace] cannot write %s\n", path.c_str()); }
      return f;
    }

#if defined(__linux__)
    // ---------------------------------------------------------------- /proc

    /// Whole file into @p out with plain read(2): no stdio, no locale.
    bool slurp(const char* path, std::string& out)
    {
      out.clear();
      const int fd = ::open(path, O_RDONLY | O_CLOEXEC);
      if (fd < 0) { return false; }
      char buf[65536];
      for (;;)
      {
        const ssize_t got = ::read(fd, buf, sizeof buf);
        if (got < 0) { ::close(fd); return false; }
        if (got == 0) { break; }
        out.append(buf, static_cast<std::size_t>(got));
      }
      ::close(fd);
      return true;
    }

    /// Value of "Key:   123 kB" lines (the number only).
    long long field(const std::string& text, const char* key)
    {
      const std::size_t klen = std::strlen(key);
      std::size_t pos = 0;
      while ((pos = text.find(key, pos)) != std::string::npos)
      {
        const char after = pos + klen < text.size() ? text[pos + klen] : '\0';
        if ((pos == 0 || text[pos - 1] == '\n') && (after == ':' || after == ' ' || after == '\t'))
        {
          std::size_t p = pos + klen;
          while (p < text.size() && (text[p] == ' ' || text[p] == '\t' || text[p] == ':')) { ++p; }
          return std::strtoll(text.c_str() + p, nullptr, 10);
        }
        pos += klen;
      }
      return -1;
    }

    double msSince(Clock::time_point t0)
    {
      return std::chrono::duration<double, std::milli>(Clock::now() - t0).count();
    }

    struct MallInfo2   // glibc >= 2.33 layout; looked up at run time
    {
      std::size_t arena, ordblks, smblks, hblks, hblkhd, usmblks, fsmblks, uordblks, fordblks, keepcost;
    };
    using MallInfo2Fn = MallInfo2 (*)();

    struct Sampler
    {
      std::thread fast, slow;
      std::mutex mutex;
      std::condition_variable cv;
      bool stop = false;
      long fast_ms = 100, numa_ms = 2000, smaps_ms = 10000, mallocinfo_ms = 2000;
      double tck = 100.0;
      MallInfo2Fn mallinfo2 = nullptr;
      std::size_t fast_samples = 0, slow_samples = 0;
      /// Per thread: state, faults, CPU times and last CPU at its previous row.
      /// A thread whose values have not changed gets no row and no further
      /// reads -- an idle OpenMP pool of 128 threads otherwise costs ~10 ms a tick.
      std::unordered_map<std::string, std::string> last_seen;

      double t() const { return std::chrono::duration<double>(Clock::now() - state().epoch).count(); }

      bool sleepUntil(Clock::time_point when)
      {
        std::unique_lock<std::mutex> lock(mutex);
        return cv.wait_until(lock, when, [this] { return stop; });
      }

      void fastLoop()
      {
        pthread_setname_np(pthread_self(), "odia-sampler");
        FILE* proc = openOut("proc.tsv");
        FILE* threads = openOut("threads.tsv");
        if (proc == nullptr || threads == nullptr)
        {
          if (proc) { std::fclose(proc); }
          if (threads) { std::fclose(threads); }
          return;
        }
        std::fputs("t_s\tphase\tsample_ms\tutime_s\tstime_s\tthreads\tn_run\tn_sleep\tn_disk\tn_other\t"
                   "rss_kb\trss_anon_kb\trss_file_kb\trss_shmem_kb\thwm_kb\tvm_kb\t"
                   "minflt\tmajflt\tvcsw\tnvcsw\tinblock\toublock\t"
                   "rchar\twchar\tsyscr\tsyscw\tread_bytes\twrite_bytes\t"
                   "mi_arena_b\tmi_mmap_b\tmi_inuse_b\tmi_free_b\tmi_keepcost_b\tmallinfo_ms\n", proc);
        std::fputs("t_s\ttid\tcomm\tstate\tcpu\tutime_s\tstime_s\trun_s\twait_s\tslices\t"
                   "vcsw\tnvcsw\tmigrations\tminflt\tmajflt\n", threads);
        auto next = Clock::now();
        for (;;)
        {
          sampleFast(proc, threads, false);
          next += std::chrono::milliseconds(fast_ms);
          const auto now_tp = Clock::now();
          if (next < now_tp) { next = now_tp; }   // fell behind: do not burst
          if (sleepUntil(next)) { break; }
        }
        sampleFast(proc, threads, true);   // the final state, every thread
        std::fclose(proc);
        std::fclose(threads);
      }

      void sampleFast(FILE* proc, FILE* threads, bool every_thread)
      {
        const auto c0 = Clock::now();
        const double ts = t();
        const char* phase = state().phase.load(std::memory_order_relaxed);
        std::string text;

        rusage ru{};
        ::getrusage(RUSAGE_SELF, &ru);

        long long rss = -1, anon = -1, file = -1, shmem = -1, hwm = -1, vm = -1;
        if (slurp("/proc/self/status", text))
        {
          rss = field(text, "VmRSS"); anon = field(text, "RssAnon"); file = field(text, "RssFile");
          shmem = field(text, "RssShmem"); hwm = field(text, "VmHWM"); vm = field(text, "VmSize");
        }
        long long rchar = -1, wchar = -1, syscr = -1, syscw = -1, rbytes = -1, wbytes = -1;
        if (slurp("/proc/self/io", text))
        {
          rchar = field(text, "rchar"); wchar = field(text, "wchar"); syscr = field(text, "syscr");
          syscw = field(text, "syscw"); rbytes = field(text, "read_bytes"); wbytes = field(text, "write_bytes");
        }

        double mi_ms = -1.0;
        MallInfo2 mi{};
        bool have_mi = false;
        if (mallinfo2 != nullptr)
        {
          const auto m0 = Clock::now();
          mi = mallinfo2();
          mi_ms = msSince(m0);
          have_mi = true;
        }

        // Per thread.
        int n_run = 0, n_sleep = 0, n_disk = 0, n_other = 0, n_threads = 0;
        std::string stat, sched, schedstat, path;
        std::unordered_map<std::string, std::string> seen;
        seen.reserve(last_seen.size() + 16);
        if (DIR* d = ::opendir("/proc/self/task"))
        {
          while (dirent* entry = ::readdir(d))
          {
            if (entry->d_name[0] < '0' || entry->d_name[0] > '9') { continue; }
            const std::string base = std::string("/proc/self/task/") + entry->d_name;
            if (!slurp((base + "/stat").c_str(), stat)) { continue; }   // exited meanwhile
            const auto open_paren = stat.find('(');
            const auto close_paren = stat.rfind(')');
            if (open_paren == std::string::npos || close_paren == std::string::npos) { continue; }
            std::string comm = sanitize(stat.substr(open_paren + 1, close_paren - open_paren - 1));
            // Fields after the command, numbered as in proc(5) from 3 (state).
            std::vector<const char*> f;
            f.reserve(50);
            for (std::size_t p = close_paren + 2; p < stat.size();)
            {
              f.push_back(stat.c_str() + p);
              const auto sp = stat.find(' ', p);
              if (sp == std::string::npos) { break; }
              p = sp + 1;
            }
            if (f.size() < 37) { continue; }
            const char state_c = *f[0];
            const long long minflt = std::strtoll(f[7], nullptr, 10);
            const long long majflt = std::strtoll(f[9], nullptr, 10);
            const double utime = static_cast<double>(std::strtoll(f[11], nullptr, 10)) / tck;
            const double stime = static_cast<double>(std::strtoll(f[12], nullptr, 10)) / tck;
            const long long cpu = std::strtoll(f[36], nullptr, 10);
            ++n_threads;
            switch (state_c)
            {
              case 'R': ++n_run; break;
              case 'S': ++n_sleep; break;
              case 'D': ++n_disk; break;
              default: ++n_other; break;
            }
            // Only the per-thread fields: num_threads, vsize and rss in the
            // same line are process-wide and change all the time.
            std::string key(1, state_c);
            for (const int k : {7, 9, 11, 12, 36})
            {
              key += ' ';
              key.append(f[static_cast<std::size_t>(k)], std::strcspn(f[static_cast<std::size_t>(k)], " "));
            }
            const auto prev = last_seen.find(entry->d_name);
            const bool changed = prev == last_seen.end() || prev->second != key;
            seen.emplace(entry->d_name, std::move(key));
            if (!changed && !every_thread) { continue; }
            double run_s = -1.0, wait_s = -1.0;
            long long slices = -1;
            if (slurp((base + "/schedstat").c_str(), schedstat))
            {
              char* end = nullptr;
              run_s = static_cast<double>(std::strtoull(schedstat.c_str(), &end, 10)) * 1e-9;
              wait_s = static_cast<double>(std::strtoull(end, &end, 10)) * 1e-9;
              slices = std::strtoll(end, nullptr, 10);
            }
            long long vcsw = -1, nvcsw = -1, migrations = -1;
            if (slurp((base + "/sched").c_str(), sched))
            {
              // "name   :   value" lines; the keys are unique prefixes.
              auto value = [&](const char* key) -> long long {
                const auto p = sched.find(key);
                if (p == std::string::npos) { return -1; }
                const auto colon = sched.find(':', p);
                return colon == std::string::npos ? -1 : std::strtoll(sched.c_str() + colon + 1, nullptr, 10);
              };
              vcsw = value("\nnr_voluntary_switches");
              nvcsw = value("\nnr_involuntary_switches");
              migrations = value("\nse.nr_migrations");
            }
            else if (slurp((base + "/status").c_str(), sched))
            {
              vcsw = field(sched, "voluntary_ctxt_switches");
              nvcsw = field(sched, "nonvoluntary_ctxt_switches");
            }
            std::fprintf(threads, "%.3f\t%s\t%s\t%c\t%lld\t%.2f\t%.2f\t%.6f\t%.6f\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\n",
                         ts, entry->d_name, comm.c_str(), state_c, cpu, utime, stime, run_s, wait_s,
                         slices, vcsw, nvcsw, migrations, minflt, majflt);
          }
          ::closedir(d);
        }
        last_seen.swap(seen);

        const double sample_ms = msSince(c0);
        std::fprintf(proc,
                     "%.3f\t%s\t%.3f\t%.3f\t%.3f\t%d\t%d\t%d\t%d\t%d\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t"
                     "%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t",
                     ts, phase, sample_ms,
                     static_cast<double>(ru.ru_utime.tv_sec) + static_cast<double>(ru.ru_utime.tv_usec) * 1e-6,
                     static_cast<double>(ru.ru_stime.tv_sec) + static_cast<double>(ru.ru_stime.tv_usec) * 1e-6,
                     n_threads, n_run, n_sleep, n_disk, n_other, rss, anon, file, shmem, hwm, vm,
                     ru.ru_minflt, ru.ru_majflt, ru.ru_nvcsw, ru.ru_nivcsw, ru.ru_inblock, ru.ru_oublock,
                     rchar, wchar, syscr, syscw, rbytes, wbytes);
        if (have_mi)
        {
          std::fprintf(proc, "%zu\t%zu\t%zu\t%zu\t%zu\t%.3f\n", mi.arena, mi.hblkhd, mi.uordblks,
                       mi.fordblks, mi.keepcost, mi_ms);
        }
        else { std::fputs("-1\t-1\t-1\t-1\t-1\t-1\n", proc); }
        ++fast_samples;
      }

      void slowLoop()
      {
        pthread_setname_np(pthread_self(), "odia-sampler2");
        FILE* numa = numa_ms > 0 ? openOut("numa.tsv") : nullptr;
        FILE* smaps = smaps_ms > 0 ? openOut("smaps.tsv") : nullptr;
        FILE* minfo = mallocinfo_ms > 0 ? openOut("mallocinfo.tsv") : nullptr;
        if (numa) { std::fputs("t_s\tphase\tcost_ms\tnode\ttotal_kb\tanon_kb\tfile_kb\thuge_kb\n", numa); }
        if (smaps)
        {
          std::fputs("t_s\tphase\tcost_ms\trss_kb\tpss_kb\tpss_anon_kb\tpss_file_kb\tpss_dirty_kb\t"
                     "anonymous_kb\tanon_huge_kb\tprivate_dirty_kb\tprivate_clean_kb\tshared_clean_kb\tswap_kb\n", smaps);
        }
        if (minfo)
        {
          std::fputs("t_s\tphase\tcost_ms\tarenas\tfast_free_b\trest_free_b\tmmap_count\tmmap_b\t"
                     "system_current_b\tsystem_max_b\taspace_total_b\taspace_mprotect_b\n", minfo);
        }
        if (!numa && !smaps && !minfo) { return; }

        // Staggered so the expensive page-table walks do not coincide.
        const auto start = Clock::now();
        auto due_numa = start + std::chrono::milliseconds(numa_ms > 0 ? numa_ms / 2 : 0);
        auto due_smaps = start + std::chrono::milliseconds(smaps_ms > 0 ? smaps_ms / 4 : 0);
        auto due_minfo = start + std::chrono::milliseconds(mallocinfo_ms > 0 ? mallocinfo_ms / 3 : 0);
        for (;;)
        {
          auto next = Clock::time_point::max();
          if (numa) { next = std::min(next, due_numa); }
          if (smaps) { next = std::min(next, due_smaps); }
          if (minfo) { next = std::min(next, due_minfo); }
          if (sleepUntil(next)) { break; }
          const auto now_tp = Clock::now();
          if (numa && now_tp >= due_numa)
          {
            sampleNuma(numa);
            due_numa = std::max(due_numa + std::chrono::milliseconds(numa_ms), Clock::now());
          }
          if (smaps && now_tp >= due_smaps)
          {
            sampleSmaps(smaps);
            due_smaps = std::max(due_smaps + std::chrono::milliseconds(smaps_ms), Clock::now());
          }
          if (minfo && now_tp >= due_minfo)
          {
            sampleMallocInfo(minfo);
            due_minfo = std::max(due_minfo + std::chrono::milliseconds(mallocinfo_ms), Clock::now());
          }
          ++slow_samples;
          for (FILE* f : {numa, smaps, minfo}) { if (f) { std::fflush(f); } }
        }
        for (FILE* f : {numa, smaps, minfo}) { if (f) { std::fclose(f); } }
      }

      void sampleNuma(FILE* out)
      {
        const auto c0 = Clock::now();
        const double ts = t();
        const char* phase = state().phase.load(std::memory_order_relaxed);
        std::string text;
        if (!slurp("/proc/self/numa_maps", text)) { return; }
        struct Node { long long total = 0, anon = 0, file = 0, huge = 0; };
        std::map<int, Node> nodes;
        std::size_t line_start = 0;
        while (line_start < text.size())
        {
          std::size_t line_end = text.find('\n', line_start);
          if (line_end == std::string::npos) { line_end = text.size(); }
          const std::string_view line(text.data() + line_start, line_end - line_start);
          long long kps = 4;
          if (const auto k = line.find(" kernelpagesize_kB="); k != std::string_view::npos)
          { kps = std::strtoll(line.data() + k + 19, nullptr, 10); }
          const bool is_file = line.find(" file=") != std::string_view::npos;
          for (std::size_t p = line.find(" N"); p != std::string_view::npos; p = line.find(" N", p + 2))
          {
            char* end = nullptr;
            const long node = std::strtol(line.data() + p + 2, &end, 10);
            if (end == line.data() + p + 2 || *end != '=') { continue; }
            const long long kb = std::strtoll(end + 1, nullptr, 10) * kps;
            auto& n = nodes[static_cast<int>(node)];
            n.total += kb;
            (is_file ? n.file : n.anon) += kb;
            if (kps > 4) { n.huge += kb; }
          }
          line_start = line_end + 1;
        }
        const double cost = msSince(c0);
        for (const auto& [node, n] : nodes)
        {
          std::fprintf(out, "%.3f\t%s\t%.3f\t%d\t%lld\t%lld\t%lld\t%lld\n", ts, phase, cost, node,
                       n.total, n.anon, n.file, n.huge);
        }
      }

      void sampleSmaps(FILE* out)
      {
        const auto c0 = Clock::now();
        const double ts = t();
        const char* phase = state().phase.load(std::memory_order_relaxed);
        std::string text;
        if (!slurp("/proc/self/smaps_rollup", text)) { return; }
        const double cost = msSince(c0);
        std::fprintf(out, "%.3f\t%s\t%.3f\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\n",
                     ts, phase, cost, field(text, "Rss"), field(text, "Pss"), field(text, "Pss_Anon"),
                     field(text, "Pss_File"), field(text, "Pss_Dirty"), field(text, "Anonymous"),
                     field(text, "AnonHugePages"), field(text, "Private_Dirty"),
                     field(text, "Private_Clean"), field(text, "Shared_Clean"), field(text, "Swap"));
      }

      void sampleMallocInfo(FILE* out)
      {
#if defined(__GLIBC__)
        const auto c0 = Clock::now();
        const double ts = t();
        const char* phase = state().phase.load(std::memory_order_relaxed);
        char* data = nullptr;
        std::size_t size = 0;
        FILE* mem = ::open_memstream(&data, &size);
        if (mem == nullptr) { return; }
        const int rc = ::malloc_info(0, mem);
        std::fclose(mem);
        const double cost = msSince(c0);
        if (rc != 0 || data == nullptr) { std::free(data); return; }
        const std::string xml(data, size);
        std::free(data);
        long long arenas = 0;
        for (std::size_t p = xml.find("<heap nr="); p != std::string::npos; p = xml.find("<heap nr=", p + 1)) { ++arenas; }
        // The process-wide totals follow the last </heap>.
        const std::size_t tail = xml.rfind("</heap>");
        const std::string totals = tail == std::string::npos ? xml : xml.substr(tail);
        auto attr = [&](const char* tag, const char* name) -> long long {
          const auto p = totals.find(tag);
          if (p == std::string::npos) { return -1; }
          const auto close = totals.find("/>", p);
          const auto q = totals.find(name, p);
          if (q == std::string::npos || (close != std::string::npos && q > close)) { return -1; }
          return std::strtoll(totals.c_str() + q + std::strlen(name), nullptr, 10);
        };
        std::fprintf(out, "%.3f\t%s\t%.3f\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\n",
                     ts, phase, cost, arenas,
                     attr("<total type=\"fast\"", "size=\""), attr("<total type=\"rest\"", "size=\""),
                     attr("<total type=\"mmap\"", "count=\""), attr("<total type=\"mmap\"", "size=\""),
                     attr("<system type=\"current\"", "size=\""), attr("<system type=\"max\"", "size=\""),
                     attr("<aspace type=\"total\"", "size=\""), attr("<aspace type=\"mprotect\"", "size=\""));
#else
        (void)out;
#endif
      }
    };

    Sampler*& sampler()
    {
      static Sampler* s = nullptr;
      return s;
    }

    /// cpu, node, package, core, L3 instance and whether the affinity mask allows it.
    void writeTopology()
    {
      FILE* out = openOut("topology.tsv");
      if (out == nullptr) { return; }
      std::fputs("cpu\tnode\tpackage\tcore_id\tl3_id\tallowed\n", out);
      cpu_set_t allowed;
      CPU_ZERO(&allowed);
      const bool have_mask = ::sched_getaffinity(0, sizeof allowed, &allowed) == 0;
      const long n = ::sysconf(_SC_NPROCESSORS_CONF);
      std::map<long, int> node_of;
      std::string text;
      for (int node = 0; node < 64; ++node)
      {
        const std::string p = "/sys/devices/system/node/node" + std::to_string(node) + "/cpulist";
        if (!slurp(p.c_str(), text)) { continue; }
        // "0-95,192-287"
        std::size_t pos = 0;
        while (pos < text.size())
        {
          char* end = nullptr;
          const long a = std::strtol(text.c_str() + pos, &end, 10);
          if (end == text.c_str() + pos) { break; }
          long b = a;
          if (*end == '-') { b = std::strtol(end + 1, &end, 10); }
          for (long c = a; c <= b; ++c) { node_of[c] = node; }
          pos = static_cast<std::size_t>(end - text.c_str());
          if (pos < text.size() && text[pos] == ',') { ++pos; } else { break; }
        }
      }
      auto readLong = [&](const std::string& p) -> long {
        return slurp(p.c_str(), text) ? std::strtol(text.c_str(), nullptr, 10) : -1;
      };
      for (long cpu = 0; cpu < n; ++cpu)
      {
        const std::string base = "/sys/devices/system/cpu/cpu" + std::to_string(cpu);
        const auto it = node_of.find(cpu);
        std::fprintf(out, "%ld\t%d\t%ld\t%ld\t%ld\t%d\n", cpu, it == node_of.end() ? -1 : it->second,
                     readLong(base + "/topology/physical_package_id"), readLong(base + "/topology/core_id"),
                     readLong(base + "/cache/index3/id"),
                     have_mask ? (CPU_ISSET(static_cast<int>(cpu), &allowed) ? 1 : 0) : -1);
      }
      std::fclose(out);
    }
#endif // __linux__

    std::string isoUtc(std::time_t t)
    {
      std::tm tm{};
#if defined(_WIN32)
      gmtime_s(&tm, &t);
#else
      gmtime_r(&t, &tm);
#endif
      char buf[32];
      std::strftime(buf, sizeof buf, "%Y-%m-%dT%H:%M:%SZ", &tm);
      return buf;
    }
  } // namespace

  namespace detail
  {
    bool readEnabled() noexcept
    {
      const char* dir = std::getenv("ODIA_TRACE");
      if (dir == nullptr || *dir == '\0') { return false; }
      try
      {
        auto& s = state();   // fixes the epoch
        s.dir = dir;
      }
      catch (...) { return false; }
      return true;
    }
  } // namespace detail

  double now() noexcept
  {
    if (!enabled()) { return 0.0; }
    return std::chrono::duration<double>(Clock::now() - state().epoch).count();
  }

  double threadCpu() noexcept
  {
#if defined(__linux__) || defined(__APPLE__)
    timespec ts{};
    if (::clock_gettime(CLOCK_THREAD_CPUTIME_ID, &ts) == 0)
    { return static_cast<double>(ts.tv_sec) + static_cast<double>(ts.tv_nsec) * 1e-9; }
#endif
    return std::numeric_limits<double>::quiet_NaN();
  }

  void Scope::begin_(const char* name, int level) noexcept
  {
    active_ = true;
    level_ = level;
    name_ = name;
    if (level == 0) { outer_phase_ = state().phase.exchange(name); }
    cpu0_ = threadCpu();
    t0_ = now();
  }

  void Scope::end_() noexcept
  {
    const double t1 = now();
    const double cpu = threadCpu() - cpu0_;
    auto& s = state();
    push({t0_, t1, cpu, currentTid(), level_,
          level_ == 0 ? name_ : s.phase.load(std::memory_order_relaxed), name_, fields_});
    if (level_ == 0) { s.phase.store(outer_phase_ != nullptr ? outer_phase_ : "-"); }
    active_ = false;
  }

  void record(const char* name, int level, double t0, double t1, double thread_cpu_s,
              const Fields& fields) noexcept
  {
    if (!enabled()) { return; }
    push({t0, t1, thread_cpu_s, currentTid(), level,
          state().phase.load(std::memory_order_relaxed), name, fields});
  }

  void nameWorker(const char* model, std::size_t session) noexcept
  {
    if (!enabled()) { return; }
#if defined(__linux__)
    char name[16];
    std::snprintf(name, sizeof name, "odia-%s-s%zu", model, session);
    pthread_setname_np(pthread_self(), name);
#else
    (void)model; (void)session;
#endif
  }

  void note(const std::string& key, const std::string& value)
  {
    if (!enabled()) { return; }
    auto& s = state();
    const std::lock_guard<std::mutex> lock(s.mutex);
    s.notes.emplace_back(key, value);
  }

  const std::string& directory()
  {
    static const std::string none;
    return enabled() ? state().dir : none;
  }

  Session::Session(int argc, const char** argv)
  {
    if (!enabled()) { return; }
    auto& s = state();
    {
      const std::lock_guard<std::mutex> lock(s.mutex);
      if (s.session_active) { return; }
      s.session_active = true;
    }
    owner_ = true;
    std::error_code ec;
    std::filesystem::create_directories(s.dir, ec);

    std::string cmd;
    for (int i = 0; i < argc; ++i) { cmd += (i ? " " : ""); cmd += argv[i]; }
    note("trace_format", "1");
    note("start_utc", isoUtc(std::time(nullptr)));
    note("start_epoch_s", std::to_string(std::chrono::duration<double>(
                            std::chrono::system_clock::now().time_since_epoch()).count()));
    note("session_start_s", std::to_string(now()));
    note("argv", sanitize(cmd));
    note("hardware_concurrency", std::to_string(std::thread::hardware_concurrency()));
    for (const char* var : {"ODIA_TRACE", "ODIA_TRACE_INTERVAL_MS", "ODIA_TRACE_NUMA_MS",
                            "ODIA_TRACE_SMAPS_MS", "ODIA_TRACE_MALLOCINFO_MS", "ODIA_ORT_PROFILE",
                            "OMP_NUM_THREADS", "MALLOC_ARENA_MAX", "GLIBC_TUNABLES"})
    {
      if (const char* v = std::getenv(var)) { note(std::string("env.") + var, sanitize(v)); }
    }
#if defined(__linux__) || defined(__APPLE__)
    note("pid", std::to_string(::getpid()));
    char host[256] = {0};
    if (::gethostname(host, sizeof host - 1) == 0) { note("host", host); }
#endif
#if defined(__GLIBC__)
    note("glibc", gnu_get_libc_version());
#endif
#if defined(__linux__)
    utsname u{};
    if (::uname(&u) == 0) { note("kernel", u.release); }
    cpu_set_t mask;
    CPU_ZERO(&mask);
    if (::sched_getaffinity(0, sizeof mask, &mask) == 0) { note("affinity_cpus", std::to_string(CPU_COUNT(&mask))); }
    // When the process started, on the trace clock: TOPP's start-up before the
    // session opened is otherwise invisible.
    {
      std::string stat;
      timespec boot{};
      if (slurp("/proc/self/stat", stat) && ::clock_gettime(CLOCK_BOOTTIME, &boot) == 0)
      {
        const auto close_paren = stat.rfind(')');
        std::size_t p = close_paren + 2;
        for (int k = 3; k < 22 && p < stat.size(); ++k) { p = stat.find(' ', p) + 1; }
        const double start_ticks = static_cast<double>(std::strtoull(stat.c_str() + p, nullptr, 10));
        const double boot_now = static_cast<double>(boot.tv_sec) + static_cast<double>(boot.tv_nsec) * 1e-9;
        const double age = boot_now - start_ticks / static_cast<double>(::sysconf(_SC_CLK_TCK));
        note("process_start_s", std::to_string(now() - age));
      }
    }
    writeTopology();

    auto* smp = new Sampler;
    smp->tck = static_cast<double>(::sysconf(_SC_CLK_TCK));
    smp->fast_ms = std::max(10L, envLong("ODIA_TRACE_INTERVAL_MS", 100));
    smp->numa_ms = envLong("ODIA_TRACE_NUMA_MS", 2000);
    smp->smaps_ms = envLong("ODIA_TRACE_SMAPS_MS", 10000);
    smp->mallocinfo_ms = envLong("ODIA_TRACE_MALLOCINFO_MS", 2000);
    smp->mallinfo2 = reinterpret_cast<MallInfo2Fn>(::dlsym(RTLD_DEFAULT, "mallinfo2"));
    note("interval_ms", std::to_string(smp->fast_ms));
    note("numa_ms", std::to_string(smp->numa_ms));
    note("smaps_ms", std::to_string(smp->smaps_ms));
    note("mallocinfo_ms", std::to_string(smp->mallocinfo_ms));
    note("mallinfo2", smp->mallinfo2 != nullptr ? "yes" : "no");
    sampler() = smp;
    smp->fast = std::thread([smp] { smp->fastLoop(); });
    smp->slow = std::thread([smp] { smp->slowLoop(); });
#endif
    std::fprintf(stderr, "[trace] writing to %s\n", s.dir.c_str());
  }

  Session::~Session()
  {
    if (!owner_) { return; }
    auto& s = state();
    const double end_s = now();
#if defined(__linux__)
    if (Sampler* smp = sampler())
    {
      {
        const std::lock_guard<std::mutex> lock(smp->mutex);
        smp->stop = true;
      }
      smp->cv.notify_all();
      if (smp->fast.joinable()) { smp->fast.join(); }
      if (smp->slow.joinable()) { smp->slow.join(); }
      note("fast_samples", std::to_string(smp->fast_samples));
      note("slow_samples", std::to_string(smp->slow_samples));
      sampler() = nullptr;
      delete smp;
    }
#endif
    buffer().flush();
    std::vector<Event> events;
    std::vector<std::pair<std::string, std::string>> notes;
    {
      const std::lock_guard<std::mutex> lock(s.mutex);
      events.swap(s.events);
      notes = s.notes;
      s.session_active = false;
    }
    std::stable_sort(events.begin(), events.end(),
                     [](const Event& a, const Event& b) { return a.t0 < b.t0; });
    if (FILE* out = openOut("events.tsv"))
    {
      std::fputs("t0_s\tt1_s\tdur_s\ttid\tcpu_s\tlevel\tphase\tname\tmodel\tcall\tidx\tsession\trows\tlen\twait_s\n", out);
      for (const auto& e : events)
      {
        const auto& f = e.fields;
        std::fprintf(out, "%.6f\t%.6f\t%.6f\t%lld\t%.6f\t%d\t%s\t%s\t%s\t%lld\t%lld\t%lld\t%lld\t%lld\t%.6f\n",
                     e.t0, e.t1, e.t1 - e.t0, static_cast<long long>(e.tid), e.cpu, e.level,
                     e.phase, e.name, f.model, static_cast<long long>(f.call),
                     static_cast<long long>(f.idx), static_cast<long long>(f.session),
                     static_cast<long long>(f.rows), static_cast<long long>(f.len), f.wait_s);
      }
      std::fclose(out);
    }
    notes.emplace_back("end_s", std::to_string(end_s));
    notes.emplace_back("events", std::to_string(events.size()));
    if (FILE* out = openOut("meta.tsv"))
    {
      std::fputs("key\tvalue\n", out);
      for (const auto& [k, v] : notes) { std::fprintf(out, "%s\t%s\n", k.c_str(), sanitize(v).c_str()); }
      std::fclose(out);
    }
    std::fprintf(stderr, "[trace] wrote %zu events to %s\n", events.size(), s.dir.c_str());
  }

} // namespace ODIA::trace
