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

    /// A number of seconds; "off" or a negative value give -1 (disabled).
    [[maybe_unused]] double envSeconds(const char* name, double fallback)
    {
      const char* v = std::getenv(name);
      if (v == nullptr || *v == '\0') { return fallback; }
      if (std::strcmp(v, "off") == 0) { return -1.0; }
      char* end = nullptr;
      const double x = std::strtod(v, &end);
      if (end == v) { return fallback; }
      return x < 0.0 ? -1.0 : x;
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

    /// Whole file into @p out with plain read(2): no stdio, no locale. Used
    /// for files of unbounded size (numa_maps) at phase boundaries only.
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

    /// A small /proc file into a caller's buffer, NUL-terminated; no allocation.
    /// Returns the length, or -1. Truncates at @p cap - 1 bytes.
    ssize_t readSmall(const char* path, char* buf, std::size_t cap)
    {
      const int fd = ::open(path, O_RDONLY | O_CLOEXEC);
      if (fd < 0) { return -1; }
      std::size_t len = 0;
      while (len + 1 < cap)
      {
        const ssize_t got = ::read(fd, buf + len, cap - 1 - len);
        if (got < 0) { ::close(fd); return -1; }
        if (got == 0) { break; }
        len += static_cast<std::size_t>(got);
      }
      ::close(fd);
      buf[len] = '\0';
      return static_cast<ssize_t>(len);
    }

    /// Value of "Key:   123 kB" or "key 123" lines (the number only), -1 if absent.
    long long field(std::string_view text, const char* key)
    {
      const std::size_t klen = std::strlen(key);
      std::size_t pos = 0;
      while ((pos = text.find(key, pos)) != std::string_view::npos)
      {
        const char after = pos + klen < text.size() ? text[pos + klen] : '\0';
        if ((pos == 0 || text[pos - 1] == '\n') && (after == ':' || after == ' ' || after == '\t'))
        {
          std::size_t p = pos + klen;
          while (p < text.size() && (text[p] == ' ' || text[p] == '\t' || text[p] == ':')) { ++p; }
          return std::strtoll(text.data() + p, nullptr, 10);
        }
        pos += klen;
      }
      return -1;
    }

    double msSince(Clock::time_point t0)
    {
      return std::chrono::duration<double, std::milli>(Clock::now() - t0).count();
    }

    double traceSeconds()
    {
      return std::chrono::duration<double>(Clock::now() - state().epoch).count();
    }

    /// The per-thread values whose change makes a thread worth a row.
    struct ThreadRec
    {
      long tid = 0;
      char state = '?';
      long long minflt = 0, majflt = 0, utime = 0, stime = 0, cpu = -1;
      bool same(const ThreadRec& o) const
      {
        return state == o.state && minflt == o.minflt && majflt == o.majflt &&
               utime == o.utime && stime == o.stime && cpu == o.cpu;
      }
    };

    /// The periodic sampler. It reads only counters that the kernel hands out
    /// without taking the allocator's arena locks or the process's mmap lock:
    /// /proc/self/status and /proc/self/io, getrusage(), /proc/vmstat, and per
    /// thread task/<tid>/stat, task/<tid>/schedstat and task/<tid>/status. In
    /// the steady state it allocates nothing: buffers are fixed or reused.
    struct Sampler
    {
      std::thread fast;
      std::mutex mutex;
      std::condition_variable cv;
      bool stop = false;
      long fast_ms = 100;
      double tck = 100.0;
      std::size_t fast_samples = 0;
      double sample_ms_max = 0.0, sample_ms_sum = 0.0;
      std::vector<ThreadRec> prev, cur;   // sorted by tid
      char buf[16384];
      char tbuf[4096];

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
        prev.reserve(4096);
        cur.reserve(4096);
        std::fputs("t_s\tphase\tsample_ms\tutime_s\tstime_s\tthreads\tn_run\tn_sleep\tn_disk\tn_other\t"
                   "rss_kb\trss_anon_kb\trss_file_kb\trss_shmem_kb\thwm_kb\tvm_kb\t"
                   "minflt\tmajflt\tvcsw\tnvcsw\tinblock\toublock\t"
                   "rchar\twchar\tsyscr\tsyscw\tread_bytes\twrite_bytes\t"
                   "sys_numa_hit\tsys_numa_miss\tsys_numa_local\tsys_numa_other\t"
                   "sys_numa_hint_faults\tsys_numa_hint_faults_local\tsys_numa_pages_migrated\n", proc);
        std::fputs("t_s\ttid\tcomm\tstate\tcpu\tutime_s\tstime_s\trun_s\twait_s\tslices\t"
                   "vcsw\tnvcsw\tminflt\tmajflt\n", threads);
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
        const double ts = traceSeconds();
        const char* phase = state().phase.load(std::memory_order_relaxed);

        rusage ru{};
        ::getrusage(RUSAGE_SELF, &ru);   // sums the thread group under the signal lock only

        long long rss = -1, anon = -1, file = -1, shmem = -1, hwm = -1, vm = -1;
        if (readSmall("/proc/self/status", buf, sizeof buf) > 0)
        {
          const std::string_view text(buf);
          rss = field(text, "VmRSS"); anon = field(text, "RssAnon"); file = field(text, "RssFile");
          shmem = field(text, "RssShmem"); hwm = field(text, "VmHWM"); vm = field(text, "VmSize");
        }
        long long rchar = -1, wchar = -1, syscr = -1, syscw = -1, rbytes = -1, wbytes = -1;
        if (readSmall("/proc/self/io", buf, sizeof buf) > 0)
        {
          const std::string_view text(buf);
          rchar = field(text, "rchar"); wchar = field(text, "wchar"); syscr = field(text, "syscr");
          syscw = field(text, "syscw"); rbytes = field(text, "read_bytes"); wbytes = field(text, "write_bytes");
        }
        // System-wide page-allocation and NUMA-balancing counters (not per
        // process): on an otherwise idle node their rate is this run's.
        long long nh = -1, nm = -1, nl = -1, no = -1, hf = -1, hfl = -1, mig = -1;
        if (readSmall("/proc/vmstat", buf, sizeof buf) > 0)
        {
          const std::string_view text(buf);
          nh = field(text, "numa_hit"); nm = field(text, "numa_miss"); nl = field(text, "numa_local");
          no = field(text, "numa_other"); hf = field(text, "numa_hint_faults");
          hfl = field(text, "numa_hint_faults_local"); mig = field(text, "numa_pages_migrated");
        }

        // Per thread. A thread whose state, faults, CPU times and last CPU are
        // unchanged since its previous row gets no row and no further reads.
        int n_run = 0, n_sleep = 0, n_disk = 0, n_other = 0, n_threads = 0;
        cur.clear();
        char path[64];
        if (DIR* d = ::opendir("/proc/self/task"))
        {
          while (dirent* entry = ::readdir(d))
          {
            if (entry->d_name[0] < '0' || entry->d_name[0] > '9') { continue; }
            const long tid = std::strtol(entry->d_name, nullptr, 10);
            std::snprintf(path, sizeof path, "/proc/self/task/%ld/stat", tid);
            const ssize_t len = readSmall(path, tbuf, sizeof tbuf);
            if (len <= 0) { continue; }   // exited meanwhile
            char* open_paren = std::strchr(tbuf, '(');
            char* close_paren = std::strrchr(tbuf, ')');
            if (open_paren == nullptr || close_paren == nullptr || close_paren + 2 >= tbuf + len) { continue; }
            // Fields after the command, numbered as in proc(5) from 3 (state).
            const char* f[40];
            std::size_t nf = 0;
            for (char* p = close_paren + 2; p < tbuf + len && nf < 40;)
            {
              f[nf++] = p;
              char* sp = std::strchr(p, ' ');
              if (sp == nullptr) { break; }
              p = sp + 1;
            }
            if (nf < 37) { continue; }
            ThreadRec r;
            r.tid = tid;
            r.state = *f[0];
            r.minflt = std::strtoll(f[7], nullptr, 10);
            r.majflt = std::strtoll(f[9], nullptr, 10);
            r.utime = std::strtoll(f[11], nullptr, 10);
            r.stime = std::strtoll(f[12], nullptr, 10);
            r.cpu = std::strtoll(f[36], nullptr, 10);
            ++n_threads;
            switch (r.state)
            {
              case 'R': ++n_run; break;
              case 'S': ++n_sleep; break;
              case 'D': ++n_disk; break;
              default: ++n_other; break;
            }
            cur.push_back(r);
            const auto it = std::lower_bound(prev.begin(), prev.end(), r.tid,
                                             [](const ThreadRec& a, long t) { return a.tid < t; });
            const bool changed = it == prev.end() || it->tid != r.tid || !it->same(r);
            if (!changed && !every_thread) { continue; }

            char comm[32];
            const std::size_t clen = std::min<std::size_t>(sizeof comm - 1,
                                                           static_cast<std::size_t>(close_paren - open_paren - 1));
            std::memcpy(comm, open_paren + 1, clen);
            comm[clen] = '\0';
            for (char* c = comm; *c; ++c) { if (*c == '\t' || *c == '\n') { *c = ' '; } }

            double run_s = -1.0, wait_s = -1.0;
            long long slices = -1;
            std::snprintf(path, sizeof path, "/proc/self/task/%ld/schedstat", tid);
            if (readSmall(path, tbuf, sizeof tbuf) > 0)
            {
              char* e = nullptr;
              run_s = static_cast<double>(std::strtoull(tbuf, &e, 10)) * 1e-9;
              wait_s = static_cast<double>(std::strtoull(e, &e, 10)) * 1e-9;
              slices = std::strtoll(e, nullptr, 10);
            }
            long long vcsw = -1, nvcsw = -1;
            std::snprintf(path, sizeof path, "/proc/self/task/%ld/status", tid);
            if (readSmall(path, tbuf, sizeof tbuf) > 0)
            {
              const std::string_view text(tbuf);
              vcsw = field(text, "voluntary_ctxt_switches");
              nvcsw = field(text, "nonvoluntary_ctxt_switches");
            }
            std::fprintf(threads, "%.3f\t%ld\t%s\t%c\t%lld\t%.2f\t%.2f\t%.6f\t%.6f\t%lld\t%lld\t%lld\t%lld\t%lld\n",
                         ts, r.tid, comm, r.state, r.cpu, static_cast<double>(r.utime) / tck,
                         static_cast<double>(r.stime) / tck, run_s, wait_s, slices, vcsw, nvcsw,
                         r.minflt, r.majflt);
          }
          ::closedir(d);
        }
        std::sort(cur.begin(), cur.end(), [](const ThreadRec& a, const ThreadRec& b) { return a.tid < b.tid; });
        prev.swap(cur);

        const double sample_ms = msSince(c0);
        sample_ms_sum += sample_ms;
        sample_ms_max = std::max(sample_ms_max, sample_ms);
        std::fprintf(proc,
                     "%.3f\t%s\t%.3f\t%.3f\t%.3f\t%d\t%d\t%d\t%d\t%d\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t"
                     "%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t"
                     "%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\n",
                     ts, phase, sample_ms,
                     static_cast<double>(ru.ru_utime.tv_sec) + static_cast<double>(ru.ru_utime.tv_usec) * 1e-6,
                     static_cast<double>(ru.ru_stime.tv_sec) + static_cast<double>(ru.ru_stime.tv_usec) * 1e-6,
                     n_threads, n_run, n_sleep, n_disk, n_other, rss, anon, file, shmem, hwm, vm,
                     ru.ru_minflt, ru.ru_majflt, ru.ru_nvcsw, ru.ru_nivcsw, ru.ru_inblock, ru.ru_oublock,
                     rchar, wchar, syscr, syscw, rbytes, wbytes, nh, nm, nl, no, hf, hfl, mig);
        ++fast_samples;
      }
    };

    Sampler*& sampler()
    {
      static Sampler* s = nullptr;
      return s;
    }

    /// The views that walk the heap or the page tables while holding a lock
    /// every other thread may need -- malloc_info (each arena's lock, and it
    /// visits every free chunk), numa_maps and smaps_rollup (the mmap lock,
    /// every page). Taken synchronously when a phase begins, i.e. while the
    /// worker threads of the previous phase are idle, and at most once per
    /// min_interval_s; never periodically. The whole snapshot is a level-0
    /// event of its own ("trace_snapshot"), so its cost is visible and kept
    /// out of the phases.
    struct Snapshots
    {
      std::mutex mutex;
      double min_interval_s = 30.0;
      double last_s = -1e300;
      FILE* numa = nullptr;
      FILE* smaps = nullptr;
      FILE* minfo = nullptr;
      std::size_t taken = 0, skipped = 0;
      double total_ms = 0.0;
      std::string text;   // reused

      void take(const char* at)
      {
        const double ts = traceSeconds();
        if (numa)
        {
          const auto c0 = Clock::now();
          if (slurp("/proc/self/numa_maps", text))
          {
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
                char* e = nullptr;
                const long node = std::strtol(line.data() + p + 2, &e, 10);
                if (e == line.data() + p + 2 || *e != '=') { continue; }
                const long long kb = std::strtoll(e + 1, nullptr, 10) * kps;
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
              std::fprintf(numa, "%.3f\t%s\t%.3f\t%d\t%lld\t%lld\t%lld\t%lld\n", ts, at, cost, node,
                           n.total, n.anon, n.file, n.huge);
            }
          }
        }
        if (smaps)
        {
          const auto c0 = Clock::now();
          if (slurp("/proc/self/smaps_rollup", text))
          {
            const double cost = msSince(c0);
            const std::string_view t(text);
            std::fprintf(smaps, "%.3f\t%s\t%.3f\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\n",
                         ts, at, cost, field(t, "Rss"), field(t, "Pss"), field(t, "Pss_Anon"),
                         field(t, "Pss_File"), field(t, "Pss_Dirty"), field(t, "Anonymous"),
                         field(t, "AnonHugePages"), field(t, "Private_Dirty"),
                         field(t, "Private_Clean"), field(t, "Shared_Clean"), field(t, "Swap"));
          }
        }
#if defined(__GLIBC__)
        if (minfo)
        {
          const auto c0 = Clock::now();
          char* data = nullptr;
          std::size_t size = 0;
          if (FILE* mem = ::open_memstream(&data, &size))
          {
            const int rc = ::malloc_info(0, mem);
            std::fclose(mem);
            const double cost = msSince(c0);
            if (rc == 0 && data != nullptr)
            {
              const std::string_view xml(data, size);
              long long arenas = 0;
              for (std::size_t p = xml.find("<heap nr="); p != std::string_view::npos; p = xml.find("<heap nr=", p + 1))
              { ++arenas; }
              // The process-wide totals follow the last </heap>.
              const std::size_t tail = xml.rfind("</heap>");
              const std::string_view totals = tail == std::string_view::npos ? xml : xml.substr(tail);
              auto attr = [&](const char* tag, const char* name) -> long long {
                const auto p = totals.find(tag);
                if (p == std::string_view::npos) { return -1; }
                const auto close = totals.find("/>", p);
                const auto q = totals.find(name, p);
                if (q == std::string_view::npos || (close != std::string_view::npos && q > close)) { return -1; }
                return std::strtoll(totals.data() + q + std::strlen(name), nullptr, 10);
              };
              const long long fast = attr("<total type=\"fast\"", "size=\"");
              const long long rest = attr("<total type=\"rest\"", "size=\"");
              const long long current = attr("<system type=\"current\"", "size=\"");
              const long long mmap_b = attr("<total type=\"mmap\"", "size=\"");
              const long long inuse = (current >= 0 && fast >= 0 && rest >= 0) ? current - fast - rest : -1;
              std::fprintf(minfo, "%.3f\t%s\t%.3f\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\n",
                           ts, at, cost, arenas, fast, rest,
                           attr("<total type=\"mmap\"", "count=\""), mmap_b, current,
                           attr("<system type=\"max\"", "size=\""),
                           attr("<aspace type=\"total\"", "size=\""), attr("<aspace type=\"mprotect\"", "size=\""),
                           inuse);
            }
          }
          std::free(data);
        }
#endif
        for (FILE* f : {numa, smaps, minfo}) { if (f) { std::fflush(f); } }
      }
    };

    Snapshots*& snapshots()
    {
      static Snapshots* s = nullptr;
      return s;
    }

    /// At the start of a phase: a snapshot if one is due (see Snapshots).
    void boundarySnapshot(const char* at, bool force) noexcept
    {
      Snapshots* snap = snapshots();
      if (snap == nullptr) { return; }
      std::unique_lock<std::mutex> lock(snap->mutex, std::try_to_lock);
      if (!lock.owns_lock()) { return; }   // another thread is taking one
      const double t0 = traceSeconds();
      if (!force && t0 - snap->last_s < snap->min_interval_s) { ++snap->skipped; return; }
      snap->last_s = t0;
      auto& s = state();
      const char* outer = s.phase.exchange("trace_snapshot");
      const double cpu0 = threadCpu();
      try { snap->take(at); }
      catch (...) {}
      const double t1 = traceSeconds();
      Fields f;
      f.idx = static_cast<std::int64_t>(snap->taken);
      push({t0, t1, threadCpu() - cpu0, currentTid(), 0, "trace_snapshot", "trace_snapshot", f});
      s.phase.store(outer != nullptr ? outer : "-");
      snap->total_ms += (t1 - t0) * 1e3;
      ++snap->taken;
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
#if defined(__linux__)
    // Before the phase opens, so that the snapshot's cost is not the phase's.
    if (level == 0) { boundarySnapshot(name, false); }
#endif
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
    note("trace_format", "2");
    note("start_utc", isoUtc(std::time(nullptr)));
    note("start_epoch_s", std::to_string(std::chrono::duration<double>(
                            std::chrono::system_clock::now().time_since_epoch()).count()));
    note("session_start_s", std::to_string(now()));
    note("argv", sanitize(cmd));
    note("hardware_concurrency", std::to_string(std::thread::hardware_concurrency()));
    for (const char* var : {"ODIA_TRACE", "ODIA_TRACE_INTERVAL_MS", "ODIA_TRACE_SNAPSHOT_MIN_S",
                            "ODIA_ORT_PROFILE", "OMP_NUM_THREADS", "MALLOC_ARENA_MAX", "GLIBC_TUNABLES",
                            "ODIA_SOLUTIONS", "ODIA_MAX_SESSIONS", "ODIA_MS2_ROWS", "ODIA_LSTM_ROWS",
                            "ODIA_MS2_BLOCK", "ODIA_LPT", "ODIA_INTERLEAVE", "ODIA_PAR_POST",
                            "ODIA_PAR_ENUM", "ODIA_PAR_WRITE", "ODIA_PAR_CREATE", "ODIA_TRIM", "ODIA_MS2_PRUNED",
                            "OPENBLAS_NUM_THREADS"})
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
    note("interval_ms", std::to_string(smp->fast_ms));
    sampler() = smp;
    smp->fast = std::thread([smp] { smp->fastLoop(); });

    const double snapshot_s = envSeconds("ODIA_TRACE_SNAPSHOT_MIN_S", 30.0);
    note("snapshot_min_s", std::to_string(snapshot_s));
    if (snapshot_s >= 0.0)
    {
      auto* snap = new Snapshots;
      snap->min_interval_s = snapshot_s;
      snap->numa = openOut("numa.tsv");
      snap->smaps = openOut("smaps.tsv");
      snap->minfo = openOut("mallocinfo.tsv");
      if (snap->numa) { std::fputs("t_s\tat\tcost_ms\tnode\ttotal_kb\tanon_kb\tfile_kb\thuge_kb\n", snap->numa); }
      if (snap->smaps)
      {
        std::fputs("t_s\tat\tcost_ms\trss_kb\tpss_kb\tpss_anon_kb\tpss_file_kb\tpss_dirty_kb\t"
                   "anonymous_kb\tanon_huge_kb\tprivate_dirty_kb\tprivate_clean_kb\tshared_clean_kb\tswap_kb\n",
                   snap->smaps);
      }
      if (snap->minfo)
      {
        std::fputs("t_s\tat\tcost_ms\tarenas\tfast_free_b\trest_free_b\tmmap_count\tmmap_b\t"
                   "system_current_b\tsystem_max_b\taspace_total_b\taspace_mprotect_b\tinuse_b\n", snap->minfo);
      }
      snapshots() = snap;
    }
#endif
    std::fprintf(stderr, "[trace] writing to %s\n", s.dir.c_str());
  }

  Session::~Session()
  {
    if (!owner_) { return; }
    auto& s = state();
#if defined(__linux__)
    // The state the run leaves behind, whatever the interval.
    boundarySnapshot("end", true);
    if (Snapshots* snap = snapshots())
    {
      note("snapshots_taken", std::to_string(snap->taken));
      note("snapshots_skipped", std::to_string(snap->skipped));
      note("snapshots_total_ms", std::to_string(snap->total_ms));
    }
#endif
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
      note("fast_samples", std::to_string(smp->fast_samples));
      note("sample_ms_mean", std::to_string(smp->fast_samples ? smp->sample_ms_sum / static_cast<double>(smp->fast_samples) : 0.0));
      note("sample_ms_max", std::to_string(smp->sample_ms_max));
      sampler() = nullptr;
      delete smp;
    }
    if (Snapshots* snap = snapshots())
    {
      {
        // Released before the delete: unlocking a freed mutex writes into
        // freed heap memory.
        const std::lock_guard<std::mutex> lock(snap->mutex);
        for (FILE* f : {snap->numa, snap->smaps, snap->minfo}) { if (f) { std::fclose(f); } }
        snapshots() = nullptr;
      }
      delete snap;
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
