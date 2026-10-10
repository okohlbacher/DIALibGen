// allocshim.so -- LD_PRELOAD allocation counter (perf/memory-cache, bid10).
// Counts malloc/calloc/realloc/memalign-family calls and requested bytes, and frees, per thread and per
// 1-second bucket since process start, plus a log2 request-size histogram. No locks on the hot path:
// each thread owns a record (mmap'd, never freed) pushed once onto a lock-free list; the dump at exit
// sums the records. Writes $ALLOCSHIM_OUT (default /tmp/allocshim.<pid>.tsv).
// Output sections (tab-separated, '#' comments):
//   bucket  t_s  mallocs  bytes  frees  threads_active   -- per 1-s bucket, all threads summed
//   size    log2_bin  calls  bytes                        -- request size in [2^bin, 2^(bin+1))
//   thread  tid  first_s  last_s  mallocs  bytes  frees   -- per thread
#define _GNU_SOURCE
#include <dlfcn.h>
#include <stdatomic.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/syscall.h>
#include <time.h>
#include <unistd.h>

#define NB 4096
#define NS 48
extern void* __libc_malloc(size_t);
extern void* __libc_calloc(size_t, size_t);
extern void* __libc_realloc(void*, size_t);
extern void __libc_free(void*);
extern void* __libc_memalign(size_t, size_t);

typedef struct Rec {
  struct Rec* next;
  long tid;
  uint32_t first, last;
  uint64_t m[NB], b[NB], f[NB];
  uint64_t sc[NS], sb[NS];
} Rec;

static _Atomic(Rec*) head = NULL;
static double t0 = -1;
static __thread Rec* me = NULL;
static __thread int busy = 0;

static inline double mono(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC_COARSE, &t); return t.tv_sec + 1e-9 * t.tv_nsec; }

static Rec* rec(void) {
  if (me) return me;
  if (busy) return NULL;
  busy = 1;
  Rec* r = mmap(NULL, sizeof(Rec), PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
  if (r == MAP_FAILED) { busy = 0; return NULL; }
  r->tid = syscall(SYS_gettid); r->first = UINT32_MAX;
  Rec* h = atomic_load(&head);
  do { r->next = h; } while (!atomic_compare_exchange_weak(&head, &h, r));
  me = r; busy = 0;
  return r;
}

static inline unsigned bucket(void) {
  if (t0 < 0) t0 = mono();
  double d = mono() - t0; unsigned k = d < 0 ? 0 : (unsigned)d;
  return k < NB ? k : NB - 1;
}
static inline void count_alloc(size_t n) {
  Rec* r = rec(); if (!r) return;
  unsigned k = bucket();
  r->m[k]++; r->b[k] += n;
  if (k < r->first) r->first = k;
  if (k > r->last) r->last = k;
  unsigned s = n ? 63 - __builtin_clzl(n) : 0; if (s >= NS) s = NS - 1;
  r->sc[s]++; r->sb[s] += n;
}
static inline void count_free(void* p) { if (!p) return; Rec* r = rec(); if (!r) return; r->f[bucket()]++; }

void* malloc(size_t n) { count_alloc(n); return __libc_malloc(n); }
void* calloc(size_t a, size_t n) { count_alloc(a * n); return __libc_calloc(a, n); }
void* realloc(void* p, size_t n) { count_alloc(n); if (p) count_free(p); return __libc_realloc(p, n); }
void free(void* p) { count_free(p); __libc_free(p); }
void* memalign(size_t al, size_t n) { count_alloc(n); return __libc_memalign(al, n); }
void* aligned_alloc(size_t al, size_t n) { count_alloc(n); return __libc_memalign(al, n); }
int posix_memalign(void** out, size_t al, size_t n) {
  count_alloc(n); void* p = __libc_memalign(al, n); if (!p) return 12; *out = p; return 0;
}
void* valloc(size_t n) { count_alloc(n); return __libc_memalign(4096, n); }

__attribute__((constructor)) static void init(void) { t0 = mono(); }

__attribute__((destructor)) static void dump(void) {
  char path[512]; const char* o = getenv("ALLOCSHIM_OUT");
  if (o && *o) snprintf(path, sizeof path, "%s", o); else snprintf(path, sizeof path, "/tmp/allocshim.%d.tsv", getpid());
  busy = 1;
  FILE* f = fopen(path, "w"); if (!f) return;
  static uint64_t M[NB], B[NB], F[NB], T[NB], SC[NS], SB[NS];
  unsigned maxk = 0;
  fprintf(f, "# allocshim: per-thread counters summed at exit; t_s = seconds since the shim constructor (CLOCK_MONOTONIC_COARSE)\n");
  fprintf(f, "#section\ttid\tfirst_s\tlast_s\tmallocs\tbytes\tfrees\n");
  for (Rec* r = atomic_load(&head); r; r = r->next) {
    uint64_t m = 0, b = 0, fr = 0;
    for (unsigned k = 0; k < NB; ++k) {
      if (r->m[k] || r->f[k]) { M[k] += r->m[k]; B[k] += r->b[k]; F[k] += r->f[k]; T[k]++; if (k > maxk) maxk = k; }
      m += r->m[k]; b += r->b[k]; fr += r->f[k];
    }
    for (unsigned s = 0; s < NS; ++s) { SC[s] += r->sc[s]; SB[s] += r->sb[s]; }
    fprintf(f, "thread\t%ld\t%u\t%u\t%lu\t%lu\t%lu\n", r->tid, r->first, r->last, m, b, fr);
  }
  fprintf(f, "#section\tt_s\tmallocs\tbytes\tfrees\tthreads_active\n");
  for (unsigned k = 0; k <= maxk; ++k) fprintf(f, "bucket\t%u\t%lu\t%lu\t%lu\t%lu\n", k, M[k], B[k], F[k], T[k]);
  fprintf(f, "#section\tlog2_bin\tcalls\tbytes\n");
  for (unsigned s = 0; s < NS; ++s) if (SC[s]) fprintf(f, "size\t%u\t%lu\t%lu\n", s, SC[s], SB[s]);
  fclose(f);
}
