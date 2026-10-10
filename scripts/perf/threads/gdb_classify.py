#!/usr/bin/env python3
"""Classify gdb all-thread backtraces (gdbsample.sh output) by phase.

usage: gdb_classify.py OUT.tsv DIAGDIR...
Output columns: diag, phase, samples (stops in that phase), category, thread_obs (thread x stop observations),
threads_per_stop (thread_obs / samples), share (of all thread observations in the phase)
Categories are decided from the stack (first matching rule, see classify()).
"""
import csv, glob, os, re, sys
from collections import defaultdict

def phases(diag):
    ev = os.path.join(diag, 'trace', 'events.tsv'); ph = []
    with open(ev) as f:
        for e in csv.DictReader(f, delimiter='\t'):
            if e['level'] == '0': ph.append((float(e['t0_s']), float(e['t1_s']), e['name']))
    meta = {}
    for line in open(os.path.join(diag, 'trace', 'meta.tsv')):
        k, _, v = line.rstrip('\n').partition('\t'); meta[k] = v
    return sorted(ph), float(meta['start_epoch_s'])

def phase_at(ph, t):
    for a, b, n in ph:
        if a <= t < b: return n
    return 'gap/startup/exit'

ALLOC = re.compile(r'\b(malloc|free|_int_malloc|_int_free|calloc|realloc|operator new|operator delete|je_|arena_|tcache|sysmalloc|malloc_consolidate|cfree|__libc_malloc|__libc_free|unlink_chunk)')
WAIT = re.compile(r'(futex|pthread_cond_wait|pthread_cond_timedwait|__poll|epoll_wait|nanosleep|clock_nanosleep|sem_wait|pthread_mutex_lock|__lll_lock_wait|sched_yield|pthread_join|__pthread_clockjoin)')
SYSC = re.compile(r'\b(__GI___libc_write|__libc_write|write|__GI___libc_read|read|mmap64|__mmap|munmap|__munmap|madvise|__madvise|openat|__open|clone3?|__clone|mprotect|__GI___mprotect|fsync|ftruncate|pwrite64)\b')

def classify(name, frames):
    s = '\n'.join(frames)
    top = frames[0] if frames else ''
    if 'blas_thread_server' in s: return 'idle pool: OpenBLAS (from CoinUtils/LAPACK)'
    if 'background_thread' in s: return 'idle pool: jemalloc background (Arrow)'
    if name.startswith('odia-sampler') or 'trace::' in s and 'sampler' in s.lower(): return 'tracer sampler'
    if 'gdb' in name: return 'other'
    # spinning ORT pool threads (allow_spinning): busy in the pool's WorkerLoop without a task
    if ('WorkerLoop' in s or 'ThreadPoolTempl' in s) and not re.search(r'Mlas|Compute|ExecuteKernel', s):
        if WAIT.search(top) or WAIT.search('\n'.join(frames[:3])): return 'ORT intra-op: blocked in pool (idle, no CPU)'
        return 'ORT intra-op: spin-wait (pool idle or end-of-section barrier; burns CPU)'
    if WAIT.search(top) or WAIT.search('\n'.join(frames[:2])):
        if re.search(r'thread::join|pthread_join|clockjoin', s): return 'waiting: join (caller waits for predictor workers)'
        if 'onnxruntime' in s: return 'waiting: inside ORT (barrier / pool sync)'
        return 'waiting: other (cond/futex/sleep)'
    if ALLOC.search('\n'.join(frames[:4])): return 'allocator (malloc/free)'
    if SYSC.search(top): return 'syscall (write/mmap/munmap/clone)'
    if re.search(r'Mlas\w*(Kernel|Gemm|Sgemm|Conv)', s): return 'ORT compute: MLAS GEMM'
    if 'onnxruntime' in s:
        if re.search(r'Compute|Mlas|ExecuteKernel', s): return 'ORT compute: other op (non-GEMM)'
        return 'ORT: session overhead (Run setup, allocation planning, copy)'
    if re.search(r'arrow|parquet', s, re.I): return 'writer: Arrow/Parquet'
    if re.search(r'PeptDeepEncoder|encode', s): return 'DIALibGen: encode (tokenise)'
    if re.search(r'OpenMS::', s):
        if 'AASequence' in s: return 'OpenMS: AASequence parse/ops'
        return 'OpenMS: other (digest, chemistry)'
    if re.search(r'ODIA::', s): return 'DIALibGen: own code (enumerate, post, write prep)'
    return 'other/unresolved'

def parse(path):
    th = []; cur = None
    for line in open(path):
        if line.startswith('Thread '):
            m = re.search(r'"([^"]*)"\)', line)
            cur = [m.group(1) if m else '?', []]; th.append(cur)
        elif line.startswith('#') and cur is not None:
            m = re.match(r'#\d+\s+(?:0x[0-9a-f]+ in )?(.*)', line.strip())
            cur[1].append(m.group(1) if m else line.strip())
    return th

def main():
    out = sys.argv[1]; rows = []
    for diag in sys.argv[2:]:
        ph, e0 = phases(diag)
        cnt = defaultdict(int); ns = defaultdict(int)
        for f in sorted(glob.glob(os.path.join(diag, 's_*.txt'))):
            t = float(os.path.basename(f)[2:-4]) - e0
            p = phase_at(ph, t); ns[p] += 1
            for name, frames in parse(f):
                cnt[(p, classify(name, frames))] += 1
        tot = defaultdict(int)
        for (p, c), n in cnt.items(): tot[p] += n
        for (p, c), n in sorted(cnt.items(), key=lambda x: (x[0][0], -x[1])):
            rows.append([os.path.basename(diag.rstrip('/')), p, ns[p], c, n, round(n / ns[p], 2), round(n / tot[p], 4)])
    with open(out, 'w') as f:
        f.write('diag\tphase\tsamples\tcategory\tthread_obs\tthreads_per_stop\tshare\n')
        for r in rows: f.write('\t'.join(map(str, r)) + '\n')

if __name__ == '__main__':
    main()
