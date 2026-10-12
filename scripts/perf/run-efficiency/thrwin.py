#!/usr/bin/env python3
"""thrwin.py RUNDIR T0 T1 [STEP]: per-thread CPU seconds in [T0,T1] (threads.tsv utime+stime deltas), and
per STEP-second bins, the number of threads that used >= 50% of a core in the bin."""
import sys, csv, collections
d, a, b = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]); step = float(sys.argv[4]) if len(sys.argv) > 4 else 10.0
last = {}; per = collections.Counter(); comm = {}
nb = int((b - a) / step) + 1; bins = [collections.Counter() for _ in range(nb)]
with open(d + '/trace/threads.tsv') as f:
    r = csv.reader(f, delimiter='\t'); h = next(r); ix = {k: i for i, k in enumerate(h)}
    for row in r:
        t = float(row[ix['t_s']]); tid = row[ix['tid']]; c = float(row[ix['utime_s']]) + float(row[ix['stime_s']])
        comm[tid] = row[ix['comm']]
        if tid in last and a <= t <= b + step:
            pt, pc = last[tid]
            dc = c - pc
            if t <= b and pt >= a - 1:
                per[tid] += dc
                k = int((t - a) / step)
                if 0 <= k < nb: bins[k][tid] += dc
        last[tid] = (t, c)
print('top threads by CPU in window:')
for tid, c in per.most_common(20): print(tid, comm[tid], '%.1f' % c)
print('bin_start\tthreads>=50%\tcpu_s')
for k, bc in enumerate(bins):
    print('%.0f\t%d\t%.1f' % (a + k * step, sum(1 for v in bc.values() if v >= 0.5 * step), sum(bc.values())))
