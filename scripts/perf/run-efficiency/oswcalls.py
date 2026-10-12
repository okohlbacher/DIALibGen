#!/usr/bin/env python3
"""oswcalls.py RUNDIR [TEAM]: per OpenSWATH call (calibration and extraction 'openswath' spans):
wall, process CPU, mean active cores, tail = time from the last moment >= 0.75*TEAM threads running to the end
(share of the call), time with <= 1.5 active cores, and seconds at >= k running threads."""
import sys, os, bisect, csv
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from effan import Run
d = sys.argv[1]; team = int(sys.argv[2]) if len(sys.argv) > 2 else 8
r = Run(d)
print('run\tphase\tidx\trows\tteam\twall_s\tcpu_s\tcores\tfull_until_s\ttail_s\ttail_share\ts_le_1p5\tidle_core_s')
for e in r.ev:
    if e['name'] != 'openswath':
        continue
    a, b = float(e['t0_s']), float(e['t1_s'])
    tm = int(e['session']) if e['session'] not in ('-1', '') else team
    i0, i1 = bisect.bisect_right(r.t, a), bisect.bisect_left(r.t, b)
    full = a
    for i in range(i0, i1):
        if r.ac[i] >= 0.75 * tm:
            full = r.t[i]
    w = b - a; c = r.cpu_at(b) - r.cpu_at(a)
    idle = tm * w - c
    print('\t'.join([r.label, e['phase'], e['idx'], e['rows'], str(tm), '%.1f' % w, '%.1f' % c, '%.2f' % (c / w), '%.1f' % (full - a), '%.1f' % (b - full), '%.2f' % ((b - full) / w), '%.1f' % r.low_seconds(a, b), '%.1f' % idle]))
