#!/usr/bin/env python3
"""bid10-threads: speed-up, efficiency and Amdahl fits per phase on sub10.

usage: scaling_amdahl.py OUTDIR SERIES_SPEC...
  SERIES_SPEC = name=run_dir:run_dir:...   (each run dir has phase_cores-like data via trace/events.tsv + time.txt)
  DIA-NN runs are given as  dn=diann_dir:...  and only contribute total elapsed (time.txt).

Writes OUTDIR/scaling_phases.tsv  series, run, threads, sessions, phase, wall_s, speedup_vs_t4, efficiency_vs_t4
       OUTDIR/amdahl_fit.tsv      series, phase, n_points, serial_s (a), parallel_s_at_1 (b), serial_fraction_at_t4,
                                  serial_fraction_at_1, max_speedup_vs_t4 (n->inf), r2
Model: T(n) = a + b/n  (n = threads given to -threads); least squares on all points of the series.
serial_fraction_at_t4 = a / T_fit(4); serial_fraction_at_1 = a / (a + b).
scaling_phases.tsv also has speed-up/efficiency vs 8 threads (the DIA-NN series starts at 8).
"""
import csv, os, re, sys
from collections import defaultdict

def elapsed(run_dir):
    for line in open(os.path.join(run_dir, 'time.txt')):
        if 'Elapsed' in line:
            v = line.split(': ', 1)[1].strip().split(':')
            return sum(float(x) * 60 ** (len(v) - 1 - i) for i, x in enumerate(v))
    return float('nan')

def threads_of(run_dir):
    meta = open(os.path.join(run_dir, 'meta.txt')).read()
    m = re.search(r'(?:-threads|--threads) (\d+)', meta)
    return int(m.group(1))

def phases(run_dir):
    p = defaultdict(float)
    ev = os.path.join(run_dir, 'trace', 'events.tsv')
    if not os.path.exists(ev): return p
    with open(ev) as f:
        r = csv.DictReader(f, delimiter='\t')
        for e in r:
            if e['level'] == '0':
                p[e['name']] += float(e['dur_s'])
    out = {}
    out['serial_pre (fasta_read+digest+enumerate+release)'] = sum(p.get(k, 0) for k in ('fasta_read', 'digest', 'enumerate', 'release'))
    out['rt'] = p.get('rt', 0); out['ms2'] = p.get('ms2', 0); out['ccs'] = p.get('ccs', 0)
    out['write'] = p.get('write', 0)
    return out

def fit(xs, ys):
    # y = a + b * (1/x)
    us = [1.0 / x for x in xs]
    n = len(us); mu = sum(us) / n; my = sum(ys) / n
    suu = sum((u - mu) ** 2 for u in us); suy = sum((u - mu) * (y - my) for u, y in zip(us, ys))
    b = suy / suu if suu else 0.0; a = my - b * mu
    ss = sum((y - (a + b * u)) ** 2 for u, y in zip(us, ys)); st = sum((y - my) ** 2 for y in ys)
    return a, b, (1 - ss / st) if st else float('nan')

def main():
    od = sys.argv[1]; os.makedirs(od, exist_ok=True)
    rows = []; fits = []
    for spec in sys.argv[2:]:
        name, dirs = spec.split('=', 1)
        pts = []
        for d in dirs.split(':'):
            n = threads_of(d); ph = phases(d); ph['TOTAL (elapsed)'] = elapsed(d)
            sess = ''
            mt = os.path.join(d, 'trace', 'meta.tsv')
            if os.path.exists(mt):
                for line in open(mt):
                    if line.startswith('predictor.MS2.sessions'): sess = line.split('\t')[1].strip()
            pts.append((n, os.path.basename(d.rstrip('/')), sess, ph))
        pts.sort()
        b4 = next((p[3] for p in pts if p[0] == 4), None); b8 = next((p[3] for p in pts if p[0] == 8), None)
        for n, run, sess, ph in pts:
            for k, v in ph.items():
                sp = b4[k] / v if b4 and v > 0 and b4.get(k) else ''
                eff = sp * 4 / n if sp != '' else ''
                sp8 = b8[k] / v if b8 and v > 0 and b8.get(k) else ''
                rows.append([name, run, n, sess, k, round(v, 3), round(sp, 3) if sp != '' else '', round(eff, 3) if eff != '' else '',
                             round(sp8, 3) if sp8 != '' else '', round(sp8 * 8 / n, 3) if sp8 != '' else ''])
        keys = pts[0][3].keys()
        for k in keys:
            xs = [p[0] for p in pts if p[3].get(k, 0) > 0]; ys = [p[3][k] for p in pts if p[3].get(k, 0) > 0]
            if len(xs) < 3: continue
            a, b, r2 = fit(xs, ys)
            t4 = a + b / 4
            fits.append([name, k, len(xs), round(a, 3), round(b, 2), round(a / t4, 4) if t4 else '', round(a / (a + b / 8), 4),
                         round(a / (a + b), 4) if a + b else '', round(t4 / a, 2) if a > 0 else 'inf', round(r2, 4)])
    with open(os.path.join(od, 'scaling_phases.tsv'), 'w') as f:
        f.write('series\trun\tthreads\tsessions\tphase\twall_s\tspeedup_vs_t4\tefficiency_vs_t4\tspeedup_vs_t8\tefficiency_vs_t8\n')
        for r in rows: f.write('\t'.join(map(str, r)) + '\n')
    with open(os.path.join(od, 'amdahl_fit.tsv'), 'w') as f:
        f.write('series\tphase\tn_points\tserial_s\tparallel_s_at_1\tserial_fraction_at_t4\tserial_fraction_at_t8\tserial_fraction_at_1\tmax_speedup_vs_t4\tr2\n')
        for r in fits: f.write('\t'.join(map(str, r)) + '\n')
    print('rows', len(rows), 'fits', len(fits))

if __name__ == '__main__':
    main()
