#!/usr/bin/env python3
"""effan.py OUTDIR RUNDIR [RUNDIR...] -- parallel-efficiency tables of traced DIALibGen -run runs.

For every RUNDIR (ODIA_TRACE output in RUNDIR/trace, /usr/bin/time -v in RUNDIR/time.txt) writes
OUTDIR/<run>.phases.tsv and OUTDIR/<run>.steps.tsv and prints them.

phases.tsv: one row per level-0 phase (+ startup, gaps, exit, TOTAL)
  wall_s            phase duration
  cpu_s             process user+system CPU over the phase (getrusage samples, interpolated)
  cores             cpu_s / wall_s (mean active cores)
  s_le_1p5          seconds inside the phase whose sampled interval had <= 1.5 active cores
  peak_rss_gb       largest VmRSS sample in the phase (GiB)
  budget            the thread budget: -threads, or for -threads 0 the physical cores of the
                    affinity mask (topology.tsv)
  eff               cores / budget (efficiency against the thread budget)
steps.tsv: level-1/2 events summed per (phase, name, level): n, wall, process CPU over their spans
  (meaningful for steps that do not overlap within one phase), cores.
"""
import bisect, csv, os, sys

def tsv(p):
    with open(p) as f:
        r = csv.reader(f, delimiter='\t')
        h = next(r)
        return [dict(zip(h, row)) for row in r if row]

class Run:
    def __init__(self, d):
        self.d = d
        self.label = os.path.basename(os.path.normpath(d))
        tr = os.path.join(d, 'trace')
        self.meta = {}
        for r in tsv(os.path.join(tr, 'meta.tsv')):
            self.meta.setdefault(r['key'], r['value'])
        self.ev = tsv(os.path.join(tr, 'events.tsv'))
        pr = tsv(os.path.join(tr, 'proc.tsv'))
        self.t = [float(r['t_s']) for r in pr]
        self.cpu = [float(r['utime_s']) + float(r['stime_s']) for r in pr]
        self.rss = [int(r['rss_kb']) for r in pr]
        self.nrun = [int(r['n_run']) - 1 for r in pr]
        self.start = float(self.meta.get('process_start_s', 0.0))
        self.end = float(self.meta.get('end_s', self.t[-1]))
        argv = self.meta['argv'].split()
        thr = int(argv[argv.index('-threads') + 1]) if '-threads' in argv else 1
        self.threads_arg = thr
        topo = tsv(os.path.join(tr, 'topology.tsv'))
        phys = {(r['package'], r['core_id']) for r in topo if r.get('allowed', '1') == '1'}
        self.budget = thr if thr > 0 else len(phys)
        # interval active cores
        self.ac = [0.0] + [(self.cpu[i] - self.cpu[i - 1]) / (self.t[i] - self.t[i - 1]) if self.t[i] > self.t[i - 1] else 0.0
                           for i in range(1, len(self.t))]

    def cpu_at(self, x):
        t, c = self.t, self.cpu
        if x <= t[0]:
            return 0.0 if x <= self.start else c[0] * (x - self.start) / (t[0] - self.start)
        if x >= t[-1]:
            return c[-1]
        i = bisect.bisect_left(t, x)
        a, b = t[i - 1], t[i]
        return c[i - 1] + (c[i] - c[i - 1]) * (x - a) / (b - a)

    def low_seconds(self, a, b, thr=1.5):
        s = 0.0
        i0 = bisect.bisect_right(self.t, a)
        for i in range(max(i0, 1), len(self.t)):
            lo, hi = max(self.t[i - 1], a), min(self.t[i], b)
            if hi > lo and self.ac[i] <= thr:
                s += hi - lo
            if self.t[i - 1] >= b:
                break
        return s

    def peak(self, a, b):
        i0, i1 = bisect.bisect_left(self.t, a), bisect.bisect_right(self.t, b)
        if i1 <= i0:
            return self.rss[min(i0, len(self.rss) - 1)]
        return max(self.rss[i0:i1])

    def phases(self):
        ph = sorted((float(e['t0_s']), float(e['t1_s']), e['name']) for e in self.ev if e['level'] == '0')
        segs, prev, first = [], self.start, True
        for a, b, n in ph:
            if a - prev > 1e-3:
                segs.append(('startup' if first else 'gap', prev, a))
            segs.append((n, a, b)); prev = b; first = False
        if self.end - prev > 1e-3:
            segs.append(('exit', prev, self.end))
        rows = []
        for n, a, b in segs:
            w = b - a
            c = self.cpu_at(b) - (0.0 if n == 'startup' else self.cpu_at(a))
            rows.append(dict(phase=n, t0=a, t1=b, wall=w, cpu=c, cores=c / w if w > 0 else 0.0,
                             low=self.low_seconds(a, b), rss=self.peak(a, b) / 1048576.0))
        return rows

    def merged_phases(self):
        out, idx = [], {}
        for r in self.phases():
            k = r['phase']
            if k not in idx:
                g = dict(r); g['n'] = 1
                idx[k] = g; out.append(g)
            else:
                g = idx[k]
                g['wall'] += r['wall']; g['cpu'] += r['cpu']; g['low'] += r['low']; g['rss'] = max(g['rss'], r['rss']); g['n'] += 1
                g['t1'] = r['t1']
        for g in out:
            g['cores'] = g['cpu'] / g['wall'] if g['wall'] > 0 else 0.0
            if g['n'] > 1:
                g['phase'] = '%s (%d)' % (g['phase'], g['n'])
        # snapshots and gaps last
        out = [g for g in out if not g['phase'].startswith(('trace_snapshot', 'gap'))] + \
              [g for g in out if g['phase'].startswith(('trace_snapshot', 'gap'))]
        tw = self.end - self.start
        out.append(dict(phase='TOTAL', t0=self.start, t1=self.end, wall=tw, cpu=self.cpu_at(self.end),
                        cores=self.cpu_at(self.end) / tw, low=self.low_seconds(self.start, self.end),
                        rss=max(self.rss) / 1048576.0))
        for r in out:
            r['eff'] = r['cores'] / self.budget
            r['share'] = r['wall'] / tw
        return out

    def steps(self):
        agg = {}
        for e in self.ev:
            if e['level'] not in ('1', '2'):
                continue
            a, b = float(e['t0_s']), float(e['t1_s'])
            k = (e['phase'], e['name'], e['level'])
            g = agg.setdefault(k, [0, 0.0, 0.0, 0.0])
            g[0] += 1; g[1] += b - a; g[2] += self.cpu_at(b) - self.cpu_at(a); g[3] += self.low_seconds(a, b)
        return [dict(phase=k[0], name=k[1], level=k[2], n=v[0], wall=v[1], cpu=v[2], cores=v[2] / v[1] if v[1] > 0 else 0, low=v[3])
                for k, v in agg.items()]


def f(v, d=1):
    return ('%.' + str(d) + 'f') % v


def main():
    out = sys.argv[1]
    os.makedirs(out, exist_ok=True)
    for d in sys.argv[2:]:
        r = Run(d)
        rows = r.merged_phases()
        with open(os.path.join(out, r.label + '.phases.tsv'), 'w') as fo:
            fo.write('run\tphase\tt0_s\tt1_s\twall_s\tcpu_s\tcores\ts_le_1p5\tpeak_rss_gb\tbudget\teff\tshare_wall\n')
            for x in rows:
                fo.write('\t'.join([r.label, x['phase'], f(x['t0'], 2), f(x['t1'], 2), f(x['wall']), f(x['cpu']), f(x['cores'], 2),
                                    f(x['low']), f(x['rss'], 2), str(r.budget), f(x['eff'], 3), f(x['share'], 3)]) + '\n')
        st = r.steps()
        with open(os.path.join(out, r.label + '.steps.tsv'), 'w') as fo:
            fo.write('run\tphase\tname\tlevel\tn\twall_s\tcpu_s\tcores\ts_le_1p5\n')
            for x in st:
                fo.write('\t'.join([r.label, x['phase'], x['name'], x['level'], str(x['n']), f(x['wall']), f(x['cpu']), f(x['cores'], 2), f(x['low'])]) + '\n')
        print('### %s (budget %d threads)' % (r.label, r.budget))
        print('| phase | wall s | CPU s | cores | s at <=1.5 | peak RSS GiB | eff |')
        for x in rows:
            print('| %s | %s | %s | %s | %s | %s | %s |' % (x['phase'], f(x['wall']), f(x['cpu']), f(x['cores'], 2), f(x['low']), f(x['rss'], 1), f(x['eff'], 3)))
        print()

if __name__ == '__main__':
    main()
