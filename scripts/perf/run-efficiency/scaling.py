#!/usr/bin/env python3
"""scaling.py OUTDIR CELL:RUNDIR[,RUNDIR...] ...  -- per-phase scaling of traced runs.

Each CELL lists its runs (any order); the thread count is read from the trace (-threads; 0 -> physical cores of the
affinity mask). Writes OUTDIR/scaling.tsv (cell, group, p, wall, cpu, cores, speed-up and efficiency against the
4-thread run, Amdahl fit T(p) = Ts + Tpar/p_eff over all p by least squares, p_eff = min(p, cap) for the phases
whose OpenSWATH team is capped at 8, with f = Ts/(Ts+Tpar) the serial fraction of a one-thread run), and
OUTDIR/scaling.png (per-cell panels: wall per phase group against threads, log-log).
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from effan import Run

GROUPS = [
    ('load+hash', {'library_load', 'run_load', 'run_hash', 'provenance', 'cache_remove', 'startup', 'exit', 'gap', 'trace_snapshot'}),
    ('prefilter', {'prefilter'}),
    ('calibration', {'calibration'}),
    ('extraction', {'extraction'}),
    ('after search', {'scoring', 'report', 'report_write', 'apply', 'tune_check', 'tune_rt', 'tune_ccs', 'rt', 'ccs', 'write'}),
]
PAL = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#4a3aa7']
# phases whose OpenMP team is capped (SearchParams::openswath_max_threads = 8): the fit uses p_eff = min(p, cap)
CAP = {'calibration': 8, 'extraction': 8}


def groups(r):
    out = {g: [0.0, 0.0] for g, _ in GROUPS}
    out['TOTAL'] = [0.0, 0.0]
    for p in r.phases():
        name = p['phase']
        for g, names in GROUPS:
            if name in names:
                out[g][0] += p['wall']; out[g][1] += p['cpu']
        out['TOTAL'][0] += p['wall']; out['TOTAL'][1] += p['cpu']
    return out


def amdahl(ps, ts):
    # least squares for T = a + b/p
    n = len(ps)
    xs = [1.0 / p for p in ps]
    mx, my = sum(xs) / n, sum(ts) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    if sxx == 0:
        return my, 0.0
    b = sum((x - mx) * (y - my) for x, y in zip(xs, ts)) / sxx
    a = my - b * mx
    return a, b


def main():
    out = sys.argv[1]
    cells = []
    for arg in sys.argv[2:]:
        name, runs = arg.split(':', 1)
        rs = []
        for d in runs.split(','):
            r = Run(d)
            rs.append((r.budget, r.threads_arg, r, groups(r)))
        rs.sort(key=lambda x: x[0])
        cells.append((name, rs))
    os.makedirs(out, exist_ok=True)
    lines = ['cell\tgroup\tthreads_arg\tp\twall_s\tcpu_s\tcores\tspeedup_vs4\teff_vs4\tamdahl_Ts\tamdahl_Tpar\tamdahl_f\tamdahl_fit_s']
    for name, rs in cells:
        base = {g: v for (p, a, r, gv) in rs if p == 4 for g, v in gv.items()}
        for g in [x for x, _ in GROUPS] + ['TOTAL']:
            ps = [min(p, CAP.get(g, p)) for p, a, r, gv in rs]
            ws = [gv[g][0] for p, a, r, gv in rs]
            A, B = amdahl(ps, ws)
            f = A / (A + B) if (A + B) > 0 else float('nan')
            for (p, a, r, gv) in rs:
                w, c = gv[g]
                s = base[g][0] / w if g in base and w > 0 else float('nan')
                e = s * 4 / p if s == s else float('nan')
                pe = min(p, CAP.get(g, p))
                lines.append('\t'.join([name, g, str(a), str(p), '%.1f' % w, '%.1f' % c, '%.2f' % (c / w if w else 0), '%.2f' % s, '%.3f' % e,
                                        '%.1f' % A, '%.1f' % B, '%.3f' % f, '%.1f' % (A + B / pe)]))
    open(os.path.join(out, 'scaling.tsv'), 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, len(cells), figsize=(4.6 * len(cells), 4.2), squeeze=False)
    for ax, (name, rs) in zip(axes[0], cells):
        ps = [p for p, a, r, gv in rs]
        for k, (g, _) in enumerate(GROUPS):
            ws = [gv[g][0] for p, a, r, gv in rs]
            if max(ws) < 5:
                continue
            ax.plot(ps, ws, marker='o', markersize=5, linewidth=2, color=PAL[k], label=g)
        tot = [gv['TOTAL'][0] for p, a, r, gv in rs]
        ax.plot(ps, tot, marker='s', markersize=5, linewidth=2, color='#0b0b0b', label='total')
        if 4 in ps:
            t4 = tot[ps.index(4)]
            ax.plot(ps, [t4 * 4 / p for p in ps], linestyle=':', color='#8a8984', linewidth=1, label='ideal (from 4)')
        ax.set_xscale('log', base=2); ax.set_yscale('log')
        ax.set_xticks(ps); ax.set_xticklabels([('%d' % p) + ('\n(-threads 0)' if a == 0 else '') for p, a, r, gv in rs], fontsize=7)
        ax.set_title(name, fontsize=10, loc='left')
        ax.set_xlabel('threads', fontsize=8); ax.set_ylabel('wall (s, log)', fontsize=8)
        ax.grid(color='#e4e3df', linewidth=0.6); ax.set_axisbelow(True)
        for s in ('top', 'right'):
            ax.spines[s].set_visible(False)
        ax.tick_params(labelsize=7)
    axes[0][0].legend(fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(out, 'scaling.png'), dpi=130, bbox_inches='tight')


if __name__ == '__main__':
    main()
