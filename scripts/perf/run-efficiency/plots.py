#!/usr/bin/env python3
"""plots.py: figures for the run-efficiency report.

  plots.py timeline OUT.png RUNDIR [TITLE]     active cores (1 s bins) and RSS, phase bands, thread budget
  plots.py cpuwall OUT.png RUNDIR [RUNDIR2]   CPU vs wall per phase (log-log), iso-core lines
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from effan import Run

# reference categorical palette (light), fixed order
PAL = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#4a3aa7', '#e34948']
GROUPS = [
    ('load and hash (serial I/O)', {'library_load', 'run_load', 'run_hash', 'provenance', 'cache_remove', 'startup'}),
    ('prefilter (MS2 model)', {'prefilter'}),
    ('calibration', {'calibration'}),
    ('extraction', {'extraction'}),
    ('scoring, report, apply', {'scoring', 'report', 'report_write', 'apply', 'tune_check'}),
    ('tuning (LibTorch)', {'tune_rt', 'tune_ccs'}),
    ('re-prediction', {'rt', 'ccs'}),
    ('write', {'write'}),
]
INK, INK2, GRID = '#0b0b0b', '#52514e', '#e4e3df'


def group_of(name):
    for i, (g, names) in enumerate(GROUPS):
        if name in names:
            return i
    return None


def style(ax):
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    for s in ('left', 'bottom'):
        ax.spines[s].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=8)
    ax.grid(axis='y', color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)


def timeline(out, d, title=None, xmax=None):
    r = Run(d)
    # 1 s bins
    t, c = r.t, r.cpu
    xs, ys, rs = [], [], []
    j = 0
    T = int(t[-1])
    for s in range(0, T):
        a, b = float(s), float(s + 1)
        xs.append((a + 0.5) / 60.0)
        ys.append(r.cpu_at(b) - r.cpu_at(a))
    import bisect
    rx = [x / 60.0 for x in t[::10]]
    ry = [v / 1048576.0 for v in r.rss[::10]]
    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(11, 5.6), sharex=True, gridspec_kw=dict(height_ratios=[3, 1.1], hspace=0.08))
    used = set()
    for p in r.phases():
        g = group_of(p['phase'])
        if g is None or p['wall'] < 0.5:
            continue
        for a_ in (ax, ax2):
            a_.axvspan(p['t0'] / 60.0, p['t1'] / 60.0, color=PAL[g], alpha=0.16, linewidth=0,
                       label=GROUPS[g][0] if (a_ is ax and g not in used) else None)
        used.add(g)
    ax.plot(xs, ys, color=INK, linewidth=0.7)
    ax.axhline(r.budget, color=INK2, linewidth=1, linestyle='--')
    ax.text(0.2, r.budget, ' thread budget %d' % r.budget, va='bottom', fontsize=8, color=INK2)
    ax.set_ylabel('active cores (1 s)', fontsize=9, color=INK2)
    ax.set_ylim(0, max(max(ys), r.budget) * 1.08)
    ax2.plot(rx, ry, color=PAL[6], linewidth=1)
    ax2.set_ylabel('RSS (GiB)', fontsize=9, color=INK2)
    ax2.set_xlabel('time since start (min)', fontsize=9, color=INK2)
    for a_ in (ax, ax2):
        style(a_)
    if xmax:
        ax.set_xlim(0, xmax)
    h, l = ax.get_legend_handles_labels()
    order = sorted(range(len(l)), key=lambda i: [g for g, _ in GROUPS].index(l[i]))
    ax.legend([h[i] for i in order], [l[i] for i in order], fontsize=7.5, ncol=4, frameon=False, loc='upper center',
              bbox_to_anchor=(0.5, 1.17))
    ax.set_title(title or r.label, fontsize=10, color=INK, pad=34, loc='left')
    fig.savefig(out, dpi=130, bbox_inches='tight')


def cpuwall(out, dirs):
    fig, ax = plt.subplots(figsize=(7.5, 6))
    markers = ['o', 's']
    for k, d in enumerate(dirs):
        r = Run(d)
        for p in r.phases():
            g = group_of(p['phase'])
            if g is None or p['wall'] < 1.0:
                continue
            ax.scatter(p['wall'], p['cpu'], s=46, marker=markers[k], color=PAL[g], edgecolor='white', linewidth=1.2, zorder=3)
            show = (k == 0 and p['wall'] > 30) or (k == 1 and p['phase'] in ('prefilter', 'run_load', 'rt', 'ccs'))
            if show:
                ax.annotate(p['phase'] + (' (t0)' if k else ''), (p['wall'], p['cpu']), fontsize=7, color=INK2,
                            xytext=(5, -10 if k == 0 else 4), textcoords='offset points')
    lo, hi = 1, 4000
    for n in (1, 4, 8, 16, 192):
        ax.plot([lo, hi], [lo * n, hi * n], color=GRID if n not in (1, 16) else '#c3c2b7', linewidth=0.9, zorder=1)
        xt = min(hi * 0.93, 150000.0 / n)
        ax.text(xt, xt * n * 0.8, '%d cores' % n, fontsize=7, color=INK2, ha='right')
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlim(1, 4000); ax.set_ylim(1, 200000)
    ax.set_xlabel('wall time per phase (s, log)', fontsize=9, color=INK2)
    ax.set_ylabel('CPU time per phase (s, log)', fontsize=9, color=INK2)
    style(ax)
    from matplotlib.lines import Line2D
    hs = [Line2D([], [], marker='o', linestyle='', color=PAL[i], label=g) for i, (g, _) in enumerate(GROUPS)]
    hs += [Line2D([], [], marker=m, linestyle='', color=INK2, label=l) for m, l in zip(markers, ['-threads 16', '-threads 0 (192)'][:len(dirs)])]
    ax.legend(handles=hs, fontsize=7, frameon=False, loc='upper left')
    ax.set_title('CPU against wall per phase (full timsTOF run)', fontsize=10, color=INK, loc='left')
    fig.savefig(out, dpi=130, bbox_inches='tight')


def phasecmp(out, phase, specs, title):
    """Active cores (1 s bins) over one phase for several runs; specs = [(rundir, label), ...]."""
    fig, ax = plt.subplots(figsize=(10, 3.8))
    for k, (d, lab) in enumerate(specs):
        r = Run(d)
        ph = [p for p in r.phases() if p['phase'] == phase][0]
        xs, ys = [], []
        s = ph['t0']
        while s < ph['t1']:
            xs.append(s - ph['t0'] + 0.5); ys.append(r.cpu_at(min(s + 1, ph['t1'])) - r.cpu_at(s)); s += 1
        ax.plot(xs, ys, color=PAL[k], linewidth=1.1, label='%s: %.0f s, %.1f cores' % (lab, ph['wall'], ph['cores']))
    ax.set_xlabel('seconds since %s started' % phase, fontsize=9, color=INK2)
    ax.set_ylabel('active cores (1 s)', fontsize=9, color=INK2)
    style(ax)
    ax.legend(fontsize=7.5, frameon=False, loc='upper right')
    ax.set_title(title, fontsize=10, color=INK, loc='left')
    fig.savefig(out, dpi=130, bbox_inches='tight')


if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'timeline':
        timeline(sys.argv[2], sys.argv[3], sys.argv[4] if len(sys.argv) > 4 else None)
    elif cmd == 'phasecmp':
        specs = [tuple(x.split('=', 1)) for x in sys.argv[5:]]
        phasecmp(sys.argv[2], sys.argv[3], [(a, b) for b, a in specs], sys.argv[4])
    elif cmd == 'cpuwall':
        cpuwall(sys.argv[2], sys.argv[3:])
