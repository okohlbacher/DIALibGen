#!/usr/bin/env python3
"""bid10-threads plots. usage: plots.py RESULTS_DIR  (reads the TSVs written by threads_analysis.py / scaling_amdahl.py
and experiments.tsv; writes PNGs next to them)."""
import os, sys
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

R = sys.argv[1]
C = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#4a3aa7', '#e34948']
PH = {'digest': '#d9d9d9', 'enumerate': '#bdbdbd', 'rt': '#cfe2f7', 'ms2': '#fde0cf', 'ccs': '#d3f0e4', 'write': '#efe3c2'}
plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False, 'axes.grid': True,
                     'grid.color': '#e6e6e6', 'grid.linewidth': 0.6})

def shade(ax, df):
    seen = set()
    for ph, g in df.groupby((df.phase != df.phase.shift()).cumsum()):
        p = g.phase.iloc[0]
        if p in PH:
            ax.axvspan(g.t_s.iloc[0] - g.dt_s.iloc[0], g.t_s.iloc[-1], color=PH[p], lw=0, zorder=0,
                       label=p if p not in seen else None)
            seen.add(p)

def timeline(runs, fname, title):
    occ = pd.read_csv(os.path.join(R, 'occupancy_timeline.tsv'), sep='\t')
    runs = [r for r in runs if r in set(occ.run)]
    if not runs: return
    fig, axes = plt.subplots(len(runs), 1, figsize=(10, 1.9 * len(runs)), sharex=True, squeeze=False)
    for ax, run in zip(axes[:, 0], runs):
        d = occ[occ.run == run].copy()
        d['c'] = d.cores.rolling(5, center=True, min_periods=1).mean()
        shade(ax, d)
        ax.plot(d.t_s, d.c, color=C[0], lw=1.2)
        ax.set_ylabel('active cores'); ax.set_title(run, loc='left', fontsize=9)
        ax.set_ylim(0, max(1, d.c.max() * 1.1))
    axes[0, 0].legend(loc='upper right', ncol=6, fontsize=7, frameon=False)
    axes[-1, 0].set_xlabel('time since trace start (s)')
    fig.suptitle(title, x=0.01, ha='left', fontsize=10)
    fig.tight_layout(); fig.savefig(os.path.join(R, fname), dpi=130); plt.close(fig)

timeline(['full_t8_pq', 'full_t16_pq', 'full_t96_pq'], 'occupancy_full.png',
         'Core occupancy (CPU-s per s, 0.5-s smoothing), full proteome, today\'s build (16-session cap)')
timeline(['sub_t96', 'u_t96', 'u_t96_r256', 'u_t96_r256_b100k', 'k_t96', 's16i6_spin1'], 'occupancy_sub10.png',
         'Core occupancy on sub10 at 96 threads: capped, uncapped, 256-row chunks, + 100k blocks, tuned (64 rows), 16 sessions x 6 intra-op')
timeline(['full_t96_pq', 'f_k_t96', 'f_k_t192'], 'occupancy_full_tuned.png',
         'Core occupancy, full proteome: today (16-session cap) vs tuned knobs (sessions = threads, 64-row chunks, 100k MS2 blocks)')

# scaling
p = os.path.join(R, 'scaling_phases.tsv')
if os.path.exists(p):
    s = pd.read_csv(p, sep='\t')
    tot = s[s.phase == 'TOTAL (elapsed)']
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
    for i, (ser, g) in enumerate(tot.groupby('series', sort=False)):
        g = g.sort_values('threads')
        ax[0].plot(g.threads, g.wall_s, marker='o', ms=4, lw=1.6, color=C[i], label=ser)
        if g.speedup_vs_t4.notna().any():
            ax[1].plot(g.threads, g.speedup_vs_t4, marker='o', ms=4, lw=1.6, color=C[i], label=ser)
    xs = [4, 8, 16, 32, 64, 96, 192]
    ax[1].plot(xs, [x / 4 for x in xs], ls='--', color='#999999', lw=1, label='ideal')
    for a in ax:
        a.set_xscale('log', base=2); a.set_xticks(xs); a.set_xticklabels(xs); a.set_xlabel('-threads')
    ax[0].set_yscale('log'); ax[0].set_ylabel('elapsed (s), sub10'); ax[1].set_ylabel('speed-up vs 4 threads')
    ax[1].set_ylim(0, 16)
    ax[0].legend(frameon=False, fontsize=7); ax[1].legend(frameon=False, fontsize=7)
    fig.tight_layout(); fig.savefig(os.path.join(R, 'scaling.png'), dpi=130); plt.close(fig)

# experiments bar chart
p = os.path.join(R, 'experiments.tsv')
if os.path.exists(p):
    e = pd.read_csv(p, sep='\t')
    e = e[e.input == 'sub10']
    fig, ax = plt.subplots(1, 2, figsize=(11, 0.32 * len(e) + 1.2), sharey=True)
    y = range(len(e))
    ax[0].barh(y, e.elapsed_s, color=C[0], height=0.7); ax[0].set_xlabel('elapsed (s)')
    ax[1].barh(y, e.cpu_s, color=C[1], height=0.7); ax[1].set_xlabel('user+sys CPU (s)')
    ax[0].set_yticks(list(y)); ax[0].set_yticklabels(e.label); ax[0].invert_yaxis()
    for i, v in enumerate(e.elapsed_s): ax[0].text(v, i, f' {v:.0f}', va='center', fontsize=7)
    for i, v in enumerate(e.cpu_s): ax[1].text(v, i, f' {v:.0f}', va='center', fontsize=7)
    fig.suptitle('sub10 experiments, one run each (traced knob build; baselines: trace build)', x=0.01, ha='left', fontsize=10)
    fig.tight_layout(); fig.savefig(os.path.join(R, 'experiments.png'), dpi=130); plt.close(fig)

# MS2 dispatch: wall vs largest chunk per call
p = os.path.join(R, 'dispatch_calls.tsv')
if os.path.exists(p):
    d = pd.read_csv(p, sep='\t'); d = d[(d.model == 'MS2')]
    runs = [r for r in ['full_t96_pq', 'sub_t96', 'u_t96', 'u_t96_r256', 'u_t96_r256_b100k'] if r in set(d.run)]
    fig, ax = plt.subplots(figsize=(6, 4))
    for i, run in enumerate(runs):
        g = d[d.run == run]
        ax.scatter(g.max_chunk_s, g.wall_s, s=9, color=C[i], label=run, alpha=0.8)
    m = d.wall_s.max(); ax.plot([0, m], [0, m], ls='--', color='#999999', lw=1)
    ax.set_xlabel('largest chunk of the call (s)'); ax.set_ylabel('dispatch wall of the call (s)')
    ax.set_title('MS2 calls: each call waits for its largest chunk', loc='left', fontsize=9)
    ax.legend(frameon=False, fontsize=7); fig.tight_layout(); fig.savefig(os.path.join(R, 'ms2_calls.png'), dpi=130); plt.close(fig)
print('plots done')

# MS2 cost per row vs concurrency
p = os.path.join(R, 'experiments.tsv')
if os.path.exists(p):
    e = pd.read_csv(p, sep='\t')
    e = e[(e.input == 'sub10') & (e.intra_op.astype(str) == '1') & e.ms2_cpu_ms_per_row.notna()]
    e = e[e.ms2_cpu_ms_per_row != '']
    e['act'] = pd.to_numeric(e.ms2_busy_workers, errors='coerce')
    e['ms'] = pd.to_numeric(e.ms2_cpu_ms_per_row, errors='coerce')
    fig, ax = plt.subplots(figsize=(6, 3.8))
    for i, (br, g) in enumerate(e.groupby('batch_rows')):
        ax.scatter(g.act, g.ms, s=28, color=C[i], label=f'{br}-row chunks', zorder=3)
        for _, r in g.iterrows(): ax.annotate(r.label, (r.act, r.ms), fontsize=6, xytext=(3, 3), textcoords='offset points')
    ax.set_xscale('log', base=2); ax.set_xlabel('mean busy MS2 workers (workers that got chunks x utilisation)')
    ax.set_ylabel('MS2 session-ms per predicted row'); ax.set_ylim(0, None)
    ax.set_title('Per-row MS2 cost rises once more than ~32 sessions run at once (sub10)', loc='left', fontsize=9)
    ax.legend(frameon=False, fontsize=7); fig.tight_layout(); fig.savefig(os.path.join(R, 'ms2_cost_per_row.png'), dpi=130); plt.close(fig)
