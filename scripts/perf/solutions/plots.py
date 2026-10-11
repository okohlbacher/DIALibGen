#!/usr/bin/env python3
"""plots.py A -- plot data and figures of perf/generate-solutions.

Writes A/results/timeline_full.tsv (plot data; 1-s bins from each run's trace/proc.tsv):
  run, system, threads, t_s (trace clock, bin start), phase (level-0 phase at the bin start),
  cores (d(utime+stime)/dt), threads_running (mean n_run), rss_gib (mean VmRSS)
and the figures
  A/results/timeline_full_t{96,16,8}.png : active cores and RSS against time, baseline / shared-data / solutions
  A/results/phases_full.png              : wall s per level-0 phase, full proteome (results/full_phases.tsv)
  A/results/changes_sub10.png            : elapsed change per switch alone, sub10, 16 and 96 threads
Run: ~/.local/bin/uv run --with matplotlib --with pandas python plots.py A
"""
import csv, os, sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

A = sys.argv[1]
R = os.path.join(A, "results")
INSTR, SD = "/scratch/kohlbach/bid10-instrument/runs", "/scratch/kohlbach/bid10-shared-data/runs"
RUNS = {96: [("baseline", f"{INSTR}/full_t96_pq"), ("shared-data", f"{SD}/sd_full_t96_pq"),
             ("solutions v4 (+pruned)", f"{A}/runs/f4_combop_t96_pq")],
        16: [("baseline", f"{INSTR}/full_t16_pq"), ("shared-data", f"{SD}/sd_full_t16_pq"),
             ("solutions v4 (+pruned)", f"{A}/runs/f4_combop_t16_pq")],
        8: [("baseline", f"{INSTR}/full_t8_pq"), ("shared-data", f"{SD}/sd_full_t8_pq"),
            ("solutions v4 (+pruned)", f"{A}/runs/f4_combop_t8_pq")]}

rows = []
for t, lst in RUNS.items():
    for system, d in lst:
        p = os.path.join(d, "trace", "proc.tsv")
        if not os.path.exists(p):
            continue
        df = pd.read_csv(p, sep="\t")
        df = df[df["sample_ms"] >= 0]
        df["cpu"] = df["utime_s"] + df["stime_s"]
        df["bin"] = df["t_s"].astype(float).floordiv(1.0)
        g = df.groupby("bin").agg(t0=("t_s", "first"), t1=("t_s", "last"), c0=("cpu", "first"), c1=("cpu", "last"),
                                  phase=("phase", "first"), run=("n_run", "mean"), rss=("rss_kb", "mean"))
        # cores from consecutive bins' first samples
        g = g.reset_index()
        g["cores"] = (g["c0"].shift(-1) - g["c0"]) / (g["t0"].shift(-1) - g["t0"])
        for _, r in g.dropna().iterrows():
            rows.append(dict(run=os.path.basename(d), system=system, threads=t, t_s=int(r["bin"]), phase=r["phase"],
                             cores=round(r["cores"], 2), threads_running=round(r["run"], 2),
                             rss_gib=round(r["rss"] / 1024 ** 2, 3)))
tl = pd.DataFrame(rows)
tl.to_csv(os.path.join(R, "timeline_full.tsv"), sep="\t", index=False)

colors = {"baseline": "#7f7f7f", "shared-data": "#1f77b4", "solutions v4 (+pruned)": "#d62728"}
for t in RUNS:
    sub = tl[tl.threads == t]
    if sub.empty:
        continue
    fig, ax = plt.subplots(2, 1, figsize=(11, 6), sharex=True)
    for system, g in sub.groupby("system"):
        ax[0].plot(g.t_s, g.cores, lw=0.8, label=system, color=colors.get(system))
        ax[1].plot(g.t_s, g.rss_gib, lw=0.8, label=system, color=colors.get(system))
    ax[0].set_ylabel("active cores")
    ax[1].set_ylabel("RSS GiB")
    ax[1].set_xlabel("s since trace start")
    ax[0].legend(loc="upper right", fontsize=8)
    ax[0].set_title(f"generate, full proteome, -threads {t} (kim, traced, one run each)")
    fig.tight_layout()
    fig.savefig(os.path.join(R, f"timeline_full_t{t}.png"), dpi=130)
    plt.close(fig)

ph = pd.read_csv(os.path.join(R, "full_phases.tsv"), sep="\t")
ph = ph[(ph.format == "parquet") & (ph.phase != "TOTAL")]
piv = ph.pivot_table(index=["threads", "system"], columns="phase", values="wall_s", aggfunc="first")
order = ["digest", "enumerate", "rt", "ms2", "ccs", "write"]
piv = piv[[c for c in order if c in piv.columns]]
fig, ax = plt.subplots(figsize=(11, 5))
piv.plot.barh(stacked=True, ax=ax)
ax.set_xlabel("wall s (trace phases)")
ax.set_title("generate, full proteome: wall time per phase (Parquet)")
fig.tight_layout()
fig.savefig(os.path.join(R, "phases_full.png"), dpi=130)
plt.close(fig)

ch = pd.read_csv(os.path.join(R, "changes_sub10.tsv"), sep="\t")
p2 = ch.pivot_table(index="change", columns="threads", values="delta_pct", aggfunc="first")
fig, ax = plt.subplots(figsize=(9, 7))
p2.sort_values(by=p2.columns[-1]).plot.barh(ax=ax)
ax.axvline(0, color="k", lw=0.6)
ax.set_xlabel("elapsed vs s_base at the same -threads, % (sub10, one run each; noise about +-2 %)")
ax.set_title("each switch alone and the combinations, sub10")
fig.tight_layout()
fig.savefig(os.path.join(R, "changes_sub10.png"), dpi=130)
plt.close(fig)
print("wrote", R)
