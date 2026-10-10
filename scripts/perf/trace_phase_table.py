#!/usr/bin/env python3
"""trace_phase_table.py RUNDIR [RUNDIR...] -- first per-phase table of a traced DIALibGen run.

Reads RUNDIR/trace/{events,proc,threads,meta}.tsv (ODIA_TRACE output, schema in docs/tracing.md) and
RUNDIR/time.txt (/usr/bin/time -v). Writes RUNDIR/phase_table.tsv, RUNDIR/predictor_table.tsv and
RUNDIR/tracer_cost.tsv and prints them as Markdown.

phase_table.tsv columns
  run          run label (directory name)
  phase        level-0 trace phase; 'startup' = process start -> first phase (TOPP/OpenMS init,
               option parsing), 'gap' = time between phases not covered by any phase (log lines),
               'exit' = last phase end -> trace end, 'TOTAL' = process start -> trace end;
               'trace_snapshot (n)' = the tracer's n phase-boundary snapshots (malloc_info, numa_maps,
               smaps_rollup), summed: tracer cost, not work of the run; 'gap (n)' likewise summed
  t0_s, t1_s   phase bounds on the trace clock (s)
  wall_s       t1 - t0
  cpu_s        process user+system CPU in the phase (getrusage, interpolated between 100 ms samples;
               includes the tracer's own sampler threads); for phases shorter than 1 s, where that
               interpolation is noise, the CPU of the thread that ran the phase plus that of the
               predictor worker threads it started (events.tsv cpu_s of the phase and of its 'worker'
               rows on other threads), see cpu_src
  cpu_src      'process' or 'threads', as above
  cores        cpu_s / wall_s = mean active cores
  run_threads  mean number of threads in state R at the samples inside the phase, the sampler excluded
  peak_rss_gb  largest VmRSS sample inside the phase (GiB, 100 ms resolution); when the process
               high-water mark (VmHWM) rose between the last sample before the phase and the first one
               after it, that VmHWM if larger, so a spike between two samples is not missed (the
               window is then up to one sample wider than the phase)
  share_wall   wall_s / TOTAL wall_s
predictor_table.tsv columns (one row per model)
  model, sessions (max worker index + 1), calls (predict calls = dispatches), chunks, rows,
  phase_wall_s (the model's level-0 phase), dispatch_wall_s (sum of dispatch spans: workers alive),
  outside_dispatch_s (phase wall - dispatch wall: serial work around the predictor calls),
  chunk_busy_s (sum of chunk spans over all workers), ort_run_s, encode_s, unpack_s (sums),
  util (chunk_busy / sum over calls of dispatch span x workers in that call),
  eff_cores (chunk_busy / phase wall)
tracer_cost.tsv columns (one row)
  run, snapshots (count), snapshot_wall_s (their summed wall time, on the phase thread, between phases),
  sampler_cpu_s (CPU of the odia-sampler thread, last threads.tsv row), fast_samples, sample_ms_mean,
  sample_ms_max (meta.tsv), sampler_share_core (sampler_cpu_s / TOTAL wall: mean fraction of one core)
"""
import bisect

SHORT = 1.0   # s; phases shorter than this take their CPU from the phase thread (events.tsv cpu_s)
import os
import sys


def read_tsv(path):
    rows = []
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            rows.append(dict(zip(header, line.rstrip("\n").split("\t"))))
    return rows


def time_v(path):
    out = {}
    if not os.path.exists(path):
        return out
    for line in open(path):
        k, _, v = line.strip().rpartition(": ")
        out[k] = v
    return out


def elapsed_s(s):
    parts = [float(x) for x in s.split(":")]
    t = 0.0
    for p in parts:
        t = t * 60 + p
    return t


def analyse(run):
    tr = os.path.join(run, "trace")
    meta = {r["key"]: r["value"] for r in read_tsv(os.path.join(tr, "meta.tsv"))}
    ev = read_tsv(os.path.join(tr, "events.tsv"))
    pr = read_tsv(os.path.join(tr, "proc.tsv"))
    ts = [float(r["t_s"]) for r in pr]
    cpu = [float(r["utime_s"]) + float(r["stime_s"]) for r in pr]
    rss = [int(r["rss_kb"]) for r in pr]
    hwm = [int(r["hwm_kb"]) for r in pr]
    nrun = [int(r["n_run"]) - 1 for r in pr]   # the sampler is running while it samples

    def cpu_at(t):
        if t <= ts[0]:   # linear from (process start, 0 CPU) to the first sample
            st = float(meta.get("process_start_s", 0.0))
            return 0.0 if t <= st else cpu[0] * (t - st) / (ts[0] - st)
        if t >= ts[-1]:
            return cpu[-1]
        i = bisect.bisect_left(ts, t)
        a, b = ts[i - 1], ts[i]
        return cpu[i - 1] + (cpu[i] - cpu[i - 1]) * (t - a) / (b - a)

    def window(t0, t1, arr):
        i0 = bisect.bisect_left(ts, t0)
        i1 = bisect.bisect_right(ts, t1)
        if i1 <= i0:   # no sample inside: nearest following one
            j = min(i0, len(ts) - 1)
            return [arr[j]]
        return arr[i0:i1]

    def peak(t0, t1):
        p = max(window(t0, t1, rss))
        i0 = bisect.bisect_right(ts, t0) - 1        # last sample at or before the start
        i1 = min(bisect.bisect_left(ts, t1), len(ts) - 1)   # first sample at or after the end
        before = hwm[i0] if i0 >= 0 else 0
        if hwm[i1] > before:
            p = max(p, hwm[i1])
        return p

    start = float(meta.get("process_start_s", 0.0))
    end = float(meta.get("end_s", ts[-1]))
    phases = sorted(((float(e["t0_s"]), float(e["t1_s"]), e["name"]) for e in ev if e["level"] == "0"))
    main_cpu = {(e["name"], float(e["t0_s"])): float(e["cpu_s"]) for e in ev if e["level"] == "0"}
    main_tid = {(e["name"], float(e["t0_s"])): e["tid"] for e in ev if e["level"] == "0"}
    workers = [(float(e["t0_s"]), e["tid"], float(e["cpu_s"])) for e in ev if e["name"] == "worker"]
    segs = []
    prev = start
    first = True
    for t0, t1, name in phases:
        if t0 - prev > 1e-3:
            segs.append(("startup" if first else "gap", prev, t0))
        segs.append((name, t0, t1))
        prev = t1
        first = False
    if end - prev > 1e-3:
        segs.append(("exit", prev, end))
    total_wall = end - start

    rows = []
    for name, t0, t1 in segs:
        wall = t1 - t0
        src = "process"
        if name == "startup":
            c = cpu_at(t1)            # CPU from process start
        elif wall < SHORT and (name, t0) in main_cpu:
            # Shorter than a few sampler ticks: the interpolated process CPU is
            # noise, so take the phase thread's own CPU and that of the
            # predictor workers it ran on other threads.
            tid = main_tid[(name, t0)]
            c = main_cpu[(name, t0)] + sum(cw for tw, wtid, cw in workers
                                           if t0 <= tw <= t1 and wtid != tid and cw == cw)
            src = "threads"
        else:
            c = cpu_at(t1) - cpu_at(t0)
        nr = window(t0, t1, nrun)
        rows.append(dict(phase=name, t0_s=t0, t1_s=t1, wall_s=wall, cpu_s=c, cpu_src=src,
                         cores=c / wall if wall > 0 else 0.0,
                         run_threads=sum(nr) / len(nr),
                         peak_rss_gb=peak(t0, t1) / 1048576.0,
                         share_wall=wall / total_wall))
    # merge the gaps, and the tracer's snapshots, into one row each, after the phases
    merged, pooled = [], {}
    for r in rows:
        if r["phase"] in ("gap", "trace_snapshot"):
            g = pooled.get(r["phase"])
            if g is None:
                g = pooled[r["phase"]] = dict(r)
                g["n"] = 1
                g["rt_sum"] = r["run_threads"] * r["wall_s"]
                g["srcs"] = {r["cpu_src"]}
            else:
                g["wall_s"] += r["wall_s"]; g["cpu_s"] += r["cpu_s"]
                g["peak_rss_gb"] = max(g["peak_rss_gb"], r["peak_rss_gb"])
                g["rt_sum"] += r["run_threads"] * r["wall_s"]; g["n"] += 1
                g["srcs"].add(r["cpu_src"])
        else:
            merged.append(r)
    for name in ("trace_snapshot", "gap"):
        g = pooled.get(name)
        if not g:
            continue
        g["cores"] = g["cpu_s"] / g["wall_s"] if g["wall_s"] > 0 else 0.0
        g["run_threads"] = g["rt_sum"] / g["wall_s"] if g["wall_s"] > 0 else 0.0
        g["share_wall"] = g["wall_s"] / total_wall
        g["phase"] = "%s (%d)" % (name, g["n"])
        g["cpu_src"] = "+".join(sorted(g["srcs"]))
        g["t0_s"] = float("nan"); g["t1_s"] = float("nan")
        merged.append(g)
    nr_all = window(start, end, nrun)
    merged.append(dict(phase="TOTAL", t0_s=start, t1_s=end, wall_s=total_wall, cpu_s=cpu_at(end), cpu_src="process",
                       cores=cpu_at(end) / total_wall, run_threads=sum(nr_all) / len(nr_all),
                       peak_rss_gb=max(max(rss), max(hwm)) / 1048576.0, share_wall=1.0))

    # predictor efficiency
    pred = []
    phase_wall = {n: t1 - t0 for t0, t1, n in phases}
    for model, ph in (("RT", "rt"), ("MS2", "ms2"), ("CCS", "ccs")):
        disp = [e for e in ev if e["name"] == "dispatch" and e["model"] == model]
        if not disp:
            continue
        chunks = [e for e in ev if e["name"] == "chunk" and e["model"] == model]
        def s(name):
            return sum(float(e["dur_s"]) for e in ev if e["name"] == name and e["model"] == model)
        busy = sum(float(e["dur_s"]) for e in chunks)
        cap = sum(float(e["dur_s"]) * int(e["session"]) for e in disp)
        dw = sum(float(e["dur_s"]) for e in disp)
        pw = phase_wall.get(ph, float("nan"))
        pred.append(dict(model=model, sessions=1 + max(int(e["session"]) for e in chunks),
                         calls=len(disp), chunks=len(chunks), rows=sum(int(e["rows"]) for e in chunks),
                         phase_wall_s=pw, dispatch_wall_s=dw, outside_dispatch_s=pw - dw,
                         chunk_busy_s=busy, ort_run_s=s("ort_run"), encode_s=s("encode"),
                         unpack_s=s("unpack"), util=busy / cap if cap > 0 else float("nan"),
                         eff_cores=busy / pw if pw > 0 else float("nan")))
    # the tracer's own cost
    snaps = [e for e in ev if e["level"] == "0" and e["name"] == "trace_snapshot"]
    sampler_cpu = float("nan")
    thp = os.path.join(tr, "threads.tsv")
    if os.path.exists(thp):
        for r in read_tsv(thp):
            if r["comm"] == "odia-sampler":
                sampler_cpu = float(r["utime_s"]) + float(r["stime_s"])
    cost = dict(snapshots=len(snaps), snapshot_wall_s=sum(float(e["dur_s"]) for e in snaps),
                sampler_cpu_s=sampler_cpu, fast_samples=meta.get("fast_samples", "-"),
                sample_ms_mean=float(meta.get("sample_ms_mean", "nan")),
                sample_ms_max=float(meta.get("sample_ms_max", "nan")),
                sampler_share_core=sampler_cpu / total_wall if total_wall > 0 else float("nan"))
    return meta, merged, pred, cost


def fmt(v):
    if isinstance(v, float):
        if v != v:
            return "-"
        return "%.3f" % v if abs(v) < 10 else "%.1f" % v
    return str(v)


def main():
    for run in sys.argv[1:]:
        label = os.path.basename(os.path.normpath(run))
        meta, rows, pred, cost = analyse(run)
        tv = time_v(os.path.join(run, "time.txt"))
        cols = ["phase", "t0_s", "t1_s", "wall_s", "cpu_s", "cpu_src", "cores", "run_threads", "peak_rss_gb", "share_wall"]
        with open(os.path.join(run, "phase_table.tsv"), "w") as f:
            f.write("run\t" + "\t".join(cols) + "\n")
            for r in rows:
                f.write(label + "\t" + "\t".join(fmt(r[c]) if c not in ("phase", "cpu_src") else r[c] for c in cols) + "\n")
        pcols = ["model", "sessions", "calls", "chunks", "rows", "phase_wall_s", "dispatch_wall_s",
                 "outside_dispatch_s", "chunk_busy_s", "ort_run_s", "encode_s", "unpack_s", "util", "eff_cores"]
        with open(os.path.join(run, "predictor_table.tsv"), "w") as f:
            f.write("run\t" + "\t".join(pcols) + "\n")
            for r in pred:
                f.write(label + "\t" + "\t".join(fmt(r[c]) for c in pcols) + "\n")
        ccols = ["snapshots", "snapshot_wall_s", "sampler_cpu_s", "fast_samples", "sample_ms_mean",
                 "sample_ms_max", "sampler_share_core"]
        with open(os.path.join(run, "tracer_cost.tsv"), "w") as f:
            f.write("run\t" + "\t".join(ccols) + "\n")
            f.write(label + "\t" + "\t".join(fmt(cost[c]) for c in ccols) + "\n")
        print(f"\n### {label}")
        if tv:
            print(f"/usr/bin/time: elapsed {tv.get('Elapsed (wall clock) time (h:mm:ss or m:ss)')}, "
                  f"user {tv.get('User time (seconds)')} s, sys {tv.get('System time (seconds)')} s, "
                  f"max RSS {int(tv.get('Maximum resident set size (kbytes)', 0)) / 1048576:.2f} GiB")
        print("| " + " | ".join(c for c in cols if c not in ("t0_s", "t1_s")) + " |")
        print("|" + "---|" * (len(cols) - 2))
        for r in rows:
            print("| " + " | ".join(fmt(r[c]) if c not in ("phase", "cpu_src") else r[c] for c in cols if c not in ("t0_s", "t1_s")) + " |")
        print()
        print("| " + " | ".join(pcols) + " |")
        print("|" + "---|" * len(pcols))
        for r in pred:
            print("| " + " | ".join(fmt(r[c]) for c in pcols) + " |")
        print()
        print("tracer: " + ", ".join("%s %s" % (c, fmt(cost[c])) for c in ccols))


if __name__ == "__main__":
    main()
