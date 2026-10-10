#!/usr/bin/env python3
"""trace_phase_table.py RUNDIR [RUNDIR...] -- first per-phase table of a traced DIALibGen run.

Reads RUNDIR/trace/{events,proc,meta}.tsv (ODIA_TRACE output, schema in docs/tracing.md) and
RUNDIR/time.txt (/usr/bin/time -v). Writes RUNDIR/phase_table.tsv and RUNDIR/predictor_table.tsv
and prints both as Markdown.

phase_table.tsv columns
  run          run label (directory name)
  phase        level-0 trace phase; 'startup' = process start -> first phase (TOPP/OpenMS init,
               option parsing), 'gap' = time between phases not covered by any phase (log lines),
               'exit' = last phase end -> trace end, 'TOTAL' = process start -> trace end
  t0_s, t1_s   phase bounds on the trace clock (s)
  wall_s       t1 - t0
  cpu_s        process user+system CPU in the phase (getrusage, interpolated between 100 ms samples;
               includes the tracer's own sampler threads); for phases shorter than 1 s the CPU of the
               thread that ran the phase (events.tsv cpu_s), see cpu_src
  cpu_src      'process' or 'thread', as above
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
            # noise, so take the phase thread's own CPU (these phases are serial).
            c = main_cpu[(name, t0)]
            src = "thread"
        else:
            c = cpu_at(t1) - cpu_at(t0)
        nr = window(t0, t1, nrun)
        rows.append(dict(phase=name, t0_s=t0, t1_s=t1, wall_s=wall, cpu_s=c, cpu_src=src,
                         cores=c / wall if wall > 0 else 0.0,
                         run_threads=sum(nr) / len(nr),
                         peak_rss_gb=peak(t0, t1) / 1048576.0,
                         share_wall=wall / total_wall))
    # merge the gaps into one row, keep order of first appearance
    merged, gap = [], None
    for r in rows:
        if r["phase"] == "gap":
            if gap is None:
                gap = dict(r)
                gap["n"] = 1
                gap["rt_sum"] = r["run_threads"] * r["wall_s"]
            else:
                gap["wall_s"] += r["wall_s"]; gap["cpu_s"] += r["cpu_s"]
                gap["peak_rss_gb"] = max(gap["peak_rss_gb"], r["peak_rss_gb"])
                gap["rt_sum"] += r["run_threads"] * r["wall_s"]; gap["n"] += 1
        else:
            merged.append(r)
    if gap:
        gap["cores"] = gap["cpu_s"] / gap["wall_s"]
        gap["run_threads"] = gap["rt_sum"] / gap["wall_s"]
        gap["share_wall"] = gap["wall_s"] / total_wall
        gap["phase"] = "gap (%d)" % gap["n"]
        gap["cpu_src"] = "process"
        gap["t0_s"] = float("nan"); gap["t1_s"] = float("nan")
        merged.append(gap)
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
    return meta, merged, pred


def fmt(v):
    if isinstance(v, float):
        if v != v:
            return "-"
        return "%.3f" % v if abs(v) < 10 else "%.1f" % v
    return str(v)


def main():
    for run in sys.argv[1:]:
        label = os.path.basename(os.path.normpath(run))
        meta, rows, pred = analyse(run)
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


if __name__ == "__main__":
    main()
