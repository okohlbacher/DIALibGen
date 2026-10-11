#!/usr/bin/env python3
"""collect.py A -- one row per run of perf/generate-solutions (and the baselines it is compared with).

Writes A/results/runs.tsv. Columns:
  label, input (sub10|full), threads, format (parquet|tsv), env (switches set), binary,
  elapsed_s, user_s, sys_s, cpu_s (user+sys), mean_cores (cpu_s/elapsed_s), maxrss_gib,
  md5, ref_md5, equal (1 = byte-identical to the reference), load1_max, other_cores_mean,
  then per level-0 phase wall seconds (from RUNDIR/phase_table.tsv): digest, enumerate, rt, ms2, ccs, write,
  and serial_s = digest+enumerate+write wall.
Reference md5: sub10 Parquet 8d8f2d22 (bid10-instrument, all thread counts), full Parquet c7d404a0,
full TSV e69f4ec9 (bid10-instrument full_t96_tsv); sub10 TSV: the s_base_tsv run at the same thread count.
"""
import csv, os, re, subprocess, sys

A = sys.argv[1]
INSTR = "/scratch/kohlbach/bid10-instrument"
REF = {("sub10", "parquet"): "8d8f2d22f588de49870c4a3206e29d35",
       ("full", "parquet"): "c7d404a0bd92af1469a62d9de32b5ad0",
       ("full", "tsv"): "e69f4ec931c9941bdfb93523d42879b9"}
PHASES = ["digest", "enumerate", "rt", "ms2", "ccs", "write"]


def meta(d):
    m = {}
    for line in open(os.path.join(d, "meta.txt")):
        k, _, v = line.rstrip("\n").partition(" ")
        m.setdefault(k, v)
    return m


def tv(d):
    t = open(os.path.join(d, "time.txt")).read()
    g = lambda pat: re.search(pat, t)
    el = g(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (\S+)").group(1).split(":")
    sec = sum(float(x) * 60 ** i for i, x in enumerate(reversed(el)))
    return (sec, float(g(r"User time \(seconds\): (\S+)").group(1)), float(g(r"System time \(seconds\): (\S+)").group(1)),
            int(g(r"Maximum resident set size \(kbytes\): (\d+)").group(1)) / 1024 ** 2)


def phases(d):
    p = os.path.join(d, "phase_table.tsv")
    if not os.path.exists(p) and os.path.isdir(os.path.join(d, "trace")):
        subprocess.run([sys.executable, os.path.join(A, "tools", "trace_phase_table.py"), d],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    out = {}
    if os.path.exists(p):
        for r in csv.DictReader(open(p), delimiter="\t"):
            out[r["phase"]] = float(r["wall_s"])
    return out


rows = []
dirs = [(os.path.join(A, "runs", x), x) for x in sorted(os.listdir(os.path.join(A, "runs"))) if x != "contaminated"]
for d, label in dirs:
    if not os.path.exists(os.path.join(d, "time.txt")) or not os.path.exists(os.path.join(d, "meta.txt")):
        continue
    if "Elapsed" not in open(os.path.join(d, "time.txt")).read():
        continue   # still running
    m = meta(d)
    cmd = m.get("cmd", "")
    inp = "sub10" if "sub10.fasta" in cmd else ("full" if "proteome.fasta" in cmd else "other")
    o = re.search(r"-out (\S+)", cmd)
    fmt = "tsv" if o and o.group(1).endswith(".tsv") else ("parquet" if o else "speclib")
    th = re.search(r"-threads (\d+)", cmd)
    env = m.get("env", "")
    env = " ".join(x for x in env.split() if not x.startswith(("ODIA_TRACE=", "ODIA_ORT_PROFILE=")) and "=" in x)
    el, us, sy, rss = tv(d)
    md5 = m.get("out_md5", "")
    ref = REF.get((inp, fmt), "")
    if not ref and inp == "sub10" and fmt == "tsv" and th:
        b = os.path.join(A, "runs", f"s_base_tsv_t{th.group(1)}", "meta.txt")
        if os.path.exists(b):
            ref = meta(os.path.dirname(b)).get("out_md5", "")
    ph = phases(d)
    row = dict(label=label, input=inp, threads=th.group(1) if th else "", format=fmt, env=env or "-",
               binary=os.path.basename(cmd.split()[0]) if cmd else "",
               elapsed_s=f"{el:.2f}", user_s=f"{us:.1f}", sys_s=f"{sy:.1f}", cpu_s=f"{us + sy:.1f}",
               mean_cores=f"{(us + sy) / el:.2f}", maxrss_gib=f"{rss:.2f}", md5=md5, ref_md5=ref,
               equal=int(bool(md5) and md5 == ref), load1_max=m.get("load1_max", ""),
               other_cores_mean=m.get("other_cores_mean", ""))
    for k in PHASES:
        row[k] = f"{ph[k]:.1f}" if k in ph else ""
    row["serial_s"] = f"{sum(ph.get(k, 0.0) for k in ('digest', 'enumerate', 'write')):.1f}"
    rows.append(row)
os.makedirs(os.path.join(A, "results"), exist_ok=True)
with open(os.path.join(A, "results", "runs.tsv"), "w") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
    w.writeheader()
    w.writerows(rows)
print("\t".join(["label", "thr", "elapsed", "cpu", "rss", "eq", "dig", "enum", "rt", "ms2", "ccs", "write", "oth"]))
for r in rows:
    print("\t".join([r["label"], r["threads"], r["elapsed_s"], r["cpu_s"], r["maxrss_gib"], str(r["equal"]),
                     r["digest"], r["enumerate"], r["rt"], r["ms2"], r["ccs"], r["write"], r["other_cores_mean"]]))
