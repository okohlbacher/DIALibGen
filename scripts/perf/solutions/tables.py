#!/usr/bin/env python3
"""tables.py A -- comparison tables of perf/generate-solutions (reads A/results/runs.tsv from collect.py).

Writes
  A/results/changes_sub10.tsv : one row per change and thread count, against s_base_t<N> (same binary family)
      change, threads, env, elapsed_s, base_elapsed_s, delta_pct, cpu_s, base_cpu_s, maxrss_gib, base_maxrss_gib,
      equal (1 = byte-identical to the reference md5), serial_s, base_serial_s
  A/results/full_compare.tsv : full proteome, every engine/build at 8/16/96 threads
      system, threads, format, elapsed_s, user_s, maxrss_gib, md5, equal_ref, source
      systems: baseline (bid10-instrument traced, e82ce88), shared-data (bid10-shared-data traced, d64c70a),
      solutions v3/v4 (this branch at 1dea09b / 039b1d3, ODIA_SOLUTIONS=1 [+ pruned MS2 graph]), DIA-NN 2.0
  A/results/full_phases.tsv : level-0 phase wall s per full-proteome run (run, phase, wall_s, cpu_s, cores, peak_rss_gb)
"""
import csv, os, re, sys

A = sys.argv[1]
INSTR, SD = "/scratch/kohlbach/bid10-instrument", "/scratch/kohlbach/bid10-shared-data"
runs = {r["label"]: r for r in csv.DictReader(open(os.path.join(A, "results", "runs.tsv")), delimiter="\t")}

out = []
for lab, r in sorted(runs.items()):
    m = re.match(r"s_(.+)_t(\d+)$", lab)
    if not m or m.group(1) == "base" or m.group(1) == "base_tsv":
        continue
    change, t = m.group(1), m.group(2)
    base = runs.get(f"s_base_tsv_t{t}" if r["format"] == "tsv" else f"s_base_t{t}")
    if not base:
        continue
    e, be = float(r["elapsed_s"]), float(base["elapsed_s"])
    out.append(dict(change=change, threads=t, env=r["env"], elapsed_s=r["elapsed_s"], base_elapsed_s=base["elapsed_s"],
                    delta_pct=f"{100 * (e - be) / be:+.1f}", cpu_s=r["cpu_s"], base_cpu_s=base["cpu_s"],
                    maxrss_gib=r["maxrss_gib"], base_maxrss_gib=base["maxrss_gib"], equal=r["equal"],
                    serial_s=r["serial_s"], base_serial_s=base["serial_s"]))
with open(os.path.join(A, "results", "changes_sub10.tsv"), "w") as f:
    w = csv.DictWriter(f, fieldnames=list(out[0]), delimiter="\t"); w.writeheader(); w.writerows(out)


def tv(d):
    t = open(os.path.join(d, "time.txt")).read()
    el = re.search(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (\S+)", t).group(1).split(":")
    return (sum(float(x) * 60 ** i for i, x in enumerate(reversed(el))),
            float(re.search(r"User time \(seconds\): (\S+)", t).group(1)),
            int(re.search(r"Maximum resident set size \(kbytes\): (\d+)", t).group(1)) / 1024 ** 2)


def md5(d):
    for line in open(os.path.join(d, "meta.txt")):
        if line.startswith("out_md5 "):
            return line.split()[1]
    return ""


REF = {"parquet": "c7d404a0bd92af1469a62d9de32b5ad0", "tsv": "e69f4ec931c9941bdfb93523d42879b9"}
full, phases = [], []
srcs = [("baseline", INSTR, {"full_t96_pq": (96, "parquet"), "full_t16_pq": (16, "parquet"), "full_t8_pq": (8, "parquet"),
                             "full_t96_tsv": (96, "tsv")}),
        ("shared-data", SD, {"sd_full_t96_pq": (96, "parquet"), "sd_full_t16_pq": (16, "parquet"),
                             "sd_full_t8_pq": (8, "parquet"), "sd_full_t96_tsv": (96, "tsv")}),
        ("solutions v3 (+pruned)", A, {"f_combop_t96_pq": (96, "parquet"), "f_combop_t16_pq": (16, "parquet"),
                                       "f_combop_t8_pq": (8, "parquet"), "f_combop_t96_tsv": (96, "tsv"),
                                       "f_combop_t16_tsv": (16, "tsv")}),
        ("solutions v3 (unpruned)", A, {"f_combo_t96_pq": (96, "parquet")}),
        ("solutions v4 (+pruned)", A, {"f4_combop_t96_pq": (96, "parquet"), "f4_combop_t16_pq": (16, "parquet"),
                                       "f4_combop_t8_pq": (8, "parquet"), "f4_combop_t96_tsv": (96, "tsv"),
                                       "f4_combop_t8_tsv": (8, "tsv")}),
        ("DIA-NN 2.0", A, {"dn_full_t96": (96, "parquet"), "dn_full_t16": (16, "parquet"), "dn_full_t8": (8, "parquet")})]
for system, root, labs in srcs:
    for lab, (t, fmt) in labs.items():
        d = os.path.join(root, "runs", lab)
        if not os.path.exists(os.path.join(d, "time.txt")):
            continue
        el, us, rss = tv(d)
        h = md5(d)
        full.append(dict(system=system, threads=t, format=fmt, elapsed_s=f"{el:.1f}", user_s=f"{us:.1f}",
                         maxrss_gib=f"{rss:.2f}", md5=h, equal_ref=("-" if system == "DIA-NN 2.0" else int(h == REF[fmt])),
                         source=d))
        p = os.path.join(d, "phase_table.tsv")
        if os.path.exists(p):
            for r in csv.DictReader(open(p), delimiter="\t"):
                if r["phase"] in ("digest", "enumerate", "rt", "ms2", "ccs", "write", "TOTAL"):
                    phases.append(dict(system=system, run=lab, threads=t, format=fmt, phase=r["phase"],
                                       wall_s=r["wall_s"], cpu_s=r["cpu_s"], cores=r["cores"], peak_rss_gb=r["peak_rss_gb"]))
for name, rows in (("full_compare.tsv", full), ("full_phases.tsv", phases)):
    with open(os.path.join(A, "results", name), "w") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)
for r in out:
    print("\t".join(str(r[k]) for k in ("change", "threads", "elapsed_s", "base_elapsed_s", "delta_pct", "cpu_s", "maxrss_gib", "equal")))
print()
for r in full:
    print("\t".join(str(r[k]) for k in ("system", "threads", "format", "elapsed_s", "user_s", "maxrss_gib", "equal_ref")))
