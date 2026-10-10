# Performance tracing

DIALibGen can record where a run spends its wall and CPU time. Tracing is off
unless `ODIA_TRACE` names a directory. When it is unset, each instrumented point
tests one boolean, which is read once per process, and does nothing else: no
allocation, no thread, no file. The library a run writes is byte-identical with
and without tracing.

```bash
ODIA_TRACE=trace-dir DIALibGen -in proteins.fasta -config cfg.json -out lib.parquet -threads 96
```

| variable | default | effect |
|---|---|---|
| `ODIA_TRACE` | unset | output directory; created if missing, existing files overwritten |
| `ODIA_TRACE_INTERVAL_MS` | 100 | period of the fast sampler (`proc.tsv`, `threads.tsv`); minimum 10 |
| `ODIA_TRACE_NUMA_MS` | 2000 | period of `numa.tsv`; 0 disables it |
| `ODIA_TRACE_SMAPS_MS` | 10000 | period of `smaps.tsv`; 0 disables it |
| `ODIA_TRACE_MALLOCINFO_MS` | 2000 | period of `mallocinfo.tsv`; 0 disables it |
| `ODIA_ORT_PROFILE` | unset | path prefix; every ONNX Runtime session writes its per-operator profile to `<prefix>_<model>_i<instance>_s<session>_<time>.json` |

`ODIA_ORT_PROFILE` works without `ODIA_TRACE`. It keeps every operator call in
memory until the predictor is destroyed, so use it on small inputs.

The samplers read `/proc` and run on Linux only. Elsewhere `events.tsv` and
`meta.tsv` are still written.

PSS and per-node page counts are not cheap. `/proc/self/smaps_rollup` and
`/proc/self/numa_maps` walk the page tables. At 15 GB resident that takes 270 ms
and 200 ms per read on the development node. That is why they have their own,
longer periods, and why each row records its own cost (`cost_ms`). RSS and its
anonymous and file parts come from the counters in `/proc/self/status` every
fast tick.

## Files

All files are tab-separated, with one header line. Times are seconds on one
monotonic clock whose zero is the first traced call (`t_s`, `t0_s`, `t1_s`).
`meta.tsv` has the wall-clock time of that zero (`start_epoch_s`) and the
process start on the same clock (`process_start_s`, negative). A value of -1
means not applicable or not available.

### events.tsv: one row per timed region

| column | meaning |
|---|---|
| `t0_s`, `t1_s`, `dur_s` | start, end and duration |
| `tid` | kernel thread id; it joins `threads.tsv` |
| `cpu_s` | CPU used by this thread inside the region (`CLOCK_THREAD_CPUTIME_ID`) |
| `level` | 0 phase, 1 step, 2 step inside a block, 3 predictor internals, 4 inside one chunk |
| `phase` | the level-0 phase that was open when the region ended, on any thread |
| `name` | region name (below) |
| `model` | `RT`, `MS2` or `CCS` for predictor regions, else `-` |
| `call` | predictor call number, shared by the dispatch, worker and chunk rows of one call |
| `idx` | block index (`block`), chunk index (`chunk`), chunk count (`chunking`, `dispatch`) or chunks done (`worker`) |
| `session` | ONNX session (= worker) index; for `dispatch`, the number of workers |
| `rows` | peptides, precursors or entries handled (`setup`: the `-threads` value; `sessions`: sessions created; `write`: transitions) |
| `len` | residues per peptide in a chunk (chunks are length-homogeneous) |
| `wait_s` | `chunk`: dispatch start to dequeue (queue wait); `worker`: dispatch start to thread start |

Regions recorded by `generate`:

| level | names |
|---|---|
| 0 | `setup`, `fasta_read`, `digest`, `enumerate`, `rt`, `irt_fit`, `ms2`, `ccs`, `decoys`, `fingerprint`, `write`, `teardown` |
| 1 | `rt`: `unique_parse`, `sessions`, `predict`, `assign`; `ms2`: `sessions`, `reserve`, `block`, `finalize`; `ccs`: `sessions`, `block`; `fingerprint`: `hash_fasta`, `hash_model`; `write`: `build_arrays`, `finish_arrays`, `write_table`, `commit` (Parquet) or `format_rows`, `close`, `commit` (TSV) |
| 2 | inside each `block`: `parse` (AASequence::fromString), `predict`, `post` (fragment enumeration, base-peak floor, ranking and cap for MS2; ion mobility for CCS) |
| 3 | `session_create`, `chunking` (length grouping), `dispatch` (workers alive), `worker` (one per thread per call), `chunk` (one batch) |
| 4 | `encode` (tokenisation and tensors), `ort_run` (ONNX Runtime `Run`), `unpack` (validation and copy-out) |

Phases follow one another on the main thread. Time between them (log output)
belongs to no phase. A phase ends when its function returns, so the part of a
phase that none of its steps covers is mostly what the function frees on the
way out: for `rt`, `ms2` and `ccs` that is the predictor, whose sessions and
arenas are released there (about 1 s per phase at 16 sessions).

### proc.tsv: process, every fast tick

| column | source and meaning |
|---|---|
| `t_s`, `phase` | sample time; level-0 phase open at that time (`-` between phases) |
| `sample_ms` | what this sample cost the sampler |
| `utime_s`, `stime_s` | `getrusage(RUSAGE_SELF)`: user and system CPU of all threads, exited ones included |
| `threads`, `n_run`, `n_sleep`, `n_disk`, `n_other` | threads in `/proc/self/task` and how many are in state R, S, D or other; `n_run` includes the sampler itself |
| `rss_kb`, `rss_anon_kb`, `rss_file_kb`, `rss_shmem_kb`, `hwm_kb`, `vm_kb` | `/proc/self/status`: VmRSS, RssAnon, RssFile, RssShmem, VmHWM, VmSize |
| `minflt`, `majflt`, `vcsw`, `nvcsw`, `inblock`, `oublock` | `getrusage`: page faults, voluntary and involuntary context switches, block I/O, cumulative |
| `rchar`, `wchar`, `syscr`, `syscw`, `read_bytes`, `write_bytes` | `/proc/self/io`, cumulative |
| `mi_arena_b`, `mi_mmap_b`, `mi_inuse_b`, `mi_free_b`, `mi_keepcost_b` | glibc `mallinfo2`: arena, hblkhd, uordblks, fordblks, keepcost (bytes; -1 when not glibc >= 2.33) |
| `mallinfo_ms` | cost of the `mallinfo2` call, which locks every arena in turn |

### threads.tsv: per thread, every fast tick

A thread gets a row at a tick when its state, CPU times, page faults or last CPU
changed since its previous row, and at its first tick and the final one. A
thread with no row at a tick was asleep and unchanged.

| column | source and meaning |
|---|---|
| `t_s`, `tid`, `comm` | sample time, thread id, thread name. Traced workers are named `odia-<model>-s<session>` and the samplers `odia-sampler` and `odia-sampler2`. Session 0 runs on the calling thread, which keeps the process name. |
| `state`, `cpu` | `/proc/self/task/<tid>/stat`: state letter and the CPU the thread last ran on (see `topology.tsv`) |
| `utime_s`, `stime_s` | the same file: user and system CPU in seconds, cumulative, at clock-tick (10 ms) resolution |
| `run_s`, `wait_s`, `slices` | `schedstat`: time on a CPU, time runnable but waiting for one (run-queue delay), timeslices; ns precision |
| `vcsw`, `nvcsw`, `migrations` | `sched`: voluntary and involuntary switches and CPU migrations (from `status` without `migrations` when `sched` is absent) |
| `minflt`, `majflt` | `stat`: page faults of this thread, cumulative |

### numa.tsv: resident pages per NUMA node, every `ODIA_TRACE_NUMA_MS`

`t_s`, `phase`, `cost_ms`, `node`, `total_kb`, `anon_kb`, `file_kb`, `huge_kb`:
one row per node and sample, summed over `/proc/self/numa_maps`. A mapping counts
as file-backed when it has a `file=` entry and as anonymous otherwise. `huge_kb`
counts pages in mappings whose page size exceeds 4 kB.

### smaps.tsv: every `ODIA_TRACE_SMAPS_MS`

`t_s`, `phase`, `cost_ms`, then `rss_kb`, `pss_kb`, `pss_anon_kb`, `pss_file_kb`,
`pss_dirty_kb`, `anonymous_kb`, `anon_huge_kb`, `private_dirty_kb`,
`private_clean_kb`, `shared_clean_kb`, `swap_kb` from `/proc/self/smaps_rollup`.

### mallocinfo.tsv: glibc heap, every `ODIA_TRACE_MALLOCINFO_MS`

`t_s`, `phase`, `cost_ms`, `arenas` (the `<heap>` count of `malloc_info`), then
the process-wide totals in bytes: `fast_free_b`, `rest_free_b`, `mmap_count`,
`mmap_b`, `system_current_b`, `system_max_b`, `aspace_total_b`,
`aspace_mprotect_b`.

### topology.tsv: once

`cpu`, `node`, `package`, `core_id`, `l3_id` (from
`/sys/devices/system/cpu/cpu*/cache/index3/id`), and `allowed` (1 if the
process's affinity mask contains the CPU). Join it with `threads.tsv` on
`cpu` to see which threads shared an L3 or crossed a NUMA node.

### meta.tsv: key, value

The command line, pid, host, kernel, glibc and ONNX Runtime versions, the
sampler periods, `ODIA_*` and allocator environment variables, the session count
of each predictor, the ONNX profile files written, the sample counts, `end_s`
and the event count.

## A first summary

`scripts/perf/trace_phase_table.py RUNDIR...` reads `RUNDIR/trace/` (and
`RUNDIR/time.txt` from `/usr/bin/time -v`, if present). It writes
`phase_table.tsv`, with wall time, CPU, mean active cores, mean running threads
and peak RSS per phase, and `predictor_table.tsv`, with session utilisation and
the serial time around the predictor calls per model. Both have their columns
documented at the top of the script. It needs only the Python standard library.
