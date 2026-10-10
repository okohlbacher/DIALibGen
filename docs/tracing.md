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
| `ODIA_TRACE_INTERVAL_MS` | 100 | period of the sampler (`proc.tsv`, `threads.tsv`); minimum 10 |
| `ODIA_TRACE_SNAPSHOT_MIN_S` | 30 | least time between two phase-boundary snapshots (`numa.tsv`, `smaps.tsv`, `mallocinfo.tsv`); 0 takes one at every phase start, `off` or a negative value none |
| `ODIA_ORT_PROFILE` | unset | path prefix; every ONNX Runtime session writes its per-operator profile to `<prefix>_<model>_i<instance>_s<session>_<time>.json` |

Measurement knobs (branch perf/threads-analysis; unset, each leaves behaviour unchanged; recorded in `meta.tsv`):

| variable | default | effect |
|---|---|---|
| `ODIA_MAX_SESSIONS` | 16 | cap on ONNX sessions per predictor (sessions = min(-threads, cap)) |
| `ODIA_BATCH_ROWS` | 2048 | rows per chunk (`MAX_BATCH_ROWS`) |
| `ODIA_MS2_BLOCK` | 20000 | precursors per MS2 block (one predictor call each) |
| `ODIA_INTRA_OP` | 1 | ONNX Runtime intra-op threads per session |
| `ODIA_SPIN` | ORT default (1) | 0 or 1: `session.intra_op.allow_spinning` and `session.inter_op.allow_spinning` |

`ODIA_BATCH_ROWS`, `ODIA_MS2_BLOCK` and `ODIA_INTRA_OP` can change the last bit of some RT/CCS/IM values; see
`scripts/perf/threads/` for the analysis scripts that use them.

`ODIA_ORT_PROFILE` works without `ODIA_TRACE`. It keeps every operator call in
memory until the predictor is destroyed, so use it on small inputs.

The samplers read `/proc` and run on Linux only. Elsewhere `events.tsv` and
`meta.tsv` are still written.

## What is read when

The sampler thread reads, every `ODIA_TRACE_INTERVAL_MS`, only counters the
kernel hands out without taking a lock that the run's own threads need:
`getrusage`, `/proc/self/status`, `/proc/self/io`, `/proc/vmstat`, and per
thread `/proc/self/task/<tid>/stat`, `schedstat` and `status`. In the steady
state it allocates nothing. Each row records what it cost (`sample_ms`).

Three views are expensive and are never read periodically:

- glibc `malloc_info` locks each arena in turn and visits every free chunk.
  After a phase has freed millions of small objects a single call takes
  seconds, and every thread that allocates from that arena waits for it.
- `/proc/self/numa_maps` and `/proc/self/smaps_rollup` hold the process's mmap
  lock while they walk every page: 200 to 270 ms at 15 GB resident, during which
  page faults and `mmap`/`munmap` in the run's threads can stall.

They are taken together as one snapshot, synchronously, on the thread that
starts a phase (level 0), just before the phase opens, so the previous phase's
worker threads are idle; at most once per `ODIA_TRACE_SNAPSHOT_MIN_S`; and once
more when the session ends. A snapshot is a level-0 event of its own
(`trace_snapshot` in `events.tsv`), so its cost is outside every phase, and each
row of the three files records the cost of its read (`cost_ms`).

The first version of the tracer called `mallinfo2` (the same walk) every
100 ms. On the full proteome, after the digest had freed its peptide map, each
call took about 2 s, the sampler held the main arena's lock almost
continuously, and a run that normally takes 17 minutes did not finish in 25.
On a tenth of the proteome the effect was invisible.

The second version (0a89832) had a use-after-free in the session's teardown:
it unlocked the snapshot mutex after deleting the object that holds it. glibc
then aborted the process at exit (`malloc(): unsorted double linked list
corrupted`), after the library was written but before `events.tsv`. ASan builds
did not show it; valgrind memcheck did. Fixed in e82ce88.

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
| 0 | `setup`, `fasta_read`, `digest`, `enumerate`, `release`, `rt`, `irt_fit`, `ms2`, `ccs`, `decoys`, `fingerprint`, `write`, `teardown`; and the tracer's own `trace_snapshot` |
| 1 | `rt`: `unique_parse`, `sessions`, `predict`, `assign`; `ms2`: `sessions`, `reserve`, `block`, `finalize`; `ccs`: `sessions`, `block`; `fingerprint`: `hash_fasta`, `hash_model`; `write`: `build_arrays`, `finish_arrays`, `write_table`, `commit` (Parquet) or `format_rows`, `close`, `commit` (TSV) |
| 2 | inside each `block`: `parse` (AASequence::fromString), `predict`, `post` (fragment enumeration, base-peak floor, ranking and cap for MS2; ion mobility for CCS) |
| 3 | `session_create`, `chunking` (length grouping), `dispatch` (workers alive), `worker` (one per thread per call), `chunk` (one batch); inside each MS2 `post`: `post_enumerate`, `post_rank`, `post_select`, `post_emit` (aggregates, below) |
| 4 | `encode` (tokenisation and tensors), `ort_run` (ONNX Runtime `Run`), `unpack` (validation and copy-out) |

The four `post_*` rows of an MS2 block are sums, not intervals: the post step
loops over precursors, so each part is timed per precursor and summed over the
block. `t0_s` is the start of the block's `post` step, `dur_s` the summed time,
`t1_s = t0_s + dur_s`, and `cpu_s` is -1. `post_enumerate` is the fragment
enumeration without the m/z cap; `post_rank` the base-peak floor and the
lookup of each fragment's predicted intensity; `post_select` the sort, the
doubly-charged quota, the cap and the minimum-fragment check; `post_emit`
appending the transitions. Precursors whose prediction failed are copied
unchanged and fall in none of the four.

Phases follow one another on the main thread. Time between them (log output)
belongs to no phase. `release` is the digest freeing its peptide map and parsed
sequences on the way out of `LibraryGenerator::generate`. `trace_snapshot` is
the tracer's phase-boundary snapshot (above); `idx` numbers them. A phase ends when its function returns, so the part of a
phase that none of its steps covers is mostly what the function frees on the
way out: for `rt`, `ms2` and `ccs` that is the predictor, whose sessions and
arenas are released there (about 1 s per phase at 16 sessions).

### proc.tsv: process, every fast tick

| column | source and meaning |
|---|---|
| `t_s`, `phase` | sample time; level-0 phase open at that time (`-` between phases, `trace_snapshot` during a snapshot) |
| `sample_ms` | what this sample cost the sampler |
| `utime_s`, `stime_s` | `getrusage(RUSAGE_SELF)`: user and system CPU of all threads, exited ones included |
| `threads`, `n_run`, `n_sleep`, `n_disk`, `n_other` | threads in `/proc/self/task` and how many are in state R, S, D or other; `n_run` includes the sampler itself |
| `rss_kb`, `rss_anon_kb`, `rss_file_kb`, `rss_shmem_kb`, `hwm_kb`, `vm_kb` | `/proc/self/status`: VmRSS, RssAnon, RssFile, RssShmem, VmHWM, VmSize. Since Linux 6.2 these are per-CPU batched counters read without summing; on 384 CPUs they can lag the exact value by up to about 1 GiB (768 pages per CPU, an upper bound). `smaps.tsv` has the exact value at the snapshots. |
| `minflt`, `majflt`, `vcsw`, `nvcsw`, `inblock`, `oublock` | `getrusage`: page faults, voluntary and involuntary context switches, block I/O, cumulative |
| `rchar`, `wchar`, `syscr`, `syscw`, `read_bytes`, `write_bytes` | `/proc/self/io`, cumulative |
| `sys_numa_hit`, `sys_numa_miss`, `sys_numa_local`, `sys_numa_other` | `/proc/vmstat`, cumulative pages allocated on the intended node or not, and by a CPU on the same node or another. System-wide, not per process: on an otherwise idle node the rate is the run's |
| `sys_numa_hint_faults`, `sys_numa_hint_faults_local`, `sys_numa_pages_migrated` | `/proc/vmstat`, cumulative: automatic NUMA balancing's hinting faults (all, and those where the page was already local) and the pages it moved. System-wide; zero growth when balancing is off (`/proc/sys/kernel/numa_balancing`) |

### threads.tsv: per thread, every fast tick

A thread gets a row at a tick when its state, CPU times, page faults or last CPU
changed since its previous row, and at its first tick and the final one. A
thread with no row at a tick was asleep and unchanged.

| column | source and meaning |
|---|---|
| `t_s`, `tid`, `comm` | sample time, thread id, thread name. Traced workers are named `odia-<model>-s<session>` and the sampler `odia-sampler`. Session 0 runs on the calling thread, which keeps the process name. |
| `state`, `cpu` | `/proc/self/task/<tid>/stat`: state letter and the CPU the thread last ran on (see `topology.tsv`) |
| `utime_s`, `stime_s` | the same file: user and system CPU in seconds, cumulative, at clock-tick (10 ms) resolution |
| `run_s`, `wait_s`, `slices` | `schedstat`: time on a CPU, time runnable but waiting for one (run-queue delay), timeslices; ns precision |
| `vcsw`, `nvcsw` | `status`: voluntary and involuntary context switches, cumulative. CPU migrations are not read (they are in `sched`, a larger file); a change of `cpu` between two rows of a thread shows at least one |
| `minflt`, `majflt` | `stat`: page faults of this thread, cumulative |

### numa.tsv, smaps.tsv, mallocinfo.tsv: at the snapshots

All three start with `t_s` (when the snapshot began), `at` (the phase about to
begin, or `end` for the one taken when the session ends) and `cost_ms` (the
cost of this read alone).

`numa.tsv`: `node`, `total_kb`, `anon_kb`, `file_kb`, `huge_kb`, one row per
node and snapshot, summed over `/proc/self/numa_maps`. A mapping counts as
file-backed when it has a `file=` entry and as anonymous otherwise. `huge_kb`
counts pages in mappings whose page size exceeds 4 kB.

`smaps.tsv`: `rss_kb`, `pss_kb`, `pss_anon_kb`, `pss_file_kb`, `pss_dirty_kb`,
`anonymous_kb`, `anon_huge_kb`, `private_dirty_kb`, `private_clean_kb`,
`shared_clean_kb`, `swap_kb` from `/proc/self/smaps_rollup`.

`mallocinfo.tsv`: glibc's `malloc_info`. `arenas` (the `<heap>` count), then
the process-wide totals in bytes: `fast_free_b` (free in fastbins),
`rest_free_b` (free in the other bins), `mmap_count` and `mmap_b` (chunks served
by `mmap`), `system_current_b` and `system_max_b` (memory the arenas hold from
the system, now and at most), `aspace_total_b`, `aspace_mprotect_b`, and
`inuse_b` = `system_current_b - fast_free_b - rest_free_b` (arena memory in
use; `mmap_b` comes on top).

### topology.tsv: once

`cpu`, `node`, `package`, `core_id`, `l3_id` (from
`/sys/devices/system/cpu/cpu*/cache/index3/id`), and `allowed` (1 if the
process's affinity mask contains the CPU). Join it with `threads.tsv` on
`cpu` to see which threads shared an L3 or crossed a NUMA node.

### meta.tsv: key, value

The command line, pid, host, kernel, glibc and ONNX Runtime versions, the
sampler period and snapshot interval, `ODIA_*` and allocator environment
variables, the session count of each predictor, the ONNX profile files written,
the sample count and the sampler's mean and largest cost per sample
(`sample_ms_mean`, `sample_ms_max`), the snapshots taken and skipped and their
summed cost (`snapshots_total_ms`), `end_s` and the event count.

## A first summary

`scripts/perf/trace_phase_table.py RUNDIR...` reads `RUNDIR/trace/` (and
`RUNDIR/time.txt` from `/usr/bin/time -v`, if present). It writes
`phase_table.tsv`, with wall time, CPU, mean active cores, mean running threads
and peak RSS per phase, `predictor_table.tsv`, with session utilisation and
the serial time around the predictor calls per model, and `tracer_cost.tsv`,
with the snapshots' wall time and the sampler's CPU. Both have their columns
documented at the top of the script. It needs only the Python standard library.
