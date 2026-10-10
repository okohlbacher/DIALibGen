#!/usr/bin/env python3
"""bid10-threads: CPU, wall and thread behaviour per phase from ODIA_TRACE run directories.

usage: threads_analysis.py OUTDIR RUNDIR...   (RUNDIR contains trace/ and optionally time.txt, meta.txt)

Writes (all tab-separated, header line; schema in OUTDIR/SCHEMA_threads.txt):
  occupancy_timeline.tsv   run, t_s, phase, dt_s, cores, n_run, threads, vcsw_per_s, nvcsw_per_s, rss_gib
  phase_cores.tsv          run, phase, wall_s, cpu_s, mean_cores, p05..p95 cores (time-weighted), frac_le_1p5,
                           frac_ge_cap (cores >= 0.9*sessions), vcsw, nvcsw, vcsw_per_s, nvcsw_per_s,
                           runq_wait_s (schedstat wait summed over threads), runq_wait_per_cpu_s, threads_max
  serial_sections.tsv      run, threads, section, kind, wall_s, cpu_s, share_wall, n
  dispatch_balance.tsv     run, model, calls, workers_mean, chunks_mean, dispatch_wall_s, chunk_busy_s,
                           capacity_core_s, idle_core_s, start_lag_core_s, tail_idle_core_s, balanced_wall_s,
                           tail_wall_s, util
  dispatch_calls.tsv       run, model, call, workers, chunks, rows, wall_s, busy_s, max_chunk_s, tail_s, util
  chunk_sizes.tsv          run, model, rows_bin, chunks, rows_sum, busy_s
  thread_classes.tsv       run, class, threads, cpu_s, vcsw, nvcsw, runq_wait_s
"""
import csv, os, sys, math, re
from collections import defaultdict

def rd(path):
    with open(path) as f:
        r = csv.reader(f, delimiter='\t')
        h = next(r)
        return [dict(zip(h, row)) for row in r]

def wq(vals, ws, q):
    """time-weighted quantile"""
    pairs = sorted(zip(vals, ws))
    tot = sum(w for _, w in pairs)
    if tot <= 0: return float('nan')
    acc = 0.0
    for v, w in pairs:
        acc += w
        if acc >= q * tot: return v
    return pairs[-1][0]

def phase_intervals(ev):
    ph = []
    for e in ev:
        if e['level'] == '0':
            ph.append((float(e['t0_s']), float(e['t1_s']), e['name']))
    ph.sort()
    return ph

def phase_at(ph, t):
    # binary search
    lo, hi = 0, len(ph) - 1
    while lo <= hi:
        m = (lo + hi) // 2
        a, b, n = ph[m]
        if t < a: hi = m - 1
        elif t >= b: lo = m + 1
        else: return n
    return 'gap'

def cls(comm):
    if comm.startswith('odia-MS2'): return 'MS2 workers (s>=1)'
    if comm.startswith('odia-RT'): return 'RT workers (s>=1)'
    if comm.startswith('odia-CCS'): return 'CCS workers (s>=1)'
    if comm.startswith('odia-sampler'): return 'tracer sampler'
    if comm.startswith('jemalloc'): return 'jemalloc background'
    return 'main + unnamed (OpenBLAS pool, ORT pools)'

def analyse(run_dir, out):
    run = os.path.basename(run_dir.rstrip('/'))
    T = os.path.join(run_dir, 'trace')
    ev = rd(os.path.join(T, 'events.tsv'))
    proc = rd(os.path.join(T, 'proc.tsv'))
    meta = {r['key']: r['value'] for r in rd(os.path.join(T, 'meta.tsv'))}
    m = re.search(r'-threads (\d+)', meta.get('argv', ''))
    threads_arg = int(m.group(1)) if m else -1
    sessions = int(meta.get('predictor.MS2.sessions', '16'))
    intra = int(meta.get('predictor.MS2.intra_op', '1') or 1)
    cap = sessions * max(1, intra)
    ph = phase_intervals(ev)

    # ---------- occupancy timeline from proc.tsv ----------
    rows = []
    prev = None
    for p in proc:
        t = float(p['t_s'])
        cpu = float(p['utime_s']) + float(p['stime_s'])
        if prev is not None:
            dt = t - prev[0]
            if dt > 0:
                ph_name = p['phase'] if p['phase'] != '-' else 'gap'
                rows.append(dict(run=run, t_s=round(t, 3), phase=ph_name, dt_s=round(dt, 4),
                                 cores=round((cpu - prev[1]) / dt, 3), n_run=int(p['n_run']), threads=int(p['threads']),
                                 vcsw_per_s=round((int(p['vcsw']) - prev[2]) / dt, 1),
                                 nvcsw_per_s=round((int(p['nvcsw']) - prev[3]) / dt, 1),
                                 rss_gib=round(int(p['rss_kb']) / 1048576, 3)))
        prev = (t, cpu, int(p['vcsw']), int(p['nvcsw']))
    out['occupancy_timeline'].extend(rows)

    # ---------- per-thread run-queue wait / csw per phase (threads.tsv deltas) ----------
    runq = defaultdict(float); runtime = defaultdict(float)
    tcls = defaultdict(lambda: dict(threads=0, cpu_s=0.0, vcsw=0, nvcsw=0, runq_wait_s=0.0))
    last = {}
    thr_path = os.path.join(T, 'threads.tsv')
    with open(thr_path) as f:
        r = csv.reader(f, delimiter='\t'); h = next(r); ix = {k: i for i, k in enumerate(h)}
        for row in r:
            tid = row[ix['tid']]; t = float(row[ix['t_s']])
            w = float(row[ix['wait_s']]); rn = float(row[ix['run_s']])
            cpu = float(row[ix['utime_s']]) + float(row[ix['stime_s']])
            v = int(row[ix['vcsw']]); nv = int(row[ix['nvcsw']])
            if tid in last:
                pw, prn = last[tid][0], last[tid][1]
                pn = phase_at(ph, t)
                runq[pn] += max(0.0, w - pw); runtime[pn] += max(0.0, rn - prn)
            else:
                c = cls(row[ix['comm']])
                tcls[c]['threads'] += 1
            last[tid] = (w, rn, cpu, v, nv, cls(row[ix['comm']]))
    for tid, (w, rn, cpu, v, nv, c) in last.items():
        tcls[c]['cpu_s'] += cpu; tcls[c]['vcsw'] += v; tcls[c]['nvcsw'] += nv; tcls[c]['runq_wait_s'] += w
    for c, d in sorted(tcls.items()):
        out['thread_classes'].append(dict(run=run, cls=c, threads=d['threads'], cpu_s=round(d['cpu_s'], 2),
                                          vcsw=d['vcsw'], nvcsw=d['nvcsw'], runq_wait_s=round(d['runq_wait_s'], 3)))

    # ---------- per-phase core distribution ----------
    byp = defaultdict(list)
    for r_ in rows: byp[r_['phase']].append(r_)
    # proc counters at phase boundaries
    pcs = [(float(p['t_s']), int(p['vcsw']), int(p['nvcsw'])) for p in proc]
    def csw_at(t):
        best = pcs[0]
        for x in pcs:
            if x[0] <= t: best = x
            else: break
        return best
    phase_wall = defaultdict(float); phase_cpu = defaultdict(float)
    for a, b, n in ph:
        phase_wall[n] += b - a
    for r_ in rows: phase_cpu[r_['phase']] += r_['cores'] * r_['dt_s']
    order = []
    for a, b, n in ph:
        if n not in order: order.append(n)
    for n in order:
        rs = byp.get(n, [])
        wall = phase_wall[n]
        v0 = n1 = 0
        vc = nvc = 0
        for a, b, nn in ph:
            if nn != n: continue
            x0 = csw_at(a); x1 = csw_at(b)
            vc += x1[1] - x0[1]; nvc += x1[2] - x0[2]
        if rs:
            vals = [x['cores'] for x in rs]; ws = [x['dt_s'] for x in rs]
            tw = sum(ws)
            cpu = sum(v * w for v, w in zip(vals, ws))
            q = {k: wq(vals, ws, k / 100) for k in (5, 25, 50, 75, 95)}
            le = sum(w for v, w in zip(vals, ws) if v <= 1.5) / tw
            ge = sum(w for v, w in zip(vals, ws) if v >= 0.9 * cap) / tw
            thmax = max(x['threads'] for x in rs)
        else:
            tw = cpu = 0; q = {k: float('nan') for k in (5, 25, 50, 75, 95)}; le = ge = float('nan'); thmax = -1
        out['phase_cores'].append(dict(run=run, threads=threads_arg, sessions=sessions, intra_op=intra, phase=n,
            wall_s=round(wall, 3), cpu_s=round(cpu, 2), mean_cores=round(cpu / tw, 3) if tw else '',
            p05=q[5], p25=q[25], p50=q[50], p75=q[75], p95=q[95], frac_le_1p5=round(le, 3), frac_ge_cap=round(ge, 3),
            vcsw=vc, nvcsw=nvc, vcsw_per_s=round(vc / wall, 1) if wall > 0 else '', nvcsw_per_s=round(nvc / wall, 1) if wall > 0 else '',
            runq_wait_s=round(runq.get(n, 0.0), 3), runq_wait_per_cpu_s=round(runq.get(n, 0.0) / runtime[n], 4) if runtime.get(n) else '',
            threads_max=thmax))

    # ---------- serial sections ----------
    total_wall = float(meta.get('end_s', ph[-1][1] if ph else 0))
    tt = os.path.join(run_dir, 'time.txt')
    if os.path.exists(tt):
        for line in open(tt):
            if 'Elapsed' in line:
                v = line.split(': ', 1)[1].strip().split(':')
                total_wall = sum(float(x) * 60 ** (len(v) - 1 - i) for i, x in enumerate(v))
    proc_start = -float(meta.get('process_start_s', '0'))
    ser = defaultdict(lambda: [0.0, 0.0, 0])
    def add(k, kind, d, c):
        ser[(k, kind)][0] += d; ser[(k, kind)][1] += max(c, 0); ser[(k, kind)][2] += 1
    add('process start to trace start (exec, dynamic loading, OpenBLAS pool, OpenMS init)', 'serial', proc_start,
        float(proc[0]['utime_s']) + float(proc[0]['stime_s']) if proc else 0)
    serial_phases = {'setup', 'fasta_read', 'digest', 'enumerate', 'release', 'irt_fit', 'decoys', 'fingerprint',
                     'write', 'teardown', 'trace_snapshot'}
    for e in ev:
        lv, n = e['level'], e['name']; d = float(e['dur_s']); c = float(e['cpu_s'])
        if lv == '0' and n in serial_phases:
            add(n, 'serial phase', d, c)
        elif lv == '1' and e['phase'] in ('rt', 'ms2', 'ccs') and n in ('unique_parse', 'sessions', 'assign', 'reserve', 'finalize'):
            add(f"{e['phase']}.{n}", 'serial step', d, c)
        elif lv == '2' and n in ('parse', 'post'):
            add(f"{e['phase']}.block.{n}", 'serial step (per block, between predictor calls)', d, c)
        elif lv == '3' and n == 'chunking':
            add(f"{e['phase']}.chunking (length grouping)", 'serial step (per call)', d, c)
    # phase remainder: rt/ms2/ccs wall minus level-1 steps (predictor teardown etc.)
    for pn in ('rt', 'ms2', 'ccs'):
        pw = phase_wall.get(pn, 0)
        covered = sum(float(e['dur_s']) for e in ev if e['level'] == '1' and e['phase'] == pn)
        if pw: add(f'{pn}.uncovered (predictor/session teardown, frees)', 'serial remainder', pw - covered, 0)
    gaps = 0.0
    for i in range(1, len(ph)):
        g = ph[i][0] - ph[i - 1][1]
        if g > 0: gaps += g
    add('between phases (logging)', 'serial', gaps, 0)
    tracer_end = total_wall - proc_start - (ph[-1][1] if ph else 0)
    add('after last phase to exit (process teardown)', 'serial', max(0.0, tracer_end), 0)
    for (k, kind), (d, c, n) in sorted(ser.items(), key=lambda x: -x[1][0]):
        out['serial_sections'].append(dict(run=run, threads=threads_arg, section=k, kind=kind, wall_s=round(d, 3),
                                           cpu_s=round(c, 3), share_wall=round(d / total_wall, 4) if total_wall else '', n=n))

    # ---------- dispatch balance ----------
    calls = defaultdict(lambda: dict(model='-', workers=0, chunks=0, rows=0, t0=0, t1=0, busy=0.0, maxc=0.0,
                                     wstart=[], wend=[]))
    for e in ev:
        if e['level'] != '3': continue
        n = e['name']
        if n not in ('dispatch', 'worker', 'chunk'): continue
        k = (e['model'], e['call'])
        c = calls[k]; c['model'] = e['model']
        if n == 'dispatch':
            c['t0'] = float(e['t0_s']); c['t1'] = float(e['t1_s']); c['workers'] = int(e['session'])
            c['chunks'] = int(e['idx']); c['rows'] = int(e['rows'])
        elif n == 'worker':
            c['wstart'].append(float(e['t0_s'])); c['wend'].append(float(e['t1_s']))
        else:
            d = float(e['dur_s']); c['busy'] += d; c['maxc'] = max(c['maxc'], d)
    agg = defaultdict(lambda: defaultdict(float))
    for (model, call), c in sorted(calls.items(), key=lambda x: (x[0][0], int(x[0][1]))):
        W = c['t1'] - c['t0']; k = max(1, c['workers'])
        if W <= 0: continue
        cap_cs = k * W
        lag = sum(s - c['t0'] for s in c['wstart'])
        tail = sum(c['t1'] - e_ for e_ in c['wend'])
        tail_wall = c['t1'] - (min(c['wend']) if c['wend'] else c['t1'])
        bal = max(c['busy'] / k, c['maxc'])
        out['dispatch_calls'].append(dict(run=run, model=model, call=call, workers=k, chunks=c['chunks'], rows=c['rows'],
                                          wall_s=round(W, 4), busy_s=round(c['busy'], 4), max_chunk_s=round(c['maxc'], 4),
                                          tail_s=round(tail_wall, 4), util=round(c['busy'] / cap_cs, 4)))
        a = agg[model]
        a['calls'] += 1; a['workers'] += k; a['chunks'] += c['chunks']; a['wall'] += W; a['busy'] += c['busy']
        a['cap'] += cap_cs; a['lag'] += lag; a['tail'] += tail; a['bal'] += bal; a['tailw'] += tail_wall
    for model, a in agg.items():
        out['dispatch_balance'].append(dict(run=run, model=model, calls=int(a['calls']),
            workers_mean=round(a['workers'] / a['calls'], 2), chunks_mean=round(a['chunks'] / a['calls'], 2),
            dispatch_wall_s=round(a['wall'], 3), chunk_busy_s=round(a['busy'], 2), capacity_core_s=round(a['cap'], 2),
            idle_core_s=round(a['cap'] - a['busy'], 2), start_lag_core_s=round(a['lag'], 3),
            tail_idle_core_s=round(a['tail'], 2), balanced_wall_s=round(a['bal'], 3), tail_wall_s=round(a['tailw'], 3),
            util=round(a['busy'] / a['cap'], 4)))

    # ---------- chunk sizes ----------
    bins = [1, 2, 8, 32, 128, 512, 1024, 2048, 10 ** 9]
    cs = defaultdict(lambda: [0, 0, 0.0])
    for e in ev:
        if e['level'] == '3' and e['name'] == 'chunk':
            rws = int(e['rows'])
            for i in range(len(bins) - 1):
                if bins[i] <= rws < bins[i + 1]:
                    lab = f"{bins[i]}-{bins[i + 1] - 1}" if bins[i + 1] < 10 ** 9 else f">={bins[i]}"
                    break
            x = cs[(e['model'], i, lab)]; x[0] += 1; x[1] += rws; x[2] += float(e['dur_s'])
    for (model, i, lab), (n, rsum, b) in sorted(cs.items()):
        out['chunk_sizes'].append(dict(run=run, model=model, rows_bin=lab, chunks=n, rows_sum=rsum, busy_s=round(b, 2)))

def main():
    od = sys.argv[1]; os.makedirs(od, exist_ok=True)
    out = defaultdict(list)
    for rdir in sys.argv[2:]:
        if not os.path.exists(os.path.join(rdir, 'trace', 'events.tsv')):
            print('skip', rdir, file=sys.stderr); continue
        analyse(rdir, out)
    for k, rows in out.items():
        if not rows: continue
        keys = list(rows[0].keys())
        with open(os.path.join(od, k + '.tsv'), 'w') as f:
            f.write('\t'.join(keys) + '\n')
            for r_ in rows:
                f.write('\t'.join('' if r_.get(x) is None else (f"{r_[x]:.4g}" if isinstance(r_[x], float) else str(r_[x])) for x in keys) + '\n')
        print('wrote', k, len(rows))

if __name__ == '__main__':
    main()
