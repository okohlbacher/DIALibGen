#!/usr/bin/env python3
"""One row per run: configuration, elapsed, CPU, memory, output md5, per-phase wall, predictor efficiency.

usage: experiments_table.py OUT.tsv DISPATCH_BALANCE.tsv RUNDIR...
Columns: label, input, binary, threads, sessions, intra_op, spin, batch_rows, ms2_block, prefix, elapsed_s, user_s, sys_s,
cpu_s, mean_cores, maxrss_gib, vcsw, nvcsw, out_md5, same_as_ref, rt_s, ms2_s, ccs_s, write_s, serial_pre_s,
ms2_cpu_ms_per_row (chunk wall summed / rows = session-ms per row), ms2_util, ms2_workers_mean (threads that got work per call),
ms2_busy_workers (workers_mean x util), ccs_util, rt_util, load1_max, other_cores_mean
same_as_ref: md5 equals the trace-build baseline of the same input (sub10 8d8f2d22..., full c7d404a0..., tiny 3998957f...).
"""
import csv, os, re, sys
REF = {'sub10': '8d8f2d22f588de49870c4a3206e29d35', 'proteome': 'c7d404a0bd92af1469a62d9de32b5ad0',
       'tiny': '3998957f3344200fe9b6647a322937ac'}
REF_TSV = 'e69f4ec931c9941bdfb93523d42879b9'   # full proteome, TSV, trace build

def kv(path):
    d = {}
    if not os.path.exists(path): return d
    for line in open(path):
        line = line.rstrip('\n')
        if line.startswith('time: '):
            k, _, v = line[6:].strip().rpartition(': '); d['time:' + k.strip()] = v.strip()
        else:
            k, _, v = line.partition(' '); d.setdefault(k, v)
    return d

def secs(s):
    v = s.split(':'); return sum(float(x) * 60 ** (len(v) - 1 - i) for i, x in enumerate(v))

def main():
    out, dbp = sys.argv[1], sys.argv[2]
    db = {}
    for r in csv.DictReader(open(dbp), delimiter='\t'):
        db[(r['run'], r['model'])] = r
    rows = []
    for d in sys.argv[3:]:
        if os.path.exists(os.path.join(d, 'SKIPPED')) or not os.path.exists(os.path.join(d, 'time.txt')): continue
        m = kv(os.path.join(d, 'meta.txt')); lab = os.path.basename(d.rstrip('/'))
        cmd = m.get('cmd', ''); env = m.get('env', '')
        inp = 'sub10' if 'sub10' in cmd else 'proteome' if 'proteome' in cmd else 'tiny' if 'tiny' in cmd else '?'
        g = lambda k, dflt: (re.search(k + r'=(\S+)', env).group(1) if re.search(k + r'=(\S+)', env) else dflt)
        th = int(re.search(r'-threads (\d+)', cmd).group(1)) if re.search(r'-threads (\d+)', cmd) else -1
        sess = g('ODIA_MAX_SESSIONS', '16'); sess = str(min(th, int(sess)))
        prefix = 'numactl' if 'numactl' in cmd else 'taskset' if 'taskset' in cmd else ''
        el = secs(m.get('time:Elapsed (wall clock) time (h:mm:ss or m:ss)', '0'))
        us = float(m.get('time:User time (seconds)', 0)); sy = float(m.get('time:System time (seconds)', 0))
        rss = float(m.get('time:Maximum resident set size (kbytes)', 0)) / 1048576
        ph = {}
        ev = os.path.join(d, 'trace', 'events.tsv')
        if os.path.exists(ev):
            for e in csv.DictReader(open(ev), delimiter='\t'):
                if e['level'] == '0': ph[e['name']] = ph.get(e['name'], 0) + float(e['dur_s'])
        rowsum = {}
        for model in ('MS2', 'CCS', 'RT'):
            rowsum[model] = db.get((lab, model))
        ms2 = rowsum['MS2']
        if 'diann' in cmd:   # DIA-NN: the DIALibGen knobs do not apply
            g = lambda k, dflt: '-'
            sess = '-'
        rows.append(dict(label=lab, input=inp, binary=os.path.basename(cmd.split()[0]) if cmd else '', threads=th, sessions=sess,
            intra_op=g('ODIA_INTRA_OP', '1'), spin=g('ODIA_SPIN', 'ort-default'), batch_rows=g('ODIA_BATCH_ROWS', '2048'),
            ms2_block=g('ODIA_MS2_BLOCK', '20000'), prefix=prefix, elapsed_s=round(el, 2), user_s=us, sys_s=sy,
            cpu_s=round(us + sy, 1), mean_cores=round((us + sy) / el, 2) if el else '', maxrss_gib=round(rss, 2),
            vcsw=m.get('time:Voluntary context switches', ''), nvcsw=m.get('time:Involuntary context switches', ''),
            out_md5=m.get('out_md5', ''), same_as_ref=int(m.get('out_md5', '') == (REF_TSV if '.tsv' in cmd.split(' -out ')[-1].split()[0] else REF.get(inp))) if ' -out ' in cmd else '',
            rt_s=round(ph.get('rt', 0), 2), ms2_s=round(ph.get('ms2', 0), 2), ccs_s=round(ph.get('ccs', 0), 2),
            write_s=round(ph.get('write', 0), 2),
            serial_pre_s=round(sum(ph.get(k, 0) for k in ('fasta_read', 'digest', 'enumerate', 'release')), 2),
            ms2_cpu_ms_per_row='', ms2_util=ms2['util'] if ms2 else '', ms2_workers_mean=ms2['workers_mean'] if ms2 else '',
            ms2_busy_workers=round(float(ms2['workers_mean']) * float(ms2['util']), 2) if ms2 else '',
            ccs_util=rowsum['CCS']['util'] if rowsum['CCS'] else '', rt_util=rowsum['RT']['util'] if rowsum['RT'] else '',
            load1_max=m.get('load1_max', ''), other_cores_mean=m.get('other_cores_mean', '')))
        if ms2 and os.path.exists(ev):
            nrows = sum(int(e['rows']) for e in csv.DictReader(open(ev), delimiter='\t') if e['name'] == 'dispatch' and e['model'] == 'MS2')
            rows[-1]['ms2_cpu_ms_per_row'] = round(1000 * float(ms2['chunk_busy_s']) / nrows, 4) if nrows else ''
    keys = list(rows[0].keys())
    with open(out, 'w') as f:
        f.write('\t'.join(keys) + '\n')
        for r in rows: f.write('\t'.join(str(r[k]) for k in keys) + '\n')
    print('rows', len(rows))

if __name__ == '__main__':
    main()
