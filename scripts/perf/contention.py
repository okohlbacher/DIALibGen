#!/usr/bin/env python3
"""contention.py CONFIGFILE OUT_TSV -- run ortbench placement/contention experiments on kim, serially.

Run the whole script under: flock /scratch/kohlbach/bid10.lock python3 tools/contention.py ...
Config lines: label|kind|rows|len|K|placement|tpi|extra
  K instances (processes) of tpi worker threads each (tpi >= 4, owner rule), each thread pinned to one CPU.
  placement: packed  -- instances fill CCDs two at a time (4 cores each, CCD = 8 cores), CCD 0,1,2,...
             spread  -- one instance per CCD, CCDs 0,1,2,... (socket 0 first: CCDs 0-11 are node 0)
             xsock   -- one instance per CCD, alternating sockets (CCD 0,12,1,13,...)
             smt     -- instance pairs share the same 4 physical cores (2nd instance on SMT siblings, +192)
             ccd8    -- tpi threads on one CCD (cores 0..tpi-1), K must be 1
             ccdspread8 -- tpi threads, one per CCD (cores 0,8,16,...), K must be 1
  extra: extra ortbench arguments (e.g. "--shared-session 1").
Output: one row per instance (its 'total' line) plus run metadata: start, uptime before/after, node busy cores
during the run (all processes, /proc/stat), own busy (sum of instances), mean scaling_cur_freq of the used CPUs
sampled at mid-run.
"""
import os, subprocess, sys, time

A = '/scratch/kohlbach/bid10-memory-cache'
MODELS = '/ceph/ibmi/abi/oliver/opt/dialibgen-0.10.1/share/DIALibGen/models'
MODEL = {'ms2': 'peptdeep_ms2_dynamic.onnx', 'ccs': 'peptdeep_ccs_dynamic.onnx', 'rt': 'peptdeep_rt_dynamic.onnx'}
SECONDS = float(os.environ.get('SECONDS_PER_RUN', '20'))


def ccd_cores(c):
    return list(range(8 * c, 8 * c + 8))


def placement(kind, K, tpi):
    insts = []
    if kind == 'packed':
        for i in range(K):
            base = 8 * (i // 2) + 4 * (i % 2)
            insts.append(list(range(base, base + 4)))
    elif kind == 'spread':
        insts = [ccd_cores(i)[:4] for i in range(K)]
    elif kind == 'xsock':
        order = [c for pair in zip(range(12), range(12, 24)) for c in pair]
        insts = [ccd_cores(order[i])[:4] for i in range(K)]
    elif kind == 'smt':
        for i in range(K):
            base = 8 * (i // 4) + 4 * ((i // 2) % 2)
            cores = list(range(base, base + 4))
            insts.append(cores if i % 2 == 0 else [c + 192 for c in cores])
    elif kind == 'ccd8':
        insts = [list(range(0, tpi))]
    elif kind == 'ccdspread8':
        insts = [[8 * j for j in range(tpi)]]
    else:
        raise SystemExit('unknown placement ' + kind)
    if tpi != 4 and kind in ('packed', 'spread', 'xsock', 'smt'):
        raise SystemExit('packed/spread/xsock/smt use 4 threads per instance')
    return insts


def stat_busy():
    with open('/proc/stat') as f:
        x = f.readline().split()
    v = list(map(int, x[1:9]))
    return sum(v), v[3] + v[4]


def freq(cpus):
    vals = []
    for c in cpus:
        try:
            with open(f'/sys/devices/system/cpu/cpu{c}/cpufreq/scaling_cur_freq') as f:
                vals.append(int(f.read()) / 1000.0)
        except OSError:
            pass
    return sum(vals) / len(vals) if vals else -1


def uptime():
    return subprocess.run(['uptime'], capture_output=True, text=True).stdout.strip()


def main():
    cfg, out = sys.argv[1], sys.argv[2]
    ncpu = os.cpu_count()
    done = set()
    if os.path.exists(out):
        for line in open(out):
            done.add(line.split('\t')[0])
    else:
        with open(out, 'w') as f:
            f.write('label\tkind\trows\tlen\tK\tplacement\ttpi\textra\tinstance\tcpus\trows_per_s\tper_thread_rows_per_s\t'
                    'cpu_s\trss_sessions_kb\trss_warm_kb\trss_end_kb\tfreq_mhz_mid\tnode_busy_cores\tother_busy_cores\t'
                    'load1_before\tload1_after\tstart_utc\n')
    for line in open(cfg):
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        label, kind, rows, ln, K, plc, tpi, extra = (line.split('|') + [''])[:8]
        if label in done:
            print('skip', label, flush=True)
            continue
        K, tpi = int(K), int(tpi)
        insts = placement(plc, K, tpi)
        start_at = time.time() + 3.0 + 0.05 * K          # sessions + warm-up finish before this
        procs = []
        lb = open('/proc/loadavg').read().split()[0]
        up_b = uptime()
        for i, cpus in enumerate(insts):
            cmd = [A + '/bin/ortbench', '--model', f'{MODELS}/{MODEL[kind]}', '--kind', kind, '--threads', str(tpi),
                   '--rows', rows, '--len', ln, '--seconds', str(SECONDS), '--cpus', ','.join(map(str, cpus)),
                   '--start-at', f'{start_at:.3f}', '--label', f'{label}.{i}'] + extra.split()
            procs.append(subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True))
        # node busy over the measured window
        time.sleep(max(0.0, start_at - time.time()) + 1.0)
        t0, i0 = stat_busy()
        time.sleep(SECONDS / 2 - 1.0)
        fm = [freq(c) for c in insts]
        time.sleep(SECONDS / 2 - 1.5)
        t1, i1 = stat_busy()
        node_busy = (1 - (i1 - i0) / (t1 - t0)) * ncpu if t1 > t0 else -1
        outs = [p.communicate()[0] for p in procs]
        la = open('/proc/loadavg').read().split()[0]
        own = K * tpi
        with open(out, 'a') as f:
            for i, o in enumerate(outs):
                tot = [l for l in o.splitlines() if '\ttotal\t' in l]
                if not tot:
                    print('FAILED', label, i, o[-500:], flush=True)
                    continue
                c = tot[0].split('\t')
                # label kind rows len threads share shared_session total -1 rows_per_s -1 wall cpu create rss0 rss_s rss_w rss_e
                f.write('\t'.join([label, kind, rows, ln, str(K), plc, str(tpi), extra or '-', str(i),
                                   ','.join(map(str, insts[i])), c[9], f'{float(c[9]) / tpi:.2f}', c[12], c[15], c[16], c[17],
                                   f'{fm[i]:.0f}', f'{node_busy:.1f}', f'{node_busy - own:.1f}', lb, la,
                                   time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime(start_at))]) + '\n')
        rates = []
        for o in outs:
            for l in o.splitlines():
                if '\ttotal\t' in l:
                    rates.append(float(l.split('\t')[9]))
        print(f'{label} K={K} {plc} rows={rows} mean/inst={sum(rates)/max(len(rates),1):.1f} node_busy={node_busy:.1f} '
              f'load {lb}->{la} | {up_b}', flush=True)
        time.sleep(2)
    print('CONTENTION_DONE', time.strftime('%FT%TZ', time.gmtime()), flush=True)


if __name__ == '__main__':
    main()
