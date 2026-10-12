#!/usr/bin/env python3
"""gantt.py PREFIX [CALL] [OUT.png]: per-batch Gantt of stock OpenSWATH extraction from ODIA_OSW_LOG output.

PREFIX.lines.tsv: t_s, tid, call, line (stock "Thread a_b will analyze N compounds and M transitions from SWATH i
(batch j out of k)", plus CALL_START/CALL_END markers); PREFIX.windows.tsv: idx, ms1, lower, upper, im bounds, spectra.
A batch runs from its line to the next line of the same thread in the same call, or to CALL_END (the last batch of a
thread also holds the window's progress update). Writes PREFIX.batches.tsv (all calls) and PREFIX.windows_calls.tsv
(per call and window: compounds, transitions, spectra, batches, summed batch seconds, first start, last end), prints
per-call summaries (wall, busy thread-seconds, tail after the team falls below 75 %, the longest window and its share),
and with CALL and OUT draws that call's Gantt (one lane per thread, bars coloured by window size).
"""
import sys, re, csv, collections

RX = re.compile(r'Thread (\d+)_(\d+) will analyze (\d+) compounds and (\d+) transitions from SWATH (\d+) \(batch (\d+) out of (\d+)\)')


def load(prefix):
    win = {}
    with open(prefix + '.windows.tsv') as f:
        for r in csv.DictReader(f, delimiter='\t'):
            win[int(r['idx'])] = r
    calls = collections.OrderedDict()
    with open(prefix + '.lines.tsv') as f:
        next(f)
        for l in f:
            t, tid, call, line = l.rstrip('\n').split('\t', 3)
            c = calls.setdefault(int(call), dict(start=None, end=None, order=None, ev=[]))
            if line.startswith('CALL_START'):
                c['start'] = float(t)
                parts = line.split()
                if 'order' in parts:
                    c['order'] = [int(x) for x in parts[parts.index('order') + 1:]]
                c['team'] = int(parts[parts.index('team') + 1])
                c['outer'] = int(parts[parts.index('outer') + 1])
            elif line.startswith('CALL_END'):
                c['end'] = float(t)
            else:
                m = RX.search(line)
                if m:
                    c['ev'].append((float(t), tid, int(m.group(5)), int(m.group(6)), int(m.group(7)), int(m.group(3)), int(m.group(4))))
    return win, calls


LASTBUSY = {}


def load_threads(path, calls):
    """Per (call, tid): the last sample time inside the call at which the thread's CPU advanced by >= 50 % of the
    sampling interval (the end of its last batch; idle team threads sleep in the barrier)."""
    spans = [(ci, c['start'], c['end']) for ci, c in calls.items() if c['start'] is not None and c['end'] is not None]
    prev = {}
    with open(path) as f:
        r = csv.reader(f, delimiter='\t'); h = next(r); ix = {k: i for i, k in enumerate(h)}
        for row in r:
            t = float(row[ix['t_s']]); tid = row[ix['tid']]
            cpu = float(row[ix['utime_s']]) + float(row[ix['stime_s']])
            if tid in prev:
                pt, pc = prev[tid]
                if cpu - pc >= 0.5 * (t - pt) and t > pt:
                    for ci, a, b in spans:
                        if a <= pt and t <= b - 0.05:
                            LASTBUSY[(ci, tid)] = max(LASTBUSY.get((ci, tid), 0.0), t)
            prev[tid] = (t, cpu)


def batches(c, ci=None):
    by = collections.defaultdict(list)
    for e in c['ev']:
        by[e[1]].append(e)
    out = []
    for tid, es in by.items():
        es.sort()
        for k, e in enumerate(es):
            t1 = es[k + 1][0] if k + 1 < len(es) else min(c['end'], max(LASTBUSY.get((ci, tid), c['end']), e[0]))
            out.append(dict(tid=tid, t0=e[0], t1=t1, pos=e[2], batch=e[3], nb=e[4], comp=e[5], trans=e[6]))
    return out


def main():
    prefix = sys.argv[1]
    win, calls = load(prefix)
    import os
    th = os.path.join(os.path.dirname(prefix), 'trace', 'threads.tsv')
    if os.path.exists(th):
        load_threads(th, calls)
    bt = open(prefix + '.batches.tsv', 'w')
    bt.write('call\ttid\tt0_s\tt1_s\tdur_s\twindow\tpos\tbatch\tbatches\tcompounds\ttransitions\tspectra\n')
    wc = open(prefix + '.windows_calls.tsv', 'w')
    wc.write('call\twindow\tpos\tcompounds\ttransitions\tspectra\tbatches\tbusy_s\tfirst_start_s\tlast_end_s\n')
    print('call\tteam\touter\twall_s\tbusy_thread_s\tmean_threads\ttail_s\ttail_share\tlongest_window\tits_s\tits_compounds\tmedian_window_s\tlast_window_end_rel')
    for ci, c in calls.items():
        order = c['order'] or sorted(win)
        bs = batches(c, ci)
        for b in bs:
            w = order[b['pos']]
            b['win'] = w
            bt.write('\t'.join(str(x) for x in [ci, b['tid'], '%.3f' % b['t0'], '%.3f' % b['t1'], '%.3f' % (b['t1'] - b['t0']), w, b['pos'], b['batch'],
                                                  b['nb'], b['comp'], b['trans'], win[w]['spectra']]) + '\n')
        perw = collections.defaultdict(lambda: dict(busy=0.0, n=0, comp=0, trans=0, t0=1e18, t1=0.0, pos=0))
        for b in bs:
            p = perw[b['win']]
            p['busy'] += b['t1'] - b['t0']; p['n'] += 1; p['comp'] = b['comp']; p['trans'] = b['trans']
            p['t0'] = min(p['t0'], b['t0']); p['t1'] = max(p['t1'], b['t1']); p['pos'] = b['pos']
        for w, p in sorted(perw.items(), key=lambda x: -x[1]['busy']):
            wc.write('\t'.join(str(x) for x in [ci, w, p['pos'], p['comp'], p['trans'], win[w]['spectra'], p['n'], '%.2f' % p['busy'],
                                                  '%.2f' % (p['t0'] - c['start']), '%.2f' % (p['t1'] - c['start'])]) + '\n')
        wall = c['end'] - c['start']
        busy = sum(b['t1'] - b['t0'] for b in bs)
        # threads busy over time -> last time >= 0.75 team
        edges = sorted([(b['t0'], 1) for b in bs] + [(b['t1'], -1) for b in bs])
        n, last_full = 0, c['start']
        for t, d in edges:
            n += d
            if n >= 0.75 * c['team']:
                last_full = t
        lw = max(perw.items(), key=lambda x: x[1]['busy'])
        ws = sorted(p['busy'] for p in perw.values())
        print('\t'.join(str(x) for x in [ci, c['team'], c['outer'], '%.1f' % wall, '%.1f' % busy, '%.2f' % (busy / wall), '%.1f' % (c['end'] - last_full),
                                          '%.2f' % ((c['end'] - last_full) / wall), lw[0], '%.1f' % lw[1]['busy'], lw[1]['comp'],
                                          '%.1f' % ws[len(ws) // 2], '%.1f' % (lw[1]['t1'] - c['start'])]))
    bt.close(); wc.close()
    if len(sys.argv) > 3:
        ci = int(sys.argv[2]); out = sys.argv[3]
        c = calls[ci]
        order = c['order'] or sorted(win)
        bs = batches(c, ci)
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib import cm, colors
        tids = sorted({b['tid'] for b in bs}, key=lambda t: min(b['t0'] for b in bs if b['tid'] == t))
        lane = {t: i for i, t in enumerate(tids)}
        comp = [b['comp'] for b in bs]
        norm = colors.Normalize(vmin=0, vmax=max(comp))
        cmap = matplotlib.colormaps['Blues']
        fig, ax = plt.subplots(figsize=(11, 0.32 * len(tids) + 1.6))
        for b in bs:
            w = order[b['pos']]
            ax.barh(lane[b['tid']], b['t1'] - b['t0'], left=b['t0'] - c['start'], height=0.78,
                    color=cmap(0.25 + 0.75 * norm(b['comp'])), edgecolor='white', linewidth=0.8)
            if b['batch'] == 0 and b['comp'] >= 900:
                ax.text(b['t0'] - c['start'] + 0.5, lane[b['tid']], 'window %d: %d compounds' % (w, b['comp']),
                        va='center', fontsize=6.5, color='#ffffff' if b['comp'] > 0.6 * max(comp) else '#0b0b0b')
        ax.set_yticks(range(len(tids))); ax.set_yticklabels(['thread %d' % (i + 1) for i in range(len(tids))], fontsize=7)
        ax.invert_yaxis()
        ax.set_xlabel('seconds since the call started', fontsize=8)
        ax.set_xlim(0, c['end'] - c['start'])
        for s in ('top', 'right'):
            ax.spines[s].set_visible(False)
        ax.tick_params(labelsize=7)
        sm = cm.ScalarMappable(norm=norm, cmap=cmap)
        cb = fig.colorbar(sm, ax=ax, pad=0.01, fraction=0.025)
        cb.set_label('compounds in the window', fontsize=7); cb.ax.tick_params(labelsize=6)
        ax.set_title('%s: extraction call %d of 3, one bar per batch of a window (%.0f s)' % (sys.argv[4] if len(sys.argv) > 4 else '', ci + 1, c['end'] - c['start']),
                     fontsize=9, loc='left')
        fig.savefig(out, dpi=130, bbox_inches='tight')


if __name__ == '__main__':
    main()
