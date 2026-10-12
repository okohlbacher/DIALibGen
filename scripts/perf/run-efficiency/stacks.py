#!/usr/bin/env python3
"""stacks.py RUNDIR [--phase NAME] [--depth N] [--threads]: aggregate gdb stack samples (RUNDIR/stacks.txt from
gdbcell.sh) per level-0 trace phase.

Each thread stack is reduced to
  state  'wait' if a frame is a known blocking/spinning primitive (futex, gomp barrier/critical/mutex,
         pthread_cond_wait, nanosleep, poll, epoll, do_spin), else 'run'
  top    the innermost frame whose function names OpenMS::, ODIA::, Ort, at::/torch, xercesc or a libgomp/libc
         primitive (the 'why')
  app    the first frames (innermost first) inside OpenMS:: or ODIA:: (the 'where'), up to --depth
Prints per phase: samples, threads per sample, then the most common (state, top, app) signatures with counts.
"""
import sys, re, os, collections, bisect

WAIT = ('futex', 'gomp_barrier_wait', 'gomp_team_barrier_wait', 'gomp_mutex_lock', 'GOMP_critical', 'do_spin', 'do_wait',
        'pthread_cond_wait', 'pthread_cond_timedwait', '__lll_lock_wait', 'nanosleep', 'clock_nanosleep', 'poll', 'epoll_wait',
        'sem_wait', 'select', '__GI___poll', 'pthread_mutex_lock', 'gomp_sem_wait', 'gomp_simple_barrier_wait', 'sched_yield',
        'gomp_thread_start')
FRAME = re.compile(r'^#(\d+)\s+(?:0x[0-9a-f]+ in )?(.+?) (?:\(|from |at )')


def fname(s):
    s = s.strip()
    # strip template args and params for brevity
    out, depth = [], 0
    for ch in s:
        if ch == '<':
            depth += 1
            if depth == 1:
                out.append('<>')
            continue
        if ch == '>':
            depth -= 1
            continue
        if depth == 0:
            out.append(ch)
    return ''.join(out)


def parse(path):
    samples = []
    cur = None
    th = None
    for line in open(path, errors='replace'):
        if line.startswith('=== sample'):
            parts = line.split()
            cur = dict(t=float(parts[3]), threads=[])
            samples.append(cur)
            th = None
        elif line.startswith('=== end'):
            cur['t_end'] = float(line.split()[3])
        elif line.startswith('Thread ') and cur is not None:
            m = re.search(r'\(LWP (\d+)\)', line)
            th = dict(lwp=m.group(1) if m else '?', name=line.split('"')[1] if '"' in line else '', frames=[], num=int(line.split()[1]))
            cur['threads'].append(th)
        elif line.startswith('#') and th is not None:
            m = FRAME.match(line)
            if m:
                th['frames'].append(fname(m.group(2)))
    return samples


def reduce(frames, depth):
    state = 'run'
    for f in frames[:6]:
        if any(w in f for w in WAIT):
            state = 'wait'
            break
    top = frames[0] if frames else '?'
    app = [f for f in frames if f.startswith('OpenMS::') or f.startswith('ODIA::') or 'Ort' in f[:4] or f.startswith('at::') or f.startswith('torch::')]
    return state, top, ' < '.join(app[:depth]) if app else '(none)'


def main():
    d = sys.argv[1]
    phase_only = None
    depth = 3
    if '--phase' in sys.argv:
        phase_only = sys.argv[sys.argv.index('--phase') + 1]
    if '--depth' in sys.argv:
        depth = int(sys.argv[sys.argv.index('--depth') + 1])
    samples = parse(os.path.join(d, 'stacks.txt'))
    meta = {}
    for l in open(os.path.join(d, 'trace', 'meta.tsv')):
        k, _, v = l.rstrip('\n').partition('\t')
        meta.setdefault(k, v)
    epoch0 = float(meta['start_epoch_s'])
    ph = []
    for l in open(os.path.join(d, 'trace', 'events.tsv')):
        f = l.split('\t')
        if len(f) > 7 and f[5] == '0':
            ph.append((float(f[0]), float(f[1]), f[7]))
    ph.sort()

    def phase_at(t):
        for a, b, n in ph:
            if a <= t <= b:
                return n
        return '-'

    agg = collections.OrderedDict()
    for s in samples:
        t = s['t'] - epoch0
        p = phase_at(t)
        if phase_only and p != phase_only:
            continue
        a = agg.setdefault(p, dict(n=0, thr=0, sig=collections.Counter(), times=[]))
        a['n'] += 1
        a['times'].append(t)
        for th in s['threads']:
            if th['name'] in ('odia-sampler',) or not th['frames']:
                continue
            if '--main' in sys.argv and th['num'] != 1:
                continue
            if '--workers' in sys.argv and th['num'] == 1:
                continue
            a['thr'] += 1
            a['sig'][reduce(th['frames'], depth)] += 1
    for p, a in agg.items():
        print('## phase %s: %d samples (t %.0f..%.0f s), %.1f threads/sample' % (p, a['n'], min(a['times']), max(a['times']), a['thr'] / a['n']))
        tot = sum(a['sig'].values())
        for (state, top, app), c in a['sig'].most_common(14):
            print('  %5.1f%%  %-4s  %-45.45s  %s' % (100.0 * c / tot, state, top, app[:220]))
        print()


if __name__ == '__main__':
    main()
