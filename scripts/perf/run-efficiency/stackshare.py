#!/usr/bin/env python3
"""stackshare.py RUNDIR PHASE PATTERN... : share of worker-or-main thread stacks (excluding sampler and idle pool threads
outside OpenSWATH) in PHASE whose frames contain each PATTERN; also the share in OpenMS::OpenSwathWorkflow at all."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import stacks as S
d, phase, pats = sys.argv[1], sys.argv[2], sys.argv[3:]
samples = S.parse(os.path.join(d, 'stacks.txt'))
meta = {}
for l in open(os.path.join(d, 'trace', 'meta.tsv')):
    k, _, v = l.rstrip('\n').partition('\t'); meta.setdefault(k, v)
e0 = float(meta['start_epoch_s'])
ph = []
for l in open(os.path.join(d, 'trace', 'events.tsv')):
    f = l.split('\t')
    if len(f) > 7 and f[5] == '0': ph.append((float(f[0]), float(f[1]), f[7]))
tot = 0; inosw = 0; cnt = {p: 0 for p in pats}; n = 0
for s in samples:
    t = s['t'] - e0
    if not any(a <= t <= b and nm == phase for a, b, nm in ph): continue
    n += 1
    for th in s['threads']:
        if th['name'] == 'odia-sampler' or not th['frames']: continue
        tot += 1
        fr = ' | '.join(th['frames'])
        if 'OpenSwath' in fr or 'MRM' in fr: inosw += 1
        for p in pats:
            if p in fr: cnt[p] += 1
print('%s %s: %d samples, %d thread stacks, %.1f%% inside OpenSWATH/MRM code' % (os.path.basename(os.path.normpath(d)), phase, n, tot, 100.0 * inosw / max(tot, 1)))
for p in pats: print('  %-28s %5.1f%% of stacks, %5.1f%% of those inside OpenSWATH' % (p, 100.0 * cnt[p] / max(tot, 1), 100.0 * cnt[p] / max(inosw, 1)))
