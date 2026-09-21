#!/usr/bin/env python3
"""-run end to end on the synthetic run (synthetic_run.h), through the tool.

  identify_e2e.py refine <DIALibGen> <identify_synth_run>
  identify_e2e.py tune   <DIALibGen> <identify_synth_run> <models-dir>

refine: -mode refine -run completes, writes the identification report and a
        refined library that consumed it; the report is byte-identical at
        -threads 1 and 4 and at two search:chunk sizes; with pyarrow, the
        report's identifications are checked against the planted truth.
tune:   -mode tune -run and the documented -mode refine -tune -no_filter -run
        complete with the RT head (the synthetic run has no ion mobility),
        write their reports and tuned libraries of the full library size.
Exits 77 (ctest SKIP) when the stock models for tune are missing.
"""
import csv
import filecmp
import json
import math
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

mode, tool, synth = sys.argv[1], sys.argv[2], sys.argv[3]
models = Path(sys.argv[4]) if len(sys.argv) > 4 else None


def fail(message):
    print('FAIL: ' + message, file=sys.stderr)
    sys.exit(1)


def run(*args, ok=True):
    r = subprocess.run([str(a) for a in args], capture_output=True, text=True, timeout=1800)
    log = r.stdout + r.stderr
    if ok and r.returncode != 0:
        print(log, file=sys.stderr)
        fail('exit %d: %s' % (r.returncode, ' '.join(str(a) for a in args)))
    return log


def rows(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter='\t'))


if mode == 'tune' and not (models and (models / 'peptdeep_rt_dynamic.onnx').is_file()):
    print('SKIP: no stock RT model in %s' % models)
    sys.exit(77)

root = Path(tempfile.mkdtemp(prefix='dialibgen-identify-e2e-'))
try:
    fixture = root / 'fixture'
    # tune needs >= 100 sequence units in each held-out cohort, hence the larger run.
    peptides, planted = ('1200', '0.4') if mode == 'refine' else ('2400', '0.4')
    print(run(synth, fixture, peptides, planted).strip())
    library, mzml = fixture / 'library.tsv', fixture / 'run.mzML'
    truth = {r['precursor_id']: float(r['apex_s']) for r in rows(fixture / 'truth.tsv')}
    n_library = len({r['Precursor.Id'] for r in rows(library)})

    def search(name, *extra, threads=1, command=('-mode', 'refine')):
        d = root / name
        d.mkdir()
        out, ids = d / 'out.tsv', d / 'ids.parquet'
        log = run(tool, *command, '-in', library, '-run', mzml, '-out', out, '-out_ids', ids,
                  '-threads', threads, *extra)
        for f in (out, Path(str(out) + '.refine.json'), ids):
            if not f.is_file() or f.stat().st_size == 0:
                print(log, file=sys.stderr)
                fail('%s: %s was not written' % (name, f.name))
        if 'search report: ' not in log or 'serves as -ids' not in log:
            print(log, file=sys.stderr)
            fail('%s: the search did not hand its report on' % name)
        leftovers = [p.name for p in d.iterdir() if p.name.startswith('.dialibgen-')]
        if leftovers:
            fail('%s: scratch left behind: %s' % (name, leftovers))
        return out, ids, json.loads(Path(str(out) + '.refine.json').read_text()), log

    if mode == 'refine':
        out, ids, prov, log = search('t1')
        s = prov['search']
        if prov['inputs']['reference'] != str(ids.resolve()) or prov['inputs']['run'] != str(mzml.resolve()):
            fail('the sidecar does not name the report as the reference and the run as input')
        if not prov['inputs'].get('run_fnv1a64'):
            fail('no run hash in the sidecar')
        cal = s['calibration']
        if cal['bootstrap'] or cal['points'] < 20 or not cal['rsq'] > 0.95:
            fail('calibration: %s' % json.dumps(cal))
        ident = s['identifications']
        print('search: %d pairs, %d peak groups, %d target precursors at q <= 0.01 (%d decoys), calibration %d points r^2 %.3f, '
              'RT window %.1f s, m/z %.1f ppm' % (s['candidates']['pairs'], s['scoring']['peak_groups'], ident['precursors'],
                                                  ident['decoy_precursors'], cal['points'], cal['rsq'], cal['rt_window_s'],
                                                  cal['mz_ppm']))
        if ident['precursors'] < 0.8 * len(truth):
            fail('%d identifications for %d planted precursors' % (ident['precursors'], len(truth)))
        refined = rows(out)
        if not refined:
            fail('the refined library is empty')
        if prov['reference']['passing'] < 0.8 * len(truth):
            fail('refine used %d passing precursors of the report' % prov['reference']['passing'])
        print('refine: %d passing report precursors, %d library precursors after' % (prov['reference']['passing'],
                                                                                        prov['library']['after']))

        # Determinism: threads and chunking change nothing in the report.
        _, ids4, _, _ = search('t4', threads=4)
        _, ids4c, _, _ = search('t4c', '-search:chunk', 600, threads=4)
        for other in (ids4, ids4c):
            if not filecmp.cmp(ids, other, shallow=False):
                fail('%s differs from %s' % (other, ids))
        print('ok   report byte-identical at -threads 1 and 4 and at search:chunk 20000 and 600')

        try:
            import pyarrow.parquet as pq
        except ImportError:
            print('pyarrow unavailable: truth checks skipped')
        else:
            table = pq.read_table(ids)
            meta = json.loads(table.schema.metadata[b'odia.identifier'])
            for key in ('chunk', 'batch_size', 'readoptions', 'cache_dir'):
                if key in meta['search']['settings']:
                    fail('the report embeds the execution setting %s' % key)
            r = table.to_pylist()
            ids_q = [x for x in r if not x['Decoy'] and x['Q.Value'] <= 0.01]
            false = [x for x in ids_q if x['Precursor.Id'] not in truth]
            # identify_extraction holds recovery to one cycle (0.5 s); here the
            # REPORTED peak group of an identification may sit up to two cycles off.
            off = [x for x in ids_q if x['Precursor.Id'] in truth and abs(60 * x['RT'] - truth[x['Precursor.Id']]) > 1.0]
            print('truth: %d identified at q <= 0.01, %d not planted (%.2f %%), %d planted but off the apex by > 1 s'
                  % (len(ids_q), len(false), 100.0 * len(false) / max(1, len(ids_q)), len(off)))
            if len(false) > 0.02 * len(ids_q) + 3:
                fail('too many identifications of unplanted precursors')
            if len(off) > 0.05 * len(ids_q) + 3:
                fail('too many identifications off the planted apex')
            if not all(math.isnan(x['IM']) if x['IM'] is not None else True for x in r):
                fail('IM must be empty for a run without ion mobility')
    else:
        common = ('-tune_heads', 'rt', '-tune_models', models, '-train:epochs', 20, '-train:warmup', 2,
                  '-stop:min_epochs', 20, '-stop:patience', 100, '-machine:device', 'cpu', '-machine:threads', 2,
                  '-machine:seed', 20260803)
        for name, command in (('tune', ('-mode', 'tune')), ('refine-tune', ('-mode', 'refine', '-tune', '-no_filter'))):
            out, ids, prov, log = search(name, *common, threads=4, command=command)
            if 'rt' not in (prov.get('tune') or {}):
                fail('%s: no RT tuning in the provenance' % name)
            if prov['tune']['rt'].get('repredicted', 0) < 1:
                fail('%s: the tuned RT model re-predicted nothing' % name)
            written = len({r['Precursor.Id'] for r in rows(out)})
            if written != n_library:
                fail('%s: %d precursors written of %d in the library' % (name, written, n_library))
            print('ok   %s: %d target precursors identified, RT tuned on them, %d library precursors re-predicted'
                  % (name, prov['search']['identifications']['precursors'], prov['tune']['rt']['repredicted']))
finally:
    shutil.rmtree(root, ignore_errors=True)
print('identify_e2e %s: PASS' % mode)
