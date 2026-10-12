#!/usr/bin/env python3
"""-run end to end on the synthetic run (synthetic_run.h), through the tool.

  identify_e2e.py refine <DIALibGen> <identify_synth_run>
  identify_e2e.py tune   <DIALibGen> <identify_synth_run> <models-dir>

Every search uses search:intensities library: the fixture plants the library's
fragments.

refine: -mode refine -run completes, writes the identification report and a
        refined library that consumed it; the report is byte-identical at
        -threads 1 and 4, at two search:chunk sizes and with
        search:ms1_in_memory on and off (in memory and cached); with pyarrow, the
        report's identifications are checked against the planted truth.
        Ion mobility: -write_im is refused on the run without it (after it is
        read, before anything is searched) and, on a synthetic diaPASEF run
        (two 1/K0 bands per window, the library's 1/K0 mis-calibrated), writes
        the observed 1/K0 of the identified precursors, within 0.01 of the
        planted one; a library with 1/K0 for too few targets is refused, one
        with 1/K0 for most searched with the loss counted; the report is
        byte-identical at -threads 1 and 4, two search:chunk sizes and
        search:ms1_in_memory on and off.
tune:   -mode tune -run and the documented -mode refine -tune -no_filter -run
        complete with the RT head (the synthetic run has no ion mobility),
        write their reports and tuned libraries of the full library size.
Exits 77 (ctest SKIP) when the stock models for tune are missing.
"""
import csv
import filecmp
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


def explain(a, b):
    """Says on stderr how two Parquet files differ: rows, columns, embedded JSON (needs pyarrow)."""
    try:
        import pyarrow.parquet as pq_diag
    except ImportError:
        print('pyarrow unavailable: no detail on the difference', file=sys.stderr)
        return
    def json_diff(x, y, path=''):
        if isinstance(x, dict) and isinstance(y, dict):
            for k in sorted(set(x) | set(y)):
                yield from json_diff(x.get(k), y.get(k), path + '/' + k)
        elif x != y:
            yield '%s: %s | %s' % (path, str(x)[:200], str(y)[:200])
    ta, tb = pq_diag.read_table(a), pq_diag.read_table(b)
    print('%s vs %s: rows %d / %d, table equal %s' % (a, b, ta.num_rows, tb.num_rows, ta.equals(tb)), file=sys.stderr)
    if ta.num_rows == tb.num_rows and not ta.equals(tb):
        for column in ta.column_names:
            if column in tb.column_names and not ta.column(column).equals(tb.column(column)):
                x, y = ta.column(column).to_pylist(), tb.column(column).to_pylist()
                rows = [i for i in range(len(x)) if x[i] != y[i]]
                print('  column %s: %d rows differ, e.g. row %d: %s | %s' % (column, len(rows), rows[0], x[rows[0]], y[rows[0]]),
                      file=sys.stderr)
    ma, mb = ta.schema.metadata or {}, tb.schema.metadata or {}
    for key in sorted(set(ma) | set(mb)):
        if ma.get(key) != mb.get(key):
            try:
                for line in list(json_diff(json.loads(ma.get(key, b'{}')), json.loads(mb.get(key, b'{}'))))[:20]:
                    print('  %s %s' % (key.decode(), line), file=sys.stderr)
            except ValueError:
                print('  %s differs (not JSON)' % key.decode(), file=sys.stderr)

mode, tool, synth = sys.argv[1], sys.argv[2], sys.argv[3]
models = Path(sys.argv[4]) if len(sys.argv) > 4 else None


def fail(message):
    print('FAIL: ' + message, file=sys.stderr)
    sys.exit(1)


def run(*args, ok=True):
    # The synthetic run plants the LIBRARY's fragments (synthetic_run.h), so
    # every search here uses them: search:intensities predicted (the default)
    # would search the model's fragments of these made-up peptides instead.
    # The predicted assays are tested in identify_decoy_exchangeability and
    # identify_pipeline_synthetic.
    if '-run' in [str(a) for a in args]:
        args = args + ('-search:intensities', 'library')
    r = subprocess.run([str(a) for a in args], capture_output=True, text=True, timeout=1800)
    log = r.stdout + r.stderr
    if ok and r.returncode != 0:
        print(log, file=sys.stderr)
        fail('exit %d: %s' % (r.returncode, ' '.join(str(a) for a in args)))
    return log


def rows(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter='\t'))


def lib_im(r):
    """A TSV library row's 1/K0: IonMobility since 0.11.2's DIA-NN column names, IM before."""
    return float(r['IonMobility'] if 'IonMobility' in r else r['IM'])


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
        # Compared resolved: macOS's temporary directory is a symlink (/var -> /private/var).
        if Path(prov['inputs']['reference']).resolve() != ids.resolve() or Path(prov['inputs']['run']).resolve() != mzml.resolve():
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
        # The report's peptide and protein q-values use the same paired
        # competition as its precursor q: refine's default gates keep the
        # identifications the search counted.
        if prov['reference']['passing'] < 0.95 * ident['precursors']:
            fail('refine passed %d of the %d identified precursors' % (prov['reference']['passing'], ident['precursors']))
        print('refine: %d passing report precursors, %d library precursors after' % (prov['reference']['passing'],
                                                                                        prov['library']['after']))

        # A failure after the search keeps the report and says how to reuse it.
        d = root / 'kept'
        d.mkdir()
        kept = d / 'ids.parquet'
        log = run(tool, '-mode', 'refine', '-in', library, '-run', mzml, '-out', d / 'out.tsv', '-out_ids', kept,
                  '-q_precursor', 1e-9, ok=False)
        if not kept.is_file() or 'was kept; rerun with -ids ' + str(kept) not in log:
            print(log, file=sys.stderr)
            fail('a failure after the search must keep the report and say so')
        run(tool, '-mode', 'refine', '-in', library, '-ids', kept, '-out', d / 'reused.tsv')
        print('ok   a failed refine keeps its identification report, and -ids reuses it')

        # The refined Parquet library embeds no timings or thread counts: the
        # same -run command writes the same bytes at any -threads.
        d = root / 'parquet'
        d.mkdir()
        parquets = []
        for threads in (1, 4):   # the same paths both times: the provenance names them
            run(tool, '-mode', 'refine', '-in', library, '-run', mzml, '-out', d / 'out.parquet', '-out_ids', d / 'ids.parquet',
                '-threads', threads)
            kept_as = d / ('out-%d.parquet' % threads)
            (d / 'out.parquet').rename(kept_as)
            (d / 'ids.parquet').rename(d / ('ids-%d.parquet' % threads))
            parquets.append(kept_as)
            (d / 'out.parquet.refine.json').unlink()
        if not filecmp.cmp(parquets[0], parquets[1], shallow=False):
            for name in ('ids', 'out'):
                explain(d / ('%s-1.parquet' % name), d / ('%s-4.parquet' % name))
            fail('the refined Parquet library differs between -threads 1 and 4')
        print('ok   refined Parquet library byte-identical at -threads 1 and 4')

        # Determinism: threads and chunking change nothing in the report.
        _, ids4, _, _ = search('t4', threads=4)
        _, ids4c, _, _ = search('t4c', '-search:chunk', 600, threads=4)
        for other in (ids4, ids4c):
            if not filecmp.cmp(ids, other, shallow=False):
                explain(ids, other)
                fail('%s differs from %s' % (other, ids))
        # -threads 0: the effective CPU budget (the physical cores), said in the log and recorded; the same report.
        _, ids0, prov0, log0 = search('t0', threads=0)
        said = re.search(r'-threads 0: (\d+) threads \(', log0)
        if not said or prov0['search']['resources']['threads'] != int(said.group(1)) or int(said.group(1)) < 1:
            print(log0, file=sys.stderr)
            fail('-threads 0 must say how many threads it means and use them: %s' % json.dumps(prov0['search']['resources']))
        if (os.cpu_count() or 1) >= 4 and int(said.group(1)) < 2:
            fail('-threads 0 means %s thread(s) on a machine with %d hardware threads' % (said.group(1), os.cpu_count()))
        if not filecmp.cmp(ids, ids0, shallow=False):
            explain(ids, ids0)
            fail('the report at -threads 0 differs from that at -threads 1')
        print('ok   -threads 0 means %s threads here, and the report is the same' % said.group(1))
        print('ok   report byte-identical at -threads 1 and 4 and at search:chunk 20000 and 600')

        # search:ms1_in_memory decides how extraction reads the MS1 map, never
        # what it reads: on and off write the same report, for a run read
        # into memory and for a cached one (whose MS1 map is otherwise read
        # from its cache file for every batch). The decision and its numbers
        # are logged and recorded in the sidecar, never in the report.
        held_by_default = prov['search']['resources']['ms1_in_memory']
        if held_by_default.get('mode') != 'auto' or not held_by_default.get('estimate_bytes'):
            fail('the sidecar does not record the default search:ms1_in_memory decision: %s' % json.dumps(held_by_default))
        cached_reference = None
        for name, extra in (('ms1-on', ('-search:ms1_in_memory', 'on')), ('ms1-off', ('-search:ms1_in_memory', 'off')),
                            ('ms1-on-cache', ('-search:ms1_in_memory', 'on', '-search:readoptions', 'cache')),
                            ('ms1-off-cache', ('-search:ms1_in_memory', 'off', '-search:readoptions', 'cache'))):
            _, other, other_prov, other_log = search(name, *extra, threads=4)
            held = other_prov['search']['resources']['ms1_in_memory']
            expect_held = 'on' in extra
            said = 'the MS1 map is held in memory' if expect_held else 'the MS1 map is read per batch'
            if held.get('held') != expect_held or said not in other_log:
                print(other_log, file=sys.stderr)
                fail('%s: the MS1 map must %sbe held in memory, and the log say so: %s'
                     % (name, '' if expect_held else 'not ', json.dumps(held)))
            if 'cache' in extra and held.get('estimate_from') != 'the MS1 cache file':
                fail('%s: a cached run estimates the MS1 map from its cache file: %s' % (name, json.dumps(held)))
            if 'cache' not in extra and held.get('estimate_from') != 'the MS1 peak count':
                fail('%s: a run in memory estimates the MS1 map from its peak count: %s' % (name, json.dumps(held)))
            # The report records the read mode, so a cached search is compared with the other cached one.
            reference = cached_reference if 'cache' in extra else ids
            if reference is None:
                cached_reference = other
            elif not filecmp.cmp(reference, other, shallow=False):
                explain(reference, other)
                fail('%s: the report differs from that of %s' % (name, reference.parent.name))
        print('ok   report byte-identical with search:ms1_in_memory on and off, in memory and cached (default: %s, %s)'
              % ('held' if held_by_default['held'] else 'streamed', held_by_default['reason']))

        # -write_im needs observed 1/K0: refused on a run without ion mobility,
        # once the run is read and before anything is searched.
        d = root / 'no-im'
        d.mkdir()
        log = run(tool, '-mode', 'refine', '-in', library, '-run', mzml, '-out', d / 'out.tsv', '-out_ids', d / 'ids.parquet',
                  '-write_im', ok=False)
        if 'needs observed 1/K0 values, and the run' not in log or 'has no ion mobility' not in log or 'search candidates' in log:
            print(log, file=sys.stderr)
            fail('-write_im on a run without ion mobility must be refused before the search')
        if (d / 'ids.parquet').exists() or (d / 'out.tsv').exists():
            fail('a refused -write_im left output behind')
        print('ok   -write_im refused on a run without ion mobility')

        # ... and accepted on a diaPASEF run, where it writes observed 1/K0.
        # A fifth of the peptides planted, not the 40 % above: an ion-mobility
        # window removes nearly all of this fixture's sparse background, and
        # with half of the searched targets present the run breaks the premise
        # of the run-level guards (search:max_target_fraction, and the label-
        # swap self-check, which then finds the decoys of present targets --
        # measured: 901 "identifications" at 40 %, 0 at 20 % and 10 %). How
        # many the swap check finds at 40 % depends on the 1/K0 window's
        # width, not on an FDR fault: 446 at the automatic 0.040, 373 at 0.1,
        # 0 at 0.3 and 0 without ion mobility, while the product's own scores
        # let no present target's decoy win its pair (0 of 885; M2 review).
        # The design document, "Scoring and FDR", says what that means for
        # rich, clean diaPASEF candidate sets.
        pasef = root / 'pasef'
        print(run(synth, pasef, peptides, '0.2', '20260921', 'im').strip())
        im_truth = {r['precursor_id']: float(r['im']) for r in rows(pasef / 'truth.tsv')}
        d = root / 'im'
        d.mkdir()
        im_out, im_ids = d / 'out.tsv', d / 'ids.parquet'
        # -q_protein 1: at a fifth planted this fixture's protein groups are
        # mostly mixed (planted and not) and none passes the protein gate;
        # this check is about 1/K0, which refine writes per precursor.
        log = run(tool, '-mode', 'refine', '-in', pasef / 'library.tsv', '-run', pasef / 'run.mzML', '-out', im_out,
                  '-out_ids', im_ids, '-write_im', '-q_protein', 1, '-threads', 4)
        im_prov = json.loads(Path(str(im_out) + '.refine.json').read_text())
        im_cal = im_prov['search']['calibration']
        mob = im_cal['detail'].get('ion_mobility') or {}
        if not im_prov['search']['run']['ion_mobility'] or not im_cal['im_window'] or im_cal['im_window'] <= 0 or not mob:
            print(log, file=sys.stderr)
            fail('the diaPASEF run was not searched with its ion mobility: %s' % json.dumps(im_cal))
        written = im_prov['library']['im_written']
        print('ion mobility: 1/K0 window %.3f, calibration %s, %d observed 1/K0 written' %
              (im_cal['im_window'], json.dumps({k: mob.get(k) for k in ('slope', 'intercept', 'inliers', 'window_assignment')}), written))
        if written < 0.7 * len(im_truth):
            fail('%d observed 1/K0 written for %d planted precursors' % (written, len(im_truth)))
        lib_im_of = {r['Precursor.Id']: lib_im(r) for r in rows(im_out) if r['Decoy'] == '0'}
        errors = sorted(abs(lib_im_of[k] - v) for k, v in im_truth.items() if k in lib_im_of and not math.isnan(lib_im_of[k]))
        within = sum(1 for e in errors if e <= 0.01)
        print('ok   -write_im on the diaPASEF run: %d of %d planted precursors within 0.01 of their true 1/K0 (median %.4f)'
              % (within, len(errors), errors[len(errors) // 2] if errors else float('nan')))
        if not errors or within < 0.9 * len(errors):
            fail('the written 1/K0 is not the planted one')

        # A library with no 1/K0 at all (no IM, no CCS) on a diaPASEF run is
        # refused once the run is read, before anything is searched, with a
        # message that names search:im_window -1 -- not by an RT calibration
        # that fails for want of seeds and points at search:allow_bootstrap.
        d = root / 'im-no-library-im'
        d.mkdir()
        with open(pasef / 'library.tsv') as f:
            reader = csv.DictReader(f, delimiter='\t')
            fields, lib_rows = reader.fieldnames, list(reader)
        with open(d / 'library.tsv', 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=[c for c in fields if c not in ('IM', 'IonMobility')], delimiter='\t', lineterminator='\n',
                                    extrasaction='ignore')
            writer.writeheader()
            writer.writerows(lib_rows)
        log = run(tool, '-mode', 'refine', '-in', d / 'library.tsv', '-run', pasef / 'run.mzML', '-out', d / 'out.tsv',
                  '-out_ids', d / 'ids.parquet', '-threads', 4, ok=False)
        if ('but only 0 of %d library targets (0.0 %%) have a 1/K0' % len({r['Precursor.Id'] for r in lib_rows}) not in log
                or 'search:im_window -1' not in log or 'search candidates' in log):
            print(log, file=sys.stderr)
            fail('a library without 1/K0 on a diaPASEF run must be refused before the search, naming search:im_window -1')
        print('ok   a library without 1/K0 is refused on a diaPASEF run before anything is searched')

        # A library with a 1/K0 for most but not all of its targets: below
        # 95 % of them (SearchParams::im_library_min_share) refused before
        # anything is searched, with the numbers; above it searched, and the
        # targets and the searched pairs that are never extracted are warned
        # about with their numbers and recorded in the provenance.
        targets = sorted({r['Precursor.Id'] for r in lib_rows if r.get('Decoy', '0') in ('0', '')})

        def partial(name, every):
            dropped = set(targets[::every])
            dd = root / name
            dd.mkdir()
            with open(dd / 'library.tsv', 'w', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=fields, delimiter='\t', lineterminator='\n')
                writer.writeheader()
                for r in lib_rows:
                    if r['Precursor.Id'] in dropped:
                        im_column = 'IonMobility' if 'IonMobility' in r else 'IM'
                        r = dict(r, **{im_column: ''}, **({'CCS': ''} if 'CCS' in r else {}))
                    writer.writerow(r)
            return dd, len(dropped)

        dd, n = partial('im-partial-refused', 12)
        log = run(tool, '-mode', 'refine', '-in', dd / 'library.tsv', '-run', pasef / 'run.mzML', '-out', dd / 'out.tsv',
                  '-out_ids', dd / 'ids.parquet', '-threads', 4, ok=False)
        if ('%d of %d library targets' % (len(targets) - n, len(targets)) not in log or 'fewer than 95 %' not in log
                or 'search:im_window -1' not in log or 'search candidates' in log):
            print(log, file=sys.stderr)
            fail('a library with 1/K0 for %d of %d targets must be refused before the search, with the numbers'
                 % (len(targets) - n, len(targets)))
        print('ok   a library with 1/K0 for %d of %d targets (%.1f %%) is refused before anything is searched'
              % (len(targets) - n, len(targets), 100.0 * (len(targets) - n) / len(targets)))
        dd, n = partial('im-partial-searched', 40)
        log = run(tool, '-mode', 'refine', '-in', dd / 'library.tsv', '-run', pasef / 'run.mzML', '-out', dd / 'out.tsv',
                  '-out_ids', dd / 'ids.parquet', '-q_protein', 1, '-threads', 4)
        p_prov = json.loads(Path(str(dd / 'out.tsv') + '.refine.json').read_text())
        p_pairs = p_prov['search']['calibration']['detail']['ion_mobility']['pairs']
        p_warn = p_prov['search']['warnings']
        library_warned = [w for w in p_warn if w.startswith('%d of %d library targets' % (n, len(targets)))]
        pairs_warned = [w for w in p_warn if w.startswith('%d of %d searched pairs' % (p_pairs['no_library_im'], p_pairs['searched']))]
        print('     %d of %d targets without 1/K0: %d of %d searched pairs never extracted; warnings %s'
              % (n, len(targets), p_pairs['no_library_im'], p_pairs['searched'], json.dumps(p_warn)))
        if p_pairs['no_library_im'] == 0 or not library_warned or not pairs_warned:
            print(log, file=sys.stderr)
            fail('a library with 1/K0 for %d of %d targets: the targets and the searched pairs without one must be warned '
                 'about with their numbers' % (len(targets) - n, len(targets)))
        print('ok   ... and one with 1/K0 for %d of %d is searched, the pairs never extracted counted and recorded'
              % (len(targets) - n, len(targets)))

        # The written 1/K0 is a MEASUREMENT of where each precursor is, not the
        # prediction it was searched with pulled a little way towards the
        # truth. On the fixture above every precursor sits at its calibrated
        # library 1/K0 (the library is off by a constant only), where any
        # estimate looks right. Here the library's 1/K0 is also off by a
        # per-precursor error (SD 0.025), and the search uses a 0.06-wide
        # 1/K0 window: the written deviation from the calibrated library 1/K0
        # must follow the planted one with slope 1. A 1/K0 read INSIDE the
        # extraction window (OpenSWATH's im_drift, which the report carried
        # until the M2 review) is pulled towards the window's centre: slope
        # well below 1, and dependent on the window width.
        scatter = root / 'pasef-scatter'
        print(run(synth, scatter, peptides, '0.2', '20260921', 'im-scatter').strip())
        s_truth = {r['precursor_id']: float(r['im']) for r in rows(scatter / 'truth.tsv')}
        s_library = {r['Precursor.Id']: lib_im(r) for r in rows(scatter / 'library.tsv')}
        d = root / 'im-scatter'
        d.mkdir()
        s_out = d / 'out.tsv'
        log = run(tool, '-mode', 'refine', '-in', scatter / 'library.tsv', '-run', scatter / 'run.mzML', '-out', s_out,
                  '-out_ids', d / 'ids.parquet', '-write_im', '-q_protein', 1, '-threads', 4, '-search:im_window', 0.06)
        s_prov = json.loads(Path(str(s_out) + '.refine.json').read_text())
        s_mob = s_prov['search']['calibration']['detail']['ion_mobility']
        a, b = s_mob['intercept'], s_mob['slope']
        # Observed values only: a precursor refine wrote no 1/K0 for keeps its
        # library's. One point per precursor (the TSV has a row per fragment).
        pts = []
        s_written = {r['Precursor.Id']: lib_im(r) for r in rows(s_out) if r['Decoy'] == '0'}
        for k, written in sorted(s_written.items()):
            if k not in s_truth or k not in s_library:
                continue
            if math.isnan(written) or abs(written - s_library[k]) < 1e-6:
                continue
            centre = a + b * s_library[k]
            pts.append((s_truth[k] - centre, written - centre, abs(written - s_truth[k])))
        if len(pts) < 100:
            print(log, file=sys.stderr)
            fail('%d observed 1/K0 written on the scattered-library fixture' % len(pts))
        mx = sum(p[0] for p in pts) / len(pts)
        my = sum(p[1] for p in pts) / len(pts)
        sxx = sum((p[0] - mx) ** 2 for p in pts)
        slope = sum((p[0] - mx) * (p[1] - my) for p in pts) / sxx
        spread = math.sqrt(sxx / len(pts))
        err = sorted(p[2] for p in pts)
        print('     1/K0 on the scattered-library fixture: %d written; planted deviation from the calibrated library 1/K0 '
              'SD %.4f; slope of the written on the planted deviation %.3f; |written - planted| median %.4f, p90 %.4f'
              % (len(pts), spread, slope, err[len(err) // 2], err[int(0.9 * (len(err) - 1))]))
        if spread < 0.012:
            fail('the fixture does not scatter the library 1/K0 (SD %.4f)' % spread)
        # Measured before the fix (the report's value was im_drift): slope
        # 0.893, |error| p90 0.0055; after it: 1.000 and 0.0003.
        if not 0.95 <= slope <= 1.05:
            fail('the written 1/K0 is pulled towards the prediction: slope %.3f against the planted deviation' % slope)
        if err[len(err) // 2] > 0.002 or err[int(0.9 * (len(err) - 1))] > 0.003:
            fail('the written 1/K0 is not the planted one (|error| median %.4f, p90 %.4f)'
                 % (err[len(err) // 2], err[int(0.9 * (len(err) - 1))]))

        # Determinism with ion mobility, on this fixture at this window, where
        # reported apexes lie outside the extraction window: -threads 1 and 4
        # and two search:chunk sizes write the same report, byte for byte --
        # the re-measured 1/K0 and its record (search.report_mobility)
        # included. (No pair of this fixture lies beyond every window;
        # identify_ion_mobility checks the nearest-window rule's extraction at
        # 1 and 2 threads and two chunk sizes.)
        s_rm = s_prov['search'].get('report_mobility') or {}
        print('     report 1/K0: %s' % json.dumps({k: s_rm.get(k) for k in ('rows', 'measured', 'outside_extraction_window')}))
        if not s_rm or s_rm['measured'] < 100 or s_rm['outside_extraction_window'] == 0:
            fail('the determinism check below would not cover the re-measurement and apexes outside the window')
        for name, threads, extra in (('im-scatter-t1', 1, ()), ('im-scatter-t4-chunk', 4, ('-search:chunk', 600)),
                                     ('im-scatter-ms1-on', 4, ('-search:ms1_in_memory', 'on')),
                                     ('im-scatter-ms1-off', 4, ('-search:ms1_in_memory', 'off'))):
            dd = root / name
            dd.mkdir()
            run(tool, '-mode', 'refine', '-in', scatter / 'library.tsv', '-run', scatter / 'run.mzML', '-out', dd / 'out.tsv',
                '-out_ids', dd / 'ids.parquet', '-write_im', '-q_protein', 1, '-threads', threads, '-search:im_window', 0.06, *extra)
            if not filecmp.cmp(d / 'ids.parquet', dd / 'ids.parquet', shallow=False):
                fail('with ion mobility the report of %s differs from that at -threads 4, search:chunk 20000' % name)
        # ... and with the run cached, where the MS1 map is otherwise read from its cache file for every batch.
        cached = []
        for name, mode in (('im-scatter-cache-on', 'on'), ('im-scatter-cache-off', 'off')):
            dd = root / name
            dd.mkdir()
            run(tool, '-mode', 'refine', '-in', scatter / 'library.tsv', '-run', scatter / 'run.mzML', '-out', dd / 'out.tsv',
                '-out_ids', dd / 'ids.parquet', '-write_im', '-q_protein', 1, '-threads', 4, '-search:im_window', 0.06,
                '-search:readoptions', 'cache', '-search:ms1_in_memory', mode)
            held = json.loads((dd / 'out.tsv.refine.json').read_text())['search']['resources']['ms1_in_memory']
            if held.get('held') != (mode == 'on') or held.get('estimate_from') != 'the MS1 cache file':
                fail('%s: %s' % (name, json.dumps(held)))
            cached.append(dd / 'ids.parquet')
        if not filecmp.cmp(cached[0], cached[1], shallow=False):
            fail('with ion mobility and a cached run the report differs between search:ms1_in_memory on and off')
        print('ok   with ion mobility the report is byte-identical at -threads 1 and 4, at search:chunk 20000 and 600, and '
              'with search:ms1_in_memory on and off (in memory and cached)')

        # 40 % planted with ion mobility, where the run-level guards' premise
        # (most candidates null) fails -- max_target_fraction and the label-
        # swap self-check would abort (above) -- so they are off here, and
        # the product's own pair competition is checked instead: no planted,
        # present target may lose its pair to its own decoy (the M2 review
        # counted 0 of 885 on this fixture).
        try:
            import pyarrow.parquet as pq_rich
        except ImportError:
            print('pyarrow unavailable: the 40 % pair-winner check is skipped')
        else:
            rich = root / 'pasef-rich'
            print(run(synth, rich, peptides, '0.4', '20260921', 'im').strip())
            rich_truth = {r['precursor_id'] for r in rows(rich / 'truth.tsv')}
            d = root / 'im-rich'
            d.mkdir()
            run(tool, '-mode', 'refine', '-in', rich / 'library.tsv', '-run', rich / 'run.mzML', '-out', d / 'out.tsv',
                '-out_ids', d / 'ids.parquet', '-q_protein', 1, '-threads', 4, '-search:max_target_fraction', 1,
                '-search:selftest', 'false')
            rep = pq_rich.read_table(d / 'ids.parquet').to_pylist()
            won = sum(1 for x in rep if not x['Decoy'] and x['Precursor.Id'] in rich_truth)
            lost = [x['Precursor.Id'] for x in rep if x['Decoy'] and x['Precursor.Id'][:-len('_decoy')] in rich_truth]
            print('     40 %% planted with ion mobility: %d of %d planted targets won their pair in the report, %d lost it to '
                  'their own decoy' % (won, len(rich_truth), len(lost)))
            if won < 0.8 * len(rich_truth):
                fail('40 %% planted: only %d of %d planted targets reported' % (won, len(rich_truth)))
            if len(lost) > max(2, 0.01 * len(rich_truth)):
                fail('40 %% planted: %d present targets lost their pair to their own decoy: %s' % (len(lost), lost[:5]))

        try:
            import pyarrow.parquet as pq
        except ImportError:
            print('pyarrow unavailable: truth checks skipped')
        else:
            table = pq.read_table(ids)
            meta = json.loads(table.schema.metadata[b'odia.identifier'])
            for key in ('chunk', 'batch_size', 'readoptions', 'cache_dir', 'ms1_in_memory'):
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
