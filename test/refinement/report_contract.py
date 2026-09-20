#!/usr/bin/env python3
"""Small report fixtures exercise refinement gates, units, joins and safe output paths."""
import csv
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile

import pyarrow as pa
import pyarrow.parquet as pq

binary = str(Path(sys.argv[1]).resolve())
columns = ['Precursor.Id', 'Modified.Sequence', 'Precursor.Charge', 'Decoy', 'RT', 'IM',
           'Precursor.Mz', 'Product.Mz', 'Relative.Intensity', 'Fragment.Type',
           'Fragment.Charge', 'Fragment.Series.Number', 'Fragment.Loss.Type', 'Protein.Group', 'CCS']

with tempfile.TemporaryDirectory(prefix='dialibgen-report-contract-') as directory:
    root = Path(directory)
    library = root / 'library.tsv'
    with library.open('w', newline='') as output:
        writer = csv.writer(output, delimiter='\t')
        writer.writerow(columns)
        for seq, charge, decoy in [('PEPTIDEK', 2, 0), ('PEPTIDEK', 2, 1), ('PEPTIDER', 1, 0), ('PEPTIDEM', 3, 0)]:
            for ordinal in (3, 4, 5):
                writer.writerow([seq + str(charge) + ('_decoy' if decoy else ''), seq, charge, decoy,
                                 50, 1.0, 500, 300 + ordinal, 0.5, 'y', 1, ordinal, 'noloss', 'P1', 400])

    def row(seq='PEPTIDEK', charge=2, **values):
        return {'Modified.Sequence': seq, 'Precursor.Charge': float(charge), 'RT': 10.0,
                'IM': 0.8, 'Q.Value': 0.001, 'Global.Q.Value': 0.001, 'PG.Q.Value': 0.001,
                'Decoy': 0.0, 'Run': 'one', 'PEP': 0.01, 'Evidence': 1.0, **values}

    def run(name, rows, *options, error=None, omit=(), out_override=None, in_override=None, warning=None):
        report, output = root / f'{name}.parquet', out_override or root / f'{name}.tsv'
        data = {key: [r.get(key) for r in rows] for key in rows[0] if key not in omit}
        pq.write_table(pa.table(data), report, row_group_size=1)
        result = subprocess.run([binary, '-mode', 'refine', '-in', str(in_override or library), '-ids', str(report),
                                 '-out', str(output), *map(str, options)], capture_output=True, text=True, timeout=60)
        log = result.stdout + result.stderr
        if warning:
            assert warning in log, (name, log)
        if error:
            assert result.returncode != 0 and error in log, (name, result.returncode, log)
            if out_override is None:
                assert not output.exists() and not Path(str(output) + '.refine.json').exists(), name
            return
        assert result.returncode == 0, (name, result.returncode, log)
        with output.open() as source:
            actual = list(csv.DictReader(source, delimiter='\t'))
        return actual, json.loads(Path(str(output) + '.refine.json').read_text())

    actual, provenance = run('basic', [row()], '-write_im')
    assert len(actual) == 6 and {r['Decoy'] for r in actual} == {'0', '1'}
    assert all(float(r['RT']) == 10 and abs(float(r['IM']) - 0.8) < 1e-6 for r in actual)
    assert all(r['CCS'] == '' or math.isnan(float(r['CCS'])) for r in actual)
    assert provenance['library']['matched_targets'] == 1 and provenance['library']['decoys_kept'] == 1
    assert provenance['residual_before']['rt']['mean'] == 40
    for field in ('Q.Value', 'Global.Q.Value', 'PG.Q.Value', 'Decoy'):
        run('missing-' + field, [row()], error='reference has no', omit=(field,))
    run('missing-rt', [row()], error='reference has no RT', omit=('RT',))
    run('missing-im', [row()], '-write_im', error='reference has no IM', omit=('IM',))
    actual, provenance = run('filter-only', [row()], '-no_write_rt', omit=('RT', 'IM'))
    assert len(actual) == 6 and {float(r['RT']) for r in actual} == {50}
    assert provenance['library']['rt_written'] == 0
    actual, provenance = run('null-rt-im', [row(RT=None, IM=None)], '-write_im')
    assert {float(r['RT']) for r in actual} == {50} and {float(r['IM']) for r in actual} == {1}
    assert provenance['library']['rt_missing'] == 1 and provenance['library']['im_missing'] == 1
    assert provenance['units']['rt'] == 'unchanged (library prediction)'
    run('mixed-missing-rt', [row(RT=10), row('PEPTIDER', 1, RT=None)],
        error='missing observed RT would mix library predictions')
    alias = row()
    for original, alternate in [('Modified.Sequence', 'FullUniModPeptideName'), ('Precursor.Charge', 'PrecursorCharge'),
                                ('RT', 'Tr_recalibrated'), ('IM', 'IonMobility'), ('Q.Value', 'QValue')]:
        alias[alternate] = alias.pop(original)
    actual, _ = run('column-aliases', [alias], '-write_im')
    assert {float(r['RT']) for r in actual} == {10}
    actual, provenance = run('empirical', [row()], '-empirical_library', omit=('Q.Value', 'Global.Q.Value', 'PG.Q.Value', 'Decoy'))
    assert set(provenance['gates_bypassed']) == {'Q.Value', 'Global.Q.Value', 'PG.Q.Value'}
    run('empirical-compact', [row(**{'Product.Mz': [303.0, 304.0]})], '-empirical_library',
        error='compact Parquet empirical references are not supported')
    for value in (None, float('nan'), float('inf'), -0.1):
        actual, provenance = run('bad-q-' + str(value), [row(**{'Q.Value': value}), row('PEPTIDER', 1)])
        assert {r['Modified.Sequence'] for r in actual} == {'PEPTIDER'}, (value, actual)
        assert provenance['reference']['q_invalid'] == 1
    for value in (None, float('nan'), 2.5, 0, 9):
        actual, provenance = run('bad-z-' + str(value), [row(charge=2, **{'Precursor.Charge': value}), row('PEPTIDER', 1)])
        assert provenance['reference']['charge_invalid'] == 1
        assert {r['Modified.Sequence'] for r in actual} == {'PEPTIDER'}
    for value in (None, ''):
        actual, provenance = run('empty-sequence-' + str(value), [row(value), row('PEPTIDER', 1)])
        assert provenance['reference']['sequence_invalid'] == 1 and provenance['reference']['precursors'] == 1
        assert {r['Modified.Sequence'] for r in actual} == {'PEPTIDER'}
    run('multi-run', [row(), row('PEPTIDER', 1, Run='two')], error='more than one run')
    run('no-join', [row('NOPEPTIDE')], error='no reference precursor matched')
    run('low-join', [row(), row('NOPEPTIDE')], '-min_match_fraction', 0.75, error='below the required fraction')
    run('all-gated', [row(**{'Q.Value': 0.5})], error='no reference observation passed')
    run('no-evidence', [row()], '-dedup', 'highest_evidence', error='needs an Evidence', omit=('Evidence',))
    run('no-fragments', [row()], '-min_fragments', 2, error='needs fragment identities')
    def fragment(seq='PEPTIDEK', charge=2, ordinal=3, **values):
        return row(seq, charge, **{'Fragment.Type': 'y', 'Fragment.Series.Number': float(ordinal),
                                  'Fragment.Charge': 1.0, **values})
    valid_fragments = [fragment('PEPTIDER', 1), fragment('PEPTIDER', 1, 4)]
    for index, invalid in enumerate([{'Fragment.Type': 'unknown'}, {'Fragment.Type': 'bogus'},
                                    {'Fragment.Type': None}, {'Fragment.Series.Number': 3.5},
                                    {'Fragment.Series.Number': 0}, {'Fragment.Series.Number': 256},
                                    {'Fragment.Series.Number': 1e100}, {'Fragment.Charge': 1.5},
                                    {'Fragment.Charge': 0}, {'Fragment.Charge': 128}, {'Fragment.Charge': 1e100}]):
        actual, provenance = run('invalid-fragment-' + str(index),
                                 [fragment(), fragment(ordinal=4, **invalid), *valid_fragments], '-min_fragments', 2)
        assert {r['Modified.Sequence'] for r in actual} == {'PEPTIDER'}
        assert provenance['reference']['fragment_invalid'] == 1 and provenance['reference']['too_few_fragments'] == 1
    actual, provenance = run('fragment-case-duplicate',
                             [fragment(), fragment(**{'Fragment.Type': 'Y'}), *valid_fragments], '-min_fragments', 2)
    assert {r['Modified.Sequence'] for r in actual} == {'PEPTIDER'}
    assert provenance['reference']['too_few_fragments'] == 1
    for name, options, expected in [('lowest-q', [], 20), ('highest-evidence', ['-dedup', 'highest_evidence'], 30)]:
        actual, _ = run(name, [row(RT=10, **{'Q.Value': 0.005}), row(RT=20, Evidence=2), row(RT=30, Evidence=3, PEP=0.02)], *options)
        assert {float(r['RT']) for r in actual} == {expected}
    actual, _ = run('dedup-tie-first-seen', [row(RT=10), row(RT=30)])
    assert {float(r['RT']) for r in actual} == {10}
    # Separate library rows may share a canonical join key; only distinct keys count toward coverage.
    with library.open() as source:
        library_rows = list(csv.DictReader(source, delimiter='\t'))
    duplicated = root / 'duplicates.tsv'
    extra = [dict(r, **{'Precursor.Id': r['Precursor.Id'] + '_duplicate'})
             for r in library_rows if r['Modified.Sequence'] == 'PEPTIDEK' and r['Decoy'] == '0']
    with duplicated.open('w', newline='') as output:
        writer = csv.DictWriter(output, columns, delimiter='\t'); writer.writeheader(); writer.writerows(library_rows + extra)
    _, provenance = run('duplicate-library-keys', [row(), row('NOPEPTIDE')], in_override=duplicated)
    assert provenance['library']['matched_targets'] == 2 and provenance['library']['match_fraction'] == 0.5
    assert provenance['reference']['unmatched'] == 1
    unknown_library = root / 'unknown-mod.tsv'
    unknown_rows = [dict(r, **{'Modified.Sequence': r['Modified.Sequence'] + '(UnresolvableTestModification)'})
                    if r['Modified.Sequence'] == 'PEPTIDER' else r for r in library_rows]
    with unknown_library.open('w', newline='') as output:
        writer = csv.DictWriter(output, columns, delimiter='\t'); writer.writeheader(); writer.writerows(unknown_rows)
    _, provenance = run('unknown-library-modification', [row()], in_override=unknown_library,
                         warning='library modification tokens could not be resolved to UniMod')
    assert provenance['library']['unknown_mod_tokens'] == 1 and provenance['reference']['unknown_mod_tokens'] == 0
    run('duplicate-library-keys-gate', [row(), row('NOPEPTIDE')], '-min_match_fraction', 0.75,
        in_override=duplicated, error='below the required fraction')
    terminal_library = root / 'terminal.tsv'
    terminal_rows = [dict(r, **{'Modified.Sequence': '.(Acetyl)' + r['Modified.Sequence']})
                     for r in library_rows if r['Modified.Sequence'] == 'PEPTIDEK']
    with terminal_library.open('w', newline='') as output:
        writer = csv.DictWriter(output, columns, delimiter='\t'); writer.writeheader(); writer.writerows(terminal_rows)
    actual, provenance = run('n-terminal-alias', [row('(UniMod:1)PEPTIDEK')], in_override=terminal_library)
    assert len(actual) == 6 and {float(r['RT']) for r in actual} == {10}
    assert provenance['library']['match_fraction'] == 1
    actual, _ = run('im-floor', [row('PEPTIDER', 1)], '-write_im')
    assert {float(r['IM']) for r in actual} == {1.0}
    actual, provenance = run('im-ramp', [row(IM=1.59)], '-write_im', '-im_ramp_top', 1.6)
    assert {float(r['IM']) for r in actual} == {1.0} and provenance['reference']['ramp_censored'] == 1
    actual, _ = run('minmax', [row(RT=10), row('PEPTIDER', 1, RT=30)], '-rt_unit', 'minmax')
    assert {float(r['RT']) for r in actual} == {0, 100}
    run('minmax-degenerate', [row()], '-rt_unit', 'minmax', error='no finite range')
    run('minmax-mixed', [row(), row('PEPTIDER', 1, RT=30), row('PEPTIDEM', 3, RT=None)], '-rt_unit', 'minmax', error='missing observed RT')
    run('minmax-unfiltered', [row()], '-rt_unit', 'minmax', '-no_filter', error='filter off')
    run('rt-unfiltered', [row()], '-no_filter', error='mix reference-run minutes')
    ccs_models = root / 'ccs-only-models'
    run('rt-unfiltered-ccs-only', [row()], '-no_filter', '-tune', '-tune_heads', 'ccs',
        '-tune_out_models', ccs_models, error='mix reference-run minutes')
    assert not ccs_models.exists(), 'CCS-only mixed RT units must be rejected before training'
    actual, provenance = run('unfiltered-predicted-rt', [row()], '-no_filter', '-no_write_rt')
    assert len(actual) == len(library_rows) and {float(r['RT']) for r in actual} == {50}
    assert provenance['library']['rt_written'] == 0
    run('intensity-unfiltered', [row()], '-write_intensity', '-no_filter', '-no_write_rt', error='needs Fragment.Info')
    # Refusing an existing file must leave the exact bytes untouched.
    protected = root / 'protected.tsv'
    protected.write_bytes(b'preserve me\n')
    run('existing', [row()], out_override=protected, error='refusing to overwrite')
    assert protected.read_bytes() == b'preserve me\n'
    collision = root / 'collision.tsv'
    run('collision', [row()], '-out_report', collision, out_override=collision, error='output paths must be distinct')
    assert not collision.exists()
    if sys.platform != 'win32':
        link = root / 'dangling.tsv'
        link.symlink_to(root / 'absent.tsv')
        run('symlink', [row()], out_override=link, error='refusing to overwrite')
        assert link.is_symlink() and not (root / 'absent.tsv').exists()
        report_link = root / 'dangling-report.tsv'
        report_link.symlink_to(root / 'absent-report.tsv')
        run('report-symlink', [row()], '-out_report', report_link, error='refusing to overwrite')
        assert report_link.is_symlink() and not (root / 'absent-report.tsv').exists()
print('PASS: report gates, duplicate ranking, joins, RT/IM units, decoys and output safety')
