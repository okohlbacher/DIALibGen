#!/usr/bin/env python3
"""Training refusal, cohort isolation and stopping contracts using real stock models."""
import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import zlib

try:
    import numpy  # synth_report.py dependency
    import pyarrow as pa
    import pyarrow.parquet as pq
except ImportError:
    print('SKIP: Python with numpy and pyarrow is required')
    sys.exit(77)

binary, models = str(Path(sys.argv[1]).resolve()), Path(sys.argv[2]).resolve()
here = Path(__file__).resolve().parent

with tempfile.TemporaryDirectory(prefix='dialibgen-training-contract-') as directory:
    root = Path(directory)
    source, library = root / 'source.parquet', root / 'library.tsv'
    fixture = subprocess.run([sys.executable, str(here / 'synth_report.py'), '-', str(source),
                              '--precursors', '1600', '--library', str(library)],
                             capture_output=True, text=True, timeout=120)
    assert fixture.returncode == 0, fixture.stdout + fixture.stderr
    table = pq.read_table(source)

    def run(name, report, *options, head='rt', error=None, training=False):
        path, output, kept = root / f'{name}.parquet', root / f'{name}.tsv', root / f'{name}-models'
        pq.write_table(report, path)
        controls = {'-filter:rt_max_minutes': 30, '-machine:threads': 2, '-train:epochs': 1,
                    '-train:warmup': 0, '-stop:min_epochs': 0, '-train:lr': '1e-300'}
        controls.update(zip(options[::2], options[1::2]))
        command = [binary, '-mode', 'tune', '-in', str(library), '-ids', str(path), '-out', str(output),
                   '-tune_models', str(models), '-tune_out_models', str(kept), '-tune_heads', head,
                   *(str(value) for pair in controls.items() for value in pair)]
        result = subprocess.run(command, capture_output=True, text=True, timeout=180)
        log = result.stdout + result.stderr
        (root / f'{name}.log').write_text(log)
        assert result.returncode != 0 and error in log, (name, result.returncode, log)
        assert not output.exists() and not Path(str(output) + '.refine.json').exists(), name
        if not training:
            assert not list(kept.glob('*')), (name, list(kept.glob('*')), log)
            return
        model = kept / f'peptdeep_{head}_dynamic.onnx'
        assert not model.exists(), (name, 'non-improving model was exported')
        provenance = json.loads(Path(str(model) + '.tune.json').read_text())
        trajectory = list(csv.DictReader(Path(str(model) + '.trajectory.tsv').open(), delimiter='\t'))
        course = provenance['course']
        assert course['updates'] > 0 and not course['exported'], (name, course)
        assert course['best_epoch'] == 0 and course['param_l2_change'] == 0, (name, course)
        assert 'output' not in provenance, (name, 'failure provenance claims an exported model')
        assert len(trajectory) == course['epochs_run'], (name, trajectory, course)
        return provenance, trajectory, log

    for column in ('Modified.Sequence', 'Precursor.Charge', 'Precursor.Mz', 'RT', 'Q.Value', 'Protein.Group'):
        run('missing-' + column, table.drop([column]), error='report lacks column ' + column)
    changed = table.set_column(table.schema.get_field_index('Run'), 'Run',
                               pa.array(['second'] + ['synthetic'] * (table.num_rows - 1)))
    run('multiple-runs', changed, error='give a single-run report')
    for column, value in (('Decoy', 1), ('Q.Value', 0.5), ('Precursor.Charge', 2.5), ('RT', float('nan'))):
        changed = table.set_column(table.schema.get_field_index(column), column, pa.array([value] * table.num_rows))
        run('no-usable-' + column, changed, error='no usable observations after filtering')
    run('small-cohorts', table.slice(0, 20), error='fewer than 100 units')

    # A positive learning rate below float32's range performs real optimizer
    # steps without changing weights, making rejection independent of convergence.
    failure = 'no checkpoint beat the stock model'
    provenance, trajectory, log = run('no-improvement', table.drop(['Decoy']), error=failure, training=True)
    assert 'no Decoy column; all report rows are treated as targets' in log
    assert provenance['course']['epochs_run'] == 1 and provenance['course']['stop_reason'] == 'horizon'
    assert trajectory[0]['val_calibrated_sd'] and trajectory[0]['stop'] == 'horizon'

    # Independent protein-group hashes and unit keys verify that conflicts across
    # CCS charge states remove both units before either can enter a cohort.
    rows = table.to_pylist()
    by_sequence = {}
    for i, row in enumerate(rows):
        by_sequence.setdefault(row['Modified.Sequence'], []).append(i)
    conflict = next(indices for indices in by_sequence.values()
                    if len({rows[i]['Precursor.Charge'] for i in indices}) > 1)
    old_group = rows[conflict[0]]['Protein.Group']
    replacement = next(f'CONFLICT{i}' for i in range(10000)
                       if (zlib.crc32(f'CONFLICT{i}'.encode()) % 5 == 0)
                       != (zlib.crc32(old_group.encode()) % 5 == 0))
    rows[conflict[-1]]['Protein.Group'] = replacement
    changed = pa.Table.from_pylist(rows, schema=table.schema)
    provenance, _, _ = run('cross-charge-conflict', changed, head='ccs', error=failure, training=True)
    groups, units = {}, {}
    for row in rows:
        seq, pg = row['Modified.Sequence'], row['Protein.Group']
        groups.setdefault(seq, set()).add(pg)
        units[(seq, row['Precursor.Charge'])] = pg
    expected = {'test': 0, 'val': 0, 'pool': 0}
    rejected = 0
    for (seq, _), pg in units.items():
        if len(groups[seq]) != 1:
            rejected += 1
            continue
        cohort = ('test' if zlib.crc32(pg.encode()) % 5 == 0 else
                  'val' if zlib.crc32(('val:' + pg).encode()) % 7 == 0 else 'pool')
        expected[cohort] += 1
    actual = provenance['cohorts']
    assert {key: actual[key] for key in expected} == expected, (actual, expected)
    assert actual['units'] == sum(expected.values()) and actual['training'] == expected['pool'], actual
    assert actual['rejected']['protein_group_conflict'] == rejected >= 2, actual
    assert not actual['val_is_test'] and not actual['full_fit'], actual

    # Budget expiry forces evaluation even between scheduled evaluations and
    # takes precedence over the minimum epoch floor.
    provenance, trajectory, _ = run('time-budget', table, '-train:epochs', 10,
        '-stop:eval_every', 10, '-stop:min_epochs', 10, '-stop:max_seconds', '1e-9',
        error=failure, training=True)
    assert provenance['course']['epochs_run'] == 1 and provenance['course']['stop_reason'] == 'max_seconds'
    assert trajectory[0]['val_calibrated_sd'] and trajectory[0]['stop'] == 'max_seconds'
    assert 'wall time' in provenance['stopping']['budget_scope']

    # No progress: anchor is epoch2, patience is satisfied at epoch4, but
    # warmup/minimum prevent stopping until the next evaluation at epoch6.
    provenance, trajectory, _ = run('patience', table, '-train:epochs', 10, '-train:warmup', 6,
        '-stop:eval_every', 2, '-stop:min_epochs', 5, '-stop:patience', 2,
        '-stop:abs_tol', 1, error=failure, training=True)
    assert provenance['course']['epochs_run'] == 6, provenance['course']
    assert provenance['course']['stop_reason'].startswith('patience (4 epochs'), provenance['course']
    assert [int(row['epoch']) for row in trajectory if row['val_calibrated_sd']] == [2, 4, 6], trajectory
    assert all(not row['stop'] for row in trajectory[:-1])
    print('PASS: training preflight, failure provenance, cross-charge cohorts, wall-time boundary and patience contracts')
