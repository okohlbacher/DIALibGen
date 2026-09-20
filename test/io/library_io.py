"""Independent Arrow fixtures exercise DIA-NN readers and persisted field semantics."""
from pathlib import Path
import csv
import subprocess
import sys
try:
    import pyarrow as pa
    import pyarrow.parquet as pq
except ImportError as e:
    print(f"SKIP: {e}", file=sys.stderr)
    sys.exit(77)

binary, work = sys.argv[1], Path(sys.argv[2])
work.mkdir(parents=True, exist_ok=True)

def run(mode, path, error=None):
    r = subprocess.run([binary, mode, str(path)], capture_output=True, text=True)
    if error is None:
        assert r.returncode == 0, (mode, path, r.stdout, r.stderr)
    else:
        assert r.returncode == 1 and error in r.stderr, (path, r.returncode, r.stderr)
    return r.stdout.strip()

run('fixtures', work)
for name in ('library.tsv', 'flat.parquet', 'compact.parquet'):
    run('roundtrip', work / name)
for name in ('LIBRARY.TSV', 'LIBRARY.PARQUET'):
    run('write', work / name)
    run('roundtrip', work / name)
run('load', work / 'library.tsv.gz', 'unsupported library extension')

flat = pq.read_table(work / 'flat.parquet')
compact = pq.read_table(work / 'compact.parquet')
assert flat.num_rows == 5 and compact.num_rows == 3
assert flat['Decoy'].to_pylist() == [0, 0, 1, 1, 0]
assert flat['Fragment.Loss.Type'].to_pylist() == ['H2O', 'noloss', 'NH3', 'noloss', 'noloss']
assert flat['RT'].to_pylist()[-1] is None and compact['RT'].to_pylist()[-1] is None
assert flat['Precursor.Mz'].to_pylist()[0] == 600.12345
assert compact.schema.metadata[b'odia.config_json'] == b'{"example":true}'
assert compact.schema.metadata[b'odia.fingerprint'] == b'abcdef:42:v4:recipe'
assert compact.schema.metadata[b'odia.fasta_fnv1a64'] == b'abcdef'
for name, table in (('flat', flat), ('compact', compact)):
    path = work / f'empty-{name}.parquet'
    pq.write_table(table.slice(0, 0), path)
    assert run('load', path) == '0 0 0 0 0'

# CRLF and numeric storage variations must preserve the same library.
tsv = (work / 'library.tsv').read_text()
(work / 'header-only.tsv').write_text(tsv.splitlines()[0] + '\n')
assert run('load', work / 'header-only.tsv') == '0 0 0 0 0'
(work / 'crlf.tsv').write_bytes(tsv.replace('\n', '\r\n').encode())
run('roundtrip', work / 'crlf.tsv')
(work / 'bom-blank-lines.tsv').write_bytes(b'\xef\xbb\xbf' + ('\r\n\r\n'.join(tsv.splitlines()) + '\r\n\r\n').encode())
run('roundtrip', work / 'bom-blank-lines.tsv')
for name, line in [('short', tsv.splitlines()[1].rsplit('\t', 1)[0]),
                   ('wide', tsv.splitlines()[1] + '\textra')]:
    path = work / f'{name}-row.tsv'
    path.write_text(tsv.splitlines()[0] + '\n' + line + '\n')
    run('load', path, 'row width')
path = work / 'duplicate-header.tsv'
path.write_text(tsv.splitlines()[0] + '\tPrecursor.Id\n')
run('load', path, 'duplicate TSV column')
for token in ('1.5garbage', 'inf', '0x1p2', '9' * 500):
    rows = list(csv.DictReader(tsv.splitlines(), delimiter='\t'))
    rows[0]['RT'] = token
    path = work / 'invalid-numeric.tsv'
    with path.open('w', newline='') as f:
        writer = csv.DictWriter(f, rows[0].keys(), delimiter='\t'); writer.writeheader(); writer.writerows(rows)
    run('load', path, 'invalid numeric value')
for field in ('Precursor.Id', 'Modified.Sequence', 'Protein.Group'):
    i = flat.schema.get_field_index(field)
    flat = flat.set_column(i, field, flat[field].cast(pa.large_string()))
flat = flat.set_column(flat.schema.get_field_index('Decoy'), 'Decoy', flat['Decoy'].cast(pa.bool_()))
pq.write_table(flat, work / 'typed.parquet', row_group_size=3)
run('roundtrip', work / 'typed.parquet')
for name, table in (('flat', flat), ('compact', compact)):
    changed = table.set_column(table.schema.get_field_index('Modified.Sequence'),
                               'Modified.Sequence', pa.array([[1]] * table.num_rows))
    path = work / f'bad-string-{name}.parquet'
    pq.write_table(changed, path)
    run('load', path, 'cannot decode Parquet string column Modified.Sequence')
# More than one decoder batch within a single precursor retains its grouping.
long = pa.concat_tables([flat.slice(0, 1)] * 65537)
pq.write_table(long, work / 'batch-boundary.parquet', row_group_size=40000)
assert run('load', work / 'batch-boundary.parquet') == '1 65537 0 0 0'

# A historical compact-v1 file may omit its all-noloss annotation.
legacy = compact.drop(['Fragment.Loss.Type'])
pq.write_table(legacy, work / 'legacy.parquet')
assert run('load', work / 'legacy.parquet') == '3 5 1 0 0'
pq.write_table(legacy.drop(['Decoy', 'Protein.Group']), work / 'optional-columns.parquet')
assert run('load', work / 'optional-columns.parquet') == '3 5 0 0 0'

(work / 'empty.tsv').write_text('')
run('load', work / 'empty.tsv', 'empty library')
(work / 'missing.tsv').write_text('Wrong\tColumns\n1\t2\n')
run('load', work / 'missing.tsv', 'missing Precursor.Id')
(work / 'corrupt.parquet').write_bytes(b'not parquet')
run('load', work / 'corrupt.parquet', 'not a Parquet file')
run('load', work / 'nonexistent.tsv', 'cannot open library')
lines = tsv.splitlines()
(work / 'ungrouped.tsv').write_text('\n'.join([lines[0], lines[1], lines[3], lines[2]]) + '\n')
run('load', work / 'ungrouped.tsv', 'not grouped by precursor')

# Malformed nested columns must report an error, never index another list's memory.
def invalid(name, column, values, dtype, message):
    changed = compact.set_column(compact.schema.get_field_index(column), column, pa.array(values, type=dtype))
    path = work / f'{name}.parquet'
    pq.write_table(changed, path)
    run('load', path, message)

invalid('wrong-child', 'Product.Mz', [['300', '400'], ['310', '410'], ['220']], pa.list_(pa.string()), 'Product.Mz')
invalid('short-list', 'Relative.Intensity', [[1.], [1., .3], [.25]], pa.list_(pa.float32()), 'list lengths')
invalid('null-list', 'Product.Mz', [None, [310., 410.], [220.]], pa.list_(pa.float64()), 'null')
invalid('null-child', 'Product.Mz', [[None, 400.], [310., 410.], [220.]], pa.list_(pa.float64()), 'null')
invalid('fractional-charge', 'Precursor.Charge', [2.5, 2, 3], pa.float64(), 'charge')
invalid('null-charge', 'Precursor.Charge', [None, 2, 3], pa.uint8(), 'charge')
invalid('null-sequence', 'Modified.Sequence', [None, 'AC(UniMod:4)DEFGK', 'PEPTIDER'], pa.string(), 'sequence')
invalid('wrong-enum', 'Fragment.Loss.Type', [[255, 0], [0, 0], [0]], pa.list_(pa.uint8()), 'loss')
invalid('wrong-type', 'Fragment.Type', [[255, 0], [0, 0], [0]], pa.list_(pa.uint8()), 'type')
invalid('zero-fragment-charge', 'Fragment.Charge', [[0, 1], [1, 1], [1]], pa.list_(pa.uint8()), 'fragment charge')
invalid('overflow-fragment-charge', 'Fragment.Charge', [[128, 1], [1, 1], [1]], pa.list_(pa.uint8()), 'fragment charge')
invalid('zero-ordinal', 'Fragment.Series.Number', [[0, 1], [1, 1], [1]], pa.list_(pa.uint8()), 'ordinal')
invalid('long-ordinal', 'Fragment.Series.Number', [[255, 1], [1, 1], [1]], pa.list_(pa.uint8()), 'ordinal')
pq.write_table(compact.drop(['Relative.Intensity']), work / 'missing-intensity.parquet')
run('load', work / 'missing-intensity.parquet', 'required columns')

# Every representation must reject a charge that would truncate or wrap.
for value in (0, -1, 2.5, 258, float('nan')):
    for kind, table in (('flat', flat), ('compact', compact)):
        path = work / f'{kind}-bad-charge-{value}.parquet'
        changed = table.set_column(table.schema.get_field_index('Precursor.Charge'),
                                   'Precursor.Charge', pa.array([value] * table.num_rows))
        pq.write_table(changed, path)
        run('load', path, 'charge')
    rows = list(csv.DictReader(tsv.splitlines(), delimiter='\t'))
    rows[0]['Precursor.Charge'] = str(value)
    path = work / f'tsv-bad-charge-{value}.tsv'
    with path.open('w', newline='') as f:
        writer = csv.DictWriter(f, rows[0].keys(), delimiter='\t'); writer.writeheader(); writer.writerows(rows)
    run('load', path, 'charge')
for kind, table in (('flat', flat), ('compact', compact)):
    path = work / f'{kind}-charge-ten.parquet'
    pq.write_table(table.set_column(table.schema.get_field_index('Precursor.Charge'),
                                   'Precursor.Charge', pa.array([10] * table.num_rows)), path)
    assert run('load', path) == '3 5 1 0 0'

# Invalid m/z is represented by the documented sentinel, not a wrapped integer.
rows = list(csv.DictReader(tsv.splitlines(), delimiter='\t'))
rows[0]['Precursor.Mz'] = '-1'
rows[0]['Product.Mz'] = '1e100'
rows[1]['Product.Mz'] = 'nan'
path = work / 'invalid-mz.tsv'
with path.open('w', newline='') as f:
    writer = csv.DictWriter(f, rows[0].keys(), delimiter='\t'); writer.writeheader(); writer.writerows(rows)
assert run('load', path) == '3 5 1 1 2'

missing = work / 'missing-values.tsv'
run('missing-values', missing)
rows = list(csv.DictReader(missing.read_text().splitlines(), delimiter='\t'))
assert rows[0]['Relative.Intensity'] == '' and all(r['IM'] == '' for r in rows)
existing = work / 'preserved.tsv'
existing.write_text('original content\n')
run('invalid-text', existing, 'tab, newline or NUL')
assert existing.read_text() == 'original content\n'

# A real OS write failure after opening the staged file must preserve every format.
if sys.platform != 'win32':
    import resource, signal
    def fail_writes():
        signal.signal(signal.SIGXFSZ, signal.SIG_IGN)
        resource.setrlimit(resource.RLIMIT_FSIZE, (1, 1))
    for name in ('quota.tsv', 'quota.parquet', 'compact-quota.parquet'):
        path = work / name; path.write_text('existing library\n')
        result = subprocess.run([binary, 'write', str(path)], capture_output=True, text=True, preexec_fn=fail_writes)
        assert result.returncode == 1, (name, result.stdout, result.stderr)
        assert path.read_text() == 'existing library\n', name
assert not list(work.glob('.dialibgen-tmp-*')), 'temporary output leaked after failure'
print('TSV/flat/compact round trips, independent Arrow values, metadata, malformed inputs and missing data passed')
