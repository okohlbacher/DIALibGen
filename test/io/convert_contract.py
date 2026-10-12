"""-mode convert rewrites a library in another format without changing its content."""
from pathlib import Path
import subprocess
import sys

library_io, dialibgen, work = sys.argv[1], sys.argv[2], Path(sys.argv[3])
# -mode convert refuses to overwrite: start from an empty directory, so a rerun passes too.
for stale in ('converted.tsv', 'converted.parquet'):
    (work / stale).unlink(missing_ok=True)
work.mkdir(parents=True, exist_ok=True)

def run(*args, ok=True):
    r = subprocess.run([str(a) for a in args], capture_output=True, text=True)
    assert (r.returncode == 0) == ok, (args, r.returncode, r.stdout, r.stderr)
    return r.stdout + r.stderr

run(library_io, 'fixtures', work)
# Parquet -> TSV -> Parquet; each output must still be the fixture library.
log = run(dialibgen, '-mode', 'convert', '-in', work / 'compact.parquet', '-out', work / 'converted.tsv')
assert 'converted 3 precursors (1 decoys) and 5 transitions' in log, log
run(library_io, 'roundtrip', work / 'converted.tsv')
run(dialibgen, '-mode', 'convert', '-in', work / 'converted.tsv', '-out', work / 'converted.parquet')
run(library_io, 'roundtrip', work / 'converted.parquet')
assert (work / 'converted.tsv').read_text() == (work / 'library.tsv').read_text(), 'convert wrote a different TSV than storeTSV'
# An existing output is never replaced, and options of other modes are refused.
sentinel = work / 'existing.tsv'
sentinel.write_text('keep\n')
assert 'refusing to overwrite' in run(dialibgen, '-mode', 'convert', '-in', work / 'compact.parquet', '-out', sentinel, ok=False)
assert sentinel.read_text() == 'keep\n'
assert 'has no effect in -mode convert' in run(dialibgen, '-mode', 'convert', '-in', work / 'compact.parquet',
                                              '-out', work / 'never.tsv', '-generation:nce', 30, ok=False)
assert 'has no effect in -mode convert' in run(dialibgen, '-mode', 'convert', '-in', work / 'compact.parquet',
                                              '-out', work / 'never.tsv', '-q_global', 0.5, ok=False)
assert not (work / 'never.tsv').exists()
print('convert: Parquet -> TSV -> Parquet keeps the library; existing outputs and foreign options refused')
