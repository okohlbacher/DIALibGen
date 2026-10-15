#!/usr/bin/env python3
"""The 0.12 command line with prediction, on the example FASTA:

- the report's 13-option Astral command and the 2-option one write the same
  library (TSV byte for byte; Parquet: the same table and fingerprint, the
  recipe differing only in its recorded sources);
- -threads: all physical cores unless the command line or the INI sets it;
- a library as -in converts without -mode, and the Parquet keeps its identity;
- an OpenSWATH output without decoys gets a hint.

Usage: cli_library.py <DIALibGen> <proteins.fasta> <model dir>
Exits 77 (ctest SKIP) without pyarrow.
"""
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET

try:
    import pyarrow.parquet as pq
except ImportError as error:
    print(f'SKIP: {error}', file=sys.stderr)
    sys.exit(77)

binary, fasta, models = (str(Path(a).resolve()) for a in sys.argv[1:4])
env = dict(os.environ, DIALIBGEN_MODEL_DIR=models, ODIA_VERBOSE='0')
env.pop('DIALIBGEN_INSTRUMENT_DIR', None)

with tempfile.TemporaryDirectory(prefix='dialibgen-cli-library-') as directory:
    root = Path(directory)

    def run(*args, ok=True):
        r = subprocess.run([binary, *map(str, args)], capture_output=True, text=True, timeout=600, env=env)
        log = r.stdout + r.stderr
        assert (r.returncode == 0) == ok, (args, r.returncode, log)
        return log

    old = ['-generation:instrument', 'Astral', '-generation:nce', 25, '-generation:precursor_mz', 380, 980,
           '-generation:fragment_mz', 150, 2000, '-generation:precursor_charges', 1, 2, 3, 4,
           '-generation:peptide_length', 7, 30, '-generation:fixed_modifications', 'Carbamidomethyl (C)',
           '-generation:variable_modifications', 'Oxidation (M)', '-generation:max_variable_modifications', 1,
           '-generation:decoys', 'none', '-generation:enzyme', 'Trypsin/P', '-generation:missed_cleavages', 1,
           '-generation:n_terminal_methionine_excision', 'true', '-threads', 4]
    new = ['-generation:instrument', 'Astral', '-generation:variable_modifications', 'Oxidation (M)']
    for ext in ('tsv', 'parquet'):
        run('-in', fasta, '-out', root / f'old.{ext}', '-machine:device', 'cpu', *old)
        log = run('-in', fasta, '-out', root / f'new.{ext}', '-machine:device', 'cpu', *new)
    assert (root / 'old.tsv').read_bytes() == (root / 'new.tsv').read_bytes(), 'Astral TSVs differ'
    a, b = pq.read_table(root / 'old.parquet'), pq.read_table(root / 'new.parquet')
    assert a.equals(b), 'Astral Parquet tables differ'
    ma, mb = a.schema.metadata, b.schema.metadata
    assert set(ma) == set(mb)
    for key in ma:
        if key != b'odia.config_json':
            assert ma[key] == mb[key], key
    ca, cb = json.loads(ma[b'odia.config_json']), json.loads(mb[b'odia.config_json'])
    differ = {k for k in set(ca) | set(cb) if ca.get(k) != cb.get(k)}
    assert differ <= {'nce_source', 'precursor_mz_source', 'fragment_mz_source', 'precursor_charges_source'}, differ
    assert cb['instrument_defaults']['instrument'] == 'astral.yaml'

    # -threads: all physical cores unless the command line or the INI says otherwise.
    log = run('-in', fasta, '-out', root / 'threads-default.tsv', '-machine:device', 'cpu')
    assert re.search(r'-threads 0: \d+ threads', log), log
    log = run('-in', fasta, '-out', root / 'threads-1.tsv', '-machine:device', 'cpu', '-threads', 1)
    assert '-threads 0:' not in log and re.search(r'\[MS2\] provider CPU device -1 sessions 1\b', log), log
    ini = root / 'settings.ini'
    run('-write_ini', ini)
    tree = ET.parse(ini)
    instance = next(n for n in tree.getroot().iter('NODE') if n.get('name') == '1')
    threads = next(n for n in instance if n.get('name') == 'threads')
    assert threads.get('value') == '0'
    for value in (1, 2):
        threads.set('value', str(value))
        tree.write(root / f'threads-{value}.ini', encoding='ISO-8859-1', xml_declaration=True)
        log = run('-ini', root / f'threads-{value}.ini', '-in', fasta, '-out', root / f'ini-{value}.tsv',
                  '-machine:device', 'cpu')
        assert '-threads 0:' not in log and re.search(rf'\[MS2\] provider CPU device -1 sessions {value}\b', log), log
    assert len({(root / n).read_bytes() for n in ('threads-default.tsv', 'threads-1.tsv', 'ini-1.tsv', 'ini-2.tsv')}) == 1

    # No -mode: a library converts, and a DIALibGen Parquet keeps its identity.
    log = run('-in', root / 'new.parquet', '-out', root / 'converted.parquet')
    assert 'no -mode: -in is a library, so -mode convert' in log, log
    c = pq.read_table(root / 'converted.parquet')
    assert c.equals(b)
    for key in (b'odia.fingerprint', b'odia.target_fingerprint', b'odia.params', b'odia.config_json'):
        assert c.schema.metadata[key] == mb[key], key
    run('-in', root / 'new.parquet', '-out', root / 'converted.tsv')
    assert (root / 'converted.tsv').read_bytes() == (root / 'new.tsv').read_bytes()

    # OpenSWATH scores against decoys: a library without them gets a hint.
    log = run('-in', fasta, '-out', root / 'lib.pqp', '-machine:device', 'cpu', '-threads', 2)
    assert 'OpenSWATH scores against the library\'s decoys and this one has none' in log, log
print('Astral 2-option command, -threads default, inferred convert, OpenSWATH hint: PASS')
