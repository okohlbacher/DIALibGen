#!/usr/bin/env python3
"""The -run command-line contract: -ids XOR -run, refusals, -out_ids never
overwritten, search: options gated to refine/tune WITH -run, INI round trips.

Needs no pyarrow and no run data: every check here is decided before the run is
read, or asserts only that the search STARTED (it names the run it reads before
reading it) and that a failed search leaves no output behind."""
import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET

binary = str(Path(sys.argv[1]).resolve())
COLUMNS = ['Precursor.Id', 'Modified.Sequence', 'Precursor.Charge', 'Decoy', 'RT', 'IM', 'Precursor.Mz',
           'Product.Mz', 'Relative.Intensity', 'Fragment.Type', 'Fragment.Charge', 'Fragment.Series.Number',
           'Fragment.Loss.Type', 'Protein.Group', 'CCS']
PEPTIDES = ['LGADEFHIVK', 'NPQSTVWYAR', 'GFHILMNDEK', 'SAVLTEDFGR', 'YWPNQHSDLK', 'ELVFGADQTR',
            'MNHIGFDSAK', 'TVSQPNYWER', 'DAEFGHLKIVR', 'QPSTNVWHYK', 'FGEDALIVTSR', 'HIMNDGSAFEK']

with tempfile.TemporaryDirectory(prefix='dialibgen-identify-cli-') as directory:
    root = Path(directory)
    library = root / 'library.tsv'
    with library.open('w', newline='') as output:
        writer = csv.writer(output, delimiter='\t')
        writer.writerow(COLUMNS)
        for n, seq in enumerate(PEPTIDES):
            for charge in (2, 3):
                for ordinal in (3, 4, 5, 6):
                    writer.writerow([seq + str(charge), seq, charge, 0, 10 + 5 * n, 0.9, 400 + 10 * n,
                                     300 + 50 * ordinal, 1.0 / ordinal, 'y', 1, ordinal, 'noloss', f'P{n // 2}', 400])
    ids = root / 'report.parquet'
    ids.write_bytes(b'not read: every check below fails first\n')
    mzml = root / 'run.mzML'
    mzml.write_text('<?xml version="1.0"?>\n')   # not a real run; see the -run cases below
    counter = iter(range(1000))

    def run(*args, ok=True):
        r = subprocess.run([binary, *map(str, args)], capture_output=True, text=True, timeout=120)
        log = r.stdout + r.stderr
        assert (r.returncode == 0) == ok, (args, r.returncode, log)
        return log

    def refused(expected, *args):
        log = run(*args, ok=False)
        assert expected in log, (expected, args, log)
        return log

    def fresh(suffix='.tsv'):
        return root / f'out-{next(counter)}{suffix}'

    # -ids XOR -run
    refused('exactly one of -ids and -run', '-mode', 'refine', '-in', library, '-ids', ids, '-run', mzml, '-out', fresh())
    refused('-in, -ids and -out are required', '-mode', 'refine', '-in', library, '-out', fresh())
    refused('exactly one of -ids and -run', '-mode', 'tune', '-in', library, '-ids', ids, '-run', mzml, '-out', fresh())
    run('-mode', 'refine', '-in', library, '-run', root / 'absent.mzML', '-out', fresh(), ok=False)

    # What -run refuses until later milestones.
    refused('-write_intensity is not available with -run', '-mode', 'refine', '-in', library, '-run', mzml,
            '-out', fresh(), '-write_intensity')
    refused('-empirical_library does not apply to -run', '-mode', 'refine', '-in', library, '-run', mzml,
            '-out', fresh(), '-empirical_library')
    refused('-min_fragments is not available with -run', '-mode', 'refine', '-in', library, '-run', mzml,
            '-out', fresh(), '-min_fragments', 3)
    refused("'mutate'", '-mode', 'refine', '-in', library, '-run', mzml, '-out', fresh(), '-search:decoys', 'mutate')
    run('-mode', 'refine', '-in', library, '-run', mzml, '-out', fresh(), '-search:report_max_q', 0.001, ok=False)
    refused('-out_ids must end in .parquet', '-mode', 'refine', '-in', library, '-run', mzml, '-out', fresh(),
            '-out_ids', root / 'ids.tsv')

    # -out_ids is never overwritten: explicit or default, and distinct from -out.
    sentinel = b'existing report\n'
    for explicit in (True, False):
        out = fresh()
        target = root / f'existing-{next(counter)}.parquet' if explicit else Path(str(out) + '.ids.parquet')
        target.write_bytes(sentinel)
        args = ['-mode', 'tune', '-in', library, '-run', mzml, '-out', out] + (['-out_ids', target] if explicit else [])
        refused('refusing to overwrite existing output', *args)
        assert target.read_bytes() == sentinel and not out.exists()
    out = fresh('.parquet')
    refused('output paths must be distinct', '-mode', 'refine', '-in', library, '-run', mzml, '-out', out, '-out_ids', out)

    # search: belongs to refine and tune, and only with -run; so does -out_ids.
    refused('has no effect in -mode generate', '-mode', 'generate', '-search:subset', 5, '-write_config', fresh('.json'))
    refused('has no effect in -mode generate', '-mode', 'generate', '-run', mzml, '-write_config', fresh('.json'))
    refused('has no effect in -mode generate', '-mode', 'generate', '-out_ids', root / 'x.parquet', '-write_config', fresh('.json'))
    for mode in ('refine', 'tune'):
        refused('-search:seed has no effect without -run', '-mode', mode, '-search:seed', 7, '-write_config', fresh('.json'))
        refused('-search:selftest has no effect without -run', '-mode', mode, '-search:selftest', '-write_config', fresh('.json'))
        refused('-out_ids has no effect without -run', '-mode', mode, '-out_ids', root / 'x.parquet', '-write_config', fresh('.json'))
        run('-mode', mode, '-run', mzml, '-search:seed', 7, '-search:selftest', '-out_ids', root / 'x.parquet',
            '-write_config', fresh('.json'))
    run('-mode', 'refine', '-tune', '-run', mzml, '-search:max_pairs', 10, '-write_config', fresh('.json'))

    # -write_config keeps its schema: the search settings are not part of it.
    plain, with_run = fresh('.json'), fresh('.json')
    run('-mode', 'refine', '-write_config', plain)
    run('-mode', 'refine', '-run', mzml, '-search:subset', 5, '-write_config', with_run)
    assert json.loads(plain.read_text()) == json.loads(with_run.read_text())

    # A full default INI (every mode's defaults, search: included) stays usable
    # in every mode; a changed search: value in it needs -run.
    ini = root / 'all.ini'
    run('-write_ini', ini)
    names = {el.get('name') for el in ET.parse(ini).iter('ITEM')}
    assert {'run', 'out_ids', 'subset', 'max_pairs', 'decoys', 'min_ids', 'report_max_q', 'selftest'} <= names, names
    for mode in ('generate', 'refine', 'tune'):
        run('-mode', mode, '-ini', ini, '-write_config', fresh('.json'))
    tree = ET.parse(ini)
    search = next(n for n in tree.iter('NODE') if n.get('name') == 'search')
    next(n for n in search if n.get('name') == 'subset').set('value', '5')
    changed = root / 'search.ini'
    tree.write(changed, encoding='utf-8', xml_declaration=True)
    refused('-search:subset has no effect without -run', '-mode', 'tune', '-ini', changed, '-write_config', fresh('.json'))
    run('-mode', 'tune', '-ini', changed, '-run', mzml, '-write_config', fresh('.json'))

    # Accepted options reach the search, which says what it reads before it
    # reads it (the run comes first: its isolation windows decide which
    # candidates can be searched). This "run" is not one, so the search fails
    # -- and leaves neither -out_ids nor -out behind.
    for mode, extra in (('tune', []), ('refine', ['-tune', '-no_filter']), ('refine', [])):
        out = fresh()
        report = root / f'report-{next(counter)}.parquet'
        log = run('-mode', mode, *extra, '-in', library, '-run', mzml, '-out', out, '-out_ids', report,
                  '-search:subset', 20, '-search:seed', 3, '-search:min_ids', 1, ok=False)
        if 'no fine-tuning stage' in log:
            continue   # a build without libtorch refuses tuning before anything else
        assert 'has no effect' not in log and 'exactly one' not in log, log
        assert 'search run: reading ' + str(mzml) in log, log
        assert 'EXPERIMENTAL' in log, log
        assert not report.exists() and not out.exists() and not Path(str(out) + '.refine.json').exists()
        leftovers = [p.name for p in root.iterdir() if p.name.startswith('.dialibgen-tmp-')]
        assert not leftovers, leftovers

print('identify CLI contract: PASS')
