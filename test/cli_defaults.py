#!/usr/bin/env python3
"""The command line without model inference: help, unknown options, mode
inference, the names removed in 0.13 (command line and INI), INIs from 0.11,
the defaults, the instrument files and the settings without an effect.

Usage: cli_defaults.py <DIALibGen> <repo data dir>
"""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET

binary = str(Path(sys.argv[1]).resolve())
data_dir = Path(sys.argv[2]).resolve()

with tempfile.TemporaryDirectory(prefix='dialibgen-cli-') as directory:
    root = Path(directory)
    home = root / 'home'
    home.mkdir()
    # Nothing from the user running the test: no instrument directory, no config dir.
    base_env = {k: v for k, v in os.environ.items() if k not in ('DIALIBGEN_INSTRUMENT_DIR', 'XDG_CONFIG_HOME', 'APPDATA')}
    base_env.update(HOME=str(home), APPDATA=str(home / 'AppData'), DIALIBGEN_DATA_DIR=str(data_dir))
    counter = iter(range(10000))

    def run(*args, ok=True, env=None):
        r = subprocess.run([binary, *map(str, args)], capture_output=True, text=True, timeout=120,
                           env=dict(base_env, **(env or {})))
        log = r.stdout + r.stderr
        assert (r.returncode == 0) == ok, (args, r.returncode, log)
        return r, log

    def fresh(suffix='.json'):
        return root / f'out-{next(counter)}{suffix}'

    def config(*args, env=None):
        out = fresh()
        _, log = run(*args, '-write_config', out, env=env)
        return json.loads(out.read_text()), log

    # ---- help ---------------------------------------------------------------
    def options(text):
        return [line.split()[0] for line in text.splitlines() if line.startswith('  -')]
    _, help_text = run('--help')
    tool = [o for o in options(help_text) if o not in ('-ini', '-threads', '-write_ini', '--help', '--helphelp')]
    visible = {'-mode', '-in', '-out', '-config', '-generation:instrument', '-generation:nce', '-generation:precursor_mz',
               '-generation:fragment_mz', '-generation:precursor_charges', '-generation:fixed_modifications',
               '-generation:variable_modifications', '-generation:enzyme', '-generation:missed_cleavages',
               '-ids', '-no_filter', '-no_write_rt', '-write_im', '-write_intensity', '-q_precursor', '-tune',
               '-tune_heads', '-tune_out_models', '-machine:device'}
    assert set(tool) == visible and len(tool) == len(visible), sorted(set(tool) ^ visible)
    assert len(options(help_text)) <= 30
    for line in help_text.splitlines():
        if line.startswith('  -'):
            key = line.split()[0].lstrip('-').split(':')[-1]
            text = line.split(None, 2)[-1] if len(line.split()) > 2 else ''
            assert not text.lower().startswith(key.replace('_', ' ').lower() + ' ('), line
    _, all_text = run('--helphelp')
    assert len(options(all_text)) > 115 and 'Deprecated' not in all_text
    assert 'share/DIALibGen/instruments' in all_text
    per_mode = {'generate': 14, 'refine': 14, 'tune': 8, 'convert': 4, 'prepare': 2}
    for mode, count in per_mode.items():
        _, text = run('--help', '-mode', mode)
        shown = options(text)
        assert len(shown) == count + 5, (mode, shown)
        assert len(text.splitlines()) <= 45, (mode, len(text.splitlines()))

    # ---- an unknown option: one line and a pointer, not the whole help -----
    r, log = run('-instrument', 'Astral', '-in', 'x.fasta', ok=False)
    assert r.returncode == 6, r.returncode
    assert "Unknown option '-instrument'" in log and 'did you mean -generation:instrument' in log, log
    assert len(log.strip().splitlines()) <= 3, log
    _, log = run('-generation:instrumnt', 'Astral', ok=False)
    assert 'did you mean -generation:instrument' in log, log

    # ---- defaults -------------------------------------------------------------
    ini = root / 'new.ini'
    run('-write_ini', ini)
    tree = ET.parse(ini)
    items = {}
    def walk(node, prefix=''):
        for child in node:
            if child.tag == 'NODE':
                walk(child, prefix + child.get('name') + ':')
            elif child.tag in ('ITEM', 'ITEMLIST'):
                items[prefix + child.get('name')] = child
    walk(tree.getroot())
    base = 'DIALibGen:1:'
    assert items[base + 'threads'].get('value') == '0', 'write_ini must record all cores'
    for old in ('tune_models', 'tune_predict_gpu', 'machine:threads', 'search:passes', 'search:allow_bootstrap'):
        assert base + old not in items, old + ' written into a new INI'
    assert items[base + 'tune_heads'].get('value') == 'auto'
    assert items[base + 'q_global'].get('value') == '1.0' and items[base + 'q_protein'].get('value') == '1.0'
    refine, _ = config('-mode', 'refine')
    assert refine['q_precursor'] == 0.01 and refine['q_global'] == 1 and refine['q_protein'] == 1
    tune, _ = config('-mode', 'tune')
    assert tune['tuning']['tune_heads'] == 'auto' and tune['tuning']['machine:train_threads'] == 0
    generate, _ = config()
    assert generate['irt_standards'] == '' and generate['instrument'] == 'QE' and generate['model_instrument'] == 'QE'

    # ---- names removed in 0.13 (deprecated in 0.12): one line with the replacement ----
    models = root / 'models'
    models.mkdir()
    for args, message in ((('-tune_models', models), '-tune_models was removed in 0.13: use -models'),
                          (('-machine:threads', 3), '-machine:threads was removed in 0.13: use -machine:train_threads'),
                          (('-tune_predict_gpu',), '-tune_predict_gpu was removed in 0.13: re-prediction follows -machine:device'),
                          (('-search:passes', 1), '-search:passes was removed in 0.13: 1 was its only value, so delete it'),
                          (('-search:allow_bootstrap',), '-search:allow_bootstrap was removed in 0.13: set ODIA_SEARCH_ALLOW_BOOTSTRAP=1')):
        out = fresh()
        r, log = run('-mode', 'tune', *args, '-write_config', out, ok=False)
        assert r.returncode == 6 and 'Error: ' + message in log and not out.exists(), (args, log)
        assert len(log.strip().splitlines()) == 1, log
    # Two dashes, and a name given twice: still one line.
    r, log = run('--tune_models', models, '-tune_models', models, ok=False)
    assert r.returncode == 6 and log.strip() == 'Error: -tune_models was removed in 0.13: use -models', log
    tuned, _ = config('-mode', 'tune', '-models', models, '-machine:train_threads', 3)
    assert tuned['tuning']['machine:train_threads'] == 3, tuned['tuning']

    # ---- an INI written by 0.11 ------------------------------------------------
    old = root / 'old.ini'
    old_tree = ET.parse(ini)
    def item(path):
        node = old_tree.getroot()
        for name in path.split(':')[:-1]:
            node = next(n for n in node if n.tag == 'NODE' and n.get('name') == name)
        return node, next((n for n in node if n.get('name') == path.split(':')[-1]), None)
    def set_value(path, value):
        node, element = item(path)
        if element is None:
            element = ET.SubElement(node, 'ITEM', name=path.split(':')[-1], type='int', value=str(value),
                                    description='', required='false', advanced='true')
        element.set('value', str(value))
    def set_list(path, values):
        node, element = item(path)
        for child in list(element):
            element.remove(child)
        for v in values:
            ET.SubElement(element, 'LISTITEM', value=str(v))
    item('DIALibGen:version')[1].set('value', '0.11.5')
    for path, value in (('DIALibGen:1:tune_heads', 'both'), ('DIALibGen:1:q_global', '0.01'),
                        ('DIALibGen:1:q_protein', '0.01'), ('DIALibGen:1:threads', '1')):
        set_value(path, value)
    set_list('DIALibGen:1:generation:precursor_mz', [350.0, 1200.0])
    set_list('DIALibGen:1:generation:fragment_mz', [200.0, 1800.0])
    old_tree.write(old, encoding='ISO-8859-1', xml_declaration=True)
    # The same with names removed in 0.13 (a full 0.11 or 0.12 INI has them).
    old_training = root / 'old-training.ini'
    set_value('DIALibGen:1:machine:threads', '2')
    set_value('DIALibGen:1:search:passes', '1')
    old_tree.write(old_training, encoding='ISO-8859-1', xml_declaration=True)
    for mode in ('generate', 'refine', 'tune'):
        effective, log = config('-mode', mode, '-ini', old)
        assert 'the INI is from DIALibGen 0.11.5 and keeps its values' in log, (mode, log)
    effective, log = config('-mode', 'tune', '-ini', old)
    assert effective['tuning']['tune_heads'] == 'both', effective['tuning']
    effective, log = config('-mode', 'refine', '-ini', old)
    assert effective['q_global'] == 0.01, 'an INI value applies as given'
    effective, _ = config('-mode', 'generate', '-ini', old, '-generation:instrument', 'Astral')
    assert effective['precursor_mz'] == [350.0, 1200.0], 'an INI window applies as given'
    for args in (('-mode', 'tune', '-write_config', fresh()), ('-write_ini', fresh('.ini'))):
        r, log = run('-ini', old_training, *args, ok=False)
        assert r.returncode == 6, log
        assert f'Error: {old_training}: machine:threads was removed in 0.13: use -machine:train_threads' in log, log
        assert f'Error: {old_training}: search:passes was removed in 0.13: 1 was its only value' in log, log
        assert len(log.strip().splitlines()) == 2 and 'Unknown (or deprecated) Parameter' not in log, log
    # The common section, a non-default instance, a TOPP JSON parameter file; the command line and the INI
    # naming the same option give one line.
    common = root / 'common.ini'
    common.write_text(old.read_text().replace('</PARAMETERS>', '<NODE name="common" description=""><ITEM name="tune_models" '
                                              'value="/m" type="string" description="" /></NODE></PARAMETERS>'))
    for args in (('-write_config', fresh()), ('-write_ini', fresh('.ini'))):
        r, log = run('-ini', common, *args, ok=False)
        assert r.returncode == 6 and f'Error: {common}: tune_models was removed in 0.13: use -models' in log, log
    second = root / 'instance2.ini'
    second.write_text(old_training.read_text().replace('<NODE name="1"', '<NODE name="2"'))
    r, log = run('-ini', second, '-instance', 2, '-write_config', fresh(), ok=False)
    assert r.returncode == 6 and f'{second}: machine:threads was removed' in log, log
    for body in ({'tune_models': '/m'}, {'DIALibGen': {'1': {'machine': {'threads': 2}}}}):
        topp_json = root / f'params-{next(counter)}.json'
        topp_json.write_text(json.dumps(body))
        r, log = run('-ini', topp_json, '-write_config', fresh(), ok=False)
        assert r.returncode == 6 and f'Error: {topp_json}: ' in log and 'was removed in 0.13' in log, (body, log)
    r, log = run('-ini', old_training, '-machine:threads', 2, '-write_config', fresh(), ok=False)
    assert log.count('machine:threads was removed') == 1, log
    upgraded = root / 'upgraded.ini'
    run('-write_ini', upgraded, '-ini', old)
    text = upgraded.read_text()
    assert 'name="threads" value="1"' in text and 'name="tune_heads" value="both"' in text, 'INI upgrade'
    assert text.count('name="threads"') == 1, text

    # ---- mode inference -----------------------------------------------------------
    library = root / 'lib.tsv'
    library.write_text('Precursor.Id\tModified.Sequence\n')
    fasta = root / 'proteins.tsv'
    fasta.write_text('\ufeff>sp|P1|X\nPEPTIDEK\n', encoding='utf-8')
    _, log = run('-in', library, '-ids', library, '-out', fresh('.tsv'), ok=False)
    assert 'give -mode refine or -mode tune' in log, log
    # -write_config keeps generate (a library as -in only names the recipe's mode).
    effective, _ = config('-in', library)
    assert 'enzyme' in effective
    # A stock INI's "generate" does not block it; a FASTA by content stays generate.
    _, log = run('-ini', ini, '-in', library, '-out', fresh('.tsv'), ok=False)
    assert 'no -mode: -in is a library, so -mode convert' in log, log
    _, log = run('-in', fasta, '-out', fresh('.tsv'), ok=False)
    assert '-mode convert' not in log, log

    # ---- instrument files -----------------------------------------------------------
    for name in ('Astral', 'astral', 'Orbitrap Astral', 'AstralZoom', 'Astral Zoom'):
        effective, log = config('-generation:instrument', name)
        assert effective['precursor_mz'] == [380.0, 980.0] and effective['fragment_mz'] == [150.0, 2000.0], name
        assert effective['nce'] == 25 and effective['model_instrument'] == 'Lumos' and effective['instrument'] == name, name
        assert effective['instrument_defaults']['instrument'] == 'astral.yaml', name
        assert 'instrument defaults: astral.yaml (Astral) over default.yaml' in log, log
    # The Orbitrap and SCIEX files (values sourced in each file's comments): every alias, as the PSI-MS/PRIDE
    # spellings carry them, resolves to its file; the model slot names themselves keep default.yaml.
    shipped = {
        'exploris.yaml': (('Exploris 480', 'Orbitrap Exploris 480', 'exploris480'),
                          dict(model_instrument='QE', nce=30, precursor_mz=[360.0, 1040.0], fragment_mz=[200.0, 1800.0],
                               precursor_charges=[1, 2, 3, 4])),
        'qexactive.yaml': (('Q Exactive HF-X', 'Q Exactive HF', 'QE HF-X', 'QE HF', 'QEHFX'),
                           dict(model_instrument='QE', nce=27, precursor_mz=[350.0, 1200.0], fragment_mz=[200.0, 1800.0],
                                precursor_charges=[1, 2, 3, 4])),
        'tribrid.yaml': (('Fusion Lumos', 'Orbitrap Fusion Lumos', 'Eclipse', 'Orbitrap Eclipse', 'Orbitrap Tribrid'),
                         dict(model_instrument='Lumos', nce=30, precursor_mz=[350.0, 1200.0], fragment_mz=[200.0, 2000.0],
                              precursor_charges=[1, 2, 3, 4])),
        'sciextof.yaml': (('TripleTOF', 'TripleTOF 6600', 'TripleTOF 5600+', 'ZenoTOF', 'ZenoTOF 7600', 'ZenoTOF 8600'),
                          dict(model_instrument='Lumos', nce=27, precursor_mz=[350.0, 1250.0], fragment_mz=[200.0, 1800.0],
                               precursor_charges=[2, 3, 4])),
    }
    for file, (names, values) in shipped.items():
        for name in names:
            effective, log = config('-generation:instrument', name)
            assert effective['instrument_defaults']['instrument'] == file, (name, effective['instrument_defaults'])
            for key, value in values.items():
                assert effective[key] == value, (name, key, effective[key], value)
            assert effective['instrument'] == name and effective['nce_source'].startswith('instrument-default:'), name
    for name, slot, nce in (('QE', 'QE', 30), ('Lumos', 'Lumos', 25), ('SciexTOF', 'SciexTOF', 30),
                            ('Exploris', 'QE', 30), ('Exploris 240', 'QE', 30)):
        effective, _ = config('-generation:instrument', name)
        assert effective['instrument_defaults']['instrument'] is None, (name, effective['instrument_defaults'])
        assert effective['model_instrument'] == slot and effective['nce'] == nce, (name, effective)
    # an explicit option still wins over the new files
    effective, _ = config('-generation:instrument', 'Eclipse', '-generation:nce', 25, '-generation:precursor_mz', 400, 1000)
    assert effective['nce'] == 25 and effective['nce_source'] == 'TOPP' and effective['precursor_mz'] == [400.0, 1000.0]
    assert effective['fragment_mz'] == [200.0, 2000.0], effective
    effective, log = config('-generation:instrument', 'timsTOF Pro 2')
    assert effective['precursor_charges'] == [2, 3, 4] and effective['model_instrument'] == 'timsTOF'
    assert 'has ion mobility: CCS and 1/K0 predicted' in log, log
    effective, log = config()
    assert effective['precursor_mz'] == [350.0, 1200.0] and effective['instrument_defaults']['instrument'] is None
    assert 'no -generation:instrument: QE assumed' in log, log
    # Precedence: the command line, then the recipe, over the instrument file.
    effective, _ = config('-generation:instrument', 'Astral', '-generation:precursor_mz', 400, 1000)
    assert effective['precursor_mz'] == [400.0, 1000.0] and effective['precursor_mz_source'] == 'TOPP'
    assert effective['fragment_mz'] == [150.0, 2000.0]
    recipe = root / 'recipe.json'
    recipe.write_text(json.dumps({'instrument': 'Astral', 'fragment_mz': [100, 1500]}))
    effective, _ = config('-config', recipe)
    assert effective['fragment_mz'] == [100.0, 1500.0] and effective['precursor_mz'] == [380.0, 980.0]
    # A user's directory: a file of the same name replaces the shipped one; a
    # new file adds an instrument; a user default.yaml replaces the defaults.
    user = root / 'user-instruments'
    user.mkdir()
    (user / 'astral.yaml').write_text('name: Astral\nmodel_instrument: Lumos\nprecursor_mz: [390, 990]\n')
    (user / 'exploris.yaml').write_text('name: Exploris 480\naliases: [E480]\nmodel_instrument: QE\nnce: 27\n'
                                        'precursor_mz: [400, 1000]\n')
    env = {'DIALIBGEN_INSTRUMENT_DIR': str(user)}
    effective, log = config('-generation:instrument', 'Astral', env=env)
    assert effective['precursor_mz'] == [390.0, 990.0] and effective['fragment_mz'] == [200.0, 1800.0], effective
    assert effective['nce'] == 25, 'Lumos model default'
    effective, _ = config('-generation:instrument', 'E480', env=env)
    assert effective['precursor_mz'] == [400.0, 1000.0] and effective['nce'] == 27 and effective['model_instrument'] == 'QE'
    # The user's config directory is searched too.
    config_dir = home / ('Library/Application Support' if sys.platform == 'darwin' else
                         'AppData' if sys.platform == 'win32' else '.config') / 'DIALibGen' / 'instruments'
    config_dir.mkdir(parents=True)
    (config_dir / 'default.yaml').write_text((data_dir / 'instruments' / 'default.yaml').read_text()
                                             .replace('instrument: QE', 'instrument: timsTOF'))
    effective, log = config()
    assert effective['instrument'] == 'timsTOF' and effective['precursor_charges'] == [2, 3, 4], effective
    (config_dir / 'default.yaml').unlink()
    # -instrument_file
    single = root / 'mine.yaml'
    single.write_text('name: My Orbitrap\nmodel_instrument: QE\nfragment_mz: [120, 1600]\nion_mobility: false\n')
    effective, log = config('-instrument_file', single)
    assert effective['fragment_mz'] == [120.0, 1600.0] and effective['instrument_defaults']['instrument'] == 'mine.yaml'
    _, log = run('-instrument_file', single, '-generation:instrument', 'Astral', '-write_config', fresh(), ok=False)
    assert 'is not the instrument of' in log, log
    effective, _ = config('-instrument_file', single, '-generation:instrument', 'my orbitrap')
    assert effective['fragment_mz'] == [120.0, 1600.0]
    # Strict validation: short errors that name the file and the key.
    for body, message in (
            ('name: X\nmodel_instrument: QE\nprecursor_window: [1, 2]\n', 'mine.yaml: precursor_window: unknown key'),
            ('name: X\nmodel_instrument: QE\nnce: high\n', 'mine.yaml: nce: must be a number'),
            ('name: X\nmodel_instrument: QE\nnce: 250\n', 'mine.yaml: nce: must be greater than 0 and at most 100'),
            ('name: X\nmodel_instrument: Orbitrap9000\n', 'mine.yaml: model_instrument:'),
            ('name: X\nmodel_instrument: QE\nprecursor_charges: [2, 2.5]\n', 'mine.yaml: precursor_charges:'),
            ('name: X\nmodel_instrument: QE\nprecursor_mz: [900, 400]\n', 'mine.yaml: precursor_mz: needs min < max'),
            ('name: X\nmodel_instrument: QE\nion_mobility: maybe\n', 'mine.yaml: ion_mobility: must be true or false'),
            ('model_instrument: QE\n', 'mine.yaml: name: missing'),
            ('name: X\nmodel_instrument: QE\nnce_by_charge: {2.5: 30}\n', 'mine.yaml: nce_by_charge: charges must be whole numbers'),
            ('name: X\nmodel_instrument: QE\nnce_by_charge: {3: 0}\n', 'mine.yaml: nce_by_charge.3: must be greater than 0'),
            ('name: X\nmodel_instrument: QE\nnce_by_charge: [30, 35]\n', 'mine.yaml: nce_by_charge: must map precursor charges'),
            ('name: X\nmodel_instrument: timsTOF\nmobility_correction:\n  intercept: 0.02\n  slope: 0.99\n', 'mine.yaml: mobility_correction.ccs_model_fnv1a64: missing'),
            ('name: X\nmodel_instrument: timsTOF\nmobility_correction:\n  ccs_model_fnv1a64: C26EB20995B991F5\n  intercept: 0.02\n  slope: 0.99\n', 'mine.yaml: mobility_correction.ccs_model_fnv1a64: must be 16 lowercase hex digits'),
            ('name: X\nmodel_instrument: timsTOF\nmobility_correction:\n  ccs_model_fnv1a64: c26eb20995b991f5\n  slope: 0.99\n', 'mine.yaml: mobility_correction: needs intercept and slope'),
            ('name: X\nmodel_instrument: timsTOF\nmobility_correction:\n  ccs_model_fnv1a64: c26eb20995b991f5\n  intercept: 0.02\n  slope: 3\n', 'mine.yaml: mobility_correction.slope: must be between'),
            ('name: X\nmodel_instrument: timsTOF\nmobility_correction:\n  ccs_model_fnv1a64: c26eb20995b991f5\n  intercept: 0.02\n  slope: 0.99\n  by_charge:\n    2.5: {intercept: 0, slope: 1}\n', 'mine.yaml: mobility_correction.by_charge: charges must be whole numbers'),
            ('name: X\nmodel_instrument: timsTOF\nmobility_correction:\n  ccs_model_fnv1a64: c26eb20995b991f5\n  intercept: 0.02\n  slope: 0.99\n  by_charge:\n    2: {intercept: 0}\n', 'mine.yaml: mobility_correction.by_charge.2.slope: must be a number'),
            ('name: X\nmodel_instrument: timsTOF\nmobility_correction:\n  ccs_model_fnv1a64: c26eb20995b991f5\n  intercept: 0.02\n  slope: 0.99\n  offset: 1\n', 'mine.yaml: mobility_correction.offset: unknown key'),
            ('[1, 2]\n', 'mine.yaml: must be a mapping'),
            ('name: [unclosed\n', 'mine.yaml: not valid YAML')):
        single.write_text(body)
        _, log = run('-instrument_file', single, '-write_config', fresh(), ok=False)
        assert message in log, (body, log)
    (user / 'broken.yaml').write_text('name: Y\nmodel_instrument: QE\ncolour: red\n')
    _, log = run('-write_config', fresh(), ok=False, env=env)
    assert 'broken.yaml: colour: unknown key' in log, log
    (user / 'broken.yaml').unlink()
    (user / 'twin.yaml').write_text('name: Exploris 480\nmodel_instrument: QE\n')
    _, log = run('-generation:instrument', 'Exploris 480', '-write_config', fresh(), ok=False, env=env)
    assert 'named by two files of one directory' in log, log
    # A custom instrument's recipe replays through its file (the recipe names it, not the model slot).
    custom = root / 'custom.yaml'
    custom.write_text('name: Custom IM\nmodel_instrument: QE\nion_mobility: true\n')
    effective, log = config('-instrument_file', custom)
    assert effective['instrument'] == 'Custom IM' and effective['model_instrument'] == 'QE'
    assert 'has ion mobility' in log, log
    replay = root / 'replay.json'
    replay.write_text(json.dumps(effective))
    _, log = config('-config', replay, '-instrument_file', custom)
    assert 'has ion mobility' in log, log
    # nce_by_charge: the instrument's per-charge NCE goes with its NCE; -generation:nce (one value) replaces both;
    # a recipe replays it (null = one NCE); a -config JSON may set it.
    by_charge = root / 'bycharge.yaml'
    by_charge.write_text('name: Per Charge\nmodel_instrument: timsTOF\nnce: 30\nnce_by_charge: {3: 25, 4: 40}\n')
    effective, log = config('-instrument_file', by_charge)
    assert effective['nce_by_charge'] == {'3': 25.0, '4': 40.0} and effective['nce_by_charge_source'] == 'instrument-default:Per Charge', effective
    assert 'by precursor charge: 3+ 25 4+ 40; other charges 30' in log, log
    effective, _ = config('-instrument_file', by_charge, '-generation:nce', 28)
    assert effective['nce'] == 28 and effective['nce_by_charge'] is None, effective
    replay = root / 'replay-bycharge.json'
    replay.write_text(json.dumps(config('-instrument_file', by_charge)[0]))
    effective, _ = config('-config', replay, '-instrument_file', by_charge)
    assert effective['nce_by_charge'] == {'3': 25.0, '4': 40.0} and effective['nce_by_charge_source'] == 'config', effective
    replay.write_text(json.dumps({'instrument': 'timsTOF', 'nce_by_charge': {'2': 22, '3': 24}}))
    effective, log = config('-config', replay)
    assert effective['nce_by_charge'] == {'2': 22.0, '3': 24.0} and 'by precursor charge: 2+ 22 3+ 24; other charges 30' in log, (effective, log)
    # an explicit -generation:nce is one NCE for every charge, over a recipe's map too
    replay.write_text(json.dumps({'instrument': 'timsTOF', 'nce': 30, 'nce_by_charge': {'3': 20}}))
    effective, _ = config('-config', replay, '-generation:nce', 28)
    assert effective['nce'] == 28 and effective['nce_by_charge'] is None and effective['nce_by_charge_source'] == 'TOPP', effective
    # a recipe written before 0.13 (nce_source, no mobility_correction) replays without the 1/K0 correction
    replay.write_text(json.dumps({'instrument': 'timsTOF', 'nce': 30.0, 'nce_source': 'instrument-default:timsTOF',
                                  'precursor_charges': [1, 2, 3, 4]}))
    effective, log = config('-config', replay)
    assert effective['mobility_correction'] == 'off' and 'predates the 1/K0 correction' in log, (effective, log)
    replay.write_text(json.dumps({'instrument': 'timsTOF'}))
    effective, _ = config('-config', replay)
    assert effective['mobility_correction'] == 'auto', effective
    for bad, message in (({'nce_by_charge': {'x': 30}}, "nce_by_charge: 'x' is not a charge"),
                         ({'nce_by_charge': {'9': 30}}, "nce_by_charge: '9' is not a charge"),
                         ({'nce_by_charge': {'2': 300}}, 'nce_by_charge: the NCE of charge 2 must be'),
                         ({'nce_by_charge': [30]}, 'nce_by_charge must be null or an object')):
        replay.write_text(json.dumps(bad))
        _, log = run('-config', replay, '-write_config', fresh(), ok=False)
        assert message in log, (bad, log)
    # "Orbitrap Exploris 480" names a file called "Exploris 480"; "Exploris" alone does not.
    (user / 'twin.yaml').unlink()
    effective, _ = config('-generation:instrument', 'Orbitrap Exploris 480', env=env)
    assert effective['nce'] == 27, effective
    effective, _ = config('-generation:instrument', 'Exploris', env=env)
    assert effective['nce'] == 30, effective
    single.write_text('name: X\nmodel_instrument: QE\nnce: 25\nnce: 30\n')
    _, log = run('-instrument_file', single, '-write_config', fresh(), ok=False)
    assert 'nce: given twice' in log, log
    # Recipe training values are range-checked like the options.
    bad = root / 'bad-tuning.json'
    bad.write_text(json.dumps({'tuning': {'cohort:train_size': -1}}))
    _, log = run('-mode', 'tune', '-config', bad, '-write_config', fresh(), ok=False)
    assert 'tuning.cohort:train_size must be between' in log, log
    # -write_ini refuses an INI with a removed name, and leaves no file.
    target = fresh('.ini')
    _, log = run('-write_ini', target, '-ini', old_training, ok=False)
    assert 'removed in 0.13' in log and not target.exists(), log
    _, log = run('-generation:instrument', 'Astrall', '-write_config', fresh(), ok=False)
    assert "unknown instrument: 'Astrall'" in log, log

    # ---- settings without an effect: refused (a WARNING in 0.12) ------------------
    for args, message in (
            (('-generation:recompute_decoy_mz', 'true'),
             '-generation:recompute_decoy_mz true has no effect with -generation:decoys none: add -generation:decoys'),
            (('-generation:max_variable_modifications', 2),
             '-generation:max_variable_modifications has no effect without -generation:variable_modifications: add'),
            (('-generation:derive_ion_mobility', 'false'),
             '-generation:derive_ion_mobility false has no effect when CCS is not predicted: remove it, or predict CCS'),
            (('-mode', 'refine', '-intensity_min_fragments', 5),
             '-intensity_min_fragments has no effect without -write_intensity: add -write_intensity or remove'),
            (('-mode', 'refine', '-allow_mixed_intensity'), '-allow_mixed_intensity has no effect without -write_intensity'),
            (('-mode', 'refine', '-im_ramp_top', 1.3), '-im_ramp_top has no effect without -write_im: add -write_im')):
        out = fresh()
        r, log = run(*args, '-write_config', out, ok=False)
        assert r.returncode == 6 and message in log and not out.exists(), (args, log)
    # ... except the two that only ever warn.
    for args, message in (
            (('-irt_standards', root / 'none.tsv'), '-irt_standards has no effect without -generation:irt_rescale true'),
            (('-machine:device', 'cpu', '-machine:coreml_units', 'cpu'), '-machine:coreml_units cpu has no effect')):
        _, log = config(*args)
        assert 'WARNING: ' + message in log, (args, log)
    # With what they need, they are taken.
    _, log = config('-mode', 'refine', '-write_intensity', '-intensity_min_fragments', 5)
    assert 'has no effect' not in log, log
    _, log = config('-generation:decoys', 'pseudo_reverse', '-generation:recompute_decoy_mz', 'true')
    assert 'has no effect' not in log, log
    _, log = config('-generation:instrument', 'timsTOF', '-generation:derive_ion_mobility', 'false')
    assert 'has no effect' not in log, log
print('command line, defaults, removed names, old INIs, instrument files, settings without an effect: PASS')
