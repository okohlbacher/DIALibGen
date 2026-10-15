#!/usr/bin/env python3
"""The 0.12 command line without model inference: help, unknown options, mode
inference, deprecated names (command line and INI), INIs from 0.11, the new
defaults, the instrument files and the warnings for settings without an effect.

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
    assert len(options(all_text)) > 120 and 'Deprecated (0.12, removed in 0.13)' in all_text
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
    assert generate['irt_standards'] == '' and generate['instrument'] == 'QE'

    # ---- deprecated names ------------------------------------------------------
    models = root / 'models'
    models.mkdir()
    tuned, log = config('-mode', 'tune', '-tune_models', models, '-machine:threads', 3)
    assert tuned['tuning']['machine:train_threads'] == 3, tuned['tuning']
    assert 'WARNING: -tune_models is deprecated (0.12, removed in 0.13): use -models' in log, log
    assert 'WARNING: -machine:threads is deprecated (0.12, removed in 0.13): use -machine:train_threads' in log, log
    same, log = config('-mode', 'tune', '-machine:threads', 3, '-machine:train_threads', 3)
    assert same['tuning']['machine:train_threads'] == 3
    _, log = run('-mode', 'tune', '-machine:threads', 3, '-machine:train_threads', 4, '-write_config', fresh(), ok=False)
    assert 'disagree' in log and 'give only -machine:train_threads' in log, log
    _, log = config('-mode', 'tune', '-tune_predict_gpu')
    assert 'WARNING: -tune_predict_gpu is deprecated' in log, log
    # A deprecated name keeps its old mode check.
    _, log = run('-mode', 'generate', '-machine:threads', 2, '-write_config', fresh(), ok=False)
    assert 'has no effect in -mode generate' in log, log

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
                        ('DIALibGen:1:q_protein', '0.01'), ('DIALibGen:1:threads', '1'),
                        ('DIALibGen:1:machine:threads', '0')):
        set_value(path, value)
    set_list('DIALibGen:1:generation:precursor_mz', [350.0, 1200.0])
    set_list('DIALibGen:1:generation:fragment_mz', [200.0, 1800.0])
    old_tree.write(old, encoding='ISO-8859-1', xml_declaration=True)
    # The same with a training setting under its old name.
    old_training = root / 'old-training.ini'
    set_value('DIALibGen:1:machine:threads', '2')
    old_tree.write(old_training, encoding='ISO-8859-1', xml_declaration=True)
    for mode in ('generate', 'refine', 'tune'):
        effective, log = config('-mode', mode, '-ini', old)
        assert 'the INI is from DIALibGen 0.11.5 and keeps its values' in log, (mode, log)
    effective, log = config('-mode', 'tune', '-ini', old_training)
    assert effective['tuning']['tune_heads'] == 'both' and effective['tuning']['machine:train_threads'] == 2
    assert 'WARNING: -machine:threads is deprecated' in log, log
    effective, log = config('-mode', 'refine', '-ini', old)
    assert effective['q_global'] == 0.01, 'an INI value applies as given'
    effective, _ = config('-mode', 'generate', '-ini', old, '-generation:instrument', 'Astral')
    assert effective['precursor_mz'] == [350.0, 1200.0], 'an INI window applies as given'
    _, log = run('-mode', 'tune', '-ini', old_training, '-machine:train_threads', 5, '-write_config', fresh(), ok=False)
    assert 'disagree' in log, log
    upgraded = root / 'upgraded.ini'
    run('-write_ini', upgraded, '-ini', old_training)
    text = upgraded.read_text()
    assert 'name="train_threads" value="2"' in text and 'name="threads" value="1"' in text, 'INI upgrade'
    assert text.count('name="threads"') == 1, 'the deprecated machine:threads must not be written'

    # ---- mode inference -----------------------------------------------------------
    library = root / 'lib.tsv'
    library.write_text('Precursor.Id\tModified.Sequence\n')
    fasta = root / 'proteins.tsv'
    fasta.write_text('﻿>sp|P1|X\nPEPTIDEK\n')
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
        assert effective['nce'] == 25 and effective['instrument'] == 'Lumos', name
        assert effective['instrument_defaults']['instrument'] == 'astral.yaml', name
        assert 'instrument defaults: astral.yaml (Astral) over default.yaml' in log, log
    effective, log = config('-generation:instrument', 'timsTOF Pro 2')
    assert effective['precursor_charges'] == [2, 3, 4] and effective['instrument'] == 'timsTOF'
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
    assert effective['precursor_mz'] == [400.0, 1000.0] and effective['nce'] == 27 and effective['instrument'] == 'QE'
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
    _, log = run('-generation:instrument', 'Astrall', '-write_config', fresh(), ok=False)
    assert "unknown instrument: 'Astrall'" in log, log

    # ---- settings without an effect: a WARNING (refused from 0.13) -----------------
    for args, message in (
            (('-irt_standards', root / 'none.tsv'), '-irt_standards has no effect without -generation:irt_rescale true'),
            (('-generation:recompute_decoy_mz', 'true'), '-generation:recompute_decoy_mz true has no effect with'),
            (('-generation:max_variable_modifications', 2), '-generation:max_variable_modifications has no effect'),
            (('-generation:derive_ion_mobility', 'false'), '-generation:derive_ion_mobility false has no effect'),
            (('-mode', 'refine', '-intensity_min_fragments', 5), '-intensity_min_fragments has no effect without -write_intensity'),
            (('-mode', 'refine', '-allow_mixed_intensity'), '-allow_mixed_intensity has no effect without -write_intensity'),
            (('-mode', 'refine', '-im_ramp_top', 1.3), '-im_ramp_top has no effect without -write_im'),
            (('-machine:device', 'cpu', '-machine:coreml_units', 'cpu'), '-machine:coreml_units cpu has no effect')):
        _, log = config(*args)
        assert 'WARNING: ' + message in log, (args, log)
    _, log = config('-mode', 'refine', '-write_intensity', '-intensity_min_fragments', 5)
    assert 'has no effect' not in log, log
print('command line, defaults, deprecated names, old INIs, instrument files: PASS')
