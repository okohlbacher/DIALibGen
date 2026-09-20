#!/usr/bin/env python3
"""Exercise TOPP/JSON precedence and pure-tuning defaults without model inference."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET

binary = str(Path(sys.argv[1]).resolve())
with tempfile.TemporaryDirectory() as directory:
    root = Path(directory)
    def run(*args, ok=True):
        r = subprocess.run([binary, *map(str, args)], capture_output=True, text=True)
        assert (r.returncode == 0) == ok, (args, r.returncode, r.stdout, r.stderr)
        return r
    def config(*args):
        out = root / 'effective.json'
        run(*args, '-write_config', out)
        return json.loads(out.read_text())
    fixture = Path(__file__).resolve().parents[1] / 'gui/src/testing/generation-defaults.json'
    expected, actual = json.loads(fixture.read_text()), config('-mode', 'generate')
    for key in ('rt_model', 'ms2_model', 'ccs_model'):
        expected.pop(key); actual.pop(key)  # Resolved paths depend on the installation.
    assert actual == expected, 'GUI default fixture differs from the native effective configuration'
    recipe = root / 'recipe.json'
    recipe.write_text(json.dumps({'missed_cleavages': 2, 'instrument': 'Lumos', 'nce': 28}))
    effective = config('-config', recipe, '-generation:missed_cleavages', 1,
                       '-generation:n_terminal_methionine_excision', 'false',
                       '-generation:precursor_charges', 2, 3, '-generation:nce', -1)
    assert effective['missed_cleavages'] == 1
    assert effective['precursor_charges'] == [2, 3]
    assert effective['n_terminal_methionine_excision'] is False
    assert effective['nce'] == 25
    ini = root / 'settings.ini'
    run('-write_ini', ini)
    # A full default INI includes every mode; untouched settings must stay usable.
    assert config('-mode', 'refine', '-ini', ini)['filter'] is True
    assert config('-mode', 'tune', '-ini', ini)['filter'] is False
    recipe.write_text('{"q_precursor":0.5}')
    assert config('-mode', 'refine', '-config', recipe)['q_precursor'] == 0.5
    assert config('-mode', 'refine', '-ini', ini, '-config', recipe)['q_precursor'] == 0.01
    assert config('-mode', 'refine', '-ini', ini, '-config', recipe, '-q_precursor', 0.2)['q_precursor'] == 0.2
    recipe.write_text(json.dumps({'missed_cleavages': 2, 'instrument': 'Lumos', 'nce': 28}))
    training_ini = root / 'training.ini'
    training_tree = ET.parse(ini)
    training = next(n for n in training_tree.iter('NODE') if n.get('name') == 'train')
    next(n for n in training if n.get('name') == 'epochs').set('value', '999')
    training_tree.write(training_ini, encoding='utf-8', xml_declaration=True)
    result = run('-mode', 'refine', '-ini', training_ini, '-write_config', root / 'invalid.json', ok=False)
    assert 'enable -tune' in result.stdout + result.stderr
    assert config('-mode', 'refine', '-ini', training_ini, '-tune')['filter'] is True
    tree = ET.parse(ini)
    generation = next(n for n in tree.iter('NODE') if n.get('name') == 'generation')
    next(n for n in generation if n.get('name') == 'instrument').set('value', 'Lumos')
    next(n for n in generation if n.get('name') == 'missed_cleavages').set('value', '3')
    tree.write(ini, encoding='utf-8', xml_declaration=True)
    effective = config('-ini', ini, '-config', recipe)
    assert effective['missed_cleavages'] == 3 and effective['nce'] == 25
    assert config('-ini', ini, '-generation:missed_cleavages', 0)['missed_cleavages'] == 0
    run('-mode', 'refine', '-ini', ini, '-write_config', root / 'invalid.json', ok=False)
    for mode in ('refine', 'tune'):
        result = run('-mode', mode, '-generation:rt_model', 'ignored.onnx',
                     '-write_config', root / 'invalid.json', ok=False)
        assert 'has no effect' in result.stdout + result.stderr
    for args in (('-tune',), ('-q_precursor', 0.5), ('-train:epochs', 1)):
        result = run('-mode', 'generate', *args, '-write_config', root / 'invalid.json', ok=False)
        assert 'has no effect' in result.stdout + result.stderr
    for args in (('-train:epochs', 999), ('-tune_heads', 'rt'), ('-machine:threads', 2), ('-filter:q_value', 0.5)):
        result = run('-mode', 'refine', *args, '-write_config', root / 'invalid.json', ok=False)
        assert 'enable -tune' in result.stdout + result.stderr
        assert config('-mode', 'refine', '-tune', *args)['filter'] is True
    for args in (('-q_precursor', 0.001), ('-q_global', 1), ('-min_fragments', 3), ('-empirical_library',)):
        result = run('-mode', 'tune', *args, '-write_config', root / 'invalid.json', ok=False)
        assert 'has no effect' in result.stdout + result.stderr
    assert not config('-mode', 'tune', '-filter:q_value', 0.001, '-no_filter', '-no_write_rt')['filter']
    for bad in ['[]', '{"peptide_length":[7.5,30]}', '{"nce":101}', '{"precursor_charges":[2.5]}', '{"missed_cleavages":-1}']:
        recipe.write_text(bad)
        run('-config', recipe, '-write_config', root / 'invalid.json', ok=False)
    tune = config('-mode', 'tune')
    assert not any(tune[k] for k in ('filter', 'write_rt', 'write_im', 'write_intensity'))
    recipe.write_text('{"filter":true}')
    run('-mode', 'tune', '-config', recipe, '-write_config', root / 'invalid.json', ok=False)
    for config_key, value in [('q_precursor', 0.001), ('q_global', 1), ('min_fragments', 3), ('require_gates', False)]:
        recipe.write_text(json.dumps({config_key: value}))
        run('-mode', 'tune', '-config', recipe, '-write_config', root / 'invalid.json', ok=False)
    recipe.write_text(json.dumps(tune))
    assert config('-mode', 'tune', '-config', recipe) == tune
print('TOPP parameters, precedence and tune isolation: PASS')
