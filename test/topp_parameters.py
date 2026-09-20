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
    tree = ET.parse(ini)
    generation = next(n for n in tree.iter('NODE') if n.get('name') == 'generation')
    next(n for n in generation if n.get('name') == 'instrument').set('value', 'Lumos')
    next(n for n in generation if n.get('name') == 'missed_cleavages').set('value', '3')
    tree.write(ini, encoding='utf-8', xml_declaration=True)
    effective = config('-ini', ini, '-config', recipe)
    assert effective['missed_cleavages'] == 3 and effective['nce'] == 25
    assert config('-ini', ini, '-generation:missed_cleavages', 0)['missed_cleavages'] == 0
    for bad in ['[]', '{"peptide_length":[7.5,30]}', '{"nce":101}', '{"missed_cleavages":-1}']:
        recipe.write_text(bad)
        run('-config', recipe, '-write_config', root / 'invalid.json', ok=False)
    tune = config('-mode', 'tune')
    assert not any(tune[k] for k in ('filter', 'write_rt', 'write_im', 'write_intensity'))
    recipe.write_text('{"filter":true}')
    run('-mode', 'tune', '-config', recipe, '-write_config', root / 'invalid.json', ok=False)
print('TOPP parameters, precedence and tune isolation: PASS')
