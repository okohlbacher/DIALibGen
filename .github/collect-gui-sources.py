#!/usr/bin/env python3
"""Accompany GUI copyleft dependencies with their checksum-pinned registry sources."""
import argparse
import base64
import json
from pathlib import Path
import re
import runpy
import shutil

helpers = runpy.run_path(str(Path(__file__).with_name('collect-runtime-licenses.py')))
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--inventory', type=Path, required=True)
parser.add_argument('--sources', type=Path, required=True)
args = parser.parse_args()
inventory = json.loads(args.inventory.read_text(encoding='utf-8'))
root = args.sources / 'gui'
root.mkdir(parents=True, exist_ok=True)
shutil.copy2(args.inventory, root / 'inventory.json')
shutil.copy2(args.inventory.with_name('THIRD_PARTY_NOTICES.txt'), root / 'THIRD_PARTY_NOTICES.txt')
components = []
for package in inventory['packages']:
    if not re.search(r'GPL|EPL|MPL|CDDL', package['license'], re.I):
        continue
    name = package['ecosystem'] + '-' + package['name'].replace('/', '-') + '-' + package['version']
    record = {'url': package['source'], 'fn': name + '.tar.gz'}
    if package['ecosystem'] == 'cargo':
        record['sha256'] = package['checksum']
    elif package['ecosystem'] == 'npm':
        algorithm, value = package['integrity'].split('-', 1)
        record[algorithm] = base64.b64decode(value).hex()
    else:
        raise RuntimeError(f'unknown GUI package registry: {package}')
    sources = helpers['collect_sources']([record], root / name, root / 'notices' / name)
    components.append({'name': name, 'license': package['license'], 'sources': sources})
(root / 'corresponding-sources.json').write_text(json.dumps(components, indent=2) + '\n', encoding='utf-8')
print(f'Accompanying GUI sources: {len(components)} components for {inventory["target"]}')
