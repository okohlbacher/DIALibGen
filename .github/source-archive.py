#!/usr/bin/env python3
"""Publish complete source archives within GitHub's per-asset size limit."""
import hashlib
import json
from pathlib import Path
import re
import sys

MAXIMUM = 2 * 1024**3
PART_SIZE = 1024**3


def package(archive, maximum=MAXIMUM, part_size=PART_SIZE):
    archive = Path(archive)
    size = archive.stat().st_size
    if not archive.name.endswith('.tar.gz') or not size:
        raise ValueError('expected a nonempty .tar.gz source archive')
    files = []
    whole = hashlib.sha256()
    split = size >= maximum
    with archive.open('rb') as source:
        if split:
            for index in range((size + part_size - 1) // part_size):
                destination = archive.with_name(f'{archive.name}.part-{index + 1:03d}')
                part_hash, written = hashlib.sha256(), 0
                with destination.open('wb') as output:
                    while written < part_size:
                        block = source.read(min(1024**2, part_size - written))
                        if not block:
                            break
                        output.write(block)
                        part_hash.update(block)
                        whole.update(block)
                        written += len(block)
                files.append({'name': destination.name, 'size': written, 'sha256': part_hash.hexdigest()})
        else:
            while block := source.read(1024**2):
                whole.update(block)
            files.append({'name': archive.name, 'size': size, 'sha256': whole.hexdigest()})
    if sum(part['size'] for part in files) != size or any(not 0 < part['size'] < maximum for part in files):
        raise ValueError('source archive parts are incomplete or exceed the asset limit')
    manifest = {'format': 1, 'archive': {'name': archive.name, 'size': size, 'sha256': whole.hexdigest()},
                'files': files,
                'instructions': 'Download every file listed here. Verify each SHA-256. If split, concatenate files in listed order into archive.name, then verify archive.sha256 and extract the .tar.gz.'}
    path = archive.with_name(archive.name.removesuffix('.tar.gz') + '.json')
    path.write_text(json.dumps(manifest, indent=2) + '\n')
    if split:
        archive.unlink()  # Only the complete, verified parts are uploaded.
    for part in files:
        print(f"{part['sha256']}  {part['name']} ({part['size']:,} bytes)")
    print(f'Source manifest: {path}')
    return manifest


def verify_assets(manifest, assets):
    """Check the release service's hashes and sizes against the build manifest."""
    archive, files = manifest['archive'], manifest['files']
    name = archive['name']
    if (manifest.get('format') != 1 or Path(name).name != name or not name.endswith('.tar.gz')
            or not re.fullmatch(r'[0-9a-f]{64}', archive['sha256']) or not files):
        raise ValueError('invalid source archive manifest')
    names = [part['name'] for part in files]
    expected = [name] if len(files) == 1 else [f'{name}.part-{index + 1:03d}' for index in range(len(files))]
    if names != expected or sum(part['size'] for part in files) != archive['size']:
        raise ValueError('source archive parts are incomplete or out of order')
    if len(files) == 1 and files[0]['sha256'] != archive['sha256']:
        raise ValueError('unsplit source archive checksum differs from manifest')
    by_name = {asset['name']: asset for asset in assets}
    for part in files:
        asset = by_name.get(part['name'])
        if (not 0 < part['size'] < MAXIMUM or not re.fullmatch(r'[0-9a-f]{64}', part['sha256'])
                or not asset or asset.get('state') != 'uploaded' or asset['size'] != part['size']
                or asset.get('digest') != 'sha256:' + part['sha256']):
            raise ValueError(f"missing, incomplete or mismatched source asset: {part['name']}")


if __name__ == '__main__':
    if len(sys.argv) >= 4 and sys.argv[1] == '--verify-assets':
        assets = json.loads(Path(sys.argv[2]).read_text(encoding='utf-8'))['assets']
        for filename in sys.argv[3:]:
            verify_assets(json.loads(Path(filename).read_text(encoding='utf-8')), assets)
            print(f'Verified source assets: {filename}')
    elif len(sys.argv) == 2:
        package(sys.argv[1])
    else:
        raise SystemExit('usage: source-archive.py ARCHIVE.tar.gz | --verify-assets RELEASE.json MANIFEST.json ...')
