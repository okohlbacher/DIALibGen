#!/usr/bin/env python3
"""Source distribution roundtrip and missing/corrupt release asset checks."""
import hashlib
import json
from pathlib import Path
import runpy
import tempfile

module = runpy.run_path(str(Path(__file__).with_name('source-archive.py')))
with tempfile.TemporaryDirectory() as temporary:
    root = Path(temporary)
    for size in (1, 63, 64, 65, 96):
        archive = root / f'sources-{size}.tar.gz'
        original = bytes(range(size))
        archive.write_bytes(original)
        manifest = module['package'](archive, maximum=64, part_size=32)
        assert json.loads(archive.with_name(archive.name.removesuffix('.tar.gz') + '.json').read_text()) == manifest
        reconstructed = b''.join((root / part['name']).read_bytes() for part in manifest['files'])
        assert reconstructed == original
        assert hashlib.sha256(reconstructed).hexdigest() == manifest['archive']['sha256']
        assert archive.exists() == (size < 64)
        assets = [{'name': part['name'], 'size': part['size'], 'digest': 'sha256:' + part['sha256'],
                   'state': 'uploaded'} for part in manifest['files']]
        module['verify_assets'](manifest, assets)
        for broken in (assets[:-1], [{**assets[0], 'digest': 'sha256:' + '0' * 64}, *assets[1:]],
                       [{**assets[0], 'size': assets[0]['size'] + 1}, *assets[1:]]):
            try:
                module['verify_assets'](manifest, broken)
            except ValueError:
                pass
            else:
                raise AssertionError('incomplete or corrupt source asset accepted')
        if len(manifest['files']) > 1:
            try:
                module['verify_assets']({**manifest, 'files': list(reversed(manifest['files']))}, assets)
            except ValueError:
                pass
            else:
                raise AssertionError('reordered parts accepted')
print('Source archives: threshold, byte-exact reconstruction, hashes and release completeness passed')
