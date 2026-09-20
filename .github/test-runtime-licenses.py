#!/usr/bin/env python3
"""Offline release inventory contracts: exact owners, notices, source hashes, failures."""
import hashlib
import io
import json
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile

script = Path(__file__).with_name('collect-runtime-licenses.py')
with tempfile.TemporaryDirectory() as temporary:
    root = Path(temporary)
    prefix, stage, cached, output = [root / name for name in ('prefix', 'stage', 'cached', 'sources')]
    for directory in (prefix / 'conda-meta', prefix / 'lib', stage / 'lib', cached / 'info/licenses', cached / 'info/recipe'):
        directory.mkdir(parents=True)
    source = root / 'source.tar.gz'
    with tarfile.open(source, 'w:gz') as archive:
        for name, content in [('example/COPYING', b'upstream license'), ('example/third-party/NOTICE.txt', b'embedded copyright'), ('example/ThirdPartyNotices.txt', b'aggregate notices'), ('example/LICENSES/LGPL-3.0-only.txt', b'LGPL terms'), ('example/main.c', b'int main() {}')]:
            info = tarfile.TarInfo(name)
            info.size = len(content)
            archive.addfile(info, io.BytesIO(content))
    checksum = hashlib.sha256(source.read_bytes()).hexdigest()
    (prefix / 'lib/libexample.so.1').write_bytes(b'original shared library')
    (stage / 'lib/libexample.so.1').write_bytes(b'relocated shared library')
    (cached / 'info/licenses/COPYING').write_text('package copyright')
    (cached / 'info/about.json').write_text(json.dumps({'home': 'https://example.invalid'}))
    (cached / 'info/recipe/meta.yaml').write_text(f'source:\n  url: {source.as_uri()}\n  sha256: {checksum}\n')
    record = {'name': 'example', 'version': '1.0', 'build': '0', 'license': 'LGPL-3.0-only',
              'files': ['lib/libexample.so.1'], 'link': {'source': str(cached)},
              'url': 'https://example.invalid/example-1.0-0.conda', 'sha256': 'a' * 64}
    (prefix / 'conda-meta/example.json').write_text(json.dumps(record))
    origins = root / 'origins.tsv'
    origins.write_text(f'lib/libexample.so.1\t{prefix}/lib/libexample.so.1\n')
    command = [sys.executable, str(script), '--stage', str(stage), '--conda-prefix', str(prefix), '--sources', str(output),
               '--source-asset', 'DIALibGen-sources-test.tar.gz', '--origins', str(origins)]
    def run(success, message=''):
        result = subprocess.run(command, text=True, capture_output=True)
        if (result.returncode == 0) != success or message not in result.stdout + result.stderr:
            raise AssertionError(result.stdout + result.stderr)
    run(True)
    inventory = json.loads((output / 'runtime-dependencies.json').read_text())
    component = inventory['components']['example-1.0-0']
    assert component['sources'][0]['sha256'] == checksum
    assert inventory['files'][0]['sha256'] == hashlib.sha256(b'relocated shared library').hexdigest()
    assert inventory['corresponding_source_asset'] == 'DIALibGen-sources-test.tar.gz'
    notices = stage / component['licenses']
    assert (notices / 'COPYING').read_text() == 'package copyright'
    assert (notices / 'source-notices/example/third-party/NOTICE.txt').read_text() == 'embedded copyright'
    assert (notices / 'source-notices/example/ThirdPartyNotices.txt').read_text() == 'aggregate notices'
    assert (notices / 'source-notices/example/LICENSES/LGPL-3.0-only.txt').read_text() == 'LGPL terms'
    assert (output / 'example-1.0-0/recipe/meta.yaml').is_file()
    (stage / 'lib/unknown.so').write_bytes(b'unknown')
    run(False, 'no licensing owner')
    (stage / 'lib/unknown.so').unlink()
    (cached / 'info/recipe/meta.yaml').write_text(f'source:\n  url: {source.as_uri()}\n  sha256: {chr(34)}{"0" * 64}{chr(34)}\n')
    run(False, 'sha256 mismatch')
    (cached / 'info/recipe/meta.yaml').write_text(f'source:\n  url: {source.as_uri()}\n  sha256: {checksum}\n')
    (cached / 'info/licenses/COPYING').unlink()
    run(False, 'no license texts')
print('runtime license inventory: ownership, copied notices, corresponding-source checksums and failure paths passed')
