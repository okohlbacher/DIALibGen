#!/usr/bin/env python3
"""Offline release inventory contracts: exact owners, notices, source hashes, failures."""
import hashlib
import io
import json
from pathlib import Path
import runpy
import subprocess
import sys
import tarfile
import tempfile
from unittest.mock import patch
import zipfile

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
    def notice_text(directory, upstream_path):
        index = json.loads((directory / 'index.json').read_text())
        matches = [filename for filename, entry in index.items() if entry['path'] == upstream_path]
        assert len(matches) == 1, matches
        return (directory / matches[0]).read_text()
    assert notice_text(notices / 'source-notices', 'example/third-party/NOTICE.txt') == 'embedded copyright'
    assert notice_text(notices / 'source-notices', 'example/ThirdPartyNotices.txt') == 'aggregate notices'
    assert notice_text(notices / 'source-notices', 'example/LICENSES/LGPL-3.0-only.txt') == 'LGPL terms'
    assert (output / 'example-1.0-0/recipe/meta.yaml').is_file()
    # Conda's libgomp runtime owns the real .so.1.0.0; _openmp_mutex owns its .so.1 symlink.
    # A recorded loader alias must resolve to the runtime package, not the mutex package.
    if sys.platform != "win32":  # This GNU/Linux alias contract does not require Windows symlink privileges.
        real_library = prefix / 'lib/libexample.so.1.0.0'
        (prefix / 'lib/libexample.so.1').rename(real_library)
        (prefix / 'lib/libexample.so.1').symlink_to(real_library.name)
        record['files'] = ['lib/libexample.so.1.0.0']
        (prefix / 'conda-meta/example.json').write_text(json.dumps(record))
        (prefix / 'conda-meta/mutex.json').write_text(json.dumps({**record, 'name': 'mutex', 'files': ['lib/libexample.so.1']}))
        run(True)
        inventory = json.loads((output / 'runtime-dependencies.json').read_text())
        assert inventory['files'][0]['component'] == 'example-1.0-0'
        (prefix / 'lib/libexample.so.1').unlink()
        real_library.rename(prefix / 'lib/libexample.so.1')
        record['files'] = ['lib/libexample.so.1']
        (prefix / 'conda-meta/example.json').write_text(json.dumps(record))
        (prefix / 'conda-meta/mutex.json').unlink()
    (stage / 'lib/unknown.so').write_bytes(b'unknown')
    run(False, 'no licensing owner')
    failure = json.loads((output / 'runtime-attribution-failure.json').read_text())
    assert failure['file'] == str((stage / 'lib/unknown.so').resolve())
    assert (output / 'runtime-origins.tsv').read_text() == origins.read_text()
    (stage / 'lib/unknown.so').unlink()
    (cached / 'info/recipe/meta.yaml').write_text(f'source:\n  url: {source.as_uri()}\n  sha256: {chr(34)}{"0" * 64}{chr(34)}\n')
    run(False, 'sha256 mismatch')
    (cached / 'info/recipe/meta.yaml').write_text(f'source:\n  url: {source.as_uri()}\n  sha256: {checksum}\n')
    (cached / 'info/licenses/COPYING').unlink()
    run(True)  # Exact upstream notices can repair missing conda license payloads.
    assert notice_text(notices / 'source-notices', 'example/COPYING') == 'upstream license'
    (cached / 'info/recipe/meta.yaml').unlink()
    (cached / 'info/recipe').rmdir()
    run(False, 'no license texts or corresponding-source recipe')
    (cached / 'info/recipe').mkdir()
    (cached / 'info/recipe/meta.yaml').write_text(f'source:\n  url: {source.as_uri()}\n  sha256: {checksum}\n')
    (cached / 'info/licenses/COPYING').write_text('package copyright')

    # Nested runtime payloads must be attributed just like top-level DLLs/SOs.
    nested = stage / 'lib/plugins/libnested.so'
    nested.parent.mkdir()
    nested.write_bytes(b'nested shared library')
    run(False, 'no licensing owner')
    original = prefix / 'lib/plugins/libnested.so'
    original.parent.mkdir()
    original.write_bytes(nested.read_bytes())
    record['files'].append('lib/plugins/libnested.so')
    (prefix / 'conda-meta/example.json').write_text(json.dumps(record))
    with origins.open('a') as stream:
        stream.write(f'lib/plugins/libnested.so\t{original}\n')
    run(True)
    inventory = json.loads((output / 'runtime-dependencies.json').read_text())
    assert any(Path(item['path']).as_posix() == 'lib/plugins/libnested.so' for item in inventory['files'])

    # Lowercase SPDX still requires sources, including notices stored in ZIPs.
    zipped = root / 'source.zip'
    with zipfile.ZipFile(zipped, 'w') as archive:
        archive.writestr('zip-project/LICENSE', 'ZIP primary license')
        archive.writestr('zip-project/vendor/NOTICE.txt', 'ZIP third-party notice')
    zip_checksum = hashlib.sha256(zipped.read_bytes()).hexdigest()
    record['license'] = 'lgpl-3.0-only'
    (prefix / 'conda-meta/example.json').write_text(json.dumps(record))
    (cached / 'info/recipe/meta.yaml').write_text(f'source:\n  url: {zipped.as_uri()}\n  sha256: {zip_checksum}\n')
    run(True)
    inventory = json.loads((output / 'runtime-dependencies.json').read_text())
    assert inventory['components']['example-1.0-0']['sources'][0]['sha256'] == zip_checksum
    assert notice_text(notices / 'source-notices', 'zip-project/LICENSE') == 'ZIP primary license'
    assert notice_text(notices / 'source-notices', 'zip-project/vendor/NOTICE.txt') == 'ZIP third-party notice'

    helpers = runpy.run_path(str(script))
    def refuses(action, diagnostic):
        try:
            action()
        except RuntimeError as error:
            assert diagnostic in str(error), str(error)
        else:
            raise AssertionError(f'accepted invalid input: expected {diagnostic}')

    # A tag cannot substitute for a source checksum; a full commit can.
    for revision in ('main', 'v1.0', 'a' * 39):
        refuses(lambda: helpers['collect_sources']([
            {'git_url': 'https://github.com/example/example.git', 'git_rev': revision}],
            output / 'git', notices), 'unpinned source checkout')
    pinned = helpers['collect_sources']([
        {'git_url': 'https://github.com/example/example.git', 'git_rev': 'a' * 40,
         'local_file': str(source)}], output / 'git', notices)
    assert pinned[0]['sha256'] == checksum
    refuses(lambda: helpers['collect_sources']([
        {'url': source.as_uri(), 'git_rev': 'a' * 40}], output / 'url', notices), 'no source checksum')
    unsafe = root / 'unsafe.zip'
    with zipfile.ZipFile(unsafe, 'w') as archive:
        archive.writestr('../LICENSE', 'must not escape')
    refuses(lambda: helpers['collect_source_notices'](unsafe, root / 'safe'), 'unsafe source notice path')
    long_zip = root / 'long-notices.zip'
    long_paths = [('project/' + 'long-subdirectory/' * 30 + f'{name}/LICENSE', name) for name in ('one', 'two')]
    with zipfile.ZipFile(long_zip, 'w') as archive:
        for path, content in long_paths:
            archive.writestr(path, content)
        archive.writestr('project/REUSE.toml', '# SPDX-FileCopyrightText: upstream author')
    short = root / 'short'
    helpers['collect_source_notices'](long_zip, short)
    for path, content in long_paths:
        assert notice_text(short, path) == content
    assert notice_text(short, 'project/REUSE.toml') == '# SPDX-FileCopyrightText: upstream author'
    assert all(len(file.name) <= 41 and file.parent == short for file in short.rglob('*'))
    flat_index = json.loads((short / 'index.json').read_text())
    assert {entry['archive_sha256'] for entry in flat_index.values()} == {hashlib.sha256(long_zip.read_bytes()).hexdigest()}

    # JSON metadata and empty license placeholders are not license notices.
    sdk_notices = root / 'sdk-notices'
    sdk_notices.mkdir()
    (sdk_notices / 'provenance.json').write_text('{}')
    providers = root / 'providers.json'
    providers.write_text(json.dumps([{'name': 'External', 'version': '1', 'license': 'MIT',
        'files': [], 'force_include': True, 'license_directory': str(sdk_notices)}]))
    command.extend(['--providers', str(providers)])
    run(False, 'missing external SDK notices')
    (sdk_notices / 'LICENSE').write_text('')
    run(False, 'missing external SDK notices')
    (sdk_notices / 'LICENSE').write_text('actual license terms')
    run(True)
    (sdk_notices / 'LICENSE').unlink()
    run(False, 'missing external SDK notices')  # Stale staged notices cannot satisfy the gate.
    (sdk_notices / 'LICENSE').write_text('actual license terms')

    torch = root / 'torch'
    torch_notices = torch / 'share/licenses/Torch'
    torch_notices.mkdir(parents=True)
    (torch_notices / 'provenance.json').write_text(json.dumps({'name': 'PyTorch', 'version': 'test', 'license': 'BSD-3-Clause'}))
    command.extend(['--torch-root', str(torch)])
    run(False, 'Torch SDK lacks a non-empty primary LICENSE')
    (torch_notices / 'LICENSE').write_text('')
    run(False, 'Torch SDK lacks a non-empty primary LICENSE')
    (torch_notices / 'LICENSE').write_text('primary Torch license')
    run(True)

    # The same extraction function used by the downloader rejects a wheel with
    # no primary license, even if an earlier SDK left one in the destination.
    extract = runpy.run_path(str(script.with_name('fetch-libtorch.py')))['extract_sdk']
    wheel = root / 'torch.whl'
    for layout in ('torch-test.dist-info/licenses/', 'torch-test.dist-info/'):
        with zipfile.ZipFile(wheel, 'w') as archive:
            archive.writestr('torch/share/cmake/Torch/TorchConfig.cmake', '# fixture')
            archive.writestr(layout + 'LICENSE', 'primary Torch license')
            archive.writestr(layout + 'NOTICE', 'Torch notice')
            archive.writestr('torch/lib/torch_python.dll', 'exclude the Python bridge')
        extract(wheel, torch)
        assert (torch_notices / 'NOTICE').read_text() == 'Torch notice'
        assert not (torch / 'lib/torch_python.dll').exists()
    for missing in (True, False):
        with zipfile.ZipFile(wheel, 'w') as archive:
            archive.writestr('torch/share/cmake/Torch/TorchConfig.cmake', '# fixture')
            archive.writestr('torch-test.dist-info/other/LICENSE', 'not the primary license')
            if not missing:
                archive.writestr('torch-test.dist-info/LICENSE', '')
        refuses(lambda: extract(wheel, torch), 'wheel lacks a non-empty primary LICENSE')

    intel = root / 'intel-openmp.whl'
    (torch / 'lib').mkdir(exist_ok=True)
    with zipfile.ZipFile(intel, 'w') as archive:
        for name in ('libiomp5md.dll', 'libiompstubs5md.dll'):
            (torch / 'lib' / name).write_bytes(name.encode())
            archive.writestr('intel_openmp-2025.3.1.data/data/Library/bin/' + name, name)
        archive.writestr('intel_openmp-2025.3.1.dist-info/LICENSE.txt', 'complete Intel EULA')
        archive.writestr('intel_openmp-2025.3.1.data/data/share/doc/compiler/licensing/openmp/third-party-programs.txt', 'complete third-party notices')
    intel_provider = runpy.run_path(str(script.with_name('fetch-libtorch.py')))['intel_openmp_provider']
    provider = intel_provider(intel, torch)
    (torch_notices / 'LICENSE').write_text('primary Torch license')
    provenance = {'name': 'PyTorch', 'version': 'test', 'license': 'BSD-3-Clause', 'native_providers': [provider]}
    (torch_notices / 'provenance.json').write_text(json.dumps(provenance))
    (stage / 'bin').mkdir(exist_ok=True)
    for name in ('libiomp5md.dll', 'libiompstubs5md.dll'):
        (stage / 'bin' / name).write_bytes(b'relocated ' + name.encode())
        with origins.open('a') as stream:
            stream.write(f'bin/{name}\t{torch}/lib/{name}\n')
    run(True)
    inventory = json.loads((output / 'runtime-dependencies.json').read_text())
    assert all(item['component'] == 'Intel-OpenMP-2025.3.1' for item in inventory['files'] if item['path'].endswith('.dll'))
    intel_notices = stage / inventory['components']['Intel-OpenMP-2025.3.1']['licenses']
    assert (intel_notices / 'LICENSE.txt').read_text() == 'complete Intel EULA'
    assert (intel_notices / 'third-party-programs.txt').read_text() == 'complete third-party notices'
    (torch / 'lib/libiomp5md.dll').write_bytes(b'wrong DLL')
    refuses(lambda: intel_provider(intel, torch), 'licensing package does not match')
    run(False, 'differs from its licensing record')
    (torch / 'lib/libiomp5md.dll').write_bytes(b'libiomp5md.dll')
    (torch / 'lib/libgfortran.so.5').write_bytes(b'GCC runtime')
    (stage / 'lib/libgfortran.so.5').write_bytes(b'GCC runtime')
    with origins.open('a') as stream:
        stream.write(f'lib/libgfortran.so.5\t{torch}/lib/libgfortran.so.5\n')
    run(False, 'no licensing owner')
    failure = json.loads((output / 'runtime-attribution-failure.json').read_text())
    assert failure['recorded_origin'] == str(torch / 'lib/libgfortran.so.5')
    assert failure['resolved_origin'] == str((torch / 'lib/libgfortran.so.5').resolve())
    (stage / 'lib/libgfortran.so.5').unlink()

    gui = root / 'gui'
    gui.mkdir()
    gui_inventory = gui / 'inventory.json'
    gui_inventory.write_text(json.dumps({'target': 'test', 'packages': [{
        'ecosystem': 'cargo', 'name': 'fixture', 'version': '1', 'license': 'mpl-2.0',
        'source': source.as_uri(), 'checksum': checksum}]}))
    (gui / 'THIRD_PARTY_NOTICES.txt').write_text('fixture notices')
    subprocess.run([sys.executable, str(script.with_name('collect-gui-sources.py')),
                    '--inventory', str(gui_inventory), '--sources', str(output)], check=True)
    assert len(json.loads((output / 'gui/corresponding-sources.json').read_text())) == 1

    # Exercise Windows provider assignment offline: curl was built separately
    # from contrib's curl archive, while Chocolatey's Eigen must match it.
    contrib, openms, provider_cache = [root / name for name in ('contrib', 'openms', 'provider-cache')]
    for directory in (contrib / 'bin', contrib / 'archives', openms, provider_cache):
        directory.mkdir(parents=True)
    for name in ('libcurl.dll', 'zlib.dll'):
        (contrib / 'bin' / name).write_bytes(b'fixture DLL')
    (openms / 'OpenMSConfigVersion.cmake').write_text('set(PACKAGE_VERSION "3.6.0")')
    eigen_header = root / 'EigenMacros.h'
    eigen_header.write_text('#define EIGEN_WORLD_VERSION 3\n#define EIGEN_MAJOR_VERSION 4\n#define EIGEN_MINOR_VERSION 0\n')
    original_rglob = Path.rglob
    def provider_rglob(path, pattern, *args, **kwargs):
        if path == Path('C:/ProgramData/chocolatey/lib/eigen'):
            return iter([eigen_header])
        return original_rglob(path, pattern, *args, **kwargs)
    def provider_manifest(eigen_version):
        lines = []
        for name in ('BZIP2', 'ZLIB', 'BOOST', 'XERCES', 'LIBSVM', 'COINOR', 'EIGEN', 'HDF5', 'ARROW', 'LIBZIP', 'CURL'):
            filename = f'{name.lower()}-{eigen_version if name == "EIGEN" else "1"}.tar.gz'
            (contrib / 'archives' / filename).write_bytes(b'fixture source')
            lines += [f'set(ARCHIVE_{name} "{filename}")', f'set(ARCHIVE_{name}_SHA256 "{"a" * 64}")']
        content = '\n'.join(lines).encode()
        with tarfile.open(provider_cache / 'contrib.tar.gz', 'w:gz') as archive:
            member = tarfile.TarInfo('contrib/CMakeLists.txt')
            member.size = len(content)
            archive.addfile(member, io.BytesIO(content))
    provider_script = script.with_name('windows-license-providers.py')
    provider_args = [str(provider_script), '--contrib', str(contrib), '--contrib-revision', 'a' * 40,
                     '--openms', str(openms), '--openms-revision', 'b' * 40, '--qt', str(root / 'qt'),
                     '--cache', str(provider_cache), '--out', str(providers)]
    with patch.object(sys, 'argv', provider_args), patch.object(Path, 'rglob', provider_rglob):
        provider_manifest('3.4.0')
        runpy.run_path(str(provider_script))
        generated = {provider['name']: provider for provider in json.loads(providers.read_text())}
        assert generated['curl']['version'] == '8.12.1'
        assert generated['curl']['files'] == [str((contrib / 'bin/libcurl.dll').resolve())]
        assert generated['OpenMS-contrib']['files'] == [str((contrib / 'bin/zlib.dll').resolve())]
        assert generated['curl']['source_recipe'] == [{
            'url': 'https://curl.se/download/curl-8.12.1.tar.gz',
            'sha256': '7b40ea64947e0b440716a4d7f0b7aa56230a5341c8377d7b609649d4aea8dbcf'}]
        provider_manifest('3.4.1')
        refuses(lambda: runpy.run_path(str(provider_script)), 'Eigen source archive does not match')
print('runtime inventory: exact ownership, nested libraries, tar/ZIP notices, source pins, Torch licenses and lowercase copyleft passed')
