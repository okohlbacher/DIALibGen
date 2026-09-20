#!/usr/bin/env python3
"""Copy exact runtime package notices and corresponding sources into release artifacts.

Requires PyYAML at packaging time only. The executable has no Python dependency.
Unknown runtime files or unavailable licensing/source evidence fail the package.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tarfile
import tempfile
import urllib.parse
import urllib.request
import zipfile

import yaml


def digest(path, algorithm='sha256'):
    value = hashlib.new(algorithm)
    with path.open('rb') as stream:
        while chunk := stream.read(1024 * 1024):
            value.update(chunk)
    return value.hexdigest()


def has_notices(directory):
    return directory.is_dir() and any(p.is_file() and p.suffix.lower() != '.json' and p.stat().st_size
                                       for p in directory.rglob('*'))


def fetch(urls, destination, checksums=None):
    urls = [urls] if isinstance(urls, str) else urls
    checksums = checksums or {}
    destination.parent.mkdir(parents=True, exist_ok=True)
    if not destination.exists():
        errors = []
        for url in urls:
            try:
                with urllib.request.urlopen(url, timeout=120) as source, destination.open('wb') as output:
                    shutil.copyfileobj(source, output, 1024 * 1024)
                break
            except Exception as error:
                destination.unlink(missing_ok=True)
                errors.append(f'{url}: {error}')
        else:
            raise RuntimeError('source download failed: ' + '; '.join(errors))
    for algorithm, expected in checksums.items():
        if digest(destination, algorithm) != expected:
            destination.unlink()
            raise RuntimeError(f'{algorithm} mismatch: {urls}')
    return {'urls': urls, 'file': destination.name, 'sha256': digest(destination)}


def conda_info(record, cache):
    directory = Path(record.get('link', {}).get('source', ''))
    if (directory / 'info' / 'about.json').is_file():
        return directory / 'info'
    name = f"{record['name']}-{record['version']}-{record['build']}"
    directory = cache / name
    if (directory / 'info' / 'about.json').is_file():
        return directory / 'info'
    url = record.get('url')
    if not url:
        raise RuntimeError(f'no original package URL for {name}')
    archive = cache / Path(urllib.parse.urlparse(url).path).name
    checksums = {k: record[k] for k in ('sha256', 'md5') if record.get(k)}
    if not checksums:
        raise RuntimeError(f'no package checksum for {name}')
    fetch(url, archive, checksums)
    directory.mkdir(parents=True, exist_ok=True)
    if archive.suffix == '.conda':
        with zipfile.ZipFile(archive) as package:
            members = [n for n in package.namelist() if n.startswith('info-') and n.endswith('.tar.zst')]
            if len(members) != 1:
                raise RuntimeError(f'expected one conda info archive: {archive}')
            compressed = cache / Path(members[0]).name
            compressed.write_bytes(package.read(members[0]))
        # bsdtar supports zstd on macOS/Windows; Linux runners provide zstd.
        tar = str(Path(os.environ['SystemRoot']) / 'System32' / 'tar.exe') if os.name == 'nt' else 'tar'
        subprocess.run([tar, '-xf', str(compressed), '-C', str(directory)], check=True)
    else:
        with tarfile.open(archive) as package:
            members = [m for m in package.getmembers() if m.name.startswith('info/')]
            package.extractall(directory, members=members, filter='data')
    if not (directory / 'info' / 'about.json').is_file():
        raise RuntimeError(f'missing package licensing metadata: {name}')
    return directory / 'info'


def source_records(recipe):
    for filename in ('rendered_recipe.yaml', 'meta.yaml', 'recipe.yaml'):
        path = recipe / filename
        if path.is_file():
            data = yaml.safe_load(path.read_text())
            data = data.get('recipe', data)
            sources = data.get('source', [])
            return [sources] if isinstance(sources, dict) else sources
    raise RuntimeError(f'no rendered source recipe: {recipe}')


def collect_source_notices(archive, target):
    """Retain upstream notices for embedded components, without extracting source trees."""
    pattern = re.compile(r'^(COPYING|LICENSE|LICENCE|NOTICE|COPYRIGHT|AUTHORS|Third[_-]?Party[_-]?Notices?)([._-].*)?$', re.I)
    archive_sha256 = digest(archive)
    manifest_file = target / 'index.json'
    manifest = json.loads(manifest_file.read_text()) if manifest_file.is_file() else {}
    def destination(name):
        path = Path(name)
        if not (pattern.match(path.name) or path.name == 'REUSE.toml' or path.suffix.lower() == '.license'
                or any(part.lower() in {'licenses', 'license'} for part in path.parts[:-1])):
            return None
        if path.is_absolute() or '..' in path.parts or '\\' in name or ':' in name:
            raise RuntimeError(f'unsafe source notice path: {name}')
        # Keep installer paths short; the manifest retains the complete upstream name.
        key = hashlib.sha256((archive_sha256 + '\0' + name).encode()).hexdigest()[:16]
        filename = key + '-' + re.sub(r'[^A-Za-z0-9._-]', '_', path.name)[:24]
        record = {'archive': archive.name, 'archive_sha256': archive_sha256, 'path': name}
        if filename in manifest and manifest[filename] != record:
            raise RuntimeError(f'source notice filename collision: {name}')
        manifest[filename] = record
        result = target / filename
        result.parent.mkdir(parents=True, exist_ok=True)
        return result

    if zipfile.is_zipfile(archive):
        with zipfile.ZipFile(archive) as source:
            for member in source.infolist():
                if not member.is_dir() and (output := destination(member.filename)):
                    with source.open(member) as stream, output.open('wb') as sink:
                        shutil.copyfileobj(stream, sink)
    else:
        with tarfile.open(archive) as source:
            for member in source:
                if member.isfile() and (output := destination(member.name)):
                    with source.extractfile(member) as stream, output.open('wb') as sink:
                        shutil.copyfileobj(stream, sink)
    if manifest:
        manifest_file.write_text(json.dumps(manifest, indent=2, sort_keys=True) + '\n')


def collect_sources(records, directory, notices):
    directory.mkdir(parents=True, exist_ok=True)
    result = []
    for source in records:
        urls = source.get('url')
        pinned_git = False
        if not urls and source.get('git_url'):
            git = source['git_url'].removesuffix('.git')
            revision = source.get('git_rev')
            if not git.startswith('https://github.com/') or not isinstance(revision, str) or not re.fullmatch(r'[a-fA-F0-9]{40}', revision):
                raise RuntimeError(f'unsupported unpinned source checkout: {source}')
            urls = f'{git}/archive/{revision}.tar.gz'
            pinned_git = True
        if not urls:
            raise RuntimeError(f'no downloadable source in recipe: {source}')
        urls = [urls] if isinstance(urls, str) else urls
        # Retain the pinned content hash while avoiding unavailable third-party Qt mirrors.
        qt_paths = [urllib.parse.urlparse(url).path for url in urls
                    if urllib.parse.urlparse(url).path.startswith(('/archive/qt/', '/official_releases/qt/'))]
        urls = list(dict.fromkeys(['https://download.qt.io' + path for path in qt_paths] + urls))
        # Debian mirrors the unmodified upstream 1.6.3 archive; preserve the recipe checksum.
        if any(url.endswith('/keyutils.git/snapshot/keyutils-1.6.3.tar.gz') for url in urls):
            urls.insert(0, 'https://deb.debian.org/debian/pool/main/k/keyutils/keyutils_1.6.3.orig.tar.gz')
        name = source.get('fn') or Path(urllib.parse.urlparse(urls[0]).path).name
        checksums = {k: source[k] for k in ('sha256', 'sha512', 'sha1', 'md5') if k in source}
        if any(not isinstance(value, str) or not re.fullmatch('[a-fA-F0-9]+', value) for value in checksums.values()):
            raise RuntimeError(f'invalid source checksum in recipe: {source}')
        if not checksums and not pinned_git:
            raise RuntimeError(f'no source checksum in recipe: {source}')
        identity = next(iter(checksums.values()), source.get('git_rev'))
        destination = directory.parent / 'archives' / f'{identity}-{name}'
        if source.get('local_file'):
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source['local_file'], destination)
        item = fetch(urls, destination, checksums)
        item['file'] = '../archives/' + destination.name
        item['patches'] = source.get('patches', [])
        # Qt and other libraries embed third parties whose notices are in their sources.
        archive = directory / item['file']
        collect_source_notices(archive, notices / 'source-notices')
        result.append(item)
    if not result:
        raise RuntimeError(f'no corresponding source records for {directory.name}')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage', type=Path, required=True)
    parser.add_argument('--conda-prefix', type=Path, required=True)
    parser.add_argument('--sources', type=Path, required=True)
    parser.add_argument('--source-asset', required=True)
    parser.add_argument('--torch-root', type=Path)
    parser.add_argument('--origins', type=Path, help='TSV of staged relative path and original absolute library path')
    parser.add_argument('--providers', type=Path, help='Additional SDK/build providers with exact files, licenses and source records')
    args = parser.parse_args()
    stage, prefix = args.stage.resolve(), args.conda_prefix.resolve()
    sources = args.sources.resolve()
    sources.mkdir(parents=True, exist_ok=True)
    notices = stage / 'share/DIALibGen/licenses/runtime'
    notices.mkdir(parents=True, exist_ok=True)
    repository_licenses = Path(__file__).resolve().parent.parent / 'licenses'
    records, owners = {}, {}
    for path in sorted((prefix / 'conda-meta').glob('*.json')):
        record = json.loads(path.read_text())
        key = f"{record['name']}-{record['version']}-{record['build']}"
        records[key] = record
        for filename in record.get('files', []):
            file = prefix / filename
            owners.setdefault(file.name.lower(), []).append((file, key))
    if not records:
        raise RuntimeError(f'no conda package inventory in {prefix}')
    providers = json.loads(args.providers.read_text()) if args.providers else []
    if args.torch_root and args.torch_root.resolve() != prefix:
        torch = args.torch_root.resolve()
        provenance = torch / 'share/licenses/Torch/provenance.json'
        if not provenance.is_file():
            raise RuntimeError('Torch SDK lacks provenance; refresh its cache using fetch-libtorch.py')
        license_file = provenance.with_name('LICENSE')
        if not license_file.is_file() or not license_file.stat().st_size:
            raise RuntimeError('Torch SDK lacks a non-empty primary LICENSE; refresh its cache using fetch-libtorch.py')
        data = json.loads(provenance.read_text())
        native = data.pop('native_providers', [])
        native_files = set()
        for provider in native:
            for filename, expected in provider['sdk_file_sha256'].items():
                if digest(torch / filename) != expected:
                    raise RuntimeError(f'Torch native dependency differs from its licensing record: {filename}')
            provider['files'] = [str((torch / filename).resolve()) for filename in provider['sdk_file_sha256']]
            native_files.update(provider['files'])
            if provider.get('license_directory'):
                provider['license_directory'] = str(torch / provider['license_directory'])
        # These are separate upstream components, not covered by PyTorch's BSD grant.
        separate = re.compile(r'^(libarm_compute|libgomp|libgfortran|libopenblas|libiomp)', re.I)
        providers.append({**data, 'files': [str(p) for p in (torch / 'lib').rglob('*')
                          if p.is_file() and str(p.resolve()) not in native_files and not separate.match(p.name)],
                          'license_directory': str(torch / 'share/licenses/Torch')})
        providers.extend(native)
    origins = dict(line.split('\t', 1) for line in args.origins.read_text().splitlines()) if args.origins else {}
    external = {}
    for provider in providers:
        for filename in provider['files']:
            external[str(Path(filename).resolve()).lower()] = provider
    bundled = [p for folder in ('lib', 'bin') for p in (stage / folder).rglob('*')
               if p.is_file() and ('.so' in p.name or p.suffix.lower() in {'.dylib', '.dll'})]
    components, files = {}, []
    with tempfile.TemporaryDirectory(prefix='dialibgen-license-packages-') as temporary:
        cache = Path(temporary)
        def conda_component(key):
            if key in components:
                return key
            record = records[key]
            info = conda_info(record, cache)
            about = json.loads((info / 'about.json').read_text())
            license_dir = info / 'licenses'
            target = notices / key
            if target.exists():
                shutil.rmtree(target)
            if has_notices(license_dir):
                shutil.copytree(license_dir, target, dirs_exist_ok=True)
            elif 'GCC-exception' in record.get('license', ''):
                target.mkdir(parents=True, exist_ok=True)
                for name in ('GPL-3.0-or-later.txt', 'GCC-Runtime-Exception-3.1.txt'):
                    shutil.copy2(repository_licenses / name, target / name)
            recipe = info / 'recipe'
            component = {'name': record['name'], 'version': record['version'], 'build': record['build'],
                         'license': record.get('license', about.get('license')),
                         'package_url': record.get('url'), 'package_sha256': record.get('sha256'),
                         'homepage': about.get('home'), 'source_repository': about.get('dev_url'),
                         'licenses': str(target.relative_to(stage)), 'source_asset': None}
            if recipe.is_dir():
                shutil.copytree(recipe, sources / key / 'recipe', dirs_exist_ok=True)
                component['source_recipe'] = source_records(recipe)
            aggregate = record['name'] in {'onnxruntime-cpp', 'libopenms', 'libparquet'} or record['name'].startswith('libarrow')
            if aggregate or re.search(r'GPL|EPL|MPL|CDDL', component['license'] or '', re.I) or not has_notices(target):
                if not recipe.is_dir():
                    raise RuntimeError(f'no license texts or corresponding-source recipe for {key}')
                selected = component['source_recipe']
                if record['name'] == 'qt6-main':
                    # Qt's cross-platform recipe also downloads a Windows software
                    # renderer. It is not a source archive or part of this CLI.
                    selected = []
                    for source in component['source_recipe']:
                        urls = source.get('url', [])
                        urls = [urls] if isinstance(urls, str) else urls
                        if any('/prebuilt/llvmpipe/windows/' in url for url in urls):
                            if any(file.name.lower() == 'opengl32sw.dll' for file in bundled):
                                raise RuntimeError('bundled Qt software renderer needs its own Mesa source provider')
                            component.setdefault('excluded_build_inputs', []).append({
                                'urls': urls, 'reason': 'Prebuilt Windows software renderer is absent from this bundle.'})
                        else:
                            selected.append(source)
                component['sources'] = collect_sources(selected, sources / key, target)
                component['source_asset'] = args.source_asset
            if not has_notices(target):
                raise RuntimeError(f'no license texts in original package or its pinned sources: {key}')
            components[key] = component
            return key

        def external_component(provider):
            key = provider['name'] + '-' + provider['version']
            if key in components:
                return key
            target = notices / key
            if target.exists():
                shutil.rmtree(target)
            license_dir = Path(provider['license_directory']) if provider.get('license_directory') else None
            if license_dir and license_dir.is_dir():
                shutil.copytree(license_dir, target, dirs_exist_ok=True)
            component = {k: v for k, v in provider.items() if k not in {'files', 'license_directory', 'force_include', 'collect_sources'}}
            component['licenses'] = str(target.relative_to(stage))
            if provider.get('collect_sources') or re.search(r'GPL|EPL|MPL|CDDL', component['license'], re.I):
                component['sources'] = collect_sources(component.get('source_recipe', []), sources / key, target)
                component['source_asset'] = args.source_asset
            if not has_notices(target):
                raise RuntimeError(f'missing external SDK notices: {key}')
            if 'source_recipe' in component:
                component['source_recipe'] = [{k: v for k, v in record.items() if k != 'local_file'}
                                              for record in component['source_recipe']]
            components[key] = component
            return key

        for file in sorted(bundled):
            matches = [(path, key) for path, key in owners.get(file.name.lower(), []) if path.is_file()]
            original = origins.get(file.relative_to(stage).as_posix())
            if original:
                original = str(Path(original).resolve()).lower()
                matches = [(path, key) for path, key in owners.get(Path(original).name.lower(), [])
                           if not path.is_symlink() and str(path.resolve()).lower() == original]
                provider = external.get(original)
            else:
                candidates = [provider for path, provider in external.items() if Path(path).name.lower() == file.name.lower()]
                if candidates and matches:
                    raise RuntimeError(f'ambiguous SDK/package origin for {file.name}; record --origins while staging')
                provider = candidates[0] if len(candidates) == 1 else None
            if provider:
                owner = external_component(provider)
            elif matches:
                keys = {key for _, key in matches}
                if len(keys) != 1:
                    raise RuntimeError(f'ambiguous package owner for {file.name}: {sorted(keys)}')
                owner = conda_component(next(iter(keys)))
            else:
                raise RuntimeError(f'no licensing owner for bundled runtime {file}; provide exact SDK metadata')
            files.append({'path': str(file.relative_to(stage)), 'sha256': digest(file), 'component': owner})
        for provider in providers:
            if provider.get('force_include'):
                external_component(provider)
        # Header-only dependencies are compiled into the program, not visible to the loader.
        for key, record in records.items():
            if record['name'] in {'nlohmann_json', 'libboost-headers', 'eigen'}:
                conda_component(key)
    inventory = {'format': 1, 'corresponding_source_asset': args.source_asset,
                 'files': files, 'components': components}
    content = json.dumps(inventory, indent=2, sort_keys=True) + '\n'
    (stage / 'share/DIALibGen/runtime-dependencies.json').write_text(content)
    (sources / 'runtime-dependencies.json').write_text(content)
    print(f'Collected {len(components)} runtime/header components and {len(files)} shared libraries')


if __name__ == '__main__':
    main()
