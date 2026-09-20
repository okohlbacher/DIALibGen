#!/usr/bin/env python3
"""Inventory the actual AppImage payload and archive its Ubuntu source packages.

Run on the Ubuntu packaging host before its dpkg database disappears. Match ELF
build IDs because linuxdeploy adjusts RPATHs; never infer ownership from a name
alone. Every payload file needs an owner or explicit generated/helper provenance.
The companion source archive is bound to the exact final AppImage SHA-256.
"""
import argparse
from collections import defaultdict
from email.parser import Parser
from functools import lru_cache
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


def run(*args):
    return subprocess.check_output(args, text=True, stderr=subprocess.PIPE, timeout=60).strip()


@lru_cache(None)
def digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


@lru_cache(None)
def is_elf(path):
    with path.open('rb') as stream:
        return stream.read(4) == b'\x7fELF'


@lru_cache(None)
def build_id(path):
    if not is_elf(path):
        return None
    result = subprocess.run(['readelf', '-n', str(path)], text=True, capture_output=True)
    match = re.search(r'Build ID: ([a-fA-F0-9]+)', result.stdout)
    return match.group(1).lower() if match else None


def same_file(first, second):
    first, second = first.resolve(), second.resolve()
    return digest(first) == digest(second) or bool(
        is_elf(first) and is_elf(second) and build_id(first) and build_id(first) == build_id(second))


def fetch(url, destination, expected=None):
    destination.parent.mkdir(parents=True, exist_ok=True)
    if not destination.exists():
        request = urllib.request.Request(url, headers={'User-Agent': 'DIALibGen-release-attribution'})
        try:
            with urllib.request.urlopen(request, timeout=120) as response, destination.open('wb') as output:
                shutil.copyfileobj(response, output, 1024 * 1024)
        except Exception:
            destination.unlink(missing_ok=True)
            raise
    actual = digest(destination)
    if expected and actual != expected:
        destination.unlink()
        digest.cache_clear()
        raise RuntimeError(f'source checksum mismatch: {url}')
    return {'url': url, 'file': destination.name, 'sha256': actual}


def get_json(url):
    with urllib.request.urlopen(url, timeout=60) as response:
        return json.load(response)


def dsc_fields(text):
    # Parse the RFC822 payload, not a verified signature. Descriptor authenticity
    # uses HTTPS to the authoritative Launchpad publication (see ubuntu_sources).
    if text.startswith('-----BEGIN PGP SIGNED MESSAGE-----'):
        text = text.split('\n\n', 1)[1].split('\n-----BEGIN PGP SIGNATURE-----', 1)[0]
        text = '\n'.join(line[2:] if line.startswith('- ') else line for line in text.splitlines())
    return Parser().parsestr(text)


def ubuntu_sources(name, version, directory):
    """Launchpad retains exact superseded Ubuntu sources absent from live apt indices."""
    directory.mkdir(parents=True, exist_ok=True)
    query = urllib.parse.urlencode({'ws.op': 'getPublishedSources', 'source_name': name,
                                   'version': version, 'exact_match': 'true'})
    listing = get_json('https://api.launchpad.net/1.0/ubuntu/+archive/primary?' + query)
    entries = [entry for entry in listing.get('entries', [])
               if entry['source_package_name'] == name and entry['source_package_version'] == version]
    if not entries:
        raise RuntimeError(f'Ubuntu has no authoritative source publication for {name}={version}')
    errors = []
    for entry in entries:
        try:
            urls = get_json(entry['self_link'] + '?ws.op=sourceFileUrls')
            by_name = {Path(urllib.parse.unquote(urllib.parse.urlparse(url).path)).name: url for url in urls}
            dsc_names = [filename for filename in by_name if filename.endswith('.dsc')]
            if len(dsc_names) != 1:
                raise RuntimeError('expected one source descriptor')
            dsc_name = dsc_names[0]
            records = [fetch(by_name[dsc_name], directory / dsc_name)]
            fields = dsc_fields((directory / dsc_name).read_text())
            if fields.get('Source') != name or fields.get('Version') != version:
                raise RuntimeError('source descriptor identity differs from installed package')
            checksums = fields.get('Checksums-Sha256', '').splitlines()
            if not any(line.strip() for line in checksums):
                raise RuntimeError('source descriptor lacks SHA-256 checksums')
            for line in checksums:
                if not line.strip():
                    continue
                checksum, size, filename = line.split()
                if Path(filename).name != filename or not re.fullmatch(r'[0-9a-f]{64}', checksum):
                    raise RuntimeError('unsafe source descriptor entry')
                if filename not in by_name:
                    raise RuntimeError(f'missing corresponding source member {filename}')
                records.append(fetch(by_name[filename], directory / filename, checksum))
                if (directory / filename).stat().st_size != int(size):
                    raise RuntimeError(f'source size mismatch: {filename}')
            return {'publication': entry['self_link'], 'files': records,
                    'descriptor_trust': 'HTTPS to Launchpad; OpenPGP signature not verified. Source members verified against descriptor SHA-256.'}
        except Exception as error:
            errors.append(str(error))
    raise RuntimeError(f'cannot retrieve exact source {name}={version}: {errors}')


def installed_packages():
    fields = '${binary:Package}\t${Version}\t${source:Package}\t${source:Version}\t${db:Status-Abbrev}\n'
    packages, basenames = {}, defaultdict(list)
    for line in run('dpkg-query', '-W', '-f=' + fields).splitlines():
        package, version, source, source_version, status = line.split('\t')
        if not status.startswith('ii'):
            continue
        packages[package] = {'name': package, 'version': version, 'source': source, 'source_version': source_version}
        file_list = Path('/var/lib/dpkg/info') / (package + '.list')
        if not file_list.exists():
            file_list = file_list.with_name(package.split(':')[0] + '.list')
        if not file_list.exists():
            raise RuntimeError(f'no installed file list for {package}')
        for filename in file_list.read_text().splitlines():
            path = Path(filename)
            basenames[path.name].append((path, package))
    return packages, basenames


def matching_owners(path, candidates):
    return [(original, owner) for original, owner in candidates
            if original.is_file() and same_file(path, original)]


def gui_asset_hashes(gui):
    # Explicit source/build outputs only; node_modules and arbitrary build files
    # do not become first-party merely by living below the GUI directory.
    assets = {}
    for folder in ('src', 'dist', 'src-tauri/icons', 'src-tauri/resources/third-party-licenses'):
        for path in (gui / folder).rglob('*'):
            if path.is_file():
                assets[digest(path)] = str(path.relative_to(gui))
    return assets


def deb_descriptor_hashes(deb, product_name):
    # These Tauri-generated descriptors must match the same-build .deb exactly.
    # Other package files are still attributed individually, never waived.
    allowed = {f'usr/share/applications/{product_name}.desktop',
               f'usr/share/metainfo/{product_name}.metainfo.xml'}
    result = {}
    with tempfile.TemporaryFile() as errors:
        with subprocess.Popen(['dpkg-deb', '--fsys-tarfile', str(deb)], stdout=subprocess.PIPE, stderr=errors) as process:
            with tarfile.open(fileobj=process.stdout, mode='r|') as archive:
                for member in archive:
                    name = member.name.removeprefix('./')
                    if member.isfile() and name in allowed:
                        with archive.extractfile(member) as content:
                            result[hashlib.file_digest(content, 'sha256').hexdigest()] = name
            if process.wait(timeout=60):
                errors.seek(0)
                raise RuntimeError('cannot read same-build Debian descriptors: ' + errors.read().decode())
    if not result:
        raise RuntimeError('same-build Debian package has no expected GUI desktop descriptor')
    return result


def generated_cache_owner(path, appdir, candidates):
    relative = path.relative_to(appdir).as_posix()
    tools = {
        'usr/share/glib-2.0/schemas/gschemas.compiled': 'glib-compile-schemas',
        'usr/lib/gtk-3.0/3.0.0/immodules.cache': 'gtk-query-immodules-3.0',
        'usr/lib/gdk-pixbuf-2.0/2.10.0/loaders.cache': 'gdk-pixbuf-query-loaders',
    }
    tool_name = tools.get(relative)
    if not tool_name:
        return None
    for tool, owner in candidates[tool_name]:
        if not tool.is_file() or not os.access(tool, os.X_OK):
            continue
        if tool_name == 'glib-compile-schemas':
            with tempfile.TemporaryDirectory(prefix='dialibgen-schema-proof-') as temporary:
                directory = Path(temporary)
                for source in path.parent.iterdir():
                    if source.name.endswith(('.gschema.xml', '.gschema.override')):
                        shutil.copy2(source, directory / source.name)
                subprocess.run([str(tool), str(directory)], check=True, capture_output=True, timeout=60)
                generated = (directory / 'gschemas.compiled').read_bytes()
        else:
            generated = subprocess.check_output([str(tool)], timeout=60)
            if tool_name == 'gtk-query-immodules-3.0':
                prefix = run('/usr/bin/pkg-config', '--variable=libdir', 'gtk+-3.0') + '/gtk-3.0/3.0.0/immodules/'
            else:
                prefix = run('/usr/bin/pkg-config', '--variable=gdk_pixbuf_moduledir', 'gdk-pixbuf-2.0') + '/'
            # Exact transformation performed by the pinned GTK deploy plugin.
            generated = generated.replace(prefix.encode(), b'')
        if hashlib.sha256(generated).hexdigest() == digest(path):
            return owner, str(tool)
    return None


def payload_manifest(directory):
    if not directory.is_dir():
        raise RuntimeError(f'AppImage payload directory does not exist: {directory}')
    result = {}
    for path in directory.rglob('*'):
        relative = path.relative_to(directory).as_posix()
        if path.is_symlink():
            result[relative] = ('symlink', os.readlink(path))
        elif path.is_file():
            result[relative] = ('file', digest(path))
        elif not path.is_dir():
            raise RuntimeError(f'unsupported AppImage payload entry: {relative}')
    return result


def validate_appimage_payload(appdir, image):
    expected = payload_manifest(appdir)
    with tempfile.TemporaryDirectory(prefix='dialibgen-appimage-proof-') as temporary:
        subprocess.run([str(image), '--appimage-extract'], cwd=temporary, check=True,
                       stdout=subprocess.DEVNULL, timeout=180)
        actual = payload_manifest(Path(temporary) / 'squashfs-root')
    if actual != expected:
        missing = sorted(expected.keys() - actual.keys())
        extra = sorted(actual.keys() - expected.keys())
        changed = sorted(path for path in expected.keys() & actual.keys() if expected[path] != actual[path])
        raise RuntimeError(f'AppImage payload differs from AppDir: missing={missing}, extra={extra}, changed={changed}')


def validate_runtime_provider(image, runtime):
    actual_offset = int(run(str(image), '--appimage-offset'))
    if runtime['size'] != actual_offset or not 0 < actual_offset < image.stat().st_size:
        raise RuntimeError('AppImage runtime offset differs from provider evidence')
    with image.open('rb') as stream:
        actual = hashlib.sha256(stream.read(actual_offset)).hexdigest()
    if actual != runtime['sha256']:
        raise RuntimeError('AppImage runtime hash differs from provider evidence')


def required_notices(used_packages, providers):
    notices = {f'packages/{name}/copyright': Path('/usr/share/doc') / name.split(':')[0] / 'copyright'
               for name in sorted(used_packages)}
    for provider in providers:
        if not provider.get('license') or not provider.get('notices'):
            raise RuntimeError(f'incomplete helper licensing: {provider["name"]}')
        for index, notice in enumerate(provider['notices']):
            notices[f'providers/{provider["name"]}/{index}-{Path(notice).name}'] = Path(notice)
    # Installed copyright files refer to this complete, small set of license texts.
    notices.update({f'common-licenses/{path.name}': path
                    for path in sorted(Path('/usr/share/common-licenses').iterdir()) if path.is_file()})
    for name, path in notices.items():
        if not path.is_file() or not path.stat().st_size:
            raise RuntimeError(f'missing or empty AppImage notice: {name}')
    return notices


def stage_notices(notices, directory):
    if directory.exists():
        shutil.rmtree(directory)
    for name, source in notices.items():
        destination = directory / name
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)


def verify_bundled_notices(notices, records):
    by_hash = defaultdict(list)
    for record in records:
        if record.get('sha256'):
            by_hash[record['sha256']].append(record['path'])
    bundled = {name: by_hash.get(digest(source), []) for name, source in notices.items()}
    missing = [name for name, paths in bundled.items() if not paths]
    if missing:
        raise RuntimeError('final AppImage lacks required notice bytes: ' + ', '.join(missing))
    return bundled


def write_failure_diagnostics(output, unknown, appdir, baseline, candidates):
    originals = set()

    def describe(path):
        result = {'path': str(path), 'resolved_path': str(path.resolve())}
        if path.is_file():
            result.update(sha256=digest(path), build_id=build_id(path))
            if is_elf(path):
                for flag, key in (('-n', 'elf_notes'), ('-d', 'elf_dynamic')):
                    inspected = subprocess.run(['readelf', flag, str(path)], text=True, capture_output=True)
                    result[key] = {'exit': inspected.returncode, 'stdout': inspected.stdout, 'stderr': inspected.stderr}
        return result

    records = []
    for record in unknown:
        path = appdir / record['path']
        matched = []
        for original, owner in [*((p, 'DIALibGen CLI baseline') for p in baseline[path.name]),
                                *candidates[path.name]]:
            matched.append({**describe(original), 'owner': owner})
            if original.is_file():
                originals.add(str(original.resolve()))
        records.append({**describe(path), 'candidates': matched})
    report = {'complete': False, 'unknown': records,
              'tools': {name: shutil.which(name) for name in ('readelf', 'patchelf', 'strip', 'glib-compile-schemas')},
              'environment': {name: os.environ.get(name) for name in ('PATH', 'LD_LIBRARY_PATH', 'PKG_CONFIG_PATH')}}
    (output / 'attribution-failure.json').write_text(json.dumps(report, indent=2) + '\n')
    (output / 'candidate-paths.txt').write_text(''.join(path + '\n' for path in sorted(originals)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--appdir', type=Path, required=True)
    parser.add_argument('--appimage', type=Path, required=True)
    parser.add_argument('--baseline', type=Path, required=True)
    parser.add_argument('--gui-binary', type=Path, required=True)
    parser.add_argument('--gui-root', type=Path, default=Path(__file__).resolve().parent.parent / 'gui')
    parser.add_argument('--deb', type=Path, required=True, help='same-build Debian package for generated descriptor proof')
    parser.add_argument('--sources', type=Path, required=True)
    parser.add_argument('--source-asset', required=True)
    parser.add_argument('--providers', type=Path)
    parser.add_argument('--prepare-notices', action='store_true',
                        help='stage exact payload notices in GUI resources for rebundling; do not collect sources')
    args = parser.parse_args()
    appdir, image = args.appdir.resolve(), args.appimage.resolve()
    validate_appimage_payload(appdir, image)
    output = args.sources.resolve() / 'appimage'
    output.mkdir(parents=True, exist_ok=True)
    providers = json.loads(args.providers.read_text()) if args.providers else []
    packages, candidates = installed_packages()
    baseline = defaultdict(list)
    for path in args.baseline.rglob('*'):
        if path.is_file():
            baseline[path.name].append(path)
    gui_assets = gui_asset_hashes(args.gui_root)
    config = json.loads((args.gui_root / 'src-tauri/tauri.conf.json').read_text())
    descriptors = deb_descriptor_hashes(args.deb, config['productName'])
    records, unknown, used_packages, used_providers = [], [], set(), set()
    for path in sorted(appdir.rglob('*')):
        rel = str(path.relative_to(appdir))
        if path.is_symlink():
            records.append({'path': rel, 'type': 'symlink', 'target': os.readlink(path)})
            continue
        if not path.is_file():
            continue
        record = {'path': rel, 'sha256': digest(path), 'build_id': build_id(path)}
        if same_file(path, args.gui_binary):
            record['owner'] = 'DIALibGen GUI (Cargo/npm inventory in resources/third-party-licenses)'
        elif any(same_file(path, original) for original in baseline[path.name]):
            record['owner'] = 'DIALibGen CLI (runtime-dependencies.json)'
        elif record['sha256'] in gui_assets:
            record.update({'owner': 'DIALibGen GUI asset', 'original': gui_assets[record['sha256']]})
        elif record['sha256'] in descriptors:
            record.update({'owner': 'DIALibGen Tauri-generated descriptor', 'original': descriptors[record['sha256']],
                           'same_build_deb_sha256': digest(args.deb)})
        else:
            matches = matching_owners(path, candidates[path.name])
            if matches:
                original, owner = matches[0]
                record.update({'owner': owner, 'original': str(original), 'original_sha256': digest(original.resolve())})
                used_packages.add(owner)
            else:
                generated = generated_cache_owner(path, appdir, candidates)
                if generated:
                    owner, tool = generated
                    record.update({'owner': owner, 'generated_by': tool, 'proof': 'identical regenerated SHA-256'})
                    used_packages.add(owner)
                    records.append(record)
                    continue
                matched = [p for p in providers if record['sha256'] in p.get('sha256', [])
                           or (record['build_id'] and record['build_id'] in p.get('build_ids', []))]
                if len(matched) == 1:
                    record['owner'] = matched[0]['name']
                    used_providers.add(matched[0]['name'])
                else:
                    unknown.append(record)
        records.append(record)
    if unknown:
        write_failure_diagnostics(output, unknown, appdir, baseline, candidates)
        raise RuntimeError('unattributed AppImage files: ' + ', '.join(item['path'] for item in unknown))
    providers = [provider for provider in providers
                 if provider['name'] in used_providers or provider.get('appimage_runtime')]
    if not any(provider.get('appimage_runtime') for provider in providers):
        raise RuntimeError('AppImage header runtime requires a corresponding-source helper provider')
    for provider in providers:
        if provider.get('appimage_runtime'):
            validate_runtime_provider(image, provider['appimage_runtime'])
    notices = required_notices(used_packages, providers)
    if args.prepare_notices:
        directory = args.gui_root / 'src-tauri/resources/third-party-licenses/appimage'
        stage_notices(notices, directory)
        print(f'Staged {len(notices)} exact AppImage notices in {directory}; rebundle before final verification')
        return
    bundled_notices = verify_bundled_notices(notices, records)
    inventory = {'appimage': image.name, 'sha256': digest(image), 'source_asset': args.source_asset,
                 'scope': 'Every regular file and symlink in the AppDir, verified against the extracted final AppImage by path, exact file SHA-256 and link target. Attribution of ELF RPATH changes may match original build IDs; data requires exact SHA-256. AppImage header runtime has separate provider proof.',
                 'files': records, 'packages': [], 'providers': [], 'unknown': unknown, 'complete': False,
                 'bundled_notice_paths': bundled_notices}
    report = output / 'inventory.json'
    report.write_text(json.dumps(inventory, indent=2) + '\n')
    # Ubuntu copyright files may refer to shared full license texts. Keep the
    # exact installed collection once; it is small and avoids fragile parsing.
    shutil.copytree('/usr/share/common-licenses', output / 'common-licenses', dirs_exist_ok=True)
    source_cache = {}
    for name in sorted(used_packages):
        package = packages[name]
        target = output / 'packages' / name
        target.mkdir(parents=True, exist_ok=True)
        copyright_file = Path('/usr/share/doc') / name.split(':')[0] / 'copyright'
        if not copyright_file.is_file():
            raise RuntimeError(f'no installed copyright notice for {name}')
        shutil.copy2(copyright_file, target / 'copyright')
        key = package['source'] + '=' + package['source_version']
        if key not in source_cache:
            source_cache[key] = ubuntu_sources(package['source'], package['source_version'], output / 'sources' / key)
        inventory['packages'].append({**package, 'copyright': str((target / 'copyright').relative_to(output)),
                                      'copyright_sha256': digest(target / 'copyright'),
                                      'bundled_copyright_paths': bundled_notices[f'packages/{name}/copyright'],
                                      'corresponding_source': source_cache[key]})
    for provider in providers:
        target = output / 'providers' / provider['name']
        target.mkdir(parents=True, exist_ok=True)
        if not provider.get('license') or not provider.get('notices'):
            raise RuntimeError(f'incomplete helper licensing: {provider["name"]}')
        notices = []
        for index, notice in enumerate(provider['notices']):
            dest = target / f'{index}-{Path(notice).name}'
            shutil.copy2(notice, dest)
            notices.append({'file': str(dest.relative_to(output)), 'sha256': digest(dest),
                            'bundled_paths': bundled_notices[f'providers/{provider["name"]}/{index}-{Path(notice).name}']})
        if re.search(r'GPL|MPL|EPL|CDDL', provider['license'], re.I) and not provider.get('sources'):
            raise RuntimeError(f'no corresponding helper source: {provider["name"]}')
        sources = [fetch(source['url'], target / source['filename'], source['sha256'])
                   for source in provider.get('sources', [])]
        inventory['providers'].append({**provider, 'notices': notices, 'sources': sources})
    (output / 'README.txt').write_text(
        'This directory accompanies the AppImage named and SHA-256 identified in inventory.json.\n'
        'Ubuntu packages retain their original copyright notices and exact source versions.\n'
        'bundled_notice_paths identifies byte-identical package copyrights, helper notices and\n'
        'common license texts inside the extracted final AppImage; absent notices fail packaging.\n'
        'Source descriptors are trusted through HTTPS to their authoritative Launchpad publication;\n'
        'their OpenPGP signatures are not verified. Source archives match descriptor SHA-256 checksums.\n'
        'sources/ contains the complete upstream archives and Debian/Ubuntu packaging/patches.\n'
        'Rebuild on matching Ubuntu: install build dependencies, run dpkg-source -x PACKAGE.dsc,\n'
        'then dpkg-buildpackage -us -uc in the extracted source directory.\n'
        'Shared-library changes made by linuxdeploy are recorded by original and bundled hashes;\n'
        'the corresponding source packages remain unmodified.\n'
        'The GUI build and AppImage commands are in the DIALibGen release source .github/workflows/ci.yml.\n')
    inventory['complete'] = True
    report.write_text(json.dumps(inventory, indent=2) + '\n')
    print(f'AppImage attribution: {len(records)} payload files/links, {len(used_packages)} Ubuntu packages, '
          f'{len(source_cache)} exact source packages; SHA256 {inventory["sha256"]}')


if __name__ == '__main__':
    main()
