#!/usr/bin/env python3
"""Inventory the pinned Windows SDKs and static OpenMS contrib inputs."""
import argparse
import json
from pathlib import Path
import re
import runpy
import tarfile

helpers = runpy.run_path(str(Path(__file__).with_name('collect-runtime-licenses.py')))
fetch = helpers['fetch']

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--contrib', type=Path, required=True)
parser.add_argument('--contrib-revision', required=True)
parser.add_argument('--openms', type=Path, required=True)
parser.add_argument('--openms-revision', required=True)
parser.add_argument('--qt', type=Path, required=True)
parser.add_argument('--cache', type=Path, required=True)
parser.add_argument('--out', type=Path, required=True)
args = parser.parse_args()
for revision in (args.contrib_revision, args.openms_revision):
    if not re.fullmatch('[0-9a-f]{40}', revision):
        raise RuntimeError('OpenMS and contrib source revisions must be full commit IDs')
args.cache.mkdir(parents=True, exist_ok=True)

# The contrib archive includes every build script and patch used by its cached build.
url = f'https://github.com/OpenMS/contrib/archive/{args.contrib_revision}.tar.gz'
archive = args.cache / 'contrib.tar.gz'
record = fetch(url, archive)
with tarfile.open(archive) as source:
    cmake = next(m for m in source if m.name.endswith('/CMakeLists.txt') and m.name.count('/') == 1)
    definitions = source.extractfile(cmake).read().decode()
# These inputs are the ALL build in the pinned contrib. GLPK is an explicit opt-in.
required = {'BZIP2', 'ZLIB', 'BOOST', 'XERCES', 'LIBSVM', 'COINOR', 'EIGEN', 'HDF5', 'ARROW', 'LIBZIP', 'CURL'}
values = dict(re.findall(r'^set\(ARCHIVE_([A-Z0-9_]+)\s+"?([^"\s)]+)"?\s*\)', definitions, re.M))
recipes = [{'url': url, 'sha256': record['sha256'], 'local_file': str(archive.resolve())}]
for name in sorted(required):
    filename, checksum = values.get(name), values.get(name + '_SHA256')
    if not filename or not checksum:
        raise RuntimeError(f'contrib source manifest lacks {name}')
    cached = args.contrib / 'archives' / filename
    if not cached.is_file():
        raise RuntimeError(f'contrib cache omitted original source archive: {cached}; rebuild its cache')
    recipes.append({'url': 'https://github.com/OpenMS/contrib-sources/releases/download/3.6.0/' + filename,
                    'sha256': checksum, 'local_file': str(cached.resolve())})
if not re.search(r'(?<![0-9])3\.4\.0(?![0-9])', values['EIGEN']):
    raise RuntimeError(f'contrib Eigen source archive does not match Chocolatey 3.4.0: {values["EIGEN"]}')

# OpenMS finds Chocolatey's Eigen headers rather than contrib's copy. Refuse an
# unrecorded version change; 3.4.0 corresponds to the source archive above.
eigen = Path('C:/ProgramData/chocolatey/lib/eigen')
headers = list(eigen.rglob('Macros.h'))
versions = set()
for header in headers:
    text = header.read_text(errors='replace')
    version = [re.search(r'^#define EIGEN_' + part + r'_VERSION\s+(\d+)', text, re.M)
               for part in ('WORLD', 'MAJOR', 'MINOR')]
    if all(version):
        versions.add('.'.join(match[1] for match in version))
if versions != {'3.4.0'}:
    raise RuntimeError(f'Chocolatey Eigen source version is not the recorded 3.4.0: {versions}')

# The pinned release/3.5.0 commit identifies itself as 3.6.0 in its source.
# Record the installed SDK's version, rather than inferring it from the tag.
openms_versions = {match[1] for path in args.openms.rglob('OpenMSConfigVersion.cmake')
                   if (match := re.search(r'set\(PACKAGE_VERSION\s+"([0-9.]+)"\)', path.read_text()))}
if len(openms_versions) != 1:
    raise RuntimeError(f'cannot identify the installed OpenMS version: {openms_versions}')

providers = [
    {'name': 'OpenMS-contrib', 'version': args.contrib_revision,
     'license': 'Multiple upstream licenses; includes EPL-1.0 and MPL-2.0',
     'scope': 'Conservative complete pinned ALL build inputs, including static/header dependencies; not all are linked.',
     'files': [str(p.resolve()) for folder in ('bin', 'lib') for p in (args.contrib / folder).glob('*.dll')
               if p.name.lower() != 'libcurl.dll'],
     'force_include': True, 'collect_sources': True, 'source_recipe': recipes},
    {'name': 'curl', 'version': '8.12.1', 'license': 'curl',
     'files': [str(p.resolve()) for folder in ('bin', 'lib') for p in (args.contrib / folder).glob('*.dll')
               if p.name.lower() == 'libcurl.dll'],
     'collect_sources': True,
     'source_recipe': [{'url': 'https://curl.se/download/curl-8.12.1.tar.gz',
                        'sha256': '7b40ea64947e0b440716a4d7f0b7aa56230a5341c8377d7b609649d4aea8dbcf'}]},
    {'name': 'OpenMS', 'version': openms_versions.pop(), 'license': 'BSD-3-Clause',
     'files': [str(p.resolve()) for p in args.openms.rglob('*.dll')],
     'collect_sources': True,
     'source_recipe': [{'git_url': 'https://github.com/OpenMS/OpenMS.git', 'git_rev': args.openms_revision}],
     'build_patch': 'windows.yml adds #include <functional> to DIAPrescoring.cpp; exact workflow is in build-instructions.'},
    {'name': 'Qt-qtbase', 'version': '6.8.3', 'license': 'LGPL-3.0-only',
     'files': [str((args.qt / 'bin' / name).resolve()) for name in ('Qt6Core.dll', 'Qt6Network.dll')],
     'source_recipe': [{'url': 'https://download.qt.io/archive/qt/6.8/6.8.3/submodules/qtbase-everywhere-src-6.8.3.tar.xz',
                        'sha256': '56001b905601bb9023d399f3ba780d7fa940f3e4861e496a7c490331f49e0b80'}]},
]
args.out.write_text(json.dumps(providers, indent=2) + '\n')
