#!/usr/bin/env python3
"""Offline fixtures for attribution gates; live package closure is checked in CI."""
import hashlib
import importlib.util
import io
import json
import shutil
from pathlib import Path
import tempfile
import tarfile
import sys
import unittest
from unittest.mock import patch

spec = importlib.util.spec_from_file_location('collector', Path(__file__).with_name('collect-appimage-licenses.py'))
collector = importlib.util.module_from_spec(spec)
spec.loader.exec_module(collector)


class AttributionTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name).resolve()
        collector.digest.cache_clear()
        collector.build_id.cache_clear()
        collector.is_elf.cache_clear()

    def test_ownership_requires_bytes_or_build_identity_not_basename(self):
        bundled, original, unrelated = [self.root / name for name in ('bundled', 'original', 'unrelated')]
        bundled.write_bytes(b'\x7fELFpatched-rpath')
        original.write_bytes(b'\x7fELForiginal-rpath')
        unrelated.write_bytes(b'\x7fELFunrelated-file')
        ids = {bundled: '012345', original: '012345', unrelated: 'abcdef'}
        with patch.object(collector, 'build_id', side_effect=lambda path: ids[path]):
            self.assertEqual(collector.matching_owners(bundled, [(unrelated, 'wrong'), (original, 'correct')]),
                             [(original, 'correct')])
        with patch.object(collector, 'build_id', return_value=None):
            self.assertEqual(collector.matching_owners(bundled, [(original, 'must-not-guess')]), [])
        identical = self.root / 'identical'
        identical.write_bytes(bundled.read_bytes())
        self.assertEqual(collector.matching_owners(bundled, [(identical, 'byte-identical')]), [(identical, 'byte-identical')])

    def test_copied_data_requires_exact_bytes_and_never_uses_build_ids(self):
        for filename in ('schema.gschema.xml', 'theme.css', 'model.onnx', 'COPYRIGHT'):
            source = self.root / ('source-' + filename)
            copied = self.root / ('copied-' + filename)
            changed = self.root / ('changed-' + filename)
            source.write_bytes(b'known original data')
            copied.write_bytes(b'known original data')
            changed.write_bytes(b'changed data with the same claimed build ID')
            with patch.object(collector, 'build_id', return_value='identical') as build_id:
                self.assertEqual(collector.matching_owners(copied, [(source, 'package')]), [(source, 'package')])
                self.assertEqual(collector.matching_owners(changed, [(source, 'package')]), [])
                build_id.assert_not_called()

    def test_gui_assets_use_explicit_source_and_build_outputs(self):
        gui = self.root / 'gui'
        files = {'src/logo.svg': b'authored logo', 'dist/index.html': b'generated frontend',
                 'src-tauri/icons/128x128.png': b'authored icon',
                 'src-tauri/resources/third-party-licenses/inventory.json': b'license inventory',
                 'node_modules/extra/data.bin': b'not automatically first party'}
        for name, contents in files.items():
            path = gui / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(contents)
        assets = collector.gui_asset_hashes(gui)
        self.assertEqual(len(assets), 4)
        self.assertNotIn(hashlib.sha256(files['node_modules/extra/data.bin']).hexdigest(), assets)

    def test_descriptors_match_only_named_same_build_debian_outputs(self):
        fixture = self.root / 'data.tar'
        desktop = b'[Desktop Entry]\nName=DIALibGen\nExec=dialibgen-gui\n'
        with tarfile.open(fixture, 'w') as archive:
            for name, data in [('usr/share/applications/DIALibGen.desktop', desktop),
                               ('usr/share/applications/unrelated.desktop', b'unrelated'),
                               ('usr/lib/extra.dat', b'arbitrary payload')]:
                info = tarfile.TarInfo('./' + name)
                info.size = len(data)
                archive.addfile(info, io.BytesIO(data))
        class Process:
            def __enter__(self):
                self.stdout = fixture.open('rb')
                return self
            def __exit__(self, *_args):
                self.stdout.close()
            def wait(self, timeout):
                return 0
        with patch.object(collector.subprocess, 'Popen', return_value=Process()):
            descriptors = collector.deb_descriptor_hashes(self.root / 'same-build.deb', 'DIALibGen')
        self.assertEqual(descriptors, {hashlib.sha256(desktop).hexdigest(): 'usr/share/applications/DIALibGen.desktop'})

    def test_generated_gtk_cache_needs_exact_regeneration(self):
        appdir = self.root / 'AppDir'
        cache = appdir / 'usr/lib/gtk-3.0/3.0.0/immodules.cache'
        cache.parent.mkdir(parents=True)
        cache.write_bytes(b'"im-cedilla.so"\n"cedilla" "Cedilla"\n')
        tool = self.root / 'gtk-query-immodules-3.0'
        tool.write_bytes(b'#!/bin/sh\n')
        tool.chmod(0o700)
        original = b'"/usr/lib/test/gtk-3.0/3.0.0/immodules/im-cedilla.so"\n"cedilla" "Cedilla"\n'
        owners = {'gtk-query-immodules-3.0': [(tool, 'libgtk-3-bin')]}
        with patch.object(collector.subprocess, 'check_output', return_value=original), \
             patch.object(collector, 'run', return_value='/usr/lib/test'):
            self.assertEqual(collector.generated_cache_owner(cache, appdir, owners), ('libgtk-3-bin', str(tool)))
            cache.write_bytes(b'arbitrary replacement cache')
            collector.digest.cache_clear()
            self.assertIsNone(collector.generated_cache_owner(cache, appdir, owners))
        extra = appdir / 'unproven.cache'
        extra.write_bytes(original)
        self.assertIsNone(collector.generated_cache_owner(extra, appdir, owners))

    def test_exact_ubuntu_source_version_and_sha256_gate(self):
        source = b'original upstream source archive'
        checksum = hashlib.sha256(source).hexdigest()
        descriptor = (f'Source: example\nVersion: 1:2.3-4ubuntu5\nChecksums-Sha256:\n'
                      f' {checksum} {len(source)} example.orig.tar.xz\n').encode()
        entry = {'source_package_name': 'example', 'source_package_version': '1:2.3-4ubuntu5',
                 'self_link': 'https://api.launchpad.net/publication/1'}
        urls = ['https://launchpad.net/example.dsc', 'https://launchpad.net/example.orig.tar.xz']
        content = dict(zip(urls, (descriptor, source)))
        def urlopen(request, **_kwargs):
            url = request.full_url if hasattr(request, 'full_url') else request
            return io.BytesIO(content[url])
        with patch.object(collector, 'get_json', side_effect=[{'entries': [entry]}, urls]), \
             patch.object(collector.urllib.request, 'urlopen', side_effect=urlopen):
            result = collector.ubuntu_sources('example', '1:2.3-4ubuntu5', self.root / 'good')
        self.assertEqual(result['files'][1]['sha256'], checksum)
        self.assertIn('OpenPGP signature not verified', result['descriptor_trust'])
        self.assertEqual((self.root / 'good/example.orig.tar.xz').read_bytes(), source)
        content[urls[1]] = b'tampered source'
        with patch.object(collector, 'get_json', side_effect=[{'entries': [entry]}, urls]), \
             patch.object(collector.urllib.request, 'urlopen', side_effect=urlopen):
            with self.assertRaisesRegex(RuntimeError, 'checksum mismatch'):
                collector.ubuntu_sources('example', '1:2.3-4ubuntu5', self.root / 'bad')
        self.assertFalse((self.root / 'bad/example.orig.tar.xz').exists())
        with patch.object(collector, 'get_json', return_value={'entries': [entry]}):
            with self.assertRaisesRegex(RuntimeError, 'no authoritative source'):
                collector.ubuntu_sources('example', '1:2.3-4ubuntu6', self.root / 'wrong-version')

    def test_clearsigned_descriptor_retains_exact_version(self):
        text = ('-----BEGIN PGP SIGNED MESSAGE-----\nHash: SHA256\n\n'
                'Source: sample\nVersion: 1:2-3\n\n-----BEGIN PGP SIGNATURE-----\nsignature\n')
        self.assertEqual(collector.dsc_fields(text)['Version'], '1:2-3')

    def test_final_appimage_payload_must_match_every_file_and_link(self):
        appdir = self.root / 'AppDir'
        appdir.mkdir()
        (appdir / 'usr').mkdir()
        (appdir / 'usr/library.so').write_bytes(b'actual library bytes')
        (appdir / 'AppRun').symlink_to('usr/library.so')
        (appdir / 'dangling').symlink_to('missing-target')
        image = self.root / 'final.AppImage'
        mutations = {
            'same': lambda tree: None,
            'missing': lambda tree: (tree / 'usr/library.so').unlink(),
            'extra': lambda tree: (tree / 'unattributed').write_bytes(b'extra payload'),
            'changed': lambda tree: (tree / 'usr/library.so').write_bytes(b'changed library bytes'),
            'link-target': lambda tree: ((tree / 'AppRun').unlink(), (tree / 'AppRun').symlink_to('different-target')),
            'file-type': lambda tree: ((tree / 'AppRun').unlink(), (tree / 'AppRun').write_text('usr/library.so')),
        }
        for name, mutate in mutations.items():
            with self.subTest(name=name):
                def extract(command, *, cwd, **kwargs):
                    self.assertEqual(command, [str(image), '--appimage-extract'])
                    self.assertTrue(kwargs['check'])
                    tree = Path(cwd) / 'squashfs-root'
                    shutil.copytree(appdir, tree, symlinks=True)
                    mutate(tree)
                with patch.object(collector.subprocess, 'run', side_effect=extract):
                    if name == 'same':
                        collector.validate_appimage_payload(appdir, image)
                    else:
                        with self.assertRaisesRegex(RuntimeError, 'payload differs from AppDir'):
                            collector.validate_appimage_payload(appdir, image)

    def test_runtime_provider_is_bound_to_exact_header(self):
        image = self.root / 'sample.AppImage'
        image.write_bytes(b'ELF header' + b'SquashFS payload')
        size = len(b'ELF header')
        provider = {'size': size, 'sha256': hashlib.sha256(b'ELF header').hexdigest()}
        with patch.object(collector, 'run', return_value=str(size)):
            collector.validate_runtime_provider(image, provider)
            with self.assertRaisesRegex(RuntimeError, 'hash differs'):
                collector.validate_runtime_provider(image, {**provider, 'sha256': '0' * 64})
            with self.assertRaisesRegex(RuntimeError, 'offset differs'):
                collector.validate_runtime_provider(image, {**provider, 'size': size + 1})
        with patch.object(collector, 'run', return_value=str(image.stat().st_size)):
            with self.assertRaisesRegex(RuntimeError, 'offset differs'):
                collector.validate_runtime_provider(image, {**provider, 'size': image.stat().st_size})

    def test_required_notices_must_survive_in_final_payload_as_exact_bytes(self):
        original = self.root / 'upstream'
        original.mkdir()
        notices = {}
        for index, name in enumerate(('packages/sample/copyright', 'providers/helper/0-LICENSE',
                                      'providers/second-helper/0-LICENSE', 'common-licenses/LGPL-2.1')):
            source = original / str(index)
            source.write_text(f'Original copyright and terms for {name}\n')
            notices[name] = source
        resources = self.root / 'resources'
        resources.mkdir()
        (resources / 'stale-notice').write_text('from an older package')
        collector.stage_notices(notices, resources)
        self.assertFalse((resources / 'stale-notice').exists())
        for name, source in notices.items():
            self.assertEqual((resources / name).read_bytes(), source.read_bytes())
        records = [{'path': 'usr/share/resources/' + name, 'sha256': collector.digest(resources / name)}
                   for name in notices]
        matches = collector.verify_bundled_notices(notices, records)
        self.assertEqual(set(matches), set(notices))
        self.assertTrue(all(len(paths) == 1 for paths in matches.values()))
        # A notice present only in sources, or changed under the right filename, cannot satisfy the gate.
        for changed in (records[:-1], [*records[:-1], {**records[-1], 'sha256': '0' * 64}]):
            with self.assertRaisesRegex(RuntimeError, 'final AppImage lacks required notice bytes.*LGPL-2.1'):
                collector.verify_bundled_notices(notices, changed)

    def test_prepare_stages_notices_without_collecting_sources_or_claiming_completion(self):
        gui = self.root / 'gui'
        (gui / 'src-tauri').mkdir(parents=True)
        (gui / 'src-tauri/tauri.conf.json').write_text('{"productName":"DIALibGen"}')
        appdir = self.root / 'AppDir'
        appdir.mkdir()
        binary = appdir / 'gui-binary'
        binary.write_bytes(b'known GUI executable')
        notice = self.root / 'original-notice'
        notice.write_text('Original helper copyright and permission\n')
        providers = self.root / 'providers.json'
        providers.write_text(json.dumps([{'name': 'runtime', 'license': 'MIT', 'notices': [str(notice)],
                                          'appimage_runtime': {'size': 1, 'sha256': '0' * 64}}]))
        argv = ['collector', '--appdir', str(appdir), '--appimage', str(self.root / 'final.AppImage'),
                '--baseline', str(self.root / 'stage'), '--gui-binary', str(binary), '--gui-root', str(gui),
                '--deb', str(self.root / 'same-build.deb'), '--sources', str(self.root / 'sources'),
                '--source-asset', 'sources.tar.gz', '--providers', str(providers), '--prepare-notices']
        with patch.object(sys, 'argv', argv), patch.object(collector, 'validate_appimage_payload'), \
             patch.object(collector, 'validate_runtime_provider'), \
             patch.object(collector, 'installed_packages', return_value=({}, {})), \
             patch.object(collector, 'deb_descriptor_hashes', return_value={}), \
             patch.object(collector, 'required_notices', return_value={'providers/runtime/LICENSE': notice}), \
             patch.object(collector, 'ubuntu_sources') as sources, patch.object(collector, 'fetch') as fetch:
            collector.main()
            sources.assert_not_called()
            fetch.assert_not_called()
        self.assertEqual((gui / 'src-tauri/resources/third-party-licenses/appimage/providers/runtime/LICENSE').read_bytes(),
                         notice.read_bytes())
        self.assertFalse((self.root / 'sources/appimage/inventory.json').exists())


if __name__ == '__main__':
    unittest.main()
