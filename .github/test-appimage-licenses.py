#!/usr/bin/env python3
"""Offline fixtures for attribution gates; live package closure is checked in CI."""
import hashlib
import importlib.util
import io
from pathlib import Path
import tempfile
import tarfile
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


if __name__ == '__main__':
    unittest.main()
