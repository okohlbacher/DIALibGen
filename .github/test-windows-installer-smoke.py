#!/usr/bin/env python3
"""Offline failures for the shared MSI/NSIS payload verifier."""
import importlib.util
from pathlib import Path
import tempfile

spec = importlib.util.spec_from_file_location('installer', Path(__file__).with_name('windows-installer-smoke.py'))
installer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(installer)

with tempfile.TemporaryDirectory() as temporary:
    root = Path(temporary)
    binary = root / 'PFiles/DIALibGen/resources/dialibgen/bin/DIALibGen.exe'
    model = binary.parents[1] / 'share/DIALibGen/models/peptdeep_rt_dynamic.onnx'
    notice = binary.parents[2] / 'third-party-licenses/LICENSE'
    gui = binary.parents[3] / 'dialibgen-gui.exe'
    for path, content in ((binary, b'CLI'), (model, b'ONNX'), (notice, b'upstream license'), (gui, b'GUI')):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(content)
    cli = installer.inventory(binary.parents[1])
    notices = installer.inventory(notice.parent)
    unpatched = root / 'build-gui.exe'
    # Include a second marker: Tauri changes only the first occurrence.
    unpatched.write_bytes(b'PE-prefix__TAURI_BUNDLE_TYPE_VAR_UNK-suffix__TAURI_BUNDLE_TYPE_VAR_UNK')
    gui_hashes = installer.gui_bundle_hashes(unpatched)
    assert gui_hashes['msi'] != gui_hashes['nsis']
    gui.write_bytes(b'PE-prefix__TAURI_BUNDLE_TYPE_VAR_MSI-suffix__TAURI_BUNDLE_TYPE_VAR_UNK')
    gui_hash = gui_hashes['msi']

    def verify():
        return installer.verify_payload(root, cli, notices, gui_hash)

    assert verify() == binary
    gui.write_bytes(b'PE-prefix__TAURI_BUNDLE_TYPE_VAR_NSS-suffix__TAURI_BUNDLE_TYPE_VAR_UNK')
    assert installer.verify_payload(root, cli, notices, gui_hashes['nsis']) == binary
    for wrong in (gui.read_bytes(), unpatched.read_bytes(),
                  b'changed-prefix__TAURI_BUNDLE_TYPE_VAR_MSI-suffix__TAURI_BUNDLE_TYPE_VAR_UNK',
                  b'PE-prefix__TAURI_BUNDLE_TYPE_VAR_MSI-suffix__TAURI_BUNDLE_TYPE_VAR_MSI'):
        gui.write_bytes(wrong)
        try:
            verify()
        except RuntimeError:
            pass
        else:
            raise AssertionError('accepted wrong installer marker or unrelated GUI byte changes')
    gui.write_bytes(b'PE-prefix__TAURI_BUNDLE_TYPE_VAR_MSI-suffix__TAURI_BUNDLE_TYPE_VAR_UNK')
    unpatched.write_bytes(b'no Tauri marker')
    try:
        installer.gui_bundle_hashes(unpatched)
    except RuntimeError:
        pass
    else:
        raise AssertionError('accepted an unrecognized Tauri build executable')
    for path in (model, notice, gui):
        original = path.read_bytes()
        for content in (None, b'changed bytes'):
            if content is None:
                path.unlink()
            else:
                path.write_bytes(content)
            try:
                verify()
            except RuntimeError:
                pass
            else:
                raise AssertionError(f'accepted missing or changed payload: {path}')
        path.write_bytes(original)
    binary.rename(binary.with_name('wrong-name.exe'))
    try:
        verify()
    except RuntimeError:
        pass
    else:
        raise AssertionError('accepted missing expected embedded CLI layout')
print('PASS: exact per-installer Tauri patch; missing or changed model, notice, GUI and CLI payloads refused')
