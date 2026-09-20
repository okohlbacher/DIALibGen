#!/usr/bin/env python3
"""Offline failures for the shared MSI/NSIS payload verifier."""
import importlib.util
import json
import ntpath
from pathlib import Path
import tempfile
from types import SimpleNamespace
from unittest.mock import patch

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

system = {'SystemDrive': 'C:', 'SystemRoot': r'C:\Windows', 'WINDIR': r'C:\Windows',
          'ProgramData': r'C:\ProgramData', 'ALLUSERSPROFILE': r'C:\ProgramData'}
with patch.dict(installer.os.environ, {**system, 'PATH': 'build-sdk', 'OPENMS_DATA_PATH': 'build-data'}, clear=True):
    environment = installer.gui_environment(Path('installed/dialibgen-gui.exe'))
assert all(environment[key] == value for key, value in system.items())
assert 'build-sdk' not in environment['PATH'] and 'OPENMS_DATA_PATH' not in environment
with patch.dict(installer.os.environ, environment, clear=True):
    assert ntpath.expandvars(r'%SystemDrive%\ProgramData\Microsoft\Windows\Caches') == r'C:\ProgramData\Microsoft\Windows\Caches'

ready = {'name': 'Reset', 'control_type': 'ControlType.Button', 'enabled': True, 'process_id': 47}
with patch.dict(installer.os.environ, system, clear=True):
    with patch.object(installer.subprocess, 'run', return_value=SimpleNamespace(stdout=json.dumps(ready))) as probe:
        assert installer.wait_for_frontend_ready(123, 57) == ready
        assert probe.call_args.kwargs['timeout'] == 57  # remaining shared startup budget
        assert 'FromHandle([IntPtr]123)' in probe.call_args.args[0][-1]
    for invalid in ({**ready, 'enabled': False}, {**ready, 'name': 'DIALibGen'},
                    {**ready, 'control_type': 'ControlType.Window'}, []):
        with patch.object(installer.subprocess, 'run', return_value=SimpleNamespace(stdout=json.dumps(invalid))):
            try:
                installer.wait_for_frontend_ready(123, 57)
            except RuntimeError:
                pass
            else:
                raise AssertionError('accepted absent, disabled, or non-button frontend control')
    for failure in (installer.subprocess.TimeoutExpired('UIA', 57),
                    installer.subprocess.CalledProcessError(1, 'UIA', stderr='UIA unavailable')):
        with patch.object(installer.subprocess, 'run', side_effect=failure):
            try:
                installer.wait_for_frontend_ready(123, 57)
            except RuntimeError:
                pass
            else:
                raise AssertionError('accepted failed or timed out native UIA probe')
    with patch.object(installer.subprocess, 'run') as probe:
        try:
            installer.wait_for_frontend_ready(123, 0)
        except RuntimeError:
            pass
        else:
            raise AssertionError('extended an exhausted startup budget')
        probe.assert_not_called()

# No native Windows API is mocked into success: these exercise the timing gate;
# the actual EnumWindows/installed-process proof runs on the Windows runner.
def window_check(find_window, exits=False):
    clock = [0.0]

    def advance(seconds):
        clock[0] += seconds

    process = SimpleNamespace(pid=42, returncode=0 if exits else None, poll=lambda: 0 if exits else None)
    with patch.object(installer.time, 'monotonic', lambda: clock[0]), patch.object(installer.time, 'sleep', advance):
        return installer.wait_for_gui_window(process, lambda pid: find_window(clock[0], pid), timeout=2, stable_seconds=1)


window = {'hwnd': 123, 'title': 'DIALibGen', 'width': 1180, 'height': 820}
sent = []
installer.post_close(window, lambda *args: sent.append(args) or True,
                     lambda: OSError(5, 'Access is denied'))
assert sent == [(123, 0x0010, 0, 0)]
try:
    installer.post_close(window, lambda *args: False, lambda: OSError(5, 'Access is denied'))
except OSError as error:
    assert error.errno == 5
else:
    raise AssertionError('accepted a failed WM_CLOSE dispatch')
assert window_check(lambda elapsed, pid: window if pid == 42 and elapsed >= 0.25 else None) == window
for finder, exits in ((lambda elapsed, pid: None, False),
                      (lambda elapsed, pid: window if elapsed < 0.5 else None, False),
                      (lambda elapsed, pid: {**window, 'hwnd': int(elapsed * 4)}, False),
                      (lambda elapsed, pid: window, True)):
    try:
        window_check(finder, exits)
    except RuntimeError:
        pass
    else:
        raise AssertionError('accepted missing, disappearing, replaced window or exited GUI process')
print('PASS: exact installer bytes; changed payloads refused; GUI liveness and close dispatch enforced')
