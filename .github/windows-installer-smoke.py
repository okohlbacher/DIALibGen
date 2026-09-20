#!/usr/bin/env python3
"""Verify Windows MSI/NSIS payloads and run their embedded CLI in a bare environment.

MSI is administratively extracted, not installed. NSIS installs into a unique
runner-temp directory and is uninstalled afterward; use a disposable runner.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time


def digest(path):
    with path.open('rb') as source:
        return hashlib.file_digest(source, 'sha256').hexdigest()


def inventory(directory):
    return {path.relative_to(directory).as_posix(): digest(path)
            for path in directory.rglob('*') if path.is_file()}


def gui_bundle_hashes(path):
    # Tauri CLI 2.11.4 patches only the first marker, then restores this build
    # binary after each installer. Compare exact per-format bytes, not a
    # normalized executable: every byte outside this fixed marker must match.
    # https://github.com/tauri-apps/tauri/blob/8909f221d1515955fc843808032bdc5d62209c96/crates/tauri-bundler/src/bundle.rs#L32-L85
    binary = path.read_bytes()
    marker = b'__TAURI_BUNDLE_TYPE_VAR_UNK'
    offset = binary.find(marker)
    if offset < 0:
        raise RuntimeError(f'unpatched Tauri bundle marker is missing: {path}')
    return {kind: hashlib.sha256(binary[:offset] + marker[:-3] + token + binary[offset + len(marker):]).hexdigest()
            for kind, token in (('msi', b'MSI'), ('nsis', b'NSS'))}


def compare_tree(expected, directory):
    actual = inventory(directory)
    missing = expected.keys() - actual.keys()
    extra = actual.keys() - expected.keys()
    changed = [name for name in expected.keys() & actual.keys() if expected[name] != actual[name]]
    if missing or extra or changed:
        raise RuntimeError(f'installer payload differs at {directory}: '
                           f'missing={sorted(missing)}, extra={sorted(extra)}, changed={sorted(changed)}')


def verify_payload(directory, cli_files, gui_notices, gui_sha256):
    candidates = list(directory.rglob('resources/dialibgen/bin/DIALibGen.exe'))
    if len(candidates) != 1:
        raise RuntimeError(f'expected one embedded CLI under resources/dialibgen/bin: {candidates}')
    binary = candidates[0]
    compare_tree(cli_files, binary.parents[1])
    compare_tree(gui_notices, binary.parents[2] / 'third-party-licenses')
    gui = binary.parents[3] / 'dialibgen-gui.exe'
    if not gui.is_file() or digest(gui) != gui_sha256:
        raise RuntimeError(f'installed GUI executable differs from the packaged binary: {gui}')
    return binary


def run(command, log, timeout):
    # A raw command string is needed only for NSIS's final unquoted /D= and _?=.
    # shell=False passes it directly to CreateProcess, never to cmd.exe.
    with log.open('w', encoding='utf-8') as output:
        output.write(f'command: {command!r}\n'); output.flush()
        result = subprocess.run(command, stdout=output, stderr=subprocess.STDOUT, timeout=timeout)
        output.write(f'\nexit: {result.returncode}\n')
    if result.returncode != 0:
        raise RuntimeError(f'command exited {result.returncode}; see {log}')


def wait_for_gui_window(process, find_window, timeout=60, stable_seconds=3):
    deadline = time.monotonic() + timeout
    first_seen, previous = None, None
    while time.monotonic() < deadline:
        if process.poll() is not None:
            raise RuntimeError(f'installed GUI exited before its window check: {process.returncode}')
        window = find_window(process.pid)
        if window is None or previous is None or window['hwnd'] != previous['hwnd']:
            first_seen = time.monotonic() if window else None
        if window and time.monotonic() - first_seen >= stable_seconds:
            return window
        previous = window
        time.sleep(0.25)
    raise RuntimeError('installed GUI did not keep a visible DIALibGen main window open')


def gui_environment(gui):
    # Windows shell/WebView cache paths expand these OS variables. Missing
    # SystemDrive creates a literal %SystemDrive% directory beneath the cwd.
    env = {key: os.environ[key] for key in ('SystemDrive', 'SystemRoot', 'WINDIR', 'ProgramData',
           'ALLUSERSPROFILE', 'TEMP', 'TMP', 'COMSPEC', 'USERPROFILE', 'APPDATA', 'LOCALAPPDATA',
           'ProgramFiles', 'ProgramFiles(x86)', 'ProgramW6432') if key in os.environ}
    env['PATH'] = str(gui.parent) + os.pathsep + str(Path(os.environ['SystemRoot']) / 'System32')
    env['OPENMS_DISABLE_UPDATE_CHECK'] = 'ON'
    return env


def exercise_gui(gui, log):
    # Native window ownership excludes another single-instance process and
    # WebView2 helper processes. This checks startup, not frontend rendering.
    import ctypes
    from ctypes import wintypes
    user32 = ctypes.WinDLL('user32', use_last_error=True)
    callback_type = ctypes.WINFUNCTYPE(wintypes.BOOL, wintypes.HWND, wintypes.LPARAM)
    user32.EnumWindows.argtypes = [callback_type, wintypes.LPARAM]
    user32.EnumWindows.restype = wintypes.BOOL
    user32.GetWindowThreadProcessId.argtypes = [wintypes.HWND, ctypes.POINTER(wintypes.DWORD)]
    user32.GetWindowThreadProcessId.restype = wintypes.DWORD
    user32.IsWindowVisible.argtypes = [wintypes.HWND]
    user32.IsWindowVisible.restype = wintypes.BOOL
    user32.GetClientRect.argtypes = [wintypes.HWND, ctypes.POINTER(wintypes.RECT)]
    user32.GetClientRect.restype = wintypes.BOOL
    user32.GetWindowTextW.argtypes = [wintypes.HWND, wintypes.LPWSTR, ctypes.c_int]
    user32.GetWindowTextW.restype = ctypes.c_int
    user32.PostMessageW.argtypes = [wintypes.HWND, wintypes.UINT, wintypes.WPARAM, wintypes.LPARAM]
    user32.PostMessageW.restype = wintypes.BOOL

    def find_window(pid):
        windows = []

        @callback_type
        def inspect(hwnd, _):
            owner = wintypes.DWORD()
            user32.GetWindowThreadProcessId(hwnd, ctypes.byref(owner))
            if owner.value == pid and user32.IsWindowVisible(hwnd):
                title = ctypes.create_unicode_buffer(256)
                rect = wintypes.RECT()
                user32.GetWindowTextW(hwnd, title, len(title))
                if (title.value == 'DIALibGen' and user32.GetClientRect(hwnd, ctypes.byref(rect))
                        and rect.right > rect.left and rect.bottom > rect.top):
                    windows.append({'hwnd': hwnd, 'title': title.value,
                                    'width': rect.right - rect.left, 'height': rect.bottom - rect.top})
            return True

        if not user32.EnumWindows(inspect, 0):
            raise ctypes.WinError(ctypes.get_last_error())
        return windows[0] if windows else None

    with log.open('w', encoding='utf-8') as output:
        process = subprocess.Popen([str(gui)], cwd=gui.parent, env=gui_environment(gui),
                                   stdout=output, stderr=subprocess.STDOUT)
        output.write(f'launched installed GUI: {gui}; pid={process.pid}\n'); output.flush()
        try:
            window = wait_for_gui_window(process, find_window)
            output.write(f'visible native main window stayed open for 3 seconds: {json.dumps(window)}\n')
            output.flush()
        finally:
            if process.poll() is None:
                # Recheck ownership immediately before closing this process's window.
                try:
                    closing = find_window(process.pid)
                    if closing:
                        user32.PostMessageW(closing['hwnd'], 0x0010, 0, 0)  # WM_CLOSE
                except OSError as error:
                    output.write(f'window enumeration failed during cleanup: {error}\n')
                try:
                    process.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    output.write('GUI did not close; forcibly cleaning its process tree\n'); output.flush()
                    try:
                        subprocess.run([str(Path(os.environ['SystemRoot']) / 'System32/taskkill.exe'),
                                        '/PID', str(process.pid), '/T', '/F'], stdout=output,
                                       stderr=subprocess.STDOUT, timeout=15, check=True)
                    finally:
                        if process.poll() is None:
                            process.kill()
                        process.wait(timeout=5)
                    raise RuntimeError(f'installed GUI required forced termination; see {log}')
            output.write(f'GUI exit: {process.returncode}\n')
        if process.returncode != 0:
            raise RuntimeError(f'installed GUI exited abnormally: {process.returncode}; see {log}')
    return {**window, 'stable_seconds': 3, 'clean_exit': True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--msi', type=Path, required=True)
    parser.add_argument('--nsis', type=Path, required=True)
    parser.add_argument('--stage', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if os.name != 'nt':
        parser.error('actual installer verification requires Windows')
    root = Path(__file__).resolve().parent.parent
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    cli_files = inventory(args.stage.resolve())
    required = {'bin/DIALibGen.exe', 'share/OpenMS/CHEMISTRY/unimod.xml',
                'share/DIALibGen/irt_standards.tsv', 'share/DIALibGen/runtime-dependencies.json'}
    required.update(f'share/DIALibGen/models/peptdeep_{head}_dynamic.onnx' for head in ('rt', 'ccs', 'ms2'))
    if not required <= cli_files.keys() or not any(name.startswith('share/DIALibGen/licenses/runtime/')
                                                   and not name.endswith('.json') for name in cli_files):
        raise RuntimeError('staged CLI is missing required models, data or runtime notices')
    gui_notices = inventory(root / 'gui/src-tauri/resources/third-party-licenses')
    if 'inventory.json' not in gui_notices or len(gui_notices) < 2:
        raise RuntimeError('GUI third-party notices are missing')
    gui_hashes = gui_bundle_hashes(root / 'gui/src-tauri/target/release/dialibgen-gui.exe')
    record = {'msi_sha256': digest(args.msi.resolve()), 'nsis_sha256': digest(args.nsis.resolve()),
              'gui_sha256_by_installer': gui_hashes,
              'cli_files': len(cli_files), 'gui_notice_files': len(gui_notices), 'checks': {}}

    def exercise(kind, directory):
        binary = verify_payload(directory, cli_files, gui_notices, gui_hashes[kind])
        record['checks'][kind] = {'embedded_cli': str(binary), 'payload_identical': True, 'all_modes_passed': False}
        run([sys.executable, str(root / '.github/release-smoke.py'), str(binary), '--bare',
             '--output', str(out / kind)], out / f'{kind}-smoke.log', 1800)
        record['checks'][kind]['all_modes_passed'] = True

    try:
        with tempfile.TemporaryDirectory(prefix='dialibgen-installer-', dir=os.environ['RUNNER_TEMP']) as temporary:
            work = Path(temporary)
            msi_root = work / 'msi'
            msi_root.mkdir()
            run([str(Path(os.environ['SystemRoot']) / 'System32/msiexec.exe'), '/a', str(args.msi.resolve()),
                 '/qn', '/norestart', f'TARGETDIR={msi_root}', '/L*v', str(out / 'msi-admin.log')],
                out / 'msi-process.log', 600)
            exercise('msi', msi_root)

            nsis_root = work / 'nsis'
            uninstaller = nsis_root / 'uninstall.exe'
            try:
                # /D must be last and unquoted even when its directory contains spaces.
                command = subprocess.list2cmdline([str(args.nsis.resolve()), '/S', '/NS']) + f' /D={nsis_root}'
                run(command, out / 'nsis-install.log', 600)
                if not uninstaller.is_file():
                    raise RuntimeError(f'NSIS did not install its uninstaller into {nsis_root}')
                exercise('nsis', nsis_root)
                record['checks']['nsis']['gui_window'] = exercise_gui(nsis_root / 'dialibgen-gui.exe', out / 'nsis-gui.log')
            finally:
                if uninstaller.is_file():
                    # Run an external copy with _?= so NSIS waits rather than
                    # spawning another process, and can delete its installed copy.
                    copied = work / 'uninstall-check.exe'
                    shutil.copy2(uninstaller, copied)
                    command = subprocess.list2cmdline([str(copied), '/S']) + f' _?={nsis_root}'
                    run(command, out / 'nsis-uninstall.log', 300)
                    remaining = [str(path.relative_to(nsis_root)) for path in nsis_root.rglob('*') if path.is_file()]
                    if remaining:
                        raise RuntimeError(f'NSIS uninstall left installed files: {remaining}')
                    record['nsis_payload_removed'] = True
    finally:
        (out / 'verification.json').write_text(json.dumps(record, indent=2) + '\n', encoding='utf-8')
    print('PASS: MSI and NSIS preserve packaged resources; embedded generate/refine/tune work; '
          'installed NSIS GUI opens a native main window; NSIS payload removed')


if __name__ == '__main__':
    main()
