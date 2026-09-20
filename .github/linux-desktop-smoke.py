#!/usr/bin/env python3
"""Install and open the Linux desktop packages on a disposable Ubuntu runner."""
import argparse
import json
import os
from pathlib import Path
import re
import signal
import subprocess
import sys
import tempfile
import time


def run(command, log, *, timeout=120, env=None):
    with log.open('w') as output:
        output.write(f'command: {command!r}\n'); output.flush()
        result = subprocess.run(command, stdout=output, stderr=subprocess.STDOUT,
                                timeout=timeout, env=env)
    if result.returncode:
        raise RuntimeError(f'command exited {result.returncode}; see {log}')


def visible_window(executable, output):
    # Each invocation has its own X server and session bus; no existing window
    # or single-instance process can satisfy this check.
    with (output / 'gui.log').open('w') as log:
        process = subprocess.Popen([str(executable)], stdout=log, stderr=subprocess.STDOUT,
                                   start_new_session=True)
        try:
            deadline = time.monotonic() + 60
            visible_since = None
            while time.monotonic() < deadline:
                if process.poll() is not None:
                    raise RuntimeError(f'GUI exited {process.returncode}; see {output / "gui.log"}')
                tree = subprocess.check_output(['xwininfo', '-root', '-tree'], text=True, timeout=10)
                (output / 'windows.log').write_text(tree)
                visible = False
                for window in re.findall(r'^\s*(0x[0-9a-fA-F]+) "DIALibGen":', tree, re.MULTILINE):
                    info = subprocess.check_output(['xwininfo', '-id', window, '-stats'], text=True, timeout=10)
                    pid = subprocess.check_output(['xprop', '-id', window, '_NET_WM_PID'], text=True, timeout=10)
                    (output / 'window.log').write_text(info + pid)
                    if ('Map State: IsViewable' in info
                            and re.search(r'=\s*' + str(process.pid) + r'\s*$', pid)):
                        visible = True
                        break
                if visible:
                    visible_since = visible_since or time.monotonic()
                    if time.monotonic() - visible_since >= 5:
                        print(f'PASS: {executable} has a visible DIALibGen window for five seconds', flush=True)
                        return
                else:
                    visible_since = None
                time.sleep(0.5)
            raise RuntimeError(f'GUI did not keep a visible DIALibGen window; see {output}')
        finally:
            try:
                os.killpg(process.pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
            try:
                process.wait(timeout=10)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait(timeout=10)


def launch(executable, output):
    output.mkdir()
    with tempfile.TemporaryDirectory(prefix='dialibgen-desktop-home-') as home:
        # Deliberately exclude conda and build-time library/model overrides.
        env = {'PATH': '/usr/bin:/bin', 'HOME': home, 'LANG': 'C.UTF-8',
               'XDG_RUNTIME_DIR': home, 'GDK_BACKEND': 'x11', 'LIBGL_ALWAYS_SOFTWARE': '1'}
        run(['xvfb-run', '-a', '--server-args=-screen 0 1280x900x24 -nolisten tcp',
             'dbus-run-session', '--', sys.executable, str(Path(__file__).resolve()),
             '--launch', str(executable), '--output', str(output)],
            output / 'session.log', timeout=120, env=env)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--deb', type=Path)
    parser.add_argument('--appdir', type=Path)
    parser.add_argument('--python', type=Path, help='Python with the existing CLI fixture dependencies')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--launch', type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if sys.platform != 'linux':
        parser.error('actual installation and startup require Linux')
    output = args.output.resolve()
    if args.launch:
        visible_window(args.launch, output)
        return
    if not all((args.deb, args.appdir, args.python)):
        parser.error('--deb, --appdir and --python are required')
    output.mkdir(parents=True, exist_ok=False)
    deb = args.deb.resolve()
    package = subprocess.check_output(['dpkg-deb', '-f', str(deb), 'Package'], text=True).strip()
    status = subprocess.run(['dpkg-query', '-W', '-f=${Status}', package], capture_output=True, text=True)
    if status.returncode == 0 and status.stdout == 'install ok installed':
        raise RuntimeError(f'{package} is already installed; use a disposable runner')
    try:
        run(['sudo', '-n', 'dpkg', '--install', str(deb)], output / 'install.log')
        files = subprocess.check_output(['dpkg-query', '-L', package], text=True).splitlines()
        (output / 'installed-files.json').write_text(json.dumps(files, indent=2) + '\n')
        gui = [Path(name) for name in files if name == '/usr/bin/dialibgen-gui']
        cli = [Path(name) for name in files if name.endswith('/resources/dialibgen/bin/DIALibGen')]
        if len(gui) != 1 or len(cli) != 1 or not all(os.access(path, os.X_OK) for path in gui + cli):
            raise RuntimeError('installed Debian package is missing its GUI or embedded CLI')
        launch(gui[0], output / 'deb-gui')
        run([str(args.python.resolve()), str(Path(__file__).with_name('release-smoke.py').resolve()),
             str(cli[0]), '--bare', '--output', str(output / 'deb-cli')],
            output / 'deb-cli.log', timeout=1800)
    finally:
        run(['sudo', '-n', 'dpkg', '--purge', package], output / 'uninstall.log')
    if any(path.exists() for path in gui + cli):
        raise RuntimeError('Debian uninstall left an installed executable behind')
    # The preceding CI collector has compared every file/link in this AppDir
    # with the final AppImage payload. AppRun exercises its real launch hooks.
    launch(args.appdir.resolve() / 'AppRun', output / 'appimage-gui')
    print('PASS: Debian install/startup/embedded CLI/uninstall and AppImage startup')


if __name__ == '__main__':
    main()
