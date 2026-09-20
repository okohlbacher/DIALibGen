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
    print('PASS: MSI and NSIS preserve packaged resources; embedded generate/refine/tune work; NSIS payload removed')


if __name__ == '__main__':
    main()
