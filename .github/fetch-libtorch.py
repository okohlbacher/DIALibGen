#!/usr/bin/env python3
"""Extract a checksum-pinned CPU wheel's C++ SDK; Python is not a runtime dependency."""
import hashlib
import json
import pathlib
from contextlib import contextmanager
import sys
import tempfile
import urllib.request
import zipfile

# Official CPU index: https://download.pytorch.org/whl/cpu/torch/
WHEELS = {
    "linux-x64": ("2.10.0", "manylinux_2_28_x86_64", "b7cb1ec66cefb90fd7b676eac72cfda3b8d4e4d0cacd7a531963bc2e0a9710ab"),
    # 2.10.0 crashes in the native ARM LSTM path. Keep this SDK outside the
    # OpenMS conda solve and require the unchanged parity and training tests.
    "linux-arm64": ("2.14.0", "manylinux_2_28_aarch64", "08894195f84541edcbd09072e6c53b791d4cdd09f650b05473f35718a848fd7b"),
    "windows-x64": ("2.10.0", "win_amd64", "17a09465bab2aab8f0f273410297133d8d8fb6dd84dccbd252ca4a4f3a111847"),
}
INTEL_OPENMP_URL = 'https://files.pythonhosted.org/packages/a9/e4/cacef8ad8570b4ae2b850f3eaab86be0b0296fc9683bdcd9ae8e10be7ecb/intel_openmp-2025.3.1-py2.py3-none-win_amd64.whl'
INTEL_OPENMP_SHA256 = 'bf325e802d9f52f95ca8d2558bb66f30d1839aed30426533063a78a3587e9035'


@contextmanager
def verified_download(url, expected):
    with tempfile.TemporaryFile() as archive:
        digest = hashlib.sha256()
        with urllib.request.urlopen(url, timeout=120) as response:
            while chunk := response.read(1024 * 1024):
                archive.write(chunk)
                digest.update(chunk)
        if digest.hexdigest() != expected:
            raise RuntimeError(f'libtorch: SHA256 mismatch for {url}')
        archive.seek(0)
        yield archive


def intel_openmp_provider(archive, out):
    """Bind both Intel notices to the exact DLLs included by the pinned Torch wheel."""
    hashes = {}
    with zipfile.ZipFile(archive) as package:
        for name in ('libiomp5md.dll', 'libiompstubs5md.dll'):
            expected = package.read('intel_openmp-2025.3.1.data/data/Library/bin/' + name)
            actual = (out / 'lib' / name).read_bytes()
            if actual != expected:
                raise RuntimeError(f'libtorch: Intel OpenMP licensing package does not match {name}')
            hashes['lib/' + name] = hashlib.sha256(actual).hexdigest()
        target = out / 'share/licenses/Intel-OpenMP'
        target.mkdir(parents=True, exist_ok=True)
        for name in ('intel_openmp-2025.3.1.dist-info/LICENSE.txt',
                     'intel_openmp-2025.3.1.data/data/share/doc/compiler/licensing/openmp/third-party-programs.txt'):
            content = package.read(name)
            if not content:
                raise RuntimeError(f'libtorch: empty Intel OpenMP notice: {name}')
            (target / pathlib.PurePosixPath(name).name).write_bytes(content)
    return {'name': 'Intel-OpenMP', 'version': '2025.3.1',
            'license': 'LicenseRef-Intel-Developer-Tools-August-2024 AND bundled component licenses',
            'package_url': INTEL_OPENMP_URL, 'package_sha256': INTEL_OPENMP_SHA256,
            'sdk_file_sha256': hashes, 'license_directory': 'share/licenses/Intel-OpenMP'}


def extract_sdk(archive, out):
    extracted_license = False
    with zipfile.ZipFile(archive) as files:
        for member in files.infolist():
            path = pathlib.PurePosixPath(member.filename)
            if member.is_dir() or ".." in path.parts:
                continue
            if str(path).startswith(("torch/include/", "torch/lib/", "torch/share/")):
                target = out.joinpath(*path.parts[1:])
            elif ".dist-info" in path.parts[0] and len(path.parts) > 2 and path.parts[1] == "licenses":
                target = out.joinpath("share/licenses/Torch", *path.parts[2:])
            elif len(path.parts) == 2 and path.parts[0].endswith('.dist-info') and path.name in {"LICENSE", "NOTICE"}:
                target = out / "share/licenses/Torch" / path.name
            else:
                continue
            # The Python bridge is never linked by Torch's C++ imported targets.
            if "torch_python" in path.name or path.name == "_C.lib":
                continue
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(files.read(member))
            if target == out / 'share/licenses/Torch/LICENSE' and member.file_size:
                extracted_license = True
    if not extracted_license:
        raise RuntimeError('libtorch: wheel lacks a non-empty primary LICENSE')
    if not (out / 'share/cmake/Torch/TorchConfig.cmake').is_file():
        raise RuntimeError('libtorch: wheel lacks TorchConfig.cmake')


def main():
    platform, destination = sys.argv[1:]
    version, wheel_platform, expected = WHEELS[platform]
    url = f"https://download.pytorch.org/whl/cpu/torch-{version}%2Bcpu-cp311-cp311-{wheel_platform}.whl"
    out = pathlib.Path(destination).resolve()
    out.mkdir(parents=True, exist_ok=True)
    with verified_download(url, expected) as archive:
        extract_sdk(archive, out)
    providers = []
    if platform == 'windows-x64':
        with verified_download(INTEL_OPENMP_URL, INTEL_OPENMP_SHA256) as archive:
            providers.append(intel_openmp_provider(archive, out))
    elif platform == 'linux-arm64':
        # Revision/version are embedded in this pinned wheel's ACL binary.
        providers.append({'name': 'Arm-Compute-Library', 'version': '53.2.0', 'license': 'MIT',
            'sdk_file_sha256': {
                'lib/libarm_compute.so': '6fa9225a9d092c7ebf6ff6acea25cd154146da036d8bdef7aa61d2c12a510bc5',
                'lib/libarm_compute_graph.so': '52d09a8bd651415cb93e68b66759201288b5ae323643aef5bdf4c00f310d4334'},
            'collect_sources': True,
            'source_recipe': [{'url': 'https://codeload.github.com/ARM-software/ComputeLibrary/tar.gz/7b256bb7965f2fd99cdee790a4b0e56dab438a8c',
                               'sha256': '6006a52b3ef375a728af688a743cddf8d031acab9e8a2b1076308116365c5eee',
                               'git_rev': '7b256bb7965f2fd99cdee790a4b0e56dab438a8c',
                               'fn': 'arm-compute-7b256bb7965f2fd99cdee790a4b0e56dab438a8c.tar.gz'}]})
    (out / "share/licenses/Torch/provenance.json").write_text(json.dumps({
        "name": "PyTorch", "version": version, "license": "BSD-3-Clause AND bundled component licenses",
        "package_url": url, "package_sha256": expected,
        "source_url": f"https://github.com/pytorch/pytorch/tree/v{version}",
        "native_providers": providers,
    }, indent=2) + "\n")
    print(f"CPU libtorch {version}: {out} (SHA256 {expected})")


if __name__ == '__main__':
    main()
