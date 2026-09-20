#!/usr/bin/env python3
"""Extract a checksum-pinned CPU wheel's C++ SDK; Python is not a runtime dependency."""
import hashlib
import pathlib
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
platform, destination = sys.argv[1:]
version, wheel_platform, expected = WHEELS[platform]
url = f"https://download.pytorch.org/whl/cpu/torch-{version}%2Bcpu-cp311-cp311-{wheel_platform}.whl"
out = pathlib.Path(destination).resolve()
out.mkdir(parents=True, exist_ok=True)
with tempfile.TemporaryFile() as archive:
    digest = hashlib.sha256()
    with urllib.request.urlopen(url, timeout=120) as response:
        while chunk := response.read(1024 * 1024):
            archive.write(chunk)
            digest.update(chunk)
    if digest.hexdigest() != expected:
        sys.exit(f"libtorch: SHA256 mismatch for {url}")
    archive.seek(0)
    with zipfile.ZipFile(archive) as files:
        for member in files.infolist():
            path = pathlib.PurePosixPath(member.filename)
            if member.is_dir() or ".." in path.parts:
                continue
            if str(path).startswith(("torch/include/", "torch/lib/", "torch/share/")):
                target = out.joinpath(*path.parts[1:])
            elif path.name == "LICENSE" and ".dist-info" in str(path):
                target = out / "share/licenses/Torch/LICENSE"
            else:
                continue
            # The Python bridge is never linked by Torch's C++ imported targets.
            if "torch_python" in path.name or path.name == "_C.lib":
                continue
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(files.read(member))
assert (out / "share/cmake/Torch/TorchConfig.cmake").is_file()
print(f"CPU libtorch {version}: {out} (SHA256 {expected})")
