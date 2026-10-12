#!/usr/bin/env python3
# Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
# SPDX-License-Identifier: BSD-3-Clause
"""Check the OpenDIAlyzer code vendored into src/odia-core/ against its manifest.

  sync-odia-core.py --check              every vendored file is listed, matches its recorded
                                         sha256 and carries the SPDX line; nothing unlisted
                                         sits in src/odia-core/
  sync-odia-core.py --check --origin DIR also verify the recorded ORIGINAL hashes against
                                         `git show <commit>:<path>` in an OpenDIAlyzer checkout,
                                         and say which originals have changed upstream since
  sync-odia-core.py --update             re-record the vendored hashes after a deliberate edit
                                         (describe the edit under "changes" by hand)

Exit status 0 when everything checks, 1 otherwise.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
MANIFEST = ROOT / "src" / "odia-core" / "MANIFEST.json"
SPDX = "SPDX-License-Identifier: BSD-3-Clause"
HEX64 = re.compile(r"^[0-9a-f]{64}$")


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def check(manifest: dict, origin: str | None) -> list[str]:
    errors = []
    for key in ("origin", "originals", "files", "changes"):
        if key not in manifest:
            errors.append(f"manifest has no '{key}'")
    if errors:
        return errors
    commit = manifest["origin"].get("commit", "")
    if not re.fullmatch(r"[0-9a-f]{40}", commit):
        errors.append(f"origin commit '{commit}' is not a full git hash")
    for path, digest in manifest["originals"].items():
        if not HEX64.match(digest):
            errors.append(f"original {path}: '{digest}' is not a sha256")
    if not manifest["changes"]:
        errors.append("manifest lists no changes; a vendored copy that differs from its origin must say how")

    listed = set()
    for entry in manifest["files"]:
        rel = entry["path"]
        listed.add(rel)
        path = ROOT / rel
        if not path.is_file():
            errors.append(f"{rel}: listed but missing")
            continue
        data = path.read_bytes()
        if sha256(data) != entry.get("sha256"):
            errors.append(f"{rel}: sha256 differs from the manifest (edited without "
                          "--update and a 'changes' entry?)")
        head = data[:600].decode("utf-8", "replace")
        if SPDX not in head:
            errors.append(f"{rel}: no '{SPDX}' line in its header")
        for src in entry.get("from", []):
            if src not in manifest["originals"]:
                errors.append(f"{rel}: derived from {src}, which has no recorded original hash")

    for path in sorted((ROOT / "src" / "odia-core").iterdir()):
        rel = path.relative_to(ROOT).as_posix()
        if path.name != "MANIFEST.json" and rel not in listed:
            errors.append(f"{rel}: in src/odia-core/ but not in the manifest")

    if origin:
        for src, digest in manifest["originals"].items():
            try:
                blob = subprocess.run(["git", "-C", origin, "show", f"{commit}:{src}"],
                                      check=True, capture_output=True).stdout
            except (OSError, subprocess.CalledProcessError) as exc:
                errors.append(f"origin {src}@{commit[:12]}: cannot read ({exc})")
                continue
            if sha256(blob) != digest:
                errors.append(f"origin {src}@{commit[:12]}: hash differs from the manifest")
            try:
                now = subprocess.run(["git", "-C", origin, "show", f"HEAD:{src}"],
                                     check=True, capture_output=True).stdout
                if sha256(now) != digest:
                    print(f"note: {src} has changed upstream since {commit[:12]}")
            except (OSError, subprocess.CalledProcessError):
                print(f"note: {src} no longer exists at the origin's HEAD")
    return errors


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--check", action="store_true")
    mode.add_argument("--update", action="store_true")
    parser.add_argument("--origin", help="an OpenDIAlyzer git checkout to verify originals against")
    args = parser.parse_args()

    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    if args.update:
        for entry in manifest["files"]:
            entry["sha256"] = sha256((ROOT / entry["path"]).read_bytes())
        MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        print(f"updated {len(manifest['files'])} hashes in {MANIFEST.relative_to(ROOT)}")
        return 0

    errors = check(manifest, args.origin)
    for error in errors:
        print(f"sync-odia-core: {error}", file=sys.stderr)
    if not errors:
        print(f"sync-odia-core: {len(manifest['files'])} vendored files match "
              f"{manifest['origin']['project']}@{manifest['origin']['commit'][:12]} plus "
              f"{len(manifest['changes'])} recorded changes")
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
