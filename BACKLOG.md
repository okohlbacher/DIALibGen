# Backlog

Work that is worth doing and is not scheduled. Anything here should be
actionable: what it is, and what it actually touches.

## Rename the tool and the repository to DIALibGen

`DIALibraryGenerator` is 21 characters to type at a prompt, and the code already
says `DIALIBGEN` wherever it had to be short — `DIALIBGEN_MODEL_DIR`,
`DIALIBGEN_DATA_DIR`, `DIALIBGEN_VERSION`. The rename makes the tool agree with
itself.

197 occurrences across 49 files, but most are textual. The ones that carry
identity, and the order they have to be done in, are what make this more than a
find-and-replace:

**Breaks on rename, plan for it**

- **Release asset names.** `DIALibraryGenerator-macos-arm64.tar.gz` and the
  twelve others are named in `ci.yml`, `windows.yml`, `release-complete`'s
  assertions, both Homebrew casks and the README. A cask pins a digest to a URL,
  so old releases keep old names for ever and only new ones move.
- **Both cask tokens**, and the tap repository itself
  (`homebrew-dialibrarygenerator`). A changed token is a different cask to
  Homebrew: existing users get no upgrade, they get an orphan. Casks support
  `old_tokens` for exactly this — use it rather than asking people to uninstall
  and reinstall.
- **The desktop app's identity.** `tauri.conf.json` carries
  `productName: DIALibraryGenerator` and `identifier: de.openms.dialibrarygenerator`.
  The identifier is what macOS and the MSI use to decide whether an install is an
  upgrade or a second copy; changing it strands the installed app rather than
  replacing it. Decide deliberately whether to change it at all — the identifier
  need not follow the name.
- **`share/DIALibraryGenerator/`**, which `dataDirCandidates()` searches. An
  installed tree from a previous version is not found by a renamed binary, which
  is fine, but the two must move together.
- **The embedded recipe and the cache fingerprint.** Both name the tool. A
  rename should invalidate old caches rather than silently match them; check
  whether `fingerprintParams` needs another version bump.

**Straightforward**

- `src/DIALibraryGenerator.cpp` and `cmake/DIALibraryGeneratorConfig.cmake.in`
  (filenames), `project()` and the target name, the exported CMake package name
  and therefore every `find_package(DIALibraryGenerator)` downstream.
- `gui/package.json`, `gui/package-lock.json`, `gui/src-tauri/Cargo.toml` and
  `Cargo.lock` (crate `dialibrarygenerator-gui`) — all five files the version
  guard already checks, so it will catch a half-done rename the same way it
  catches a half-done bump.
- The doc URL in the tool description, which points at
  `TOPP_DIALibraryGenerator.html` upstream.
- `README.md`, `CHANGELOG.md` (historical entries stay as they were), the tap's
  README, and the `software-signing` skill's per-repo table.

**Sequence.** Rename the repo last: GitHub redirects the old path, so the tap and
any clone keep working while the rest lands. Cut a release under the new names
before pointing the casks at them, or the updater will compute digests for URLs
that do not exist yet.
