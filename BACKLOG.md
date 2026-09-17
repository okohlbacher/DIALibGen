# Backlog

Work that is worth doing and is not scheduled. Anything here should be
actionable: what it is, and what it actually touches.

## The macOS CLI takes 5½ minutes on its first run

Measured on an M-series Mac against the published v0.9.0 tarball: **335 s on the
first invocation, 1 s on every one after.** The unsigned v0.2.3 tarball started in
0.4 s. It is not a hang — a `sample` of the process shows dyld parked in
`Loader::mapSegments`, validating a code signature, while syspolicyd logs one XPC
round trip roughly every second. 145 Developer-ID-signed dylibs, evaluated one at
a time, verdict cached afterwards.

So signing made the first run two orders of magnitude slower, and most people will
conclude it has hung and kill it. This is now the worst thing about the macOS
build, and it reaches Homebrew too: the CLI cask installs that same tarball.

**A `.tar.gz` cannot carry a stapled ticket**, so each dylib is checked with Apple
individually. The fix is to ship the CLI in a container that can: a signed,
notarized, **stapled** `.dmg` or `.pkg`, whose ticket covers the whole payload and
is consulted once, offline. FASTag already does exactly this — see its
"Package the CLI as a signed, notarized disk image" step, which is the same problem
solved.

Options, cheapest first:
- Ship `DIALibGen-macos-<arch>.dmg` alongside the tarball and point the
  CLI cask at it. The tarball stays for people who want a plain archive.
- Or a `.pkg`, which installs to a fixed prefix and is what `brew install --cask`
  handles most cleanly.
- Measure before and after. The claim to verify is first-run wall time from a
  fresh download, not that the ticket exists.

Until then, say so in the README: the first run is slow, it is not stuck.
