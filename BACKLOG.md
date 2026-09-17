# Backlog

Work that is worth doing and is not scheduled. Anything here should be
actionable: what it is, and what it actually touches.

## The macOS CLI is slow on its first run — measure before packaging

Measured on an M-series Mac against the published v0.9.0 **tarball**: 335 s for
the first invocation, 1 s afterwards. The unsigned v0.2.3 tarball started in
0.4 s. It is not a hang — `sample` shows dyld parked in `Loader::mapSegments`
validating a code signature while syspolicyd logs about one XPC round trip a
second, once per Developer-ID-signed dylib, cached per machine afterwards.

**But the v0.10.0 cask install of the same kind of payload took 33 s, not 335.**
Ten times faster, from a different delivery route, and nobody has explained why.
That gap is the most informative thing here and it has not been chased.

### The obvious fix does not work, and this entry used to recommend it

Shipping the CLI in a signed, notarized, **stapled** `.dmg` or `.pkg` does *not*
help. A stapled ticket travels with the container, not with files copied out of
it: drag an app out of a stapled disk image and the copy has no ticket, which is
the documented behaviour and a common packaging mistake. Our payload is worse
than an app — it is loose files, and `stapler` cannot target a bare Mach-O at
all, only a bundle, a `.dmg` or a `.pkg`. A `.pkg` has the same problem for its
installed payload.

So the only container whose ticket survives the copy is a **bundle**. Wrapping
the whole tree in a `DIALibGen.app` and stapling that would keep the ticket, at
the cost of an odd shape for a command-line tool and a wrapper to reach
`Contents/MacOS/DIALibGen`.

### Do this first, before building anything

1. From a machine that has never run the build in question (Gatekeeper caches
   its verdict per code-directory hash, machine-wide, so a second measurement on
   the same Mac is meaningless), time the first run of 0.10.1 via **both**
   routes: `curl` + `tar`, and `brew install --cask`. That says whether the
   problem is the tarball specifically, and reproduces or kills the 33 s figure.
2. Only then decide. If the cask route is genuinely ~30 s, this is a much smaller
   problem than the headline number suggests and may not be worth an `.app`
   wrapper at all.

Until it is measured, the README says the first run is slow and not stuck.

