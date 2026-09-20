# AppImage helper sources and runtime replacement

The companion source archive contains the complete upstream sources, notices,
build scripts and patches for the separately built AppRun helper and AppImage
runtime. `appimage-helper-pins.json` records the exact revisions and SHA-256
checksums. The package inventory also records the final runtime prefix hash.

AppRun is the MIT-licensed `src/AppRun.c` from AppImageKit revision
`5735cc5bed206497cddfbd2a75e1982c2606c35d`. Its original binary is byte-identical
to Tauri's `apprun-old` mirror. It uses the host C library dynamically. AppRun's
source is self-contained and can be rebuilt with a C compiler, for example:

```sh
cc -std=gnu99 -O2 -D_FILE_OFFSET_BITS=64 src/AppRun.c -o AppRun
```

The runtime is revision `75849dce7cc37e4319b633df1f116ca895c71a12` of
[AppImage/type2-runtime](https://github.com/AppImage/type2-runtime/tree/75849dce7cc37e4319b633df1f116ca895c71a12).
The original [build log](https://github.com/AppImage/type2-runtime/actions/runs/28063784345)
identifies the Alpine image digest and package versions recorded in the pins.
The runtime links libfuse 3.15.0, squashfuse 0.5.2, musl 1.2.5-r11,
zstd 1.5.6-r2, zlib 1.3.2-r0 and mimalloc2 2.1.7-r0 statically.
libfuse's library files are LGPL-2.1; its separate tools use GPL-2.0.
Both original license texts accompany the complete libfuse source.

To rebuild or modify the runtime, extract `type2-runtime.tar.gz` and the other
source archives. Follow its `BUILD.md`, `scripts/docker/Dockerfile`,
`scripts/common/install-dependencies.sh` and `src/runtime/Makefile`. Use the
recorded Alpine image digest. The accompanying `alpine-aports.tar.gz` contains
the original package recipes and distribution patches under `main/musl`,
`main/zstd`, `main/zlib` and `community/mimalloc2`; their `APKBUILD` files give
the compiler settings and source checksums used for those static libraries.
The runtime source includes the applied libfuse patch at
`patches/libfuse/mount.c.diff`. All upstream archives referenced by these four
package recipes and by the runtime dependency script are included alongside
their recipes, so the original source can also be used when distribution
package revisions leave public binary mirrors.

For a modified libfuse, apply that patch to the supplied libfuse source, make
your changes, then build and install its static library in the runtime build
environment with Meson/Ninja as shown in `install-dependencies.sh`. Rebuild
the runtime with `scripts/build-runtime.sh`. It invokes the supplied Makefile
and produces `out/runtime-x86_64` or `out/runtime-aarch64`; the complete runtime source is included,
so replacing libfuse and relinking does not depend on a proprietary object.

You can run the GUI without the runtime by extracting the original AppImage:

```sh
./DIALibGen.AppImage --appimage-extract
./squashfs-root/AppRun
```

To use a rebuilt runtime, keep the extracted payload and repackage it with
AppImage's `appimagetool --runtime-file /path/to/out/runtime-x86_64 squashfs-root rebuilt.AppImage`.
The replacement image has a different checksum and any original signature
no longer applies. DIALibGen imposes no signature check on this replacement.
The AppRun wrapper and GTK environment hook are attributed to their exact linuxdeploy source templates; generated paths are validated separately. The original runtime code is verified during packaging; only its documented
image-specific checksum/update/signature data sections may differ from the
pinned upstream binary.
