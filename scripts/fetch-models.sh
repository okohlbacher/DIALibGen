#!/usr/bin/env bash
# Download the three AlphaPeptDeep ONNX models and put them where the tool looks.
#
# Installed as `dialibgen-fetch-models`; in a checkout it is scripts/fetch-models.sh.
#
#   dialibgen-fetch-models                  # install beside the DIALibraryGenerator on PATH
#   dialibgen-fetch-models --dir DIR        # install into DIR instead
#   dialibgen-fetch-models --prefix DIR     # install into DIR/share/DIALibraryGenerator/models
#   dialibgen-fetch-models --check          # verify what is already installed, download nothing
#
# The models are NOT redistributed with this project -- see THIRD-PARTY-NOTICES.md
# and BACKLOG.md -- so this fetches them from OpenMS's archive, which is where a
# source build of OpenMS gets them too.
#
# Every file is checked against a pinned SHA256. A mismatch is deleted rather than
# installed: these drive every prediction the tool makes, and a silently different
# model is a silently different library.
set -euo pipefail

BASE=https://archive.openms.de/openms/models

# name:sha256. Pinned, not "latest": a model that changes underneath a cached
# library would make two runs of the same config disagree with no visible cause.
MODELS="
peptdeep_rt_dynamic.onnx  3f3b847ea37a9333ac163409f3d4b343d021471548091c61aed6cfa398ce62fd
peptdeep_ms2_dynamic.onnx 4141ea233d663b6f5b98ad80b8923ee685dae5cfdbe12243e2e30cd9c44ecf00
peptdeep_ccs_dynamic.onnx 41816f325cdd838897b7e7dcaf291da5d4d3903a1532fc5dc9db7e0725015801
"

DIR=""
PREFIX=""
CHECK=0
while [ $# -gt 0 ]; do
  case "$1" in
    --dir)    DIR=${2:?--dir needs a directory}; shift 2 ;;
    --prefix) PREFIX=${2:?--prefix needs a directory}; shift 2 ;;
    --check)  CHECK=1; shift ;;
    -h|--help) sed -n '2,14p' "$0"; exit 0 ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
done

# shasum on macOS, sha256sum on Linux. Both ship in a base install; neither is
# everywhere, so pick whichever is present rather than requiring one.
if command -v sha256sum > /dev/null 2>&1; then
  digest() { sha256sum "$1" | cut -d' ' -f1; }
elif command -v shasum > /dev/null 2>&1; then
  digest() { shasum -a 256 "$1" | cut -d' ' -f1; }
else
  echo "need sha256sum or shasum to verify the downloads" >&2; exit 1
fi

# Where the tool itself looks, in its own order: $DIALIBGEN_MODEL_DIR first, then
# <the binary>/../share/DIALibraryGenerator/models. Installing into the second
# means nothing has to be set in the environment afterwards.
# The tool resolves share/DIALibraryGenerator RELATIVE TO ITS REAL EXECUTABLE, so
# the prefix has to come from the real path, not from the name that was invoked.
# A Homebrew cask makes this bite: $(brew --prefix)/bin/DIALibraryGenerator is a
# symlink to a wrapper inside the Caskroom, and taking dirname without resolving
# it gives /opt/homebrew/share -- which is writable, so nothing would fail. The
# models would land there, the tool would look in the Caskroom, and this script
# would report success.
#
# readlink -f is not portable to older macOS, so the links are walked by hand.
resolve_link() {
  local p=$1 n=0
  while [ -L "$p" ]; do
    n=$((n + 1))
    [ "$n" -gt 32 ] && { echo "$p"; return; }   # a loop is not worth dying over
    local t
    t=$(readlink "$p")
    case "$t" in
      /*) p=$t ;;
      *)  p=$(dirname "$p")/$t ;;
    esac
  done
  echo "$p"
}

prefix_of() { cd -P "$(dirname "$(resolve_link "$1")")/.." && pwd; }

# BIN is the binary this install is FOR -- the one under --prefix when that was
# given, not whatever happens to be on PATH. Getting that wrong made the script
# tell you to set DIALIBGEN_MODEL_DIR for an installation that needs nothing set.
BIN=""
if [ -z "$DIR" ]; then
  if [ -n "$PREFIX" ]; then
    DIR="$PREFIX/share/DIALibraryGenerator/models"
    BIN="$PREFIX/bin/DIALibraryGenerator"
  elif BIN=$(command -v DIALibraryGenerator 2> /dev/null); then
    DIR="$(prefix_of "$BIN")/share/DIALibraryGenerator/models"
  else
    echo "no DIALibraryGenerator on PATH -- pass --dir or --prefix" >&2
    exit 2
  fi
else
  BIN=$(command -v DIALibraryGenerator 2> /dev/null || true)
fi

# A Homebrew cask stages into the Caskroom, which brew owns and an upgrade
# replaces wholesale. Writing there would work today and vanish on the next
# upgrade, so fall back to the user's own data directory and say so.
if ! mkdir -p "$DIR" 2> /dev/null || ! [ -w "$DIR" ]; then
  FALLBACK="${XDG_DATA_HOME:-$HOME/.local/share}/dialibgen/models"
  echo "cannot write to $DIR -- using $FALLBACK instead"
  DIR=$FALLBACK
  mkdir -p "$DIR"
fi

echo "models: $DIR"
missing=0
got=0

while read -r name want; do
  [ -n "$name" ] || continue
  target="$DIR/$name"

  if [ -f "$target" ]; then
    have=$(digest "$target")
    if [ "$have" = "$want" ]; then
      printf '  %-26s ok\n' "$name"
      continue
    fi
    printf '  %-26s WRONG DIGEST (have %.12s..., want %.12s...)\n' "$name" "$have" "$want"
    [ "$CHECK" = 1 ] && { missing=$((missing + 1)); continue; }
    rm -f "$target"
  elif [ "$CHECK" = 1 ]; then
    printf '  %-26s missing\n' "$name"
    missing=$((missing + 1))
    continue
  fi

  printf '  %-26s downloading...' "$name"
  # To a temporary name first: an interrupted download must not leave something
  # that looks installed. --retry because this is one server and one chance.
  tmp="$target.part"
  if ! curl -fsSL --retry 3 --retry-delay 2 -o "$tmp" "$BASE/$name"; then
    echo " FAILED"
    rm -f "$tmp"
    echo "could not download $BASE/$name" >&2
    exit 1
  fi
  have=$(digest "$tmp")
  if [ "$have" != "$want" ]; then
    echo " BAD DIGEST"
    rm -f "$tmp"
    echo "$name: got $have, expected $want" >&2
    echo "the archive served something other than the pinned model; not installing it" >&2
    exit 1
  fi
  mv "$tmp" "$target"
  echo " ok"
  got=$((got + 1))
done <<EOF
$MODELS
EOF

if [ "$CHECK" = 1 ]; then
  [ "$missing" -eq 0 ] || { echo; echo "$missing model(s) missing or corrupt"; exit 1; }
  echo; echo "all three present and verified"
  exit 0
fi

echo
[ "$got" -gt 0 ] && echo "downloaded $got, $(( 3 - got )) already present" || echo "nothing to do"

# Whether the tool will actually FIND them, which is the only claim worth making.
# It is not the same question as whether the files are on disk: the directory has
# to be one the binary searches.
found=0
if [ -n "$BIN" ] && [ -x "$BIN" ]; then
  want="$(prefix_of "$BIN")/share/DIALibraryGenerator/models"
  [ "$DIR" = "$want" ] && found=1
fi
if [ "$found" = 1 ]; then
  echo "$BIN finds these on its own; nothing to set."
else
  echo "Set this so the tool finds them:"
  echo
  echo "  export DIALIBGEN_MODEL_DIR=$DIR"
fi
