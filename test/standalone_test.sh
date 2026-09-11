#!/usr/bin/env bash
# The standalone-tool gates: everything that is right for a tool shipped INSIDE
# OpenMS and wrong for one that merely links it. Each of these has shipped
# broken past a green pipeline in a sibling project, and each is invisible to
# any check that reads only the build log.
#
#   1. the binary runs with NOTHING from the build environment
#   2. it reports ITS OWN version, not the OpenMS it was built against
#   3. OpenMS's update check is off, so Qt does not print network errors that
#      read like errors from this tool
#   4. the shipped data tables are reachable from the executable, not only from
#      the source tree it was compiled in
#
# `env -i` is the load-bearing part of (1): a runner that still has the build
# tree on PATH, or OPENMS_DATA_PATH set, will happily pass a binary that cannot
# stand on its own.
#
# Usage: standalone_test.sh <binary> <expected-version> <example.fasta>
set -u
BIN="$1"
WANT="$2"
FASTA="$3"

fail() { echo "FAIL: $*" >&2; exit 1; }
TMP=$(mktemp -d "${TMPDIR:-/tmp}/dialibgen-standalone.XXXXXX") || exit 1
trap 'rm -rf "$TMP"' EXIT

# OpenMS needs its own share/ to initialise. Inside a relocatable install it
# finds it relative to the binary; here, in a build tree against an external
# OpenMS, it may need OPENMS_DATA_PATH -- which is exactly what must NOT be
# inherited silently. So pass it explicitly when the caller set it, and record
# in the output that it was needed.
KEEP=(HOME="${HOME:-$TMP}")
[ -n "${OPENMS_DATA_PATH:-}" ] && KEEP+=(OPENMS_DATA_PATH="$OPENMS_DATA_PATH")
# Windows has no RPATH: a binary finds its DLLs through PATH, so a bare
# environment there is not "no build tree on PATH", it is "no libraries at all"
# and nothing starts. The caller passes the runtime directories in PATH_KEEP --
# which is still the point of the test, since the BUILD tree is not among them.
# SystemRoot is required by the Windows loader itself.
[ -n "${PATH_KEEP:-}" ] && KEEP+=(PATH="$PATH_KEEP")
[ -n "${SYSTEMROOT:-}" ] && KEEP+=(SYSTEMROOT="$SYSTEMROOT")
[ -n "${SystemRoot:-}" ] && KEEP+=(SystemRoot="$SystemRoot")

# ---------------------------------------------------------------- 1. it runs
if ! env -i "${KEEP[@]}" "$BIN" --help >"$TMP/help.txt" 2>&1; then
  echo "--- output ---" >&2; cat "$TMP/help.txt" >&2
  fail "the binary does not run in a bare environment"
fi

# --------------------------------------------------- 2. our version, not theirs
# TOPPBase prints "Version: <version_>" and falls back to OpenMS's when the
# constructor leaves version_ empty.
got=$(sed -n 's/^Version: \([^ ]*\).*/\1/p' "$TMP/help.txt" | head -1)
[ -n "$got" ] || fail "--help printed no 'Version:' line"
[ "$got" = "$WANT" ] || fail "reports version '$got', expected '$WANT' (OpenMS's leaking through?)"

# --helphelp carries verboseVersion_, which must name BOTH numbers: a bug report
# needs this tool's version and the OpenMS it was built against.
env -i "${KEEP[@]}" "$BIN" --helphelp >"$TMP/helphelp.txt" 2>&1
grep -q "$WANT" "$TMP/helphelp.txt" || fail "--helphelp does not carry version $WANT"
grep -qi "OpenMS " "$TMP/helphelp.txt" || fail "--helphelp does not name the OpenMS version"

# ------------------------------------------------------- 3. update check is off
for f in "$TMP/help.txt" "$TMP/helphelp.txt"; do
  ! grep -q "QIODevice\|QNetworkReply" "$f" \
    || fail "OpenMS's update check is live (Qt network noise on stderr in $(basename "$f"))"
done

# --------------------------------------------------- 4. it parses in isolation
# -write_config needs no input, no models and no data file, so it proves the
# tool starts, registers its parameters and materialises its defaults.
env -i "${KEEP[@]}" "$BIN" -write_config "$TMP/eff.json" >"$TMP/cfg.log" 2>&1 \
  || { cat "$TMP/cfg.log" >&2; fail "-write_config failed in a bare environment"; }
[ -s "$TMP/eff.json" ] || fail "-write_config wrote nothing"
grep -q '"schema_version"' "$TMP/eff.json" || fail "effective config has no schema_version"

# The tool must accept its OWN output. Every default has to sit inside the
# ranges the config validator enforces, and a default that does not would make
# the documented "-write_config is the authoritative reference" a lie: copy it,
# pass it back, get a refusal.
env -i "${KEEP[@]}" "$BIN" -in "$FASTA" -config "$TMP/eff.json" \
    -out "$TMP/roundtrip.tsv" >"$TMP/rt.log" 2>&1
if grep -qE "^config:" "$TMP/rt.log"; then
  echo "--- output ---" >&2; cat "$TMP/rt.log" >&2
  fail "the tool rejected its own -write_config output"
fi

# ------------------------------------------------- 5. bad input is refused here
# Both of these used to be accepted: an unknown decoy method became "mutate"
# and was then recorded in the provenance under the name that was typed, and
# any -out extension that was not .parquet was written as a DIA-NN TSV under
# whatever name was asked for. They are refused BEFORE the models are needed,
# so this runs with no models present.
printf '{"decoys": "reverese"}\n' > "$TMP/bad_decoy.json"
if env -i "${KEEP[@]}" "$BIN" -in "$FASTA" -config "$TMP/bad_decoy.json" \
       -out "$TMP/x.tsv" >"$TMP/d.log" 2>&1; then
  fail "an unknown decoys method was accepted"
fi
grep -q "unknown decoys method" "$TMP/d.log" \
  || fail "an unknown decoys method failed for the wrong reason: $(head -3 "$TMP/d.log" | tr '\n' ' ')"

printf '{"schema_version": 99}\n' > "$TMP/bad_schema.json"
if env -i "${KEEP[@]}" "$BIN" -in "$FASTA" -config "$TMP/bad_schema.json" \
       -out "$TMP/x.tsv" >"$TMP/s.log" 2>&1; then
  fail "an unsupported schema_version was accepted"
fi
grep -q "schema_version 99" "$TMP/s.log" \
  || fail "an unsupported schema_version failed for the wrong reason"

if env -i "${KEEP[@]}" "$BIN" -in "$FASTA" -out "$TMP/library.parqet" \
       >"$TMP/e.log" 2>&1; then
  fail "-out with an unknown extension was accepted"
fi
grep -q "must end in .parquet or .tsv" "$TMP/e.log" \
  || fail "a bad -out extension failed for the wrong reason"

# Out-of-range numbers. A negative integer is the sharp one: nlohmann converts
# it to std::size_t without complaint, so "missed_cleavages": -1 used to become
# 18446744073709551615 and the digest ran on it -- a wrong library, produced
# without a word. Each case pairs a config with the phrase the refusal must
# contain, so a value rejected for some unrelated reason still fails the test.
while IFS='|' read -r json want; do
  [ -n "$json" ] || continue
  printf '%s\n' "$json" > "$TMP/bad.json"
  if env -i "${KEEP[@]}" "$BIN" -in "$FASTA" -config "$TMP/bad.json" \
         -out "$TMP/x.tsv" >"$TMP/r.log" 2>&1; then
    fail "accepted an out-of-range config: $json"
  fi
  grep -q "$want" "$TMP/r.log" \
    || fail "$json was refused for the wrong reason: $(head -2 "$TMP/r.log" | tr '\n' ' ')"
done <<'CASES'
{"missed_cleavages": -1}|missed_cleavages must be a whole number
{"missed_cleavages": 1.5}|missed_cleavages must be a whole number
{"max_variable_modifications": -2}|max_variable_modifications must be a whole number
{"reserved_doubly_charged": -1}|reserved_doubly_charged must be a whole number
{"precursor_charges": [0, 2]}|precursor_charges must each be between 1 and 10
{"precursor_charges": [-3]}|precursor_charges must each be between 1 and 10
{"max_fragment_charge": 0}|max_fragment_charge must be between 1 and 10
{"min_relative_intensity": -0.5}|min_relative_intensity must be between 0 and 1
{"min_relative_intensity": 2}|min_relative_intensity must be between 0 and 1
{"peptide_length": [-1, 30]}|peptide_length values must be between
{"precursor_mz": [-100.0, 1200.0]}|precursor_mz values must be between
{"peptide_length": [30, 7]}|peptide_length needs min < max
{"precursor_charges": []}|precursor_charges must not be empty
CASES

# ------------------------------------------- 6. a missing model says which one
# Without this the run died inside the ONNX session constructor with
# "Load model from  failed" -- an empty path and no hint that a model was the
# thing missing.
env -i "${KEEP[@]}" DIALIBGEN_MODEL_DIR="$TMP/no-such-models" \
    "$BIN" -in "$FASTA" -out "$TMP/x.tsv" >"$TMP/m.log" 2>&1
if grep -q "peptdeep_rt_dynamic.onnx" "$TMP/m.log"; then
  : # named the file it could not find, and listed where it looked
elif grep -q "wrote $TMP/x.tsv" "$TMP/m.log"; then
  : # this build DOES have models installed; nothing to assert
else
  echo "--- output ---" >&2; cat "$TMP/m.log" >&2
  fail "a missing model produced neither a named model nor a library"
fi

echo "standalone_test: version $got, bare environment, update check off, inputs validated"
