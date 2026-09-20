#!/usr/bin/env bash
# Bounded end-to-end fine-tune on a synthetic run, both heads, through the
# MERGED tool: a library and a reference report go in, a refined library comes
# out, and the tuned models are kept so the write-back can still be checked.
#
#   tune_e2e.sh <DIALibGen> <models-dir> <proteins.fasta> [python]
# Exits 77 (ctest SKIP) when python/pyarrow/numpy are unavailable.
#
# The fixture covers MORE precursors than the report identifies. Those extras
# are the whole reason the stage exists -- refinement cannot touch them, only a
# tuned model can -- so a fixture where everything matched would test nothing.
set -u
BIN="$1"; MODELS="$2"; FASTA="$3"; PY="${4:-python3}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fail() { echo "FAIL: $*" >&2; exit 1; }
TMP=$(mktemp -d "${TMPDIR:-/tmp}/dlr-tune-e2e.XXXXXX") || exit 1
trap 'rm -rf "$TMP"' EXIT

"$PY" -c "import pyarrow, numpy" 2>/dev/null || { echo "SKIP: python with pyarrow+numpy needed" >&2; exit 77; }
"$PY" "$HERE/synth_report.py" "$FASTA" "$TMP/report.parquet" --precursors 1600 \
      --library "$TMP/library.tsv" || fail "synthetic fixture"
for head in rt ccs; do
  [ -s "$MODELS/peptdeep_${head}_dynamic.onnx" ] || fail "no stock model for $head in $MODELS"
done

# Foreign mass-only modifications must fail before either head trains.
"$PY" - "$TMP/library.tsv" "$TMP/unencodable.tsv" <<'PY_ENCODING' || exit 1
import csv, sys
with open(sys.argv[1]) as source:
    reader = csv.DictReader(source, delimiter="\t"); fields = reader.fieldnames; rows = list(reader)
first = rows[0]["Precursor.Id"]
for row in rows:
    if row["Precursor.Id"] == first:
        row["Modified.Sequence"] += "[+123.456789]"
with open(sys.argv[2], "w", newline="") as target:
    writer = csv.DictWriter(target, fields, delimiter="\t"); writer.writeheader(); writer.writerows(rows)
PY_ENCODING
if "$BIN" -mode tune -in "$TMP/unencodable.tsv" -ids "$TMP/report.parquet" -out "$TMP/not-written.tsv" \
   -tune_models "$MODELS" -tune_out_models "$TMP/encoding-models" \
   > "$TMP/encoding.log" 2>&1; then fail "unencodable library was accepted"; fi
grep -q 'tuning cannot encode library precursor' "$TMP/encoding.log" || { cat "$TMP/encoding.log"; fail "wrong library encoding error"; }
if find "$TMP/encoding-models" -type f 2>/dev/null | grep -q .; then fail "trained before library encoding validation"; fi

# Missing CCS data and refinement gates must fail before spending any epochs
# on RT or writing kept model artifacts.
"$PY" - "$TMP/report.parquet" "$TMP/no-im.parquet" <<'PY_PREFLIGHT' || exit 1
import sys
import pyarrow.parquet as pq
pq.write_table(pq.read_table(sys.argv[1]).drop(["IM"]), sys.argv[2])
PY_PREFLIGHT
if "$BIN" -mode tune -in "$TMP/library.tsv" -ids "$TMP/no-im.parquet" -out "$TMP/not-written.tsv" \
   -tune_models "$MODELS" -tune_out_models "$TMP/preflight-models" \
   > "$TMP/preflight.log" 2>&1; then fail "missing CCS IM was accepted"; fi
grep -q 'IM' "$TMP/preflight.log" || { cat "$TMP/preflight.log"; fail "wrong preflight error"; }
if find "$TMP/preflight-models" -type f 2>/dev/null | grep -q .; then fail "RT artifacts written before CCS validation"; fi
if "$BIN" -mode refine -in "$TMP/library.tsv" -ids "$TMP/report.parquet" -out "$TMP/not-written.tsv" \
   -tune -tune_models "$MODELS" -tune_out_models "$TMP/gate-models" \
   > "$TMP/gate.log" 2>&1; then fail "missing refinement gates accepted"; fi
if find "$TMP/gate-models" -type f 2>/dev/null | grep -q .; then fail "models trained before refinement validation"; fi

# One invocation: tune both heads, re-predict, refine, write a library.
# -q_global/-q_protein 1 disable the gates whose columns a synthetic report does
# not carry; the precursor gate still applies.
"$BIN" -mode refine -in "$TMP/library.tsv" -ids "$TMP/report.parquet" -out "$TMP/refined.tsv" \
   -q_global 1 -q_protein 1 -no_filter \
   -tune -tune_models "$MODELS" -tune_out_models "$TMP/tuned" \
   -filter:rt_max_minutes 30 \
   -train:epochs 20 -train:warmup 2 -stop:min_epochs 20 -stop:patience 100 \
   -machine:device cpu -machine:threads 2 -machine:seed 20260803 \
   > "$TMP/run.log" 2>&1 || { cat "$TMP/run.log" >&2; fail "DIALibGen -mode refine -tune exited non-zero"; }

# The deliverable is a LIBRARY. This is the claim the merge added.
[ -s "$TMP/refined.tsv" ] || fail "no refined library written"
[ "$(wc -l < "$TMP/refined.tsv")" -gt 1 ] || fail "the refined library has no rows"
[ -s "$TMP/refined.tsv.refine.json" ] || fail "no refine provenance sidecar"
grep -q '"tool"' "$TMP/refined.tsv.refine.json" || fail "the sidecar names no tool"

# Unfiltered observed RT is safe only after all predictions are in the run's minutes.
"$PY" - "$TMP/library.tsv" "$TMP/refined.tsv" "$TMP/report.parquet" <<'PY_RT_UNITS' || exit 1
import csv, json, math, sys
import pyarrow.parquet as pq
def rows(path):
    with open(path) as source:
        return list(csv.DictReader(source, delimiter="\t"))
before, after = rows(sys.argv[1]), rows(sys.argv[2])
observed = {(r["Modified.Sequence"], int(r["Precursor.Charge"])): r["RT"]
            for r in pq.read_table(sys.argv[3]).to_pylist()}
assert len(before) == len(after), "unfiltered RT tuning dropped library transitions"
matched = unseen = 0
for old, new in zip(before, after):
    assert old["Precursor.Id"] == new["Precursor.Id"], "unfiltered RT tuning changed precursor order"
    key = new["Modified.Sequence"], int(new["Precursor.Charge"])
    rt = float(new["RT"])
    assert math.isfinite(rt)
    if key in observed:
        assert math.isclose(rt, observed[key], rel_tol=2e-7, abs_tol=1e-6), "matched RT is not observed minutes"
        matched += 1
    else:
        unseen += not math.isclose(rt, float(old["RT"]), abs_tol=1e-5)
assert matched > 0 and unseen > 0, "fixture must exercise observations and unmatched re-predictions"
provenance = json.load(open(sys.argv[2] + ".refine.json"))
assert provenance["library"]["before"] == provenance["library"]["after"]
assert provenance["library"]["rt_repredicted_in_reference_minutes"] is True
assert provenance["tune"]["rt"]["repredicted"] == provenance["library"]["before"]
assert "reference-run minutes" in provenance["units"]["rt"]
print("ok   unfiltered RT refinement retains the whole library in reference-run minutes")
PY_RT_UNITS

# The write-back guarantees, unchanged, against the kept models.
for head in rt ccs; do
  stock="$MODELS/peptdeep_${head}_dynamic.onnx"
  tuned="$TMP/tuned/peptdeep_${head}_dynamic.onnx"
  [ -s "$tuned" ] || fail "$head: no tuned ONNX kept"
  [ -s "$tuned.tune.json" ] || fail "$head: no provenance sidecar"
  cmp -s "$stock" "$tuned" && fail "$head: the tuned ONNX is byte-identical to the stock one"
  "$PY" - "$tuned.tune.json" "$head" "$tuned" "$stock" <<'PYEOF' || exit 1
import hashlib, json, pathlib, sys
p = json.load(open(sys.argv[1])); head = sys.argv[2]
tuned_bytes = pathlib.Path(sys.argv[3]).read_bytes()
stock_sha256 = hashlib.sha256(pathlib.Path(sys.argv[4]).read_bytes()).hexdigest()
tuned_sha256 = hashlib.sha256(tuned_bytes).hexdigest()
if p["inputs"]["model_in_sha256"] != stock_sha256:
    raise AssertionError("stock model SHA-256 differs from hashlib")
if p["output"]["model_out_sha256"] != tuned_sha256:
    raise AssertionError("tuned model SHA-256 differs from hashlib")
# Read only the protobuf fields needed for metadata: no Python ONNX dependency.
def fields(data):
    i = 0
    def varint():
        nonlocal i
        value = shift = 0
        while True:
            b = data[i]; i += 1; value |= (b & 127) << shift
            if not b & 128: return value
            shift += 7
    while i < len(data):
        key = varint(); number, wire = key >> 3, key & 7
        if wire == 0: value = varint()
        elif wire in (1, 5):
            size = 8 if wire == 1 else 4; value = data[i:i+size]; i += size
        elif wire == 2:
            size = varint(); value = data[i:i+size]; i += size
        else: raise ValueError(f"unsupported protobuf wire {wire}")
        yield number, value
metadata = {}
for field, value in fields(tuned_bytes):
    if field == 14:
        entry = dict(fields(value)); metadata[entry[1].decode()] = entry[2].decode()
embedded = json.loads(metadata["org.openms.dialibgen.training"])
expected = dict(p); expected.pop("output", None)
if embedded != expected: raise AssertionError("embedded training provenance differs from the sidecar")
if embedded.get("tool") != "DIALibGen": raise AssertionError("embedded provenance names a different tool")
c = p["course"]; ev = p["evaluation"]
def die(m): print("FAIL: " + head + ": " + m, file=sys.stderr); sys.exit(1)
if c["updates"] < 1: die("no optimizer updates")
if not c["param_l2_change"] > 0: die("parameters did not move")
key = "calibrated_sd"
s, t = ev["stock"]["val"][key], ev["tuned"]["val"][key]
if not (t < s): die(f"validation {key} did not improve: stock {s:.5f} tuned {t:.5f}")
import csv, math
trajectory = list(csv.DictReader(open(sys.argv[3] + ".trajectory.tsv"), delimiter="\t"))
selected = [r for r in trajectory if r["val_calibrated_sd"]]
best_row = min(selected, key=lambda r: float(r["val_calibrated_sd"]))
if int(best_row["epoch"]) != c["best_epoch"]: die("best epoch does not match the trajectory minimum")
if not math.isclose(float(best_row["val_calibrated_sd"]), t, rel_tol=1e-5, abs_tol=1e-6): die("exported best-checkpoint metric differs from the trajectory")
if p["cohorts"]["test"] < 100 or p["cohorts"]["val"] < 100: die("cohorts too small: " + str(p["cohorts"]))
print(f"ok   {head}: val {key} {s:.4f} -> {t:.4f} (best epoch {c['best_epoch']}) in {c['epochs_run']} epochs, {c['updates']} updates, TEST {ev['stock']['test'][key]:.4f} -> {ev['tuned']['test'][key]:.4f}")
PYEOF
done

# And that the tuning actually REACHED the library. Not that the stage ran --
# that it re-predicted precursors. predictRetentionTimes returns what it could
# NOT do, so a caller reading it as a success count reports 0 when everything
# worked, and every assertion above still passes.
"$PY" - "$TMP/refined.tsv.refine.json" "$MODELS" "$TMP/tuned" <<'PYEOF2' || exit 1
import hashlib, json, pathlib, sys
p = json.load(open(sys.argv[1]))
t = p.get("tune") or {}
def die(m): print("FAIL: " + m, file=sys.stderr); sys.exit(1)
if not t: die("the refine sidecar has no tune section")
for head in ("rt", "ccs"):
    h = t.get(head) or die(f"no {head} section in the tune provenance")
    if h.get("repredicted", 0) < 1:
        die(f"{head}: re-predicted {h.get('repredicted')} precursors -- the stage ran and changed nothing")
    if not h.get("model_sha256") or h["model_sha256"] == h.get("stock_sha256"):
        die(f"{head}: the tuned model hash equals the stock one")
    training = json.loads((pathlib.Path(sys.argv[3]) / f"peptdeep_{head}_dynamic.onnx.tune.json").read_text())
    if h.get("training") != training: die(f"{head}: library lost the complete training recipe")
    for key, directory in (("stock_sha256", sys.argv[2]), ("model_sha256", sys.argv[3])):
        digest = hashlib.sha256((pathlib.Path(directory) / f"peptdeep_{head}_dynamic.onnx").read_bytes()).hexdigest()
        if h[key] != digest: die(f"{head}: refine {key} differs from hashlib")
print("ok   re-predicted rt=%d ccs=%d precursors" % (t["rt"]["repredicted"], t["ccs"]["repredicted"]))
PYEOF2

# Pure tuning must preserve every precursor and fragment, including the library
# precursors absent from the report. RT-only adaptation must preserve mobility.
"$BIN" -mode tune -in "$TMP/library.tsv" -ids "$TMP/report.parquet" -out "$TMP/pure-tuned.tsv" \
   -tune_heads rt -tune_models "$MODELS" -tune_out_models "$TMP/pure-models" \
   -filter:rt_max_minutes 30 \
   -train:epochs 20 -train:warmup 2 -stop:min_epochs 20 -stop:patience 100 \
   -machine:device cpu -machine:threads 2 -machine:seed 20260803 \
   > "$TMP/pure.log" 2>&1 || { cat "$TMP/pure.log" >&2; fail "pure tuning exited non-zero"; }
"$PY" - "$TMP/library.tsv" "$TMP/pure-tuned.tsv" "$TMP/report.parquet" <<'PYEOF3' || exit 1
import csv, json, math, sys
import pyarrow.parquet as pq

def read(path):
    with open(path) as source:
        return list(csv.DictReader(source, delimiter="\t"))
before, after = read(sys.argv[1]), read(sys.argv[2])
assert len(before) == len(after), "pure tuning changed transition count"
# The TSV writer uses canonical numeric formatting and omits Stripped.Sequence.
keys = ("Precursor.Id", "Modified.Sequence", "Protein.Group", "Fragment.Type")
numeric = ("Precursor.Charge", "Precursor.Mz", "Product.Mz", "Relative.Intensity", "Fragment.Charge", "Fragment.Series.Number", "IM", "CCS", "Decoy")
identified = set(pq.read_table(sys.argv[3], columns=["Modified.Sequence"]).column(0).to_pylist())
unseen_changed = 0
for old, new in zip(before, after):
    assert all(old[k] == new[k] for k in keys), "pure tuning changed precursor/fragment identity or order"
    assert all(math.isclose(float(old[k]), float(new[k]), rel_tol=2e-7, abs_tol=1e-7) for k in numeric), "RT-only tuning changed fragment or mobility values"
    assert math.isfinite(float(new["RT"]))
    if old["Modified.Sequence"] not in identified and not math.isclose(float(old["RT"]), float(new["RT"]), abs_tol=1e-5):
        unseen_changed += 1
assert unseen_changed > 0, "unidentified precursors were not re-predicted"
p = json.load(open(sys.argv[2] + ".refine.json"))
assert p["mode"] == "tune" and p["reference"] is None
assert p["library"]["before"] == p["library"]["after"]
assert "ccs" not in p["tune"] and p["library"]["rt_written"] == 0
print("ok   pure RT tuning preserves all keys, fragments, CCS and IM; changes unseen RT")
PYEOF3

# The RT head trained twice from stock with one seed on this CPU build. Compare
# its real checkpoints and trajectory, excluding only timing and output paths.
"$PY" - "$TMP/tuned/peptdeep_rt_dynamic.onnx" "$TMP/pure-models/peptdeep_rt_dynamic.onnx" <<'PY_REPEAT' || exit 1
import csv, json, pathlib, sys
paths = list(map(pathlib.Path, sys.argv[1:]))
provenance = [json.loads(path.with_suffix(path.suffix + '.tune.json').read_text()) for path in paths]
for p in provenance:
    assert p['device'] == 'cpu' and p['recipe']['seed'] == 20260803
    assert p['course']['updates'] > 0 and p['course']['param_l2_change'] > 0
    p.pop('output')
    p['course'].pop('train_seconds'); p['course'].pop('eval_seconds')
assert provenance[0] == provenance[1], 'same CPU seed changed cohorts, training course or evaluation'
trajectories = []
for path in paths:
    with open(str(path) + '.trajectory.tsv') as source:
        rows = list(csv.DictReader(source, delimiter='\t'))
    for row in rows:
        row.pop('train_s'); row.pop('eval_s')
    trajectories.append(rows)
assert trajectories[0] == trajectories[1], 'same CPU seed changed the numeric epoch trajectory'
# ONNX metadata includes timing, so compare the complete non-metadata wire
# fields, including every initializer byte, rather than the whole-file hash.
def without_metadata(data):
    i = 0
    def varint():
        nonlocal i
        value = shift = 0
        while True:
            b = data[i]; i += 1; value |= (b & 127) << shift
            if not b & 128: return value
            shift += 7
    result = []
    while i < len(data):
        start = i; key = varint(); wire = key & 7
        if wire == 0: varint()
        elif wire in (1, 5): i += 8 if wire == 1 else 4
        elif wire == 2:
            length = varint(); i += length
        else: raise ValueError(f'unsupported protobuf wire {wire}')
        assert i <= len(data)
        if key >> 3 != 14: result.append(data[start:i])
    return b''.join(result)
assert without_metadata(paths[0].read_bytes()) == without_metadata(paths[1].read_bytes()), 'same CPU seed changed the trained ONNX weights'
print('ok   same CPU seed reproduces cohorts, trajectory, evaluation and trained ONNX weights')
PY_REPEAT

echo "PASSED"
