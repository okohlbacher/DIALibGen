#!/usr/bin/env python3
"""Exercise all three modes of an installed or relocated binary, including CPU training."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import xml.etree.ElementTree as ET

p = argparse.ArgumentParser()
p.add_argument("binary", type=Path)
p.add_argument("--bare", action="store_true", help="run the tool without the build environment")
p.add_argument("--output", type=Path, required=True)
a = p.parse_args()
binary = a.binary.resolve()
out = a.output.resolve()
out.mkdir(parents=True, exist_ok=True)
root = Path(__file__).resolve().parent.parent
env = os.environ.copy()
if a.bare:
    env = {k: os.environ[k] for k in ("HOME", "SystemRoot", "WINDIR", "TEMP", "TMP") if k in os.environ}
    if os.name == "nt":
        env["PATH"] = str(binary.parent) + os.pathsep + str(Path(os.environ["SystemRoot"]) / "System32")
    else:
        env["PATH"] = "/usr/bin:/bin"
env["OPENMS_DISABLE_UPDATE_CHECK"] = "ON"


def run(name, args, allowed=(0,)):
    result = subprocess.run([str(binary), *map(str, args)], env=env, capture_output=True, text=True, timeout=900)
    (out / f"{name}.log").write_text(result.stdout + result.stderr)
    if result.returncode not in allowed:
        sys.exit(f"{name}: exit {result.returncode}\n{result.stdout}{result.stderr}")
    return result


run("help", ["--helphelp"])
run("ini", ["-write_ini", out / "DIALibGen.ini"])
items = {el.get("name"): el for el in ET.parse(out / "DIALibGen.ini").iter("ITEM")}
assert all(mode in items["mode"].get("restrictions", "") for mode in ("generate", "refine", "tune")), "INI lacks the three modes"
assert "tune_out_models" in items, "training options missing from the release"
for mode in ("generate", "refine", "tune"):
    run(f"config-{mode}", ["-mode", mode, "-write_config", out / f"{mode}.json"])
run("generate", ["-mode", "generate", "-in", root / "example/proteins.fasta", "-out", out / "generated.tsv", "-threads", "2"])
assert len((out / "generated.tsv").read_text().splitlines()) > 1
# Python creates the fixture and checks provenance; the executable receives a
# bare environment and never calls Python.
subprocess.run([sys.executable, str(root / "test/refinement/synth_report.py"), "-", str(out / "report.parquet"), "--precursors", "1600", "--library", str(out / "library.tsv")], check=True)
common = ["-in", out / "library.tsv", "-ids", out / "report.parquet"]
run("refine", ["-mode", "refine", *common, "-q_global", "1", "-q_protein", "1", "-out", out / "refined.tsv", "-write_im"])
assert len((out / "refined.tsv").read_text().splitlines()) > 1
assert json.loads((out / "refined.tsv.refine.json").read_text())["tool"] == "DIALibGen"
for head in ("rt", "ccs"):
    models = out / head
    result = run(f"tune-{head}", ["-mode", "tune", *common, "-out", out / f"tuned-{head}.tsv", "-tune_heads", head, "-tune_out_models", models,
        "-filter:rt_max_minutes", "30", "-train:epochs", "20", "-train:warmup", "2", "-stop:min_epochs", "20", "-machine:threads", "2", "-machine:device", "cpu"])
    model = models / f"peptdeep_{head}_dynamic.onnx"
    provenance = json.loads(Path(str(model) + ".tune.json").read_text())
    course = provenance["course"]
    assert course["updates"] > 0 and course["epochs_run"] == 20, f"{head}: no optimizer work"
    assert provenance["device"] == "cpu"
    assert course["exported"] and model.stat().st_size > 0, f"{head}: no tuned model"
    assert len((out / f"tuned-{head}.tsv").read_text().splitlines()) > 1
print("PASS: generate, refine and both CPU training heads work" + (" without the build environment" if a.bare else ""))
