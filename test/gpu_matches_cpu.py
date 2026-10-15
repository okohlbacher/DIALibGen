#!/usr/bin/env python3
"""GPU generation: deterministic across CUDA sessions, and equal to the CPU.

usage: gpu_matches_cpu.py DIALibGen FASTA WORKDIR

Skips (77) when no CUDA session can be created. With TF32 off (the default)
the GPU library must select exactly the CPU's precursors and fragments and
agree on every value to float rounding; 1 and 4 CUDA sessions must write the
same bytes.
"""
import csv
import os
import shutil
import subprocess
import sys
from pathlib import Path

tool, fasta, work = sys.argv[1], sys.argv[2], Path(sys.argv[3])
shutil.rmtree(work, ignore_errors=True)
work.mkdir(parents=True)


def generate(name, *args, cpu=False):
    env = dict(os.environ)
    if cpu:
        env["CUDA_VISIBLE_DEVICES"] = ""
    out = work / name
    result = subprocess.run([tool, "-in", fasta, "-out", str(out), "-threads", "4",
                             "-generation:predict_ion_mobility", "true", *args],
                            env=env, capture_output=True, text=True)
    return out, result


gpu1, result = generate("gpu_s1.tsv", "-machine:require_gpu", "-machine:gpu_sessions", "1")
if result.returncode != 0:
    if "the GPU is required" in result.stdout + result.stderr:
        print("no usable CUDA device; skipped")
        sys.exit(77)
    sys.exit(result.stdout + result.stderr)
gpu4, result = generate("gpu_s4.tsv", "-machine:require_gpu", "-machine:gpu_sessions", "4")
assert result.returncode == 0, result.stdout + result.stderr
assert gpu1.read_bytes() == gpu4.read_bytes(), "1 and 4 CUDA sessions wrote different libraries"
cpu, result = generate("cpu.tsv", cpu=True)
assert result.returncode == 0, result.stdout + result.stderr


def rows(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


a, b = rows(cpu), rows(gpu1)
assert len(a) == len(b), f"{len(a)} CPU rows, {len(b)} GPU rows"
# The TSV's DIA-NN column names (iRT, IonMobility, LibraryIntensity).
tolerance = {"iRT": 1e-4, "IonMobility": 1e-4, "CCS": 0.01, "LibraryIntensity": 1e-4}
for x, y in zip(a, b):
    for key in x:
        if key in tolerance and x[key] != "":
            assert abs(float(x[key]) - float(y[key])) <= tolerance[key], (key, x, y)
        elif key not in tolerance:
            assert x[key] == y[key], (key, x, y)
print(f"{len(a)} transitions: GPU equals CPU within float rounding; sessions 1 == 4")
