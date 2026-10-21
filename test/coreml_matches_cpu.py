#!/usr/bin/env python3
"""Apple GPU generation: the CoreML library equals the CPU library.

usage: coreml_matches_cpu.py DIALibGen FASTA WORKDIR

Skips (77) when the ONNX Runtime has no CoreML provider (every platform but
macOS) or the MS2 model cannot run on it. On CoreML the MS2 model runs on the
GPU in float32 with static shapes and padded batches, so its spectra agree
with the CPU's to float rounding rather than exactly: measured on
example/proteins.fasta, intensities within 1.3e-6 and no fragment set
differs. RT and CCS stay on the CPU (CoreML has no LSTM operator) and must
be bit-identical. Two CoreML runs, with 1 and 4 concurrent batches, must
write the same bytes, and -machine:require_gpu must refuse the RT fallback.
RT and CCS run before MS2 by default and, with ODIA_OVERLAP=1, on a thread of
their own while the GPU predicts MS2: the same bytes either way.
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


OVERLAP_LINE = "while the GPU predicts MS2"


def generate(name, *args, overlap=None):
    env = dict(os.environ)
    env["CUDA_VISIBLE_DEVICES"] = ""
    env.pop("ODIA_OVERLAP", None)
    env.pop("ODIA_STREAM", None)  # the overlap needs the streaming writer
    if overlap is not None:
        env["ODIA_OVERLAP"] = overlap
    out = work / name
    result = subprocess.run([tool, "-in", fasta, "-out", str(out), "-threads", "4",
                             "-generation:predict_ion_mobility", "true", *args],
                            env=env, capture_output=True, text=True)
    return out, result


coreml4, result = generate("coreml_c4.tsv", "-machine:device", "coreml", "-machine:gpu_sessions", "4")
log = result.stdout + result.stderr
if result.returncode != 0:
    # A negative code is the signal that killed it: Core ML took the process
    # down in spite of the probe child, which is the bug to report.
    sys.exit(f"CoreML generate failed with return code {result.returncode}\n" + log)
if "[MS2] provider CoreML" not in log:
    # The probe child or the placement check kept MS2 on the CPU (a Mac or a
    # virtual machine whose Core ML cannot run the model): nothing to compare.
    reason = [line for line in log.splitlines() if "MS2" in line and "WARNING" in line]
    print("MS2 does not run on CoreML; skipped:", reason[0] if reason else "no CoreML provider line")
    sys.exit(77)
assert "WARNING: RT: " in log and "LSTM" in log, "the RT fallback (no LSTM in CoreML) was not reported"
assert OVERLAP_LINE not in log, "RT and CCS overlapped MS2 on CoreML without ODIA_OVERLAP=1"
overlapped, result = generate("coreml_overlap.tsv", "-machine:device", "coreml", "-machine:gpu_sessions", "4",
                              overlap="1")
assert result.returncode == 0, result.stdout + result.stderr
assert OVERLAP_LINE in result.stdout + result.stderr, "ODIA_OVERLAP=1 did not overlap RT and CCS with MS2"
assert overlapped.read_bytes() == coreml4.read_bytes(), "ODIA_OVERLAP=1 wrote a different library"
coreml1, result = generate("coreml_c1.tsv", "-machine:device", "coreml", "-machine:gpu_sessions", "1")
assert result.returncode == 0, result.stdout + result.stderr
assert coreml1.read_bytes() == coreml4.read_bytes(), "1 and 4 concurrent CoreML batches wrote different libraries"
_, result = generate("require.tsv", "-machine:device", "coreml", "-machine:require_gpu")
assert result.returncode != 0 and "the GPU is required" in result.stdout + result.stderr, \
    "-machine:require_gpu accepted the RT fallback"
cpu, result = generate("cpu.tsv")
assert result.returncode == 0, result.stdout + result.stderr


def rows(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


a, b = rows(cpu), rows(coreml4)
assert len(a) == len(b), f"{len(a)} CPU rows, {len(b)} CoreML rows"
# The TSV's DIA-NN column names. RT, mobility and CCS come from the CPU path
# on both sides; the intensities from CoreML's float32 GPU math.
tolerance = {"iRT": 0.0, "IonMobility": 0.0, "CCS": 0.0, "LibraryIntensity": 1e-4}
for x, y in zip(a, b):
    for key in x:
        if key in tolerance and x[key] != "":
            assert abs(float(x[key]) - float(y[key])) <= tolerance[key], (key, x, y)
        elif key not in tolerance:
            assert x[key] == y[key], (key, x, y)
print(f"{len(a)} transitions: CoreML equals CPU within float rounding; concurrent 1 == 4; overlap == in order")
