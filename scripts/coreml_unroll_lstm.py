#!/usr/bin/env python3
"""Experiment: LSTM-free PeptDeep RT/CCS graphs for ONNX Runtime's CoreML provider.

CoreML has no LSTM operator, so DIALibGen predicts RT and CCS on the CPU
even with -machine:device coreml. This script rewrites a PeptDeep RT or CCS
model for ONE sequence length: every LSTM node is unrolled into per-step
MatMul/Add/Sigmoid/Tanh/Mul nodes (ONNX's LSTM semantics, gate order
i, o, f, c, no peepholes, no clipping), which the provider can compile. It
then measures, for that length and a batch of 2048 rows:

  1. the CPU difference between the stock graph (ONNX Runtime's LSTM kernel)
     and the unrolled graph, both on the CPU provider;
  2. the difference between the unrolled graph on CoreML and the stock graph
     on the CPU (what the predictor would ship);
  3. the throughput of the stock graph on the CPU (one thread, as DIALibGen
     runs its sessions) and of the unrolled graph on CoreML.

usage: coreml_unroll_lstm.py MODEL.onnx LENGTH [--rows 2048] [--out DIR] [--cache DIR]

Result on an M4 Max (2026-10-13, ONNX Runtime 1.31, 2048 rows, random
inputs; docs/usage.md#apple-gpu and the 0.11.4 changelog): the unrolled
graphs are exact enough -- on the CPU they agree with ONNX Runtime's LSTM
kernel to 4.5e-7 (RT, normalised) and 3.1e-4 A^2 (CCS), on CoreML to
6e-7 and 4.3e-4 A^2, both within the tolerances the MS2 path is held to --
but not fast: CoreML runs the unrolled RT graph at 16.2 k rows/s at 16
tokens and 9.5 k at 32 (CCS 18.5 k and 10.9 k), against 18.9 k and 10.2 k
rows/s for ONE CPU thread with the stock LSTM graph (CCS 18.9 k and
10.1 k), and DIALibGen's CPU path runs one session per core (16 on this
machine). An LSTM is a chain of small dependent steps, which the GPU
cannot batch any better than a core can, and the sessions' Tile nodes
stay on the CPU besides ("Data-dependent shapes were disabled"). So RT
and CCS stay on the CPU, overlapped with the GPU's MS2 phase (0.11.4)
instead.
"""
import argparse
import os
import sys
import time
from pathlib import Path

import numpy as np
import onnx
from onnx import helper, numpy_helper


def unroll(model, length):
    graph = model.graph
    inits = {i.name: i for i in graph.initializer}
    new_nodes = []
    counter = [0]

    def name(base):
        counter[0] += 1
        return f"unrolled/{base}_{counter[0]}"

    def add_init(base, array):
        n = name(base)
        graph.initializer.append(numpy_helper.from_array(np.ascontiguousarray(array), n))
        return n

    for node in graph.node:
        if node.op_type != "LSTM":
            new_nodes.append(node)
            continue
        attrs = {a.name: helper.get_attribute_value(a) for a in node.attribute}
        hidden = int(attrs["hidden_size"])
        direction = attrs.get("direction", b"forward").decode() if isinstance(attrs.get("direction"), bytes) else attrs.get("direction", "forward")
        dirs = 2 if direction == "bidirectional" else 1
        X, Wn, Rn, Bn = node.input[0], node.input[1], node.input[2], node.input[3]
        # PeptDeep's LSTMs start from learned states (inputs 5 and 6,
        # [dirs, batch, H]), not from zeros: one [batch, H] slice per direction.
        H0 = node.input[5] if len(node.input) > 5 and node.input[5] else None
        C0 = node.input[6] if len(node.input) > 6 and node.input[6] else None
        W = numpy_helper.to_array(inits[Wn])   # [dirs, 4H, I]
        R = numpy_helper.to_array(inits[Rn])   # [dirs, 4H, H]
        B = numpy_helper.to_array(inits[Bn])   # [dirs, 8H]
        Y, Y_h, Y_c = (list(node.output) + ["", "", ""])[:3]
        # X is [seq, batch, input]: one slice per step.
        steps = [name("x_t") for _ in range(length)]
        split_x = add_init("split_x", np.ones(length, dtype=np.int64))
        new_nodes.append(helper.make_node("Split", [X, split_x], steps, axis=0, name=name("split_x")))
        step_rows = []   # squeezed [batch, input]
        axes0 = add_init("axes0", np.array([0], dtype=np.int64))
        for s in steps:
            o = name("x_row")
            new_nodes.append(helper.make_node("Squeeze", [s, axes0], [o], name=name("squeeze")))
            step_rows.append(o)
        split_dirs = add_init("split_dirs", np.ones(dirs, dtype=np.int64))
        h0s, c0s = [None] * dirs, [None] * dirs
        if H0 is not None:
            parts = [name("h0_d") for _ in range(dirs)]
            new_nodes.append(helper.make_node("Split", [H0, split_dirs], parts, axis=0, name=name("split_h0")))
            for d, part in enumerate(parts):
                h0s[d] = name("h0")
                new_nodes.append(helper.make_node("Squeeze", [part, axes0], [h0s[d]], name=name("squeeze")))
        if C0 is not None:
            parts = [name("c0_d") for _ in range(dirs)]
            new_nodes.append(helper.make_node("Split", [C0, split_dirs], parts, axis=0, name=name("split_c0")))
            for d, part in enumerate(parts):
                c0s[d] = name("c0")
                new_nodes.append(helper.make_node("Squeeze", [part, axes0], [c0s[d]], name=name("squeeze")))
        y_dirs, h_last, c_last = [], [], []
        for d in range(dirs):
            wt = add_init("WT", W[d].T.astype(np.float32))           # [I, 4H]
            rt = add_init("RT", R[d].T.astype(np.float32))           # [H, 4H]
            bias = add_init("b", (B[d][:4 * hidden] + B[d][4 * hidden:]).astype(np.float32))   # Wb + Rb, [4H]
            split_g = add_init("split_g", np.full(4, hidden, dtype=np.int64))
            order = range(length) if d == 0 else range(length - 1, -1, -1)
            h, c = h0s[d], c0s[d]
            outs = [None] * length
            for t in order:
                xw = name("xw")
                new_nodes.append(helper.make_node("MatMul", [step_rows[t], wt], [xw], name=name("mm_x")))
                pre = name("pre")
                if h is None:
                    new_nodes.append(helper.make_node("Add", [xw, bias], [pre], name=name("add_b")))
                else:
                    hr = name("hr")
                    new_nodes.append(helper.make_node("MatMul", [h, rt], [hr], name=name("mm_h")))
                    s1 = name("s1")
                    new_nodes.append(helper.make_node("Add", [xw, hr], [s1], name=name("add_h")))
                    new_nodes.append(helper.make_node("Add", [s1, bias], [pre], name=name("add_b")))
                gi, go, gf, gc = (name(g) for g in ("gi", "go", "gf", "gc"))
                new_nodes.append(helper.make_node("Split", [pre, split_g], [gi, go, gf, gc], axis=1, name=name("split_g")))
                i_, o_, f_, c_ = (name(g) for g in ("i", "o", "f", "cc"))
                new_nodes.append(helper.make_node("Sigmoid", [gi], [i_], name=name("sig")))
                new_nodes.append(helper.make_node("Sigmoid", [go], [o_], name=name("sig")))
                new_nodes.append(helper.make_node("Sigmoid", [gf], [f_], name=name("sig")))
                new_nodes.append(helper.make_node("Tanh", [gc], [c_], name=name("tanh")))
                ic = name("ic")
                new_nodes.append(helper.make_node("Mul", [i_, c_], [ic], name=name("mul")))
                if c is None:
                    c_new = ic
                else:
                    fc = name("fc")
                    new_nodes.append(helper.make_node("Mul", [f_, c], [fc], name=name("mul")))
                    c_new = name("c")
                    new_nodes.append(helper.make_node("Add", [fc, ic], [c_new], name=name("add_c")))
                tc = name("tc")
                new_nodes.append(helper.make_node("Tanh", [c_new], [tc], name=name("tanh")))
                h_new = name("h")
                new_nodes.append(helper.make_node("Mul", [o_, tc], [h_new], name=name("mul")))
                h, c = h_new, c_new
                outs[t] = h
            # Y[:, d] is [seq, batch, H]: unsqueeze each step and concat on axis 0.
            unsq = []
            for t in range(length):
                u = name("y_t")
                new_nodes.append(helper.make_node("Unsqueeze", [outs[t], axes0], [u], name=name("unsq")))
                unsq.append(u)
            yd = name("y_dir")
            new_nodes.append(helper.make_node("Concat", unsq, [yd], axis=0, name=name("concat_t")))
            yd2 = name("y_dir4")
            axes1 = add_init("axes1", np.array([1], dtype=np.int64))
            new_nodes.append(helper.make_node("Unsqueeze", [yd, axes1], [yd2], name=name("unsq")))   # [seq, 1, batch, H]
            y_dirs.append(yd2)
            hl = name("h_last"); cl = name("c_last")
            new_nodes.append(helper.make_node("Unsqueeze", [h, axes0], [hl], name=name("unsq")))
            new_nodes.append(helper.make_node("Unsqueeze", [c, axes0], [cl], name=name("unsq")))
            h_last.append(hl); c_last.append(cl)
        if Y:
            new_nodes.append(helper.make_node("Concat", y_dirs, [Y], axis=1, name=name("concat_dir")))
        if Y_h:
            new_nodes.append(helper.make_node("Concat", h_last, [Y_h], axis=0, name=name("concat_h")))
        if Y_c:
            new_nodes.append(helper.make_node("Concat", c_last, [Y_c], axis=0, name=name("concat_c")))
    del graph.node[:]
    graph.node.extend(new_nodes)
    # Fix the free dimensions: the unrolled graph is for this length only.
    for v in list(graph.input) + list(graph.output):
        for dim in v.type.tensor_type.shape.dim:
            if dim.dim_param == "seq_length":
                dim.dim_value = length
                dim.dim_param = ""
    del graph.value_info[:]
    return model


def inputs_for(model, rows, length, rng):
    feeds = {}
    for inp in model.graph.input:
        shape = [d.dim_value if d.dim_value else (rows if d.dim_param == "batch_size" else length)
                 for d in inp.type.tensor_type.shape.dim]
        if inp.name in ("input_sequences", "aa_indices"):
            feeds[inp.name] = rng.integers(1, 21, size=shape, dtype=np.int64)
        elif inp.name == "mod_x":
            feeds[inp.name] = np.zeros(shape, dtype=np.float32)
        elif inp.name == "charges":
            feeds[inp.name] = (rng.integers(2, 4, size=shape) * 0.1).astype(np.float32)
        else:
            feeds[inp.name] = np.zeros(shape, dtype=np.float32)
    return feeds


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("model"); ap.add_argument("length", type=int)
    ap.add_argument("--rows", type=int, default=2048)
    ap.add_argument("--out", default=None, help="write the unrolled model here")
    ap.add_argument("--cache", default=None, help="CoreML compiled-model cache directory")
    ap.add_argument("--repeat", type=int, default=5)
    a = ap.parse_args()
    import onnxruntime as ort
    model = onnx.load(a.model)
    unrolled = unroll(onnx.load(a.model), a.length)
    onnx.checker.check_model(unrolled)
    out = Path(a.out) if a.out else Path(a.model).with_name(Path(a.model).stem + f"_unrolled_L{a.length}.onnx")
    onnx.save(unrolled, out)
    print(f"unrolled graph: {out} ({len(unrolled.graph.node)} nodes, {out.stat().st_size / 1e6:.1f} MB)")

    rng = np.random.default_rng(20260803)
    feeds = inputs_for(model, a.rows, a.length, rng)
    so = ort.SessionOptions(); so.intra_op_num_threads = 1; so.inter_op_num_threads = 1
    stock = ort.InferenceSession(a.model, so, providers=["CPUExecutionProvider"])
    cpu_un = ort.InferenceSession(str(out), so, providers=["CPUExecutionProvider"])
    ref = stock.run(None, feeds)[0]
    got = cpu_un.run(None, feeds)[0]
    print(f"CPU stock vs CPU unrolled: max |diff| {np.abs(ref - got).max():.3g} (output mean {np.abs(ref).mean():.3g})")
    t0 = time.perf_counter()
    for _ in range(a.repeat):
        stock.run(None, feeds)
    cpu_s = (time.perf_counter() - t0) / a.repeat
    print(f"CPU stock, 1 thread: {cpu_s * 1e3:.0f} ms per {a.rows} rows = {a.rows / cpu_s:.0f} rows/s per thread")

    if "CoreMLExecutionProvider" not in ort.get_available_providers():
        print("no CoreML provider in this onnxruntime; stopping after the CPU check")
        return
    opts = {"ModelFormat": "MLProgram", "MLComputeUnits": "CPUAndGPU", "RequireStaticInputShapes": "1",
            "EnableOnSubgraphs": "0"}
    if a.cache:
        os.makedirs(a.cache, exist_ok=True)
        opts["ModelCacheDirectory"] = a.cache
    so2 = ort.SessionOptions(); so2.log_severity_level = 3
    so2.add_free_dimension_override_by_name("batch_size", a.rows)
    t0 = time.perf_counter()
    cml = ort.InferenceSession(str(out), so2, providers=[("CoreMLExecutionProvider", opts), "CPUExecutionProvider"])
    print(f"CoreML session built in {time.perf_counter() - t0:.1f} s")
    got2 = cml.run(None, feeds)[0]
    print(f"CPU stock vs CoreML unrolled: max |diff| {np.abs(ref - got2).max():.3g}")
    t0 = time.perf_counter()
    for _ in range(a.repeat):
        cml.run(None, feeds)
    cml_s = (time.perf_counter() - t0) / a.repeat
    print(f"CoreML unrolled: {cml_s * 1e3:.0f} ms per {a.rows} rows = {a.rows / cml_s:.0f} rows/s")


if __name__ == "__main__":
    main()
