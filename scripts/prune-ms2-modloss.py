#!/usr/bin/env python3
"""Derive peptdeep_ms2_pruned.onnx from peptdeep_ms2_dynamic.onnx.

The MS2 model ends in two heads concatenated on the channel axis: b/y ions
without a loss (channels 0-3) and with a modification's neutral loss (4-7).
DIALibGen reads only channels 0-3. This removes the modloss branch and puts
exact zeros in channels 4-7 (ConstantOfShape of the no-loss head's shape), so
the output stays [batch, nAA - 1, 8] and DIALibGen's shape checks and Spectrum
layout are unchanged. About a fifth fewer FLOPs per prediction.

usage: prune-ms2-modloss.py <peptdeep_ms2_dynamic.onnx> <peptdeep_ms2_pruned.onnx>

Before writing it verifies with ONNX Runtime on the CPU that channels 0-3 are
bit-identical to the full model and 4-7 are zero, for nAA 7..30 at batch 1, 2
and 64. The shipped file was made with onnx 1.23.0, onnxruntime 1.30.0 and
numpy 2.5.3 from the full model with SHA256
4141ea233d663b6f5b98ad80b8923ee685dae5cfdbe12243e2e30cd9c44ecf00
(scripts/fetch-models.sh pins both).
"""
import sys

import numpy as np
import onnx
import onnxruntime as ort
from onnx import helper, TensorProto

src, dst = sys.argv[1], sys.argv[2]
m = onnx.load(src)
g = m.graph
out_name = g.output[0].name
prod = {o: n for n in g.node for o in n.output}
node = prod[out_name]
while node.op_type != 'Concat':
    node = prod[node.input[0]]
concat = node
axis = [a.i for a in concat.attribute if a.name == 'axis'][0]
keep = [i for i in concat.input if 'modloss' not in i and 'modloss' not in prod[i].name]
drop = [i for i in concat.input if i not in keep]
assert len(keep) == 1 and len(drop) == 1, list(concat.input)
print('concat', concat.name, 'axis', axis, 'keep', keep[0], 'drop', drop[0])
shape_out, zeros_out = keep[0] + '__shape', keep[0] + '__zeros'
new_nodes = [
    helper.make_node('Shape', [keep[0]], [shape_out], name='prune_shape'),
    helper.make_node('ConstantOfShape', [shape_out], [zeros_out], name='prune_zeros',
                     value=helper.make_tensor('zero', TensorProto.FLOAT, [1], [0.0])),
]
for k, name in enumerate(concat.input):
    if name == drop[0]:
        concat.input[k] = zeros_out
# insert the new nodes just before the concat (topological order)
nodes = list(g.node)
idx = nodes.index(concat)
nodes[idx:idx] = new_nodes
prod = {o: n for n in nodes for o in n.output}
needed, stack = set(), [o.name for o in g.output]
while stack:
    t = stack.pop()
    n = prod.get(t)
    if n is not None and id(n) not in needed:
        needed.add(id(n))
        stack.extend(i for i in n.input if i)
kept = [n for n in nodes if id(n) in needed]
del g.node[:]
g.node.extend(kept)
used = {i for n in g.node for i in n.input}
inits = [i for i in g.initializer if i.name in used]
before = sum(int(np.prod(i.dims)) for i in m.graph.initializer)
del g.initializer[:]
g.initializer.extend(inits)
after = sum(int(np.prod(i.dims)) for i in g.initializer)
onnx.checker.check_model(m)
print('nodes', len(nodes), '->', len(kept), 'params', before, '->', after)

so = ort.SessionOptions()
so.intra_op_num_threads = 1
so.graph_optimization_level = ort.GraphOptimizationLevel.ORT_ENABLE_ALL
a = ort.InferenceSession(src, so, providers=['CPUExecutionProvider'])
b = ort.InferenceSession(m.SerializeToString(), so, providers=['CPUExecutionProvider'])
rng = np.random.default_rng(7)
names = [i.name for i in a.get_inputs()]
for nAA in range(7, 31):
    for B in (1, 2, 64):
        L = nAA + 2
        aa = rng.integers(1, 21, (B, L)).astype(np.int64)
        aa[:, 0] = 0
        aa[:, -1] = 0
        mod = np.zeros((B, L, 109), np.float32)
        feed = {names[0]: aa, names[1]: mod,
                names[2]: rng.integers(2, 5, (B, 1)).astype(np.float32) / 10,
                names[3]: np.full((B, 1), 0.3, np.float32),
                names[4]: np.zeros((B,), np.int64)}
        ya = a.run(None, feed)[0]
        yb = b.run(None, feed)[0]
        assert ya.shape == yb.shape, (ya.shape, yb.shape)
        assert np.array_equal(ya[..., :4], yb[..., :4]), (nAA, B, np.max(np.abs(ya[..., :4] - yb[..., :4])))
        assert not np.any(yb[..., 4:]), (nAA, B)
print('verified: channels 0-3 bit-identical, 4-7 zero, nAA 7-30 x batch 1/2/64')
onnx.save(m, dst)
print('wrote', dst)
