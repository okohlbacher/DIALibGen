#!/usr/bin/env python
"""A reference PeptDeep encoder and RT predictor, written from the spec.

This exists to cross-check ODIA's C++ implementation. It is written from
`doc/04-peptdeep-encoding.md` only -- not from the C++ -- so agreement between
the two is evidence that the spec was implemented, rather than evidence that one
implementation copied the other.

That distinction has mattered repeatedly here: OpenMS's own PeptDeep binding
takes its constants from the same upstream files, so agreeing with it would only
show a shared reading. It is also unusable from outside OpenMS -- the headers are
compiled into libOpenMS.so but never installed.

Usage:
  peptdeep_reference.py encode <modified-sequence>      # dump the tensors
  peptdeep_reference.py rt <model.onnx> <seq> [seq...]  # predict iRT
"""

import json
import os
import sys

import numpy as np

# --- the 109 mod_elements, in yaml order; index IS the feature position -------
# Read from the extracted data file rather than typed: a hand-typed copy came
# out with 112 entries and diverged at index 38.
_ELEMENTS_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                              "..", "data", "peptdeep_mod_elements.txt")
MOD_ELEMENTS = [ln.strip() for ln in open(_ELEMENTS_FILE)
                if ln.strip() and not ln.startswith("#")]
assert len(MOD_ELEMENTS) == 109, f"expected 109 elements, got {len(MOD_ELEMENTS)}"

# Compositions for the modifications the fixtures use, as UniMod gives them.
COMPOSITION = {
    "Acetyl": {"H": 2, "C": 2, "O": 1},
    "Carbamidomethyl": {"H": 3, "C": 2, "N": 1, "O": 1},
    "Oxidation": {"O": 1},
    "GG": {"H": 6, "C": 4, "N": 2, "O": 2},
    "Deamidated": {"H": -1, "N": -1, "O": 1},          # counts are signed
    "Phospho": {"H": 1, "O": 3, "P": 1},
    "Amidated": {"H": 1, "N": 1, "O": -1},             # a C-terminal modification
    "Dehydrated": {"H": -2, "O": -1},
    "2": {"H": 1, "N": 1, "O": -1},
    "23": {"H": -2, "O": -1},
    "7": {"H": -1, "N": -1, "O": 1},
    "1": {"H": 2, "C": 2, "O": 1},
    "4": {"H": 3, "C": 2, "N": 1, "O": 1},
    "35": {"O": 1},
    "21": {"H": 1, "O": 3, "P": 1},
    "121": {"H": 6, "C": 4, "N": 2, "O": 2},
}

# Applied by the CCS/MS2 models, which this reference does not yet drive.
CHARGE_SCALE = 0.1
NCE_SCALE = 0.01


def mod_vector(name):
    """Composition -> 109-vector. Unknown elements accumulate in the '?' slot."""
    v = np.zeros(len(MOD_ELEMENTS), dtype=np.float32)
    for element, count in COMPOSITION[name].items():
        if element in MOD_ELEMENTS:
            v[MOD_ELEMENTS.index(element)] = count      # assigned, not added
        else:
            v[-1] += count                              # accumulated
    return v


def parse(seq):
    """(residues, [(site, mod_name)]).

    Sites follow AlphaPeptDeep: 0 is the peptide N-terminus, 1..nAA are
    residues, and -1 is the peptide C-terminus (row nAA+1).

    OpenMS's AASequence::toString marks a terminal modification with a '.'
    immediately against the bracket -- ".(Acetyl)PEPTIDEK" and
    "PEPTIDER.(Amidated)" -- which is the only thing distinguishing a C-terminal
    modification from one on the last residue. Ignoring the marker put an
    Amidated C-terminus on residue 8 instead of row 9, and the resulting iRT
    error (0.0495) was larger than the modification's entire true effect
    (0.0429), silently.
    """
    residues, mods, i = [], [], 0
    while i < len(seq):
        c = seq[i]
        if c in "([":
            close = ")" if c == "(" else "]"
            j = seq.index(close, i)
            name = seq[i + 1:j]
            if name.startswith("UniMod:"):
                name = name.split(":")[1]
            terminal = i > 0 and seq[i - 1] == "."
            if terminal and residues:
                mods.append((-1, name))                 # C-terminus
            else:
                mods.append((len(residues), name))      # N-terminus, or a residue
            i = j + 1
            continue
        if c.isalpha():
            residues.append(c)
        i += 1
    return residues, mods


def encode(seq):
    """aa_indices [n+2] and mod_x [n+2, 109]."""
    residues, mods = parse(seq)
    n = len(residues)
    aa = np.zeros(n + 2, dtype=np.int64)
    for k, r in enumerate(residues):
        aa[k + 1] = ord(r) - ord("A") + 1                # A->1 .. Z->26, 0 pads
    if n == 0:
        raise ValueError("empty peptide sequence")
    if any(ord(r) - ord("A") + 1 > 26 or ord(r) - ord("A") + 1 < 1 for r in residues):
        # Out-of-range indices are one-hot encoded to an all-off row, i.e. they
        # become indistinguishable from padding -- and padding is not inert.
        raise ValueError(f"residue outside A-Z in {seq!r}")
    mod_x = np.zeros((n + 2, len(MOD_ELEMENTS)), dtype=np.float32)
    for site, name in mods:
        # site -1 indexes the last row, n+1, exactly as numpy negative indexing
        # does in AlphaPeptDeep's own parse_mod_feature.
        mod_x[site] += mod_vector(name)                  # accumulates
    return aa, mod_x


def predict_rt(model_path, sequences):
    import onnxruntime as ort

    options = ort.SessionOptions()
    # Matches ODIA_ORT_THREADS on the C++ side; reduction order depends on it.
    threads = int(os.environ.get("ODIA_ORT_THREADS", "0"))
    if threads > 0:
        options.intra_op_num_threads = threads
    session = ort.InferenceSession(model_path, options,
                                   providers=["CPUExecutionProvider"])
    names = [i.name for i in session.get_inputs()]
    out = [None] * len(sequences)
    # One batch per encoded length: padding is not inert, because index 0 is
    # one-hot encoded and no model applies a padding mask.
    by_length = {}
    for k, s in enumerate(sequences):
        by_length.setdefault(len(parse(s)[0]), []).append(k)
    for _, group in sorted(by_length.items()):
        aas, mods = zip(*(encode(sequences[k]) for k in group))
        feed = {names[0]: np.stack(aas), names[1]: np.stack(mods)}
        values = session.run(None, feed)[0].reshape(-1)
        for k, v in zip(group, values):
            out[k] = float(v)          # by position: duplicates must not collapse
    return out


# Instrument index, from featurize.py. Anything unrecognised maps to
# max_instrument_num - 1 = 7, not to 0 -- mapping an unknown instrument onto QE
# would silently predict for the wrong one.
INSTRUMENTS = {"QE": 0, "LUMOS": 1, "TIMSTOF": 2, "SCIEXTOF": 3, "THERMOTOF": 4}
UNKNOWN_INSTRUMENT = 7

# MS2 output channels, from alphabase's sort_charged_frag_types:
# sorted(no_loss) + sorted(loss) over [b, y, b_modloss, y_modloss] x charge 1-2.
MS2_CHANNELS = ["b_z1", "b_z2", "y_z1", "y_z2",
                "b_modloss_z1", "b_modloss_z2", "y_modloss_z1", "y_modloss_z2"]


def predict_ms2(model_path, sequences, charges, nce=30.0, instrument="QE"):
    """Predict fragment intensities: [n_peptides][nAA-1][8].

    charges is per peptide. The scale factors are not cosmetic: passing raw NCE
    gives a spectrum with cosine 0.0028 against the correct one, and raw charge
    0.6377 -- unrelated output, with no error.
    """
    import onnxruntime as ort

    options = ort.SessionOptions()
    threads = int(os.environ.get("ODIA_ORT_THREADS", "0"))
    if threads > 0:
        options.intra_op_num_threads = threads
    session = ort.InferenceSession(model_path, options,
                                   providers=["CPUExecutionProvider"])
    names = [i.name for i in session.get_inputs()]
    instrument_index = INSTRUMENTS.get(instrument.upper(), UNKNOWN_INSTRUMENT)

    out = [None] * len(sequences)
    by_length = {}
    for k, s in enumerate(sequences):
        by_length.setdefault(len(parse(s)[0]), []).append(k)

    for _, group in sorted(by_length.items()):
        aas, mods = zip(*(encode(sequences[k]) for k in group))
        feed = {
            names[0]: np.stack(aas),
            names[1]: np.stack(mods),
            names[2]: np.array([[charges[k] * CHARGE_SCALE] for k in group],
                               dtype=np.float32),
            names[3]: np.full((len(group), 1), nce * NCE_SCALE, dtype=np.float32),
            # Rank 1, unlike every other meta input.
            names[4]: np.full((len(group),), instrument_index, dtype=np.int64),
        }
        values = session.run(None, feed)[0]
        for i, k in enumerate(group):
            out[k] = values[i]
    return out


def predict_ccs(model_path, sequences, charges):
    """Predict collision cross-section, one value per peptide.

    Three inputs, not two and not five, and a rank-1 output like the RT model's.
    Charge carries the same 0.1 scale as in MS2.
    """
    import onnxruntime as ort

    options = ort.SessionOptions()
    threads = int(os.environ.get("ODIA_ORT_THREADS", "0"))
    if threads > 0:
        options.intra_op_num_threads = threads
    session = ort.InferenceSession(model_path, options,
                                   providers=["CPUExecutionProvider"])
    names = [i.name for i in session.get_inputs()]

    out = [float("nan")] * len(sequences)
    by_length = {}
    for k, s in enumerate(sequences):
        by_length.setdefault(len(parse(s)[0]), []).append(k)

    for _, group in sorted(by_length.items()):
        aas, mods = zip(*(encode(sequences[k]) for k in group))
        values = session.run(None, {
            names[0]: np.stack(aas),
            names[1]: np.stack(mods),
            names[2]: np.array([[charges[k] * CHARGE_SCALE] for k in group],
                               dtype=np.float32),
        })[0]
        for i, k in enumerate(group):
            out[k] = float(values[i])
    return out


if __name__ == "__main__":
    if sys.argv[1] == "encode":
        aa, mod_x = encode(sys.argv[2])
        print(json.dumps({
            "aa_indices": aa.tolist(),
            "mod_x_nonzero": {str(i): mod_x[i][mod_x[i] != 0].tolist()
                              for i in range(mod_x.shape[0]) if mod_x[i].any()},
            "mod_x_nonzero_idx": {str(i): np.nonzero(mod_x[i])[0].tolist()
                                  for i in range(mod_x.shape[0]) if mod_x[i].any()},
        }))
    elif sys.argv[1] == "ms2":
        model, seq, charge = sys.argv[2], sys.argv[3], int(sys.argv[4])
        spectrum = predict_ms2(model, [seq], [charge])[0]
        print(json.dumps({
            "shape": list(spectrum.shape),
            "channels": MS2_CHANNELS,
            "values": [[round(float(v), 6) for v in row] for row in spectrum],
        }))
    elif sys.argv[1] == "rt":
        seqs = sys.argv[3:]
        for seq, rt in zip(seqs, predict_rt(sys.argv[2], seqs)):
            print(f"{seq}\t{rt:.6f}")
