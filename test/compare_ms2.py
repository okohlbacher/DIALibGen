#!/usr/bin/env python
"""Compare ODIA's predicted fragment spectra against the Python reference.

Covers the inputs the RT model does not have -- charge, NCE and the rank-1
instrument index -- each of which fails silently rather than loudly when wrong:
passing raw NCE gives a spectrum with cosine 0.0028 against the correct one,
and raw charge 0.6377.

Two things this file learned the hard way.

*The cases must reach every channel.* An earlier version varied charge, NCE and
instrument but carried no modification that can lose a neutral fragment, so
channels 4-7 -- the four modloss channels -- were below 5.2e-8 in every case,
25x under the tolerance. Zeroing all four, or swapping them pairwise, passed.
Methionine oxidation is in essentially every real library and puts 36% of the
intensity there, with the base peak in y_modloss_z1. Every channel is now
asserted to be *reached* by at least one case, so a channel cannot quietly
become decoration again.

*The batch must have more than one peptide.* Predicting one at a time leaves
the row striding, the length grouping and the per-peptide charge unexercised:
using charges[0] for every row, or reading every spectrum from output row 0,
passed a suite built on single-peptide calls.

Usage: compare_ms2.py <odia_predict_ms2> <model.onnx>
"""
import json
import math
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import peptdeep_reference as ref

# Sequence, charge, NCE, instrument. Charge, NCE and instrument are varied
# deliberately: holding them fixed would let a defect that ignores one of them
# pass, which is how the scale factors were missed in the first place.
CASES = [
    ("ELVISLIVESK", 2, 30.0, "QE"),
    ("ELVISLIVESK", 3, 30.0, "QE"),
    ("PEPTIDEK", 2, 27.0, "Lumos"),
    ("PEPC(Carbamidomethyl)TIDEK", 2, 35.0, "timsTOF"),
    (".(Acetyl)PEPTIDEK", 2, 30.0, "QE"),
    ("PEPTIDER.(Amidated)", 2, 30.0, "SciexTOF"),
    ("AVVPASLSGQDVGSFAYLTIK", 3, 28.0, "NotAnInstrument"),
    # The two that reach the modloss channels at all. Oxidised methionine puts
    # 36% of the intensity there and the base peak in y_modloss_z1;
    # phosphoserine, 51%.
    ("PEPTIDEM(Oxidation)K", 2, 30.0, "QE"),
    ("AAAAAAAS(Phospho)K", 2, 30.0, "QE"),
    ("S(Phospho)AAAAAAAK", 3, 32.0, "Lumos"),
    # The doubly-charged b channels, 1 and 5, need more than a modification:
    # a b ion only carries two charges when the N-terminal fragment holds two
    # basic residues, so these have missed cleavages and a high precursor
    # charge. Without them both channels stayed near 1e-6 -- a hundredth of the
    # comparison tolerance -- even with oxidation and phosphorylation present.
    ("HKRS(Phospho)AAAAAAAKR", 4, 30.0, "QE"),
    ("RHKM(Oxidation)PEPTIDEKR", 4, 33.0, "Lumos"),
]

# Batches, run as a single call each, so the striding and grouping are covered.
# Mixed lengths force the grouping to permute the input and un-permute the
# output; distinct charges within one length group are what catches a row that
# takes the first peptide's charge.
BATCHES = [
    ([("ELVISLIVESK", 2), ("SAMPLERSAMK", 3), ("GGGGGGGGGGK", 6)], 30.0, "QE"),
    ([("PEPTIDEK", 2), ("ELVISLIVESK", 4), ("AAAK", 1),
      ("PEPTIDEM(Oxidation)K", 3), ("AVVPASLSGQDVGSFAYLTIK", 2),
      ("PEPC(Carbamidomethyl)TIDEK", 1)], 28.0, "timsTOF"),
]

TOLERANCE = 1e-5
CHANNELS = ["b_z1", "b_z2", "y_z1", "y_z2",
            "b_modloss_z1", "b_modloss_z2", "y_modloss_z1", "y_modloss_z2"]
# A channel counts as reached when the reference itself puts something in it
# well above the comparison tolerance. 100x is arbitrary but the gap in
# practice is either ~1e-8 (nothing there) or ~1e-1 (a real peak).
REACHED = 100 * TOLERANCE


def run(tool, model, nce, instrument, pairs):
    args = [tool, model, str(nce), instrument] + [f"{s}:{z}" for s, z in pairs]
    return subprocess.run(args, capture_output=True, text=True,
                          env=dict(os.environ, ODIA_ORT_THREADS="1"))


def compare(label, got, want, reached):
    """Returns a message on disagreement, None on agreement."""
    if list(got["shape"]) != list(want.shape):
        return f"shape {got['shape']} vs {list(want.shape)}"
    worst = 0.0
    for p, row in enumerate(got["values"]):
        for c, value in enumerate(row):
            if not math.isfinite(value):
                return f"non-finite intensity at position {p}, channel {c}"
            reached[c] = max(reached[c], abs(float(want[p][c])))
            worst = max(worst, abs(value - float(want[p][c])))
    if worst > TOLERANCE:
        return f"worst {worst:.2e}"
    return None


def check_error_paths(tool, model):
    """Inputs that must be refused, and one that must not poison its batch."""
    bad = 0

    # A charge that is not a charge. The model gives no sign: charge 0 returns a
    # base-peak-normalised spectrum like any other, and the DIA-NN reader has a
    # live path to it, because a null Precursor.Charge parses to 0.
    for pairs, why in (([("PEPTIDEK", 0)], "charge 0"),
                       ([("PEPTIDEK", -3)], "charge -3"),
                       ([("PEPTIDEK", 99)], "charge 99")):
        proc = run(tool, model, 30.0, "QE", pairs)
        if proc.returncode == 0:
            print(f"  FAIL {why} was accepted")
            bad += 1

    # NCE is one value for the whole call, so it is refused outright rather
    # than reported per peptide.
    for nce, why in ((0.0, "nce 0"), (-30.0, "nce -30"), (1e9, "nce 1e9")):
        proc = run(tool, model, nce, "QE", [("PEPTIDEK", 2)])
        if proc.returncode == 0:
            print(f"  FAIL {why} was accepted")
            bad += 1

    # One unencodable peptide must cost one peptide, not its whole chunk. The
    # earlier version reported every peptide in the chunk as failed, each with
    # a reason naming another peptide's modification.
    # C[999] is a bare mass shift with no elemental composition, which is what
    # PeptDeep cannot encode. C[+57.0] would NOT do: OpenMS resolves it to
    # Carbamidomethyl, which has one.
    #
    # All three are eight residues, which is the point. Peptides are batched by
    # length, so a bad peptide of a length nothing else shares sits alone in its
    # chunk and discarding the whole chunk costs exactly one peptide -- the
    # defect this is here to catch would pass. It must share a chunk with good
    # ones.
    pairs = [("PEPTIDEK", 2), ("PEPTIDC[999]K", 2), ("SAMPLERK", 3)]
    proc = subprocess.run(
        [tool, model, "30.0", "QE"] + [f"{s}:{z}" for s, z in pairs],
        capture_output=True, text=True,
        env=dict(os.environ, ODIA_ORT_THREADS="1", ODIA_MS2_TOLERATE_FAILURES="1"))
    if proc.returncode != 0:
        print(f"  FAIL tolerated batch died: {proc.stderr.strip()}")
        bad += 1
    else:
        out = ref.load_json(proc, "tolerated batch")
        if out["failed"] != [1]:
            print(f"  FAIL one bad peptide reported as {out['failed']}, expected [1]")
            bad += 1
        empty = [k for k, sp in enumerate(out["spectra"]) if sp["shape"][0] == 0]
        if empty != [1]:
            print(f"  FAIL one bad peptide emptied spectra {empty}, expected [1]")
            bad += 1
    if not bad:
        print("  ok   error paths (bad charge, bad NCE, one bad peptide in a batch)")
    return bad


def main(tool, model):
    failures = 0
    reached = [0.0] * len(CHANNELS)
    os.environ["ODIA_ORT_THREADS"] = "1"

    for seq, charge, nce, instrument in CASES:
        proc = run(tool, model, nce, instrument, [(seq, charge)])
        if proc.returncode != 0:
            print(f"  FAIL {seq} z={charge}: {proc.stderr.strip()}")
            failures += 1
            continue
        want = ref.predict_ms2(model, [seq], [charge], nce, instrument)[0]
        bad = compare(seq, ref.load_json(proc, seq)["spectra"][0], want, reached)
        if bad:
            print(f"  FAIL {seq} z={charge} nce={nce} {instrument}: {bad}")
            failures += 1
        else:
            print(f"  ok   {seq} z={charge} nce={nce} {instrument}")

    for pairs, nce, instrument in BATCHES:
        proc = run(tool, model, nce, instrument, pairs)
        if proc.returncode != 0:
            print(f"  FAIL batch of {len(pairs)}: {proc.stderr.strip()}")
            failures += 1
            continue
        got = ref.load_json(proc, "batch")["spectra"]
        seqs = [s for s, _ in pairs]
        charges = [z for _, z in pairs]
        want = ref.predict_ms2(model, seqs, charges, nce, instrument)
        if len(got) != len(pairs):
            print(f"  FAIL batch of {len(pairs)}: got {len(got)} spectra")
            failures += 1
            continue
        bad_any = False
        for k, (seq, charge) in enumerate(pairs):
            # The tool echoes what it was given, so a batch that silently
            # dropped or reordered a peptide cannot look like agreement.
            if got[k]["peptide"] != seq or got[k]["charge"] != charge:
                print(f"  FAIL batch position {k}: got {got[k]['peptide']} "
                      f"z={got[k]['charge']}, sent {seq} z={charge}")
                failures += 1
                bad_any = True
                continue
            bad = compare(seq, got[k], want[k], reached)
            if bad:
                print(f"  FAIL batch position {k} ({seq} z={charge}): {bad}")
                failures += 1
                bad_any = True
        if not bad_any:
            print(f"  ok   batch of {len(pairs)}, nce={nce} {instrument}")

    failures += check_error_paths(tool, model)

    # A channel no case reaches is a channel this comparison does not test.
    for c, name in enumerate(CHANNELS):
        if reached[c] < REACHED:
            print(f"  FAIL channel {c} ({name}) peaks at {reached[c]:.2e} across "
                  f"every case, under the {REACHED:.0e} needed for the comparison "
                  f"to see it. Add a case that produces this ion.")
            failures += 1

    print(f"MS2 comparison: {failures} failures over {len(CASES)} cases "
          f"and {len(BATCHES)} batches")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1], sys.argv[2]))
