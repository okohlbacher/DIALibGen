# Historical benchmark: 0.10.1

This summary covers the completed 20 September 2026 benchmark of
**DIALibGen 0.10.1**, **DIALibRefine 0.3.0-dev** and **DIA-NN 2.0** on K562
from a timsTOF diaPASEF acquisition. It is not a benchmark of 0.11.0.
In particular, 0.11.0 restores Met-excised missed-cleavage peptides, changing
unrestricted generation outputs.

The protocol used a development run to select configurations and two
previously exposed technical replicates for comparison. All primary
libraries were restricted to the same 5,754,771 precursor keys. Counts below
are for run 2 / run 3 at a matched 1% combined-estimate entrapment FDP,
an upper-bound estimate.

| Workflow | Precursor identifications | Protein groups |
|---|---:|---:|
| DIA-NN, untuned | 112,391 / 114,586 | 7,213 / 7,228 |
| DIALibGen, untuned | 112,333 / 114,210 | 7,242 / 7,223 |
| DIA-NN, tuned | 118,310 / 119,510 | 7,397 / 7,378 |
| DIALibGen + DIALibRefine, tuned | 117,380 / 119,102 | 7,418 / 7,386 |

The selected DIALibGen configuration was timsTOF/NCE 30 followed by RT/CCS
training at learning rate 0.0001. DIA-NN selected learning rate 0.0001 with
fragment-model tuning. Tuned DIALibGen was 0.79% / 0.34% lower in whole-run
precursor yield, within the registered 1% margin. The protein-group margins
were also met. These categories do not establish equivalence.

The qualification is material: DIALibGen was 1.56% / 1.11% lower on the held-out
protein cohort, and one automatic-window sensitivity gave a 1.11% whole-run
deficit. DIALibGen had tighter RT residuals on unseen sequences in the common
identification set; DIA-NN had higher spectral-angle agreement. Those are
search-derived descriptive measurements and do not isolate a cause for the
identification difference.

This demonstrates cross-run transfer within the tested acquisition method.
It does not establish general superiority, transfer to another instrument or
method, or quantitative equivalence. Two exposed technical replicates from
one sample are not independent confirmation. The tools had different tuning
candidate budgets and only DIA-NN tuned fragment prediction.

Observed-value write-in reversed direction between runs and carried an
entrapment-asymmetry warning. Filtered empirical libraries could not be ranked
against the full libraries at the registered operating point. Neither is an
established improvement from this experiment.

## Provenance

Official scoring version 2 completed on 20 September 2026 after a disclosed
post-unblinding fix to recognize an entrapment identifier after a pipe.
Independent recomputation confirmed that the fix changed guard/status fields
for two reports and did not change primary counts. Earlier cache-marker and
registered G0/G1 alias plumbing amendments were also archived. G1 selected the
same NCE as G0 and supplies no independent evidence.

- Frozen library record SHA256:
  `e8d4c38a86dda9fdd1bd7f13f44f818e58a1de175ab4d5474bfa135e3a177729`.
- Final script-manifest SHA256:
  `5f2ec18622192e8881a283729c3ad28ef0e58dfdcb4b8f2678d813d76aea9ed2`.

The full experiment archive is maintained separately from the product source;
this summary does not distribute raw data or claim a fresh reproduction.
