# Design: built-in identification for `refine` and `tune`

Status: **accepted, in implementation** (milestones M1 and M3). The option
stays *experimental* until the honesty gate (a) below passes.

## Goal

`refine` and `tune` currently need `-ids report.parquet`, the report of an
external DIA-NN search of one run. With `-run run.mzML` DIALibGen finds its own
confident identifications in that run by targeted trace extraction, scoring and
target-decoy FDR, and feeds them to the existing refinement and tuning code. The
external route stays available and unchanged.

```bash
DIALibGen -mode tune -tune_heads rt -in predicted.parquet -run run.mzML -out tuned.tsv -threads 8
DIALibGen -mode refine -tune -tune_heads rt -no_filter -in predicted.parquet -run run.mzML -out adapted.tsv
DIALibGen -mode refine -in predicted.parquet -ids report.parquet -out refined.tsv   # unchanged
```

(M1 measures no 1/K0, so the CCS head and `-write_im` are refused with `-run`
until M2.)

The approach follows OpenDIAlyzer (ODIA): OpenSWATH's extraction and
sub-scoring from OpenMS, with the harness around it (candidate selection,
calibration, semi-supervised LDA, FDR) owned by the tool. Unlike ODIA it uses
**stock OpenMS 3.5.0** only, the version every DIALibGen release links, so it
builds on all five platforms without patches.

## Interface

- `-run <file>` (mzML). Exactly one of `-run` and `-ids` is required.
- `-out_ids <file>`: the identification report, default `<out>.ids.parquet`,
  never overwritten. It is written in the columns the existing readers expect
  (`Run`, `Precursor.Id`, `Modified.Sequence`, `Precursor.Charge`,
  `Precursor.Mz`, `Protein.Group`, `Decoy`, `RT` in minutes, `IM` as 1/K0,
  `Q.Value`, `Global.Q.Value`, `PG.Q.Value`, `PEP`), and then serves as `-ids`
  for the rest of the invocation. The report readers, gates and provenance code
  do not change.
- TOPP subsection `search:` with, among others: `candidates evidence|random`
  (default evidence since M3), `prefilter_depth 5`, `prefilter_top_peaks 1000`,
  `prefilter_ppm 10`, `subset` (random only), `max_pairs 200000`,
  `decoys shuffle|pseudo_reverse`, `intensities predicted|library` (default
  predicted since the M3 review), `instrument auto`, `nce -1`, `seed`,
  `passes`, `rt_window`, `mz_ppm`,
  `im_window` (0 = automatic), `ms1`, `rt_im_scores`,
  `calibration_min_rsq 0.70`, `calibration_min_coverage 0.30`, `readoptions`,
  `min_ids 200`, `entrapment_tag`, `selftest true`, `chunk`, `batch_size`.
- Refused with `-run` until later milestones: `-write_intensity`,
  `-empirical_library`, `-min_fragments`, `-write_im` and `-tune_heads
  ccs|both` (M1 measures no 1/K0). Everything that would fail after the search
  (tuning recipe, models, a library without protein groups, a directory such
  as a Bruker `.d` given as `-run`, the MS2 model that `intensities
  predicted` needs, a library whose RT has no range) is checked before the
  run is read; a failure after the search keeps the report and says how to
  reuse it. The MS2 model comes from where `-tune` finds its models
  (`-tune_models`, allowed with `-run` in `refine` for this, then
  `$DIALIBGEN_MODEL_DIR`, then the bundled models).
- The provenance sidecar gains `inputs.run` and a `search` block: settings,
  calibration, candidate counts, targets and decoys at each gate, resources.
  The refined Parquet library embeds it without the resources (timings,
  threads), so the same command writes the same bytes.
- Every guard aborts with counts in the message. Nothing falls back silently.
  `search:min_ids` is the one guard that writes the report first (its
  embedded metadata records the abort, and the message names the file), so
  a search that is honest but too small can be inspected; every other guard
  fires before a report exists, and before `min_ids` when both apply.
- A run stock OpenMS cannot parse fails with the reason OpenMS logged (stock
  `MzMLFile` rethrows a `ParseError` that names only XMLHandler.cpp@108; the
  reason goes to its log stream, which the loader now copies), e.g.
  `Required attribute 'defaultDataProcessingRef' not present!` for the mzML
  that mzpeak-convert 0.12.5 writes, with a hint.

## Architecture

New static target `odia_search` (`include/odia/search/`, `src/search/`). It links
only OpenMS and OpenSwathAlgo, which `odia_library` already links, so release
bundles do not grow.

1. **Library.** Loaded as today and never mutated; a side index orders
   precursors by m/z. Decoys for the search are built in memory (shuffle by
   default); decoys already in the file are ignored for the search and left
   untouched. The DIA-NN-derived `mutate` method is not offered, and neither
   is `reverse`, which moves the tryptic C-terminus and makes decoys separable
   by construction. A target whose interior cannot be rearranged has no decoy
   (no residue substitution as a fallback), and an arrangement that
   reproduces the target's fragment masses (I/L) or puts one outside the
   m/z range of the library's target fragments is re-drawn. Each rule reads
   the target only, so a failure removes a whole pair.
   **Both termini keep TWO residues** (the generator's own decoys keep one;
   `SEARCH_DECOY_KEEP_NTERM`). The second residue from each end fixes the
   short ions -- y1, y2, b1, b2 and, with the composition, y(n-1), y(n-2),
   b(n-1), b(n-2) -- and a whole proteome shares those: every peptide ending
   in "PK" has the same y2, and the run's spectra are full of it because
   cleavage N-terminal to proline makes it intense in the peptides that ARE
   there. A decoy that moves the second residue trades its target's short ion
   for another composition, so a null pair is decided by which member drew
   the commoner one. Measured on the Astral entrapment library, 2.69 M
   known-absent pairs with predicted assays: keeping one residue, the deeper
   member of a null pair was its decoy 821,653 times against the target's
   790,097 (sign test z -25), and null targets whose second-to-last residue
   is P beat their decoys 1.73 : 1 while the rest lost; keeping two,
   753,786 against 753,635 (z -0.1) with every positional split flat, and the
   present targets' excess at the full prefilter depth untouched (10,671
   against their decoys' 4,578). Cost: a peptide of six residues or fewer has
   no decoy (3,415 of 6.16 M eligible targets against 357).
   **Fragments and intensities of both members (`search:intensities
   predicted`, the default since the M3 review).** The target and its decoy
   are each predicted by the same PeptDeep MS2 model (the instrument and NCE
   the library recorded when DIALibGen generated it, else `search:instrument`
   / timsTOF for an ion-mobility run / QE) from their OWN sequences, and each
   takes its own most intense b/y fragments (charge 1 to min(2, precursor
   charge), no loss, within the library's target fragment m/z range, the
   generator's rank order), the same number for both: min(the target's
   library transition count, max(3, the fewer of the two members' fragments
   above 1e-4 of their base peak)). A pair the model cannot predict, whose
   decoy's chosen fragments all lie within 10 ppm of the target's, or with
   fewer than 3 fragments leaves whole. BOTH members are predicted, always,
   whatever the library holds: where the library is DIALibGen's own with the
   same model, the target's assay comes OUT equal to its library assay cut to
   that count (a test checks it with the real model), but it is never read
   out of the library -- see "The library-assay shortcut, removed".
   Why: a decoy in its target's fragment slots with its target's intensities
   (the M1-M3 rule, still `search:intensities library` for comparison, which
   now warns loudly on every run that selects it) is weaker than a null
   target. The slots were chosen by a predictor for the TARGET's sequence --
   the compositions real peptides fragment into, Pro-directed y2/y3 ions
   above all -- and the decoy puts random compositions into them. On the M3
   acceptance runs, known-absent entrapment targets reached the prefilter
   depth 1.18x (timsTOF) and 1.11x (Astral) as often as their own decoys
   (sign tests z 81 and 79), and won 240 : 168 of the entrapment pairs
   identified at q <= 0.01 (z 3.6; 2,594 : 2,078 at q <= 0.1): q-values about
   1.4x too optimistic, the whole entrapment excess over nominal. The same
   slots chosen at random or from the bottom of the ranking shrank the
   imbalance; only each member's own prediction removes it by construction.
   Cost: both members of every eligible pair are predicted (about 7,000
   peptides per second at 16 CPU sessions), and the searched set once more --
   24 minutes of a 3 h timsTOF run. Predicting a target whose assay the
   library already holds is exactly that much work for the same answer, and
   taking the library's assay instead is not the same answer: it puts the
   library's own fragment cap into the pair's count rule ("The library-assay
   shortcut, removed").
2. **Run.** `SwathFile::loadMzML` (`normal` or `cache`). Ion-mobility window
   limits are normalised (lower/upper swapped where reversed, as in current
   mzpeak-convert output). diaPASEF is detected from window limits *and* a
   per-peak 1/K0 array. The per-window cache files hold every peak
   uncompressed: 1.45x the mzML on an Orbitrap Astral run and 2.56x on a
   diaPASEF run (whose 1/K0 array is cached too); the log quotes that range
   before the read and the measured size after it. Later (M4): a
   DIALibGen-owned streaming store that keeps only the peaks inside active
   assay boxes, so memory stays bounded.
3. **Candidates.** Eligibility first, for both selections: a target whose
   precursor m/z lies in none of the run's isolation windows is ineligible
   before any choice and the cap (its decoy shares the m/z; 19.5 % of a
   random draw on an Astral run), and so are targets with too few fragments,
   no RT, or a (sequence, charge) held twice.
   - `search:candidates evidence` (default since M3), `EvidencePrefilter`:
     every eligible target gets its search decoy computed (`searchDecoy`, the
     one function `appendSearchDecoys` is built on, so the prefilter scores
     exactly the decoy that is later searched; a test checks the m/z bit for
     bit), and both members are predicted (above). Both members enter one
     m/z-sorted index per isolation window with their own 6 most intense
     predicted fragments (with `intensities library`: the target's by
     library intensity, the decoy in the same slots). One sweep over the run's MS2 spectra matches each spectrum's
     `prefilter_top_peaks` most intense peaks at +-`prefilter_ppm`; per member
     it records the depth (distinct indexed fragments matched in one
     spectrum), the spectra reaching `prefilter_depth`, and the run time of the
     best spectrum (deepest, then most matched intensity, then earliest). A
     pair is kept when its target OR its decoy reaches the depth (pair-union).
     The kept set is checked to hold as many decoys as targets (a ratio guard
     that catches a selection keeping members instead of pairs), then capped
     to `max_pairs`: strata are isolation window x library-RT decile x charge,
     each gets its largest-remainder share, and within a stratum pairs rank by
     the better member's (depth, spectra), then by draw key. Every decision
     about a pair reads its two members through a symmetric function, never
     the label; indexing every target fragment as its decoy's (a test hook)
     gives the exact swap of the evidence and the identical pair set.
   - `search:candidates random` (M1): a deterministic, paired random subset
     (the draw hashes the key target and decoy share, so it is label-blind);
     with predicted intensities, the first `max_pairs` drawn pairs that the
     model can predict.
4. **Assays.** Per chunk, an `OpenSwath::LightTargetedExperiment` built directly
   from the library arrays. Never the whole library. Peptide sequences are left
   empty for targets and decoys, which makes OpenSWATH's scoring sequence-blind
   and therefore symmetric: DIALibGen decoys carry their target's sequence.
5. **Calibration.** `OpenSwathCalibrationWorkflow::performRTNormalization` on seed
   assays, with thresholds suited to predicted libraries. It returns the RT
   transformation, m/z and 1/K0 windows and the ion-mobility correction; its
   points are refit robustly and validated (>= 20 points, r^2, slope, span).
   Failure aborts. Seeds: every iRT/CiRT kit precursor the library holds,
   plus, with evidence candidates (M3), the prefilter's seeds -- targets that
   reach the depth THEMSELVES and beat their own decoy clearly (at least 3
   spectra at the depth and 3x their decoy's), one per peptide, the best
   (depth, spectra) of each of 100 bins, at most 2,000; a bin without such a
   target stays empty -- and with random candidates (M1) a stock-sampled
   share of the searched targets. Every seed, kit precursors included, comes
   from the central library-RT range (the 0.1 % to 99.9 % quantiles of the
   targets' library RT), and the 100 bins span that range: before the M3
   review they spanned the raw minimum to maximum, so ONE library row with an
   absurd RT cut the evidence seeds 10-fold and, if it was a seed itself,
   failed the stock binned-coverage check and the whole search (the abort
   message now names the seeds' and the library's RT ranges). Random seeds are what
   limited M1: on a whole-proteome library of which DIA-NN identifies 0.23 %
   of the precursors in the Astral run, 2,191 seeds gave 16 points and the
   calibration failed. The prefilter's best-spectrum RT of each seed is
   compared with the point the calibration picked for it and recorded (a
   diagnostic, not a filter). A LOWESS second fit (stock
   `TransformationModelLowess`, span chosen by cross-validation, linear
   interpolation, the line through the end points outside them, its own
   iterative 3-sigma outlier cut on the stock points, at most 10 fits)
   replaces the line when at least 200 points remain, the map increases in
   both directions, and it predicts held-out points better: 5-fold
   cross-validation of both models on ONE point set (the union of both
   inlier sets), residuals capped at 3 robust line scales, LOWESS's mean
   error below 0.97x the line's. Until the M3 review the rule compared the
   two RT windows (0.99 quantile of each model's own inliers): on the Astral
   run it kept the line, which ran a median 89 s late in the top library-RT
   decile and missed 26 % of DIA-NN's identifications there -- a bias the
   global quantile cannot see, not "the predictor's scatter". The median
   residual per library-RT decile of the chosen model is recorded. The RT
   window (0.99 quantile of the residuals x 1.3, doubled, at least 30 s) is
   never narrower than 2 x 1.3 x 2.576 robust SDs of the residuals, so a
   0.99 quantile of 20 residuals (their maximum) cannot size it alone, and a
   calibration on fewer than 100 points warns. Measured before the review:
   on Astral the evidence seeds gave 944 points (M1's random seeds 16); on
   timsTOF LOWESS narrowed the window (661 s against 756 s). Calibration memory grows with seeds x fragments x spectra per
   window (M1: 3.5 GiB on the 6.4 GB Astral run, 1.5 GiB above extraction).
6. **Extraction.** `OpenSwathWorkflow::performExtraction` (stock 13-argument form)
   per chunk, with an inactive OSW writer and in-memory features. After each
   chunk the features become compact score rows and are freed. Stock OpenSWATH
   runs on at most 8 OpenMP threads: its feature scoring serialises on the
   process-wide MetaInfoRegistry lock, and more threads only spin (extraction
   was twice as slow at 32 threads as at 8). It parallelises over isolation
   windows, so since M3 a chunk is not m/z-contiguous: the pairs, in m/z
   order, are cut into pieces of one OpenSWATH batch (`batch_size`
   precursors) dealt round-robin over the chunks, so every chunk holds pieces
   of many windows while each window still sees whole batches (M1's
   contiguous 20,000-precursor chunks covered 1-2 of a diaPASEF run's
   25-Th windows: 43 ms per precursor against 11-17 ms over about 8).
7. **Scoring and FDR** (below), then the report.

Shared code with ODIA: the dependency-free LDA and FDR headers are vendored into
`src/odia-core/` with a manifest of origin commit and hashes, and a CI check.
ODIA's DIA-NN-shaped network classifier, its Bruker/mzPeak readers and its
patched-OpenMS prefilter are not ported.

## Scoring and FDR

- Features: about twelve non-collinear OpenSWATH sub-scores, plus MS1 and ion-
  mobility scores when available; elution-model and ion-series scores off. The
  MS1-MS2 co-elution scores are the `_contrast` variants: stock 3.5.0 computes
  the plain ones only from two or more precursor isotope traces, and only the
  monoisotopic trace is extracted.
- Classifier: ODIA's semi-supervised LDA (mProphet/pyProphet lineage), three
  folds, three iterations, with four changes: folds assigned per target-decoy
  pair; fold scores rescaled label-blind (median/MAD of all held-out scores);
  RT and 1/K0 deviation scores hidden from the calibration discriminant; GBT and
  network code removed.
- Precursor level: concatenated target-decoy competition. Each pair keeps its
  better member (ties to the decoy); q = (D + 1)/T, monotonised, no pi0.
  Pairing is structural by (modified sequence, charge).
- Peptide and protein level: roll-up by canonical sequence and by the library's
  `Protein.Group` string, picked competition at both levels with the same
  estimator. The existing `refine`/`tune` gates stay the control.
- Why paired, not pooled. Concatenated competition needs only that a NULL
  target and its own decoy are equally likely to win their pair. With both
  termini fixed, a decoy shares b(n-1) and y(n-1) with its target (11 % of all
  decoy fragments lie within 10 ppm of one of their target's on an Astral
  run), so the decoy of a PRESENT target lights up with it -- and loses its
  pair, which is correct. A pooled estimator ranks those losing decoys too and
  counts them as false discoveries: it is conservative, by a lot (365 vs 825
  precursors at 1 % on that run; 375 of 825 passed a pooled peptide gate
  before the peptide level became picked). The pooled count stays in the
  provenance as a diagnostic; "pooled below paired" is the expected direction,
  and only pooled above twice the paired count (targets losing pairs they
  should win) is a warning.
- What tests exchangeability, and what does not.
  - The **null-pair balance** is a BULK diagnostic, not the exchangeability
    test. Among complete pairs whose winner is in the lowest quarter of
    winning scores, targets and decoys must win about equally often; logged,
    recorded, a warning at |z| > 3. It assumes most candidates are null
    (depth 6 on the Astral run, 54 % of the targets identified, puts present
    targets into that quarter: z 3.2 without a fault). And it cannot see an
    asymmetry that grows with the score, which is where the 1 % cut falls:
    the M3 Astral entrapment search passed it (23,676 : 24,082, z -1.9) while
    its entrapment pairs were won 439 : 347 by the target at q <= 0.1
    (z 3.3); on timsTOF the ratio grew into the tail (1.25 at q <= 0.1, 1.43
    at 0.01, 2.0 at 0.001).
  - The **entrapment winner test** is. With `search:entrapment_tag`, every
    pair whose target is an entrapment entity is a known null pair: target
    and decoy must win it equally often. The winners at q <= 0.01 and
    q <= 0.1 are counted at precursor, peptide and protein-group level
    (picked entities), recorded with their binomial z, and z > 3 at any of
    them warns that decoys are weaker than null targets. It needs no
    database ratio and no isomer screening, and it is a standing acceptance
    metric of gate (a). For production runs without entrapment, a hidden
    probe set (shuffled copies of a few percent of the targets, competing
    with their own decoys and left out of the report) would give the same
    check; it is not built yet.
- Entrapment (`search:entrapment_tag`, a validation aid, not a gate): a
  protein is entrapment when its id starts with the tag or carries it right
  after a '|' (`sp|ENTRAP_P12345|...`, what `generate` writes; before the M3
  review only the start counted, and the estimate silently read "no
  estimate"), the same rule for every member of a group (all tagged: an
  entrapment group; some: shared, left out). A tag that matches no protein
  warns. The combined estimator weights each entrapment identification by 1 + 1/r, r the
  database ratio of entrapment to real targets. With random candidates r is
  counted over the searched pairs; with evidence candidates over the eligible
  pairs BEFORE the prefilter, which keeps present targets far more often
  than absent ones (entrapment and null targets pass by the same rule). The
  searched ratio was 0.083 against 1.00 in the library in a diaPASEF search
  at depth 6, and the estimate read 7.9 % where the library ratio gives
  1.2 %.
- Run-level guards, each aborting: decoy:target ratio out of band, more than
  half of the targets at q <= 0.01 (the same premise as the null-pair
  balance), too few identifications (`search:min_ids`, which writes the report
  first), failed calibration, the candidate set's decoy:target count ratio
  (whole pairs: exactly 1), and the self-checks (`search:selftest`, on by default):
  scoring with swapped labels and with each pair's labels exchanged on a coin
  must identify (almost) nothing. They catch a classifier that leaks labels.
  They do NOT catch decoys built weaker than null targets -- random pair labels
  symmetrise the construction away (reverse decoys, whose null pairs won 1.66:1
  by the target, passed both with 0 and 0). Only known nulls do: the
  entrapment winner test above.
- Decoy construction and what is known about it: targets and decoys have
  identical precursor m/z, RT, charge and protein group. In M1 (Astral run)
  extraction found peak groups for 199,999 targets and 200,000 decoys of
  200,000 pairs, and the null-pair balance was 25,016 : 24,983; but null
  targets beat their shuffle decoys on the MS2 fragment sub-scores in
  51.2-51.5 % of the pairs of an entrapment search (z 7-9 over 105,433
  pairs), and the M3 review traced this, and the entrapment excess, to the
  decoy's inherited fragment slots and intensities (Architecture, 1). Since
  the review both members are predicted (`search:intensities predicted`).

## Resources

Measured ODIA figures show the cost is dominated by materialising the library as
OpenSWATH objects and by holding the run in memory. This design never does the
first, bounds the candidate set (`search:max_pairs`), extracts once by default,
and (from M4) streams the run. Projected on 8-12 laptop cores: about 15 min for a
6 GB Orbitrap Astral mzML and 20-30 min for a 29 GB diaPASEF mzML, with 4-10 GB
peak memory. These figures are derived, not measured; M1 and M4 measure them.

M1, measured on the 6.4 GB Astral mzML (a 254,767-precursor library from a
FASTA subset of which DIA-NN identified 2.9 % in that run; 8 threads on a
shared node; `-mode refine`, whole invocation):

| `search:max_pairs` | identified at 1 % | decoys | pass refine's gates | tune val / test units | in DIA-NN's report | wall | peak RSS |
|---|---|---|---|---|---|---|---|
| 40,000 | 1,045 | 9 | 992 | 93 / 191 | 85.6 % | 9.3 min* | 3.5 GB |
| 100,000 | 2,804 | 27 | 2,659 | 250 / 511 | 84.0 % | 9.4 min | 3.5 GB |
| 200,000 (default) | 5,821 | 57 | 5,541 | 454 / 1,030 | 84.5 % | 20.4 min | 4.5 GB |

(*a slow read of the run: 280 s instead of 73 s.) Median |dRT| to DIA-NN is
0.4 s at every size. Identifications scale with the share of the library that
is present: a whole-proteome library, where that share is several times
smaller, needs proportionally more pairs for tuning's 100-unit cohorts, which
is what the evidence prefilter (M3) is for.

### M3: the evidence prefilter on whole-proteome libraries

Measured with the build of `feature/builtin-identification` that introduced it,
on the shared 384-core node, 16 threads (OpenSWATH 8), 1-4 other jobs of this
measurement running at the same time. "DIA-NN in the set" counts the DIA-NN
2.0 report's precursors (targets, q <= 0.01) among the searched targets;
"recovered" is the share of those we identify at q <= 0.01. The searched set
of each setting comes from `identify_prefilter_probe` (test/search), which
runs the product's prefilter code; for every run below all of our
identifications lie inside it.

**Orbitrap Astral, whole-proteome library** (human reviewed FASTA, 4,020,619
precursors; DIA-NN identifies 9,261 = 0.23 % of them in the run; 3,079,208
eligible pairs, 189 targets without a decoy). M1 stopped at calibration on
this library (16 points from 2,191 random seeds).

| | M1, random 200,000 pairs | M3 defaults |
|---|---|---|
| searched pairs | 200,000 | 100,728 (every pair that passed; the cap did not bind) |
| DIA-NN in the set | 573 | 7,698 (83.1 %) |
| calibration | failed: 16 points | 944 points from 1,927 seeds (1,734 evidence, 193 kit), r^2 0.982, RT window 317 s |
| target precursors at q <= 0.01 (decoys) | - | 4,755 (46); 4,143 peptides, 300 protein groups; 3,954 pass refine's three gates |
| in DIA-NN's report / recovered | - | 94.5 % / 58.4 % |
| median / p99 \|dRT\| to DIA-NN | - | 0.006 / 0.066 min |
| null-pair balance / self-checks | - | 12,713 : 12,469 (z 1.5) / 0 and 0 |
| `-mode tune -tune_heads rt` | aborted | completed: test 887 / val 329 units; TEST calibrated SD 0.910 -> 0.500 min (raw 2.081 -> 0.500) |
| search: load, prefilter, calibration, extraction | 72 s read, then the abort | 72 s, 12 s, 102 s, 380 s (9:31 in all) |
| whole invocation, peak RSS | 3:49, 6.0 GB | 11:56, 16.7 GB (the tuner's whole-library re-prediction; the search alone peaks at 5.6-7.2 GB in refine runs) |

The `-out_ids` report and the tuned library are byte-identical at 8 and 16
threads (8 threads: 11:52, 17.1 GB), and the report is byte-identical again
from the documented workflow `-mode refine -tune -tune_heads rt -no_filter`
at `search:chunk 7000` (29 chunks instead of 11; 12:36, 17.0 GB), which
wrote 3,954 observed RTs with a residual SD of 0.353 min after tuning.

The prefilter itself took 12.4 s: 7.1 s for the decoys of 3.08 M targets,
0.3 s for the 308 MB index (36.95 M fragments of both classes in 150
windows), 4.4 s for the sweep of 303,701 spectra. Members by depth, 0 to 6:
targets 0 / 7,568 / 806,611 / 1,697,513 / 505,967 / 52,969 / 8,580, decoys
0 / 10,189 / 876,800 / 1,676,802 / 468,627 / 44,854 / 1,936. Chance
co-occurrence is common up to depth 4, and the targets' excess (signal)
lies at depths 5 and 6. Of 974 calibration points matched to an evidence
seed, 968 lie within 30 s of the seed's best prefilter spectrum (median
0.37 s). LOWESS did not narrow the RT window (321 s against the line's
317 s, span 0.2), and the rule of the time kept the line. That was read as
"per-peptide scatter, not curvature"; the M3 review showed otherwise: the
line runs a median 89 s late in the top library-RT decile and misses 26 % of
DIA-NN's identifications there, and LOWESS on the same points halves that
bias (the rule now cross-validates, Architecture, 5).

Choosing the defaults (first the probe, 4 min per grid, then full searches in
`-mode refine`):

| top peaks | +-ppm | depth | pairs searched | DIA-NN in the set | at q <= 0.01 (decoys) | recovered | wall |
|---|---|---|---|---|---|---|---|
| 1,000 | 10 | 4 | 200,000 of 872,036 | 7,779 (84.0 %) | 4,404 (43) | 53.6 % | 16:17 |
| **1,000** | **10** | **5** | **100,728** | **7,698 (83.1 %)** | **4,755 (46)** | **58.4 %** | 11:56 with tuning |
| 1,000 | 10 | 6 | 10,314 | 5,899 (63.7 %) | 5,567 (54), aborted by `max_target_fraction` (54 % of targets) | 88.8 % | 4:14 |
| 300 | 10 | 4 | 123,865 | 7,240 (78.2 %) | 4,460 (43) | 58.4 % | 12:47 |
| 2,000 | 7 | 5 | 137,130 | 8,058 (87.0 %) | 4,550 (44) | 53.4 % | 12:47 |

Probe only (pairs / DIA-NN in the set): 1,000 peaks at 20 ppm, depth 5:
200,000 of 508,424 / 71.2 %; 100 peaks at 10 ppm, depth 5: 6,153 / 48.9 %;
2,000 at 5 ppm, depth 4: 200,000 of 649,345 / 84.2 %; 500 at 7 ppm, depth 4:
166,142 / 83.6 %. Wider tolerances and more peaks raise chance matches faster
than signal; fewer peaks lose real fragments. (A probe sweep counts qualifying
spectra at the lowest depth it evaluates, so its capped rows at higher depths
rank slightly differently from the product; the full-search rows used probes
at the searched depth and match the product's set exactly.)

**timsTOF diaPASEF, whole-proteome library with entrapment** (5,754,771
precursors, half of them shuffled-twin entrapment, library ratio r = 1.00;
DIA-NN 2.0 identifies 116,267 = 2.0 % in the run at q <= 0.01, combined
entrapment FDP 1.06 %; 5,135,763 eligible pairs, 650 targets without a decoy).
Ion mobility is still ignored, as in M1: each 25-Th window is split into two
1/K0 maps, and every compound is extracted in both.

| | M1: random 200,000, `chunk 100000` | random 200,000, `chunk 20000` (M3 chunking) | M3 defaults |
|---|---|---|---|
| searched pairs / DIA-NN in the set | 200,000 / 4,554 | 200,000 / 4,554 | 200,000 of 226,175 / 38,364 (33.0 %) |
| calibration | 32 points from 2,160 seeds, window 483 s | the same | 1,077 points from 1,917 seeds (1,742 evidence), LOWESS (span 0.15): window 661 s against the line's 756 s |
| extraction | 4 chunks, 6,174 s | 20 chunks over 52-56 of 64 maps, 5,109 s (231-277 s each) | 20 chunks over 40-50 maps, 5,976 s |
| target precursors at q <= 0.01 (decoys) | 1,031 (9) | 1,031 (9): the report's rows identical to M1's | 32,193 (320); 29,068 peptides, 4,820 protein groups |
| in DIA-NN's report / recovered | 97.8 % / 22.1 % | the same | 96.1 % / 80.6 % |
| median / p99 \|dRT\| | 0.018 / 0.100 min | the same | 0.019 / 0.172 min |
| null-pair balance / self-checks | 25,111 : 24,889 / 0, 0 | the same | 25,574 : 24,426 (z 5.1, warns) / 0, 0 |
| entrapment, combined (lower) | 1.55 % (0.78 %), 8 of 1,031 | the same | 1.49 % (0.75 %), 240 of 32,193 |
| tune (RT): cohorts; TEST calibrated SD | 212 / 122 units; 1.71 -> 1.26 min | the same | 5,936 / 3,554 units; 2.16 -> 1.17 min |
| whole invocation, peak RSS | 1:54:38, 21.0 GB | 1:36:41, 18.8 GB | 1:56:22, 18.1 GB |

With M1's contiguous chunks at the default `chunk 20000`, chunk 1 took 857 s
over 1-2 maps (projected 4.8 h in all); the dealt-out chunks take 231-330 s
over 40-56 maps, as fast per precursor as M1's 100,000-precursor chunks at a
fifth of their size. Extraction dominates: 8 OpenSWATH threads at 12-15 ms
per precursor, twice over for the two 1/K0 halves (M2). The prefilter took
25.8 s (decoys 6.9 s, index 0.5 s and 514 MB for 61.6 M fragments, sweep
17.0 s over 236,052 spectra); the cache took 78.1 GB, 2.56x the mzML.
The entrapment numbers are a quick look with the entrapment design's
peptide classes (T, E; shared excluded) and r = 1. The acceptance read them
as "not shown to be anti-conservative" (lower bound minus 2 SE below 1 %)
and credited isomeric twins with part of the excess (leaving out the 47
entrapment identifications whose target twin co-elutes gave 1.20 %). The
M3 review corrected both: the lower bound assumes no false null-target
discoveries and cannot show anti-conservativeness, and the direct test
does -- entrapment pairs were won by the target 240 : 168 at precursor
q <= 0.01 (z 3.6), 223 : 150 at peptide level and 44 : 25 at protein level;
leaving co-eluting twins out of BOTH sides (the decoys of those pairs
co-elute with the twin too) still gives 196 : 150 (z 2.5). With r = 1
confirmed in the searched set (80,524 entrapment against 80,963 real
DIA-NN-absent pairs kept), the excess is real: gate (a) failed on timsTOF
by a direct exchangeability test, independent of isomer screening, and had
E wins equalled decoy wins the combined FDP would have been 1.05 %. The
cause was the decoys' inherited fragment slots (Architecture, 1). The null-pair balance warning comes with 38,364 of DIA-NN's
precursors among 200,000 searched targets: the lowest quarter of pair
winners is no longer only null pairs.

The probe grid on this run agrees with Astral on the defaults. At 1,000 peaks
and 10 ppm the cap keeps 38,800 DIA-NN precursors at depth 5 (226,175 union),
37,097 at depth 4 (2,270,482) and 25,490 at depth 6 (31,442, no cap). At
20 ppm it keeps 32,009 at depth 5. At 3,000 and 10,000 peaks it keeps 6-18 %,
no better than chance: frame-merged diaPASEF spectra are dense. Depth 6 with
`search:max_target_fraction 1` and `search:entrapment_tag` gave 25,314
identifications (252 decoys), 96.2 % in DIA-NN's report and 95.5 % of the
DIA-NN precursors in its set recovered, in 33:33 at 9.3 GB. Its entrapment
estimate is 1.19 % (lower 0.60 %). Its null-pair balance, 5,188 : 2,672
(z 28), shows the premise failing rather than the error rate.

**The depth, and the guards.** On the Astral run, identifications RISE as
the searched set gets smaller and purer, because every null pair searched is
one more chance for a decoy to outscore a real target at the 1 % cut. Depth 6
on Astral, run again with `search:max_target_fraction 1`, gave 5,567
identifications (54 decoys). 94.1 % of them are in DIA-NN's report, and
88.8 % of the DIA-NN precursors in its set were recovered. It took 4:33 at
5.7 GB.

On the timsTOF run, depth 6 gave fewer identifications than depth 5
(25,314 against 32,193), because the searched set lost too much. Its
entrapment estimate (1.19 %) is no worse than depth 5's (1.49 %), so the
selective prefilter shows no sign of loosening the error rate.

Depth 6 is still not the default. It breaks the premise of two safeguards,
both of which assume that most candidates are null:

- more than half of its scored targets pass q <= 0.01 (54 % on Astral, 80 %
  on timsTOF), so the `max_target_fraction` guard aborts;
- the lowest quarter of pair winners is no longer null, so the null-pair
  balance warns (z 3.2 and z 28).

Neither observation shows that the FDR is wrong, since the selection is
label-symmetric. But neither guard can then tell a real problem from a rich
candidate set, and even the default warns on the rich timsTOF run (z 5.1).
Depth 5 is the default: on Astral it gives the most identifications while
every guard keeps its meaning, and on timsTOF it keeps the most DIA-NN
precursors.

Making depth 6, or a depth adapted to the run, the default needs two things:

- guards stated relative to the prefilter's own evidence. For example, the
  excess of passing targets over passing decoys estimates the present
  targets, which no honest search can exceed by much. It is 6,644 of 10,314
  pairs at depth 6 on Astral and 14,759 of 100,728 at depth 5; on timsTOF it
  is 27,124 of 31,442 at depth 6;
- the entrapment campaign of M5.

### After the M3 review: predicted decoys (`search:intensities predicted`)

Measured with the review-fix build on the same node (16 threads; up to six
jobs of this measurement at once, so wall times are inflated). Report
analysis: `tools/ana.py`, `tools/dn.py` in `/scratch/kohlbach/bid3r-m3r-fix/`,
and the acceptance's `entrap_score.py` and `cmp.py` unchanged. "E pairs" are
pairs whose target is an entrapment peptide of the design table (known
absent); their winners are counted target : decoy.

**Astral, entrapment library** (8,045,757 precursors; M3 defaults otherwise):

| | M3 acceptance (library rule) | this build, `intensities library` | this build, `predicted` |
|---|---|---|---|
| pairs searched | 191,035 | 191,035 | 200,000 of 210,144 (decoys now pass as often as targets: 109,455 against 112,486 at depth 5) |
| calibration | 744 points, linear | 916 points, LOWESS (CV error 2.20 against 2.58) | 908 points, LOWESS |
| target precursors at q <= 0.01 (decoys) | 3,772 (36) | 3,964 (38) | 3,119 (30) |
| combined FDP, precursors / peptides / protein groups | 1.43 % / 1.58 % / 2.64 % | 1.56 % / 1.66 % / 2.62 % | **0.90 % / 0.65 % / 1.97 %** (2 entrapment of 203 groups) |
| E-pair winners, precursor q <= 0.01 / 0.1 | 27 : 18 / 439 : 347 (z 3.3) | 31 : 19 / 467 : 376 (z 3.1; the product's winner test warns) | **14 : 18 / 280 : 311** (z -0.7 / -1.3) |
| E-pair winners, peptide / protein group, q <= 0.1 | 375 : 301 / 44 : 36 | 383 : 319 / 44 : 37 | 240 : 278 / 25 : 35 |
| prefilter, E pairs: depth >= 5 target / decoy; sign test | 50,955 / 45,874 (1.11); z 78.8 | as the acceptance | 50,926 / 54,035 (0.94); z -25.0 |
| null-pair balance | z -1.9 | 23,869 : 23,889 | 25,130 : 24,870 |
| E among the calibration seeds | 437 of 1,754 | 446 of 1,964 | 448 of 1,964 |
| wall / peak RSS | 15:53 / 9.9 GB | 19:03 / 10.1 GB | 30:59 / 13.0 GB (14 min predicting 6.2 M decoys) |

The winner test, the direct exchangeability check, is balanced at every
level with predicted decoys, and the combined FDP is 0.90 % at nominal 1 %.
With the library rule the same build reproduces the acceptance's target
favour (467 : 376, z 3.1), which the new in-product winner test flags while
the null-pair balance does not. The prefilter's E-pair sign test is no longer
target-favoured but now decoy-favoured (a 6 % excess of decoys at depth 5):
conservative, not anti-conservative, and not yet explained; it costs
identifications (3,119 against the library rule's 3,964, most of them the
library rule's excess). The entrapment seeds are twins of present peptides
(isomers that share their fragments): the robust calibration fit drops their
points, but the "beat your own decoy" rule does not exclude them.

**Astral, whole-proteome library** (`-mode tune -tune_heads rt`):

| | M3 acceptance | predicted |
|---|---|---|
| library check | - | 2,000 of 2,000 sampled targets are the model's own prediction: decoys only predicted (the shortcut, since removed; 0 of 99,931 of this library's pairs change assay without it, so the counts above stand) |
| pairs searched / DIA-NN in the set | 100,728 / 7,698 (83.1 %) | 110,530 / 7,704 (83.2 %) |
| calibration | 944 points from 1,927 seeds, linear, 317 s | 1,141 points from 2,143 seeds (99 of 100 bins), LOWESS (CV 2.16 against 2.54), 325 s |
| target precursors at q <= 0.01 (decoys) | 4,755 (46) | 4,423 (43) |
| in DIA-NN's report / recovered overall / within the searched set | 94.5 % / 48.5 % / 58.4 % | 94.3 % / 45.0 % / 54.1 % |
| median \|dRT\| | 0.006 min | 0.006 min |
| tuning (RT): TEST calibrated SD | 0.910 -> 0.500 min | 1.016 -> 0.486 min |
| prefilter / extraction / whole invocation, peak RSS | 12 s / 373 s / 12:39, 17.1 GB | 451 s (418 s predicting 3.08 M decoys) / 439 s / 20:45, 17.9 GB |

**timsTOF (diaPASEF, RT only), whole-proteome library with shuffled-twin
entrapment** (5,754,771 precursors, `-mode tune -tune_heads rt`). The library
is another predictor's, so both members are predicted:

| | M3 acceptance | predicted |
|---|---|---|
| prefilter, depth >= 5: targets / decoys | 144,970 / 89,754 | 129,212 / 91,046 |
| pairs searched | 200,000 of 226,175 | 200,000 of 212,958 |
| calibration | 1,077 points, LOWESS, window 661 s | 1,238 points from 2,134 seeds (100 of 100 bins), LOWESS (CV error 1.37 against the line's 1.52), window 675 s |
| target precursors at q <= 0.01 (decoys) | 32,193 (320) | 30,018 (299) |
| peptides / protein groups | 29,068 / 4,820 | 27,063 / 4,168 |
| in DIA-NN's report | 96.1 % | 97.0 % (29,102 of 30,018) |
| combined FDP, precursors / peptides / protein groups | 1.49 % / 1.54 % / 1.83 % | **0.91 % / 0.89 % / 0.82 %** (DIA-NN on the same design: 1.06 / 1.20 / 1.28 %) |
| E-pair winners, precursor q <= 0.01 / 0.1 | 240 : 168 (z 3.6) / 2,594 : 2,078 (z 7.6) | **137 : 158 (z -1.2) / 1,933 : 1,916 (z 0.3)** |
| E-pair winners, peptide / protein group, q <= 0.01 | 223 : 150 / 44 : 25 | 120 : 142 / 17 : 25 |
| targets at a matched combined FDP of 1 % | q* 0.0067: 30,895 | q* 0.0111: 30,204 |
| prefilter, E pairs: depth >= 5 target / decoy; sign test | 1.18; z 81 | 47,617 / 48,245 (0.99); z -13.1 |
| null-pair balance | 25,574 : 24,426 (z 5.1, warns) | 25,286 : 24,714 (z 2.6) |
| E among the calibration seeds | 70 of 1,771 | 54 of 1,991 |
| tuning (RT): TEST calibrated SD | 2.161 -> 1.167 min | 2.149 -> 1.126 min |
| wall / peak RSS | 2:04:42 / 17.6 GB | 3:15:40 / 25.6 GB (24 min predicting both members of 5.1 M pairs; extraction 9,140 s against 6,467 s at 2 to 4 times the node load, same 400,000 precursors) |

The run that settled gate (a) is now the run that passes it: the entrapment
winner test is flat at every level and every threshold (the largest |z| over
precursors, peptides and protein groups at q <= 0.01, 0.02, 0.05 and 0.1 is
2.2), and the three combined FDP figures fall from 1.49 / 1.54 / 1.83 % to
0.91 / 0.89 / 0.82 %, below DIA-NN's on the same design. It costs
identifications: 30,018 against 32,193 at q <= 0.01, and 30,204 against
30,895 where the entrapment estimate itself reads 1 % -- 2.2 % of them, the
price of q-values that are now slightly conservative rather than 1.5x
optimistic. The timsTOF targets also lose their library's own predicted
intensities (the built-in model's replace them), which the Astral runs do
not.

### The library-assay shortcut, removed

A commit on this branch (`d0ec486`) gave `search:intensities predicted` a
shortcut: a sample of 2,000 library targets was re-predicted, and when they
all carried exactly the model's own top fragments -- a library DIALibGen
generated with the same model, instrument and NCE -- the targets KEPT their
library assays and only the decoys were predicted, half the model's work.
It is gone. Both members are predicted, always, by the same model, from
their own sequences, under the same rule.

**What it broke.** The pair's count is min(the target's library transition
count, max(3, the fewer of the two members' fragments above 1e-4 of their
base peak)). Reading the target's assay out of the library puts the LIBRARY's
fragment cap where the TARGET's above-floor count belongs: the count is then
set by the decoy's above-floor count alone, and wherever the target's own
prediction is the shorter of the two, the pair gets more fragments than the
rule allows: the target's extra ones are library fragments its own
prediction puts BELOW its floor, while the decoy's are all above its own.
Why that ends up favouring the target is not settled; that it does is
measured, below. Assays that agree fragment for fragment are not enough --
the rule that CHOOSES how many has to read both members the same way.

**The probe** (`probes/tims_counts.json` under `/scratch/kohlbach/bid4-label/`,
99,921 pairs of the timsTOF entrapment library, targets and decoys predicted
in one batch): 20,737 targets get a different assay under the shortcut than
from their own prediction, and in EVERY one of them the shortcut's is the
longer (count bigger 20,715 times, smaller 0). That library holds 12
fragments per precursor (mean 12.0) where the model puts a mean 16.9 above
the floor -- but fewer than 12 for 41,081 of the pairs, which is where the
counts part: mean 10.05 fragments per pair with the shortcut against 9.31
without. On the Astral whole-proteome library the same probe finds 0 of
99,931 (mean 11.6 library fragments against 20.6 above the floor, never
fewer), which is why no Astral run ever showed this.

**The measurement.** Four `-mode tune -tune_heads rt` searches of the same
timsTOF run against the same shuffled-twin entrapment library (200,000 pairs,
`-search:entrapment_tag ENTRAP_`), two by two: the shortcut against both
members predicted, crossed with the decoy rule that keeps one terminal
residue against the two the branch settled on. Column three is the
"predicted" column of the timsTOF table above. Every figure is the run's own
in-product entrapment block, so the four are comparable with each other; the
acceptance script quoted above differs in the last decimal (0.93 % here
against its 0.91 % for the same run). The two new runs (columns one and
four) ran side by side on one node, so their wall times are inflated alike.

| | **both predicted, keep 2** (the rule) | shortcut, keep 2 | both predicted, keep 1 | shortcut, keep 1 |
|---|---|---|---|---|
| target precursors at q <= 0.01 (decoys) | 28,843 (287) | 25,930 (258) | 30,018 (299) | 27,210 (271) |
| peptides / protein groups | 25,922 / 4,175 | 23,463 / 3,646 | 27,063 / 4,168 | 24,631 / 3,338 |
| combined entrapment FDP | **0.74 %** | 1.30 % | 0.93 % | **3.04 %** |
| E-pair winners, precursor q <= 0.01 | **107 : 126** (z -1.24) | 168 : 114 (z 3.22) | 140 : 159 (z -1.10) | **414 : 136** (z 11.85) |
| E-pair winners, peptide / protein group, q <= 0.01 | 96 : 113 / 15 : 26 | 157 : 99 / 28 : 23 | 123 : 146 / 17 : 25 | 380 : 123 / 45 : 17 |
| E-pair winners, precursor q <= 0.1 | 1,891 : 1,802 (z 1.46) | - | - | 3,663 : 2,119 (z 20.3) |
| wall / peak RSS | 3:06 h / 19.0 GB | - | - | 2:50 h / 18.9 GB |

Read along each row: with both members predicted the winner test is flat
(|z| <= 1.3 at q <= 0.01, at every level) and the FDP sits below the nominal
1 %, under either decoy rule. With the shortcut the targets win (z 3.2
keeping two terminal residues, z 11.9 keeping one) and the FDP rises to
1.30 % and 3.04 %. The shortcut is the cause; the two-residue decoy rule is
not -- the winner test is flat with it and without it, and it has its own
measurement (Architecture, decoys). Nor did the shortcut buy
identifications: 25,930 against 28,843 at q <= 0.01 keeping two residues.
What it bought was 16 minutes of a three-hour run (2:50 h against 3:06 h,
the decoy rule differing too), for the entrapment gate.

**The 2,000-target sample check was not evidence.** On this very library it
reported 2,000 of 2,000 sampled targets carrying exactly the model's
prediction, while 21 % of that library's pairs come out with a different
assay. The check compares the library's fragment LIST against the model's
top fragments; it passes whenever the library holds a prefix of the model's
ranking, and says nothing about where the target's floor sits, which is what
the count reads. A check that could see this would have to compare the
target's above-floor count with its library count for every precursor -- and
computing that count IS predicting the target, so there is nothing left to
save. And a sample of 2,000 out of millions can only ever speak for the
precursors in it: what it reports about a library is not a property the
search may then rely on for every pair.

What stays from that commit: the fixed-point fragment-range test (a double
bound from `fromFixed` could sit a rounding above the library's own extreme
fragment and drop it), the prediction progress log, and releasing the model
once the searched set has its assays. `search:intensities library`, the old
asymmetric rule, stays as the comparison the exchangeability test measures
`predicted` against -- and only as that: every run that selects it is warned,
before the run is read, in the log and in the report's `search.warnings`.

## Formats

First release: centroided DIA **mzML**. timsTOF diaPASEF needs a frame-merged
mzML with a per-peak 1/K0 array (`mzpeak-convert --to mzml`; stock OpenMS 3.5.0
has no Bruker reader). A file with window limits but no per-peak array is usable
for RT only, with a warning. Later: native Bruker `.d` (opentims-based reader or
an OpenMS upgrade), and `.mzpeak` on POSIX builds.

## Milestones

- **M1** End-to-end slice on Orbitrap/Astral data: CLI, paired random subset,
  in-memory decoys, sequence-blind assays, stock calibration, chunked
  extraction, LDA/FDR port with the fixes, report hand-off; synthetic-run
  tests. Cluster check against a DIA-NN report of the same run.
- **M2** diaPASEF through mzML: ion-mobility window, scores and calibration.
- **M3** Evidence prefilter and evidence-seeded calibration. Done: the
  prefilter (default), evidence seeds with a LOWESS second fit, chunks spread
  over the windows, OpenMS's own parse errors, the report kept on
  `search:min_ids`, the cache-size log; measured on Astral and timsTOF
  (Resources). Its review added predicted decoys (`search:intensities
  predicted`), the entrapment winner test and the tag rule for UniProt ids,
  seeds from the central RT range that beat their decoy, and LOWESS chosen
  by cross-validation (Resources, "After the M3 review"). The library-assay
  shortcut that followed it was measured and removed (Resources, "The
  library-assay shortcut, removed").
- **M4** Streaming store with a memory budget; 16 GB / 8-core gate.
- **M5** Honesty and purpose campaign: entrapment, concordance, tuning transfer.
- **M6** Desktop app, README figure (identification becomes a DIALibGen step),
  release smoke on all platforms.
- **M7** Second pass, fragment intensities (`-write_intensity`), GBT opt-in.
- **M8** Native Bruker `.d`, `.mzpeak`.

## Validation gates

- **(a) Honesty:** entrapment on the timsTOF test runs and on an Astral run:
  combined FDP <= 1.5 % at nominal 1 % for precursors, peptides and proteins,
  reported with its lower bound; AND the entrapment winner test (pairs whose
  target is entrapment won by target and decoy equally often: binomial
  |z| < 3 at q <= 0.01 and q <= 0.1, at precursor, peptide and protein level),
  which needs no database ratio and no isomer screening; and the E-pair
  depth sign test of the prefilter (|z| < 3). Identification counts and the
  internal decoy rate never decide this. The entrapment peptides must be screened against the
  sample's full proteome -- with I = L, isomers (same composition and termini)
  and missed-cleavage joins -- or come from a foreign proteome: shuffled twins
  that keep each tryptic piece's composition are isomeric with real peptides,
  share their b(n-1)/y(n-1) ions and precursor m/z, and some are real peptides,
  which inflates the combined FDP (a first measurement found 20 of 24 clean
  trap identifications had an identified isomeric parent).
- **(b) Concordance** with DIA-NN 2.0 on the same library and run: >= 60 % of
  its precursors recovered within the searched set; median |dRT| <= 0.1 min and
  |d1/K0| <= 0.01 on shared precursors.
- **(c) Purpose:** models tuned on built-in identifications are within 10 % of
  models tuned on a DIA-NN report, on held-out proteins and on sibling runs.
- **(d) Resources:** 16 GB memory limit, 8 cores: no out-of-memory, <= 60 min.
- **(e) Self-checks** on every run; **(f)** determinism across threads and chunk
  sizes, and the `-ids` route byte-identical to 0.11.0.

## Decisions taken (defaults, revisable)

1. mzML first; timsTOF by conversion. Native `.d` decided after M4 numbers.
2. Hand-off through the report file, not an in-memory interface.
3. Refine with a capped candidate set keeps only confirmed candidates; the
   documented built-in workflow is `-mode refine -tune -no_filter`.
4. Shuffle decoys built in memory, sequence-blind scoring, no `mutate`.
5. `-run` ships as experimental and becomes the documented default route only
   after gates (a)-(d) pass on an Orbitrap and a diaPASEF run.
6. Candidates by fragment evidence (M3): the 6 most intense predicted
   fragments, a spectrum's 1,000 most intense peaks, +-10 ppm, depth 5, cap
   200,000 pairs stratified by window x library-RT decile x charge and ranked
   by the better member's evidence. Chosen on the Astral run (the most
   identifications with every guard meaningful) and the best of the timsTOF
   probe grid too; depth 6 finds more on Astral and waits for guards stated
   relative to the prefilter's evidence and for M5 (Resources, M3).
7. LOWESS replaces the linear RT map when it predicts held-out points
   better (5-fold cross-validation on one point set, error below 0.97x the
   line's) with at least 200 points; until the M3 review, when it narrowed
   the window by 10 %. Stock `TransformationModelLowess`; ODIA's PAVA/LOESS
   code is not ported.
8. Both members of every search pair predicted by one MS2 model, each with
   its own most intense fragments (`search:intensities predicted`, M3
   review) -- always, including targets whose library assay already IS that
   prediction: keeping it instead reads the library's fragment cap into the
   pair's count rule ("The library-assay shortcut, removed"). The decoy in
   its target's slots (`library`) stays as the comparison the exchangeability
   test measures against, and warns loudly on every run that selects it.

## Main risk

The stock OpenSWATH + LDA stack has not yet produced a validated identification
set on diaPASEF data. Yield (the tuner needs at least 100 units per held-out
cohort, ideally thousands of identifications) and the honesty of the anchors are
unknown until M2 and M5; on the Astral run above, the M1 default clears the
cohort floor 4.5-fold with a FASTA-subset library, and with M3 a whole-proteome
library clears it too (Astral: test 887 / val 329 units, 58 % of DIA-NN's
precursors in the searched set recovered where gate (b) asks 60 %, 54 % with
predicted decoys; timsTOF, ion mobility still ignored: 5,936 / 3,554 units,
81 % recovered). The prefilter is the one selection step that reads the run
before the decoys are scored: it is label-symmetric by construction and by
test -- but label symmetry of the SELECTION is not enough: the M3 review
found the decoys themselves weaker than null targets (their inherited
fragment slots), which paired competition cannot survive. Predicted decoys
fixed the entrapment winner test on Astral; the prefilter's evidence is now
slightly decoy-favoured, which is conservative and not yet understood. The
entrapment campaign of M5 must confirm the error rate on screened
entrapment. A biased anchor set would silently degrade every later
search of the tuned library, which is why every selection step is label-
symmetric, every failure aborts, and the option stays experimental until
entrapment confirms the error rate.
