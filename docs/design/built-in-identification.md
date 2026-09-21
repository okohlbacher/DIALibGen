# Design: built-in identification for `refine` and `tune`

Status: **accepted, in implementation** (milestone M1). The option stays
*experimental* until the honesty gate (a) below passes.

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
- TOPP subsection `search:` with, among others: `candidates random|evidence`,
  `subset`, `max_pairs 200000`, `decoys shuffle|pseudo_reverse`, `seed`,
  `passes`, `rt_window`, `mz_ppm`, `im_window` (0 = automatic), `ms1`,
  `rt_im_scores`, `calibration_min_rsq 0.70`, `calibration_min_coverage 0.30`,
  `readoptions`, `memory_gb`, `min_ids 200`, `entrapment_tag`, `selftest true`.
- Refused with `-run` until later milestones: `-write_intensity`,
  `-empirical_library`, `-min_fragments`, `-write_im` and `-tune_heads
  ccs|both` (M1 measures no 1/K0). Everything that would fail after the search
  (tuning recipe, models, a library without protein groups, a directory such
  as a Bruker `.d` given as `-run`) is checked before the run is read; a
  failure after the search keeps the report and says how to reuse it.
- The provenance sidecar gains `inputs.run` and a `search` block: settings,
  calibration, candidate counts, targets and decoys at each gate, resources.
  The refined Parquet library embeds it without the resources (timings,
  threads), so the same command writes the same bytes.
- Every guard aborts with counts in the message. Nothing falls back silently.

## Architecture

New static target `odia_search` (`include/odia/search/`, `src/search/`). It links
only OpenMS and OpenSwathAlgo, which `odia_library` already links, so release
bundles do not grow.

1. **Library.** Loaded as today and never mutated; a side index orders
   precursors by m/z. Decoys for the search are built in memory (shuffle by
   default); decoys already in the file are ignored for the search and left
   untouched. The DIA-NN-derived `mutate` method is not offered, and neither
   is `reverse`, which moves the tryptic C-terminus and makes decoys separable
   by construction. A decoy must differ from its target in fragment m/z only:
   a target whose interior cannot be rearranged has no decoy (no residue
   substitution as a fallback), every decoy fragment lies in the m/z range of
   the library's target fragments, and an arrangement that reproduces the
   target's fragment masses (I/L) is re-drawn. Each rule reads the target
   only, so a failure removes a whole pair.
2. **Run.** `SwathFile::loadMzML` (`normal` or `cache`). Ion-mobility window
   limits are normalised (lower/upper swapped where reversed, as in current
   mzpeak-convert output). diaPASEF is detected from window limits *and* a
   per-peak 1/K0 array. Later (M4): a DIALibGen-owned streaming store that keeps
   only the peaks inside active assay boxes, so memory stays bounded.
3. **Candidates.** M1: a deterministic, paired random subset (the draw hashes the
   key target and decoy share, so it is label-blind), selected after the run is
   loaded: a target whose precursor m/z lies in none of the run's isolation
   windows is ineligible before the draw and the cap (its decoy shares the m/z;
   19.5 % of a random draw on an Astral run). M3: an evidence prefilter
   (top-6 predicted fragments co-occurring in one spectrum), applied to targets
   and decoys with the same rule and keeping a pair when either member passes.
4. **Assays.** Per chunk, an `OpenSwath::LightTargetedExperiment` built directly
   from the library arrays. Never the whole library. Peptide sequences are left
   empty for targets and decoys, which makes OpenSWATH's scoring sequence-blind
   and therefore symmetric: DIALibGen decoys carry their target's sequence.
5. **Calibration.** `OpenSwathCalibrationWorkflow::performRTNormalization` on seed
   assays (iRT/CiRT kit peptides plus sampled library precursors in M1,
   evidence seeds from M3), with thresholds suited to predicted libraries. It
   returns the RT transformation, m/z and 1/K0 windows and the ion-mobility
   correction. Failure aborts. What limits it in M1 is the library, not the
   cutoffs: in the Astral reference run DIA-NN identified 2.9 % of a
   FASTA-derived library's precursors, so ~60 of 2,122 random seeds are
   findable at all, and the stock quality cutoffs kept 29-51 of them in five
   runs. The linear fit through them leaves a median residual of ~30 s, and
   the RT window (0.99 quantile x 1.3, doubled) is 240-320 s. Evidence-seeded calibration (M3)
   and a non-linear second pass are the improvements; more random seeds would
   add points only in proportion to calibration memory, which grows with seeds
   x fragments x spectra per window and is the run's peak (3.5 GiB on the
   6.4 GB run, 1.5 GiB above extraction).
6. **Extraction.** `OpenSwathWorkflow::performExtraction` (stock 13-argument form)
   per chunk, with an inactive OSW writer and in-memory features. After each
   chunk the features become compact score rows and are freed. Stock OpenSWATH
   runs on at most 8 OpenMP threads: its feature scoring serialises on the
   process-wide MetaInfoRegistry lock, and more threads only spin (extraction
   was twice as slow at 32 threads as at 8).
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
- What tests exchangeability: the **null-pair balance**. Among complete pairs
  whose winner is in the lowest quarter of winning scores (pairs where nothing
  is present), targets and decoys must win about equally often; logged,
  recorded, and a warning at |z| > 3.
- Run-level guards, each aborting: decoy:target ratio out of band, more than
  half of the targets at q <= 0.01, too few identifications (`search:min_ids`),
  failed calibration, and the self-checks (`search:selftest`, on by default):
  scoring with swapped labels and with each pair's labels exchanged on a coin
  must identify (almost) nothing. They catch a classifier that leaks labels.
  They do NOT catch decoys built weaker than null targets -- random pair labels
  symmetrise the construction away (reverse decoys, whose null pairs won 1.66:1
  by the target, passed both with 0 and 0). That is the null-pair balance's job.
- Decoy construction and what is known about it (Astral run): targets and
  decoys have identical precursor m/z, RT, charge, protein group, fragment
  slots and intensities; extraction finds peak groups for 199,999 targets and
  200,000 decoys of 200,000 pairs, and the null-pair balance is 25,016 : 24,983. A review measured a small residual asymmetry on
  null pairs of an entrapment search, before classification: null targets beat
  their shuffle decoys on the MS2 fragment sub-scores in 51.2-51.5 % of pairs
  (z 7-9 over 105,433 pairs), persisting with every decoy fragment in range.
  Decoys copy their target's predicted-intensity slot choice; re-predicting
  intensities for the decoy sequence is the candidate fix and is not done in
  M1. Until the balance on null pairs is near 0.5, `-run` stays experimental.

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
- **M3** Evidence prefilter and evidence-seeded calibration.
- **M4** Streaming store with a memory budget; 16 GB / 8-core gate.
- **M5** Honesty and purpose campaign: entrapment, concordance, tuning transfer.
- **M6** Desktop app, README figure (identification becomes a DIALibGen step),
  release smoke on all platforms.
- **M7** Second pass, fragment intensities (`-write_intensity`), GBT opt-in.
- **M8** Native Bruker `.d`, `.mzpeak`.

## Validation gates

- **(a) Honesty:** entrapment on the timsTOF test runs and on an Astral run:
  combined FDP <= 1.5 % at nominal 1 % for precursors, peptides and proteins,
  reported with its lower bound. Identification counts and the internal decoy
  rate never decide this. The entrapment peptides must be screened against the
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

## Main risk

The stock OpenSWATH + LDA stack has not yet produced a validated identification
set on diaPASEF data. Yield (the tuner needs at least 100 units per held-out
cohort, ideally thousands of identifications) and the honesty of the anchors are
unknown until M2 and M5; on the Astral run above, the M1 default clears the
cohort floor 4.5-fold with a FASTA-subset library, and a whole-proteome library
would not without M3. A biased anchor set would silently degrade every later
search of the tuned library, which is why every selection step is label-
symmetric, every failure aborts, and the option stays experimental until
entrapment confirms the error rate.
