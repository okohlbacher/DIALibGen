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
DIALibGen -mode tune -in predicted.parquet -run run.mzML -out tuned.tsv -threads 8
DIALibGen -mode refine -tune -no_filter -in predicted.parquet -run run.mzML -out adapted.tsv -write_im
DIALibGen -mode refine -in predicted.parquet -ids report.parquet -out refined.tsv   # unchanged
```

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
  `subset`, `max_pairs`, `decoys shuffle|pseudo_reverse|reverse`, `seed`,
  `passes`, `rt_window`, `mz_ppm`, `im_window` (0 = automatic), `ms1`,
  `rt_im_scores`, `calibration_min_rsq 0.70`, `calibration_min_coverage 0.30`,
  `readoptions`, `memory_gb`, `min_ids 200`, `entrapment_tag`, `selftest`.
- Refused with `-run` until later milestones: `-write_intensity`,
  `-empirical_library`, `-min_fragments`.
- The provenance sidecar gains `inputs.run` and a `search` block: settings,
  calibration, candidate counts, targets and decoys at each gate, resources.
- Every guard aborts with counts in the message. Nothing falls back silently.

## Architecture

New static target `odia_search` (`include/odia/search/`, `src/search/`). It links
only OpenMS and OpenSwathAlgo, which `odia_library` already links, so release
bundles do not grow.

1. **Library.** Loaded as today and never mutated; a side index orders
   precursors by m/z. Decoys for the search are built in memory (shuffle by
   default); decoys already in the file are ignored for the search and left
   untouched. The DIA-NN-derived `mutate` method is not offered.
2. **Run.** `SwathFile::loadMzML` (`normal` or `cache`). Ion-mobility window
   limits are normalised (lower/upper swapped where reversed, as in current
   mzpeak-convert output). diaPASEF is detected from window limits *and* a
   per-peak 1/K0 array. Later (M4): a DIALibGen-owned streaming store that keeps
   only the peaks inside active assay boxes, so memory stays bounded.
3. **Candidates.** M1: a deterministic, paired random subset (the draw hashes the
   key target and decoy share, so it is label-blind). M3: an evidence prefilter
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
   correction. Failure aborts.
6. **Extraction.** `OpenSwathWorkflow::performExtraction` (stock 13-argument form)
   per chunk, with an inactive OSW writer and in-memory features. After each
   chunk the features become compact score rows and are freed.
7. **Scoring and FDR** (below), then the report.

Shared code with ODIA: the dependency-free LDA and FDR headers are vendored into
`src/odia-core/` with a manifest of origin commit and hashes, and a CI check.
ODIA's DIA-NN-shaped network classifier, its Bruker/mzPeak readers and its
patched-OpenMS prefilter are not ported.

## Scoring and FDR

- Features: about twelve non-collinear OpenSWATH sub-scores, plus MS1 and ion-
  mobility scores when available; elution-model and ion-series scores off.
- Classifier: ODIA's semi-supervised LDA (mProphet/pyProphet lineage), three
  folds, three iterations, with four changes: folds assigned per target-decoy
  pair; fold scores rescaled label-blind (median/MAD of all held-out scores);
  RT and 1/K0 deviation scores hidden from the calibration discriminant; GBT and
  network code removed.
- Precursor level: concatenated target-decoy competition. Each pair keeps its
  better member (ties to the decoy); q = (D + 1)/T, monotonised, no pi0.
  Pairing is structural by (modified sequence, charge).
- Peptide and protein level: roll-up by canonical sequence and by the library's
  `Protein.Group` string, picked competition at protein level. The existing
  `refine`/`tune` gates stay the control.
- Run-level guards, each aborting: decoy:target ratio out of band, more than
  half of the targets at q <= 0.01, too few identifications (`search:min_ids`),
  failed calibration. `search:selftest` runs a label-swap and a pure-null check.

## Resources

Measured ODIA figures show the cost is dominated by materialising the library as
OpenSWATH objects and by holding the run in memory. This design never does the
first, bounds the candidate set (`search:max_pairs`), extracts once by default,
and (from M4) streams the run. Projected on 8-12 laptop cores: about 15 min for a
6 GB Orbitrap Astral mzML and 20-30 min for a 29 GB diaPASEF mzML, with 4-10 GB
peak memory. These figures are derived, not measured; M1 and M4 measure them.

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

- **(a) Honesty:** entrapment with shuffled twin proteins on the timsTOF test
  runs and on an Astral run: combined FDP <= 1.5 % at nominal 1 % for
  precursors, peptides and proteins. Identification counts and the internal
  decoy rate never decide this.
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
unknown until M2 and M5. A biased anchor set would silently degrade every later
search of the tuned library, which is why every selection step is label-
symmetric, every failure aborts, and the option stays experimental until
entrapment confirms the error rate.
