# Parameter reference

Generated from DIALibGen 0.11.0 with `scripts/update-parameters.py`. Do not edit this table by hand.

The schema lists all modes. `-mode generate` is the default; refinement and training options apply to their respective modes. Generation booleans take `true`/`false`; list items are separate arguments. See [usage](usage.md).

| Option | Type | Default | Constraints | Description |
|---|---|---|---|---|
| `-mode` | string | `generate` | generate,refine,tune | generate: FASTA to library; refine: apply observed values; tune: train RT/CCS models and re-predict the whole library. |
| `-in` | input-file | (empty) |  | Protein FASTA (generate) or spectral library (refine/tune). |
| `-out` | output-file | (empty) |  | DIA-NN spectral library; Parquet embeds provenance. |
| `-config` | input-file | (empty) |  | Optional JSON configuration for the selected mode. Explicit CLI/INI generation settings override JSON. |
| `-write_config` | output-file | (empty) |  | Write the effective mode configuration and exit. |
| `-irt_standards` | input-file | (empty) |  | iRT calibration standards; defaults to the bundled table. |
| `-ids` | input-file | (empty) |  | DIA-NN report.parquet, or a pre-filtered library with -empirical_library. Modification names are canonicalised. |
| `-run` | input-file | (empty) |  | EXPERIMENTAL, in development: a centroided DIA run (mzML) that the built-in identification step searches instead of reading -ids. Give exactly one of -ids and -run. Settings: search:. |
| `-out_ids` | output-file | (empty) |  | With -run: the identification report (Parquet, DIA-NN column names), which then serves as -ids. Default: <out>.ids.parquet. Never overwritten. |
| `-out_report` | output-file | (empty) |  | Per-axis residual report (TSV), measured BEFORE the overwrite. |
| `-no_filter` | bool | `false` |  | Keep precursors the reference did not identify; requires -no_write_rt or successful RT tuning so RT units stay consistent. |
| `-empirical_library` | bool | `false` |  | Declare -ids a pre-filtered empirical library rather than a report: gates whose columns are absent are BYPASSED and each bypass is recorded. Without this, a missing gate column is an error. |
| `-no_write_rt` | bool | `false` |  | Keep predicted RT rather than replacing it with observed RT. |
| `-write_im` | bool | `false` |  | Also overwrite 1/K0 with the observed value for charges >= -im_min_charge; off by default. |
| `-write_intensity` | bool | `false` |  | Replace predicted fragment intensities with the reference run's observed ones. Needs Fragment.Info/Fragment.Quant.Raw (1.9) or Fr.N.Id/Quantity (2.x), exported with --report-lib-info or --export-quant respectively. Every match is cross-checked on fragment m/z, and a run that replaces nothing is an error. |
| `-intensity_min_correlation` | double | `0.0` |  | Require fragment quality greater than this: correlation in DIA-NN 1.9, Score in DIA-NN 2.x. -1 disables the quality gate. |
| `-intensity_no_restrict` | bool | `false` |  | Replace a precursor only when EVERY one of its transitions is trusted, else keep its predictions whole. Neutral-loss transitions cannot match report fragments, so precursors carrying them keep their predictions. Transition counts cannot change. |
| `-intensity_no_rerank` | bool | `false` |  | Keep a replaced precursor's transitions in their original order. |
| `-intensity_min_fragments` | int | `3` | 0: | A replaced precursor keeps at least this many transitions or keeps its predictions whole. |
| `-intensity_norm` | string | `library_max` | library_max,base_peak,sum,raw | How an observed area becomes a library intensity. library_max scales to the maximum that precursor already held, which preserves the file's own convention -- base peak = 1 is not an invariant of these libraries. |
| `-intensity_min_relative` | double | `1.0e-04` |  | Drop observed intensities below this fraction of the kept maximum. |
| `-intensity_mz_tol_ppm` | double | `20.0` |  | A fragment matches only when the report's m/z agrees with the library's to within this. |
| `-intensity_max_mz_mismatch` | double | `0.01` |  | Refuse when more than this fraction of identity matches fail the m/z check. 1 surveys a suspect pairing instead of failing. |
| `-allow_mixed_intensity` | bool | `false` |  | Permit -write_intensity with -no_filter, which leaves observed and predicted intensities in one library. Recorded in the provenance. |
| `-q_precursor` | double | `0.01` |  | Precursor q-value gate (>= 1 disables). |
| `-q_global` | double | `0.01` |  | Global/peptide q-value gate (>= 1 disables). |
| `-q_protein` | double | `0.01` |  | Protein q-value gate (>= 1 disables). |
| `-min_fragments` | int | `0` | 0: | Minimum DISTINCT reference fragments per precursor; needs fragment identities in -ids. 0 = off. |
| `-rt_unit` | string | `observed` | observed,minmax | observed = the reference's own units; minmax = rescaled to 0..100 over the matched set (refused with -no_filter). |
| `-dedup` | string | `lowest_q` | lowest_q,highest_evidence | Which observation wins when a precursor was seen more than once: lowest_q (ties: lower PEP, then first seen) or highest_evidence (needs an Evidence column). |
| `-im_min_charge` | int | `2` |  | Charges below this never receive an observed 1/K0 (z1 is censored at the ramp top on timsTOF diaPASEF). |
| `-im_ramp_top` | double | `0.0` |  | The instrument's mobility ramp top, if known; observations within -im_ramp_margin of it are treated as censored. 0 = unknown. |
| `-im_ramp_margin` | double | `0.02` |  | See -im_ramp_top. |
| `-min_match_fraction` | double | `0.0` |  | Refuse when fewer than this fraction of passing reference precursors match the library. 0 = refuse only when nothing matches. |
| `-tune` | bool | `false` |  | Fine-tune the RT and CCS models on -ids and re-predict the whole library through them BEFORE refining, so precursors the reference never identified are corrected too. Needs a build with the fine-tuning stage. Turns a seconds-long refinement into a training run plus whole-library inference. |
| `-tune_models` | string | (empty) |  | Directory holding the stock peptdeep_{rt,ccs,ms2}_dynamic.onnx: -tune starts from the RT and CCS models, and -run predicts both members of every searched target-decoy pair with the MS2 model. Default: $DIALIBGEN_MODEL_DIR or the bundled models. |
| `-tune_heads` | string | `both` | rt,ccs,both | Which models to tune. |
| `-tune_out_models` | string | (empty) |  | Keep the tuned ONNX files and their .tune.json and .trajectory.tsv sidecars here. Default: a scratch directory, removed on exit -- the deliverable is the library. |
| `-tune_predict_gpu` | bool | `false` |  | Use the GPU for the re-prediction pass (the ONNX one, not training). |
| `-tune_predict_sessions` | int | `0` | 0: | ONNX Runtime sessions for the re-prediction pass; 0 = default. |
| `-tune_keep_free_cysteine_offset` | bool | `false` |  | Keep DIALibGen's free-cysteine RT offset when re-predicting with a TUNED RT model. Off by default: that offset was fitted against the STOCK model, and a model tuned on this run's own identifications has had the chance to learn the effect itself -- applying both counts it twice. |
| `-log` | string | (empty) |  | Name of log file (created only when specified) |
| `-debug` | int | `0` |  | Sets the debug level |
| `-threads` | int | `1` |  | Sets the number of threads allowed to be used by the TOPP tool |
| `-no_progress` | bool | `false` |  | Disables progress logging to command line |
| `-force` | bool | `false` |  | Overrides tool-specific checks |
| `-test` | bool | `false` |  | Enables the test mode (needed for internal use only) |
| `-generation:ccs_model` | string | (empty) |  | ccs model |
| `-generation:decoys` | string | `none` | none,mutate,pseudo_reverse,reverse,shuffle | decoys |
| `-generation:derive_ion_mobility` | string | `true` | true,false | derive ion mobility |
| `-generation:enzyme` | string | `Trypsin/P` |  | enzyme |
| `-generation:fixed_modifications` | string list | `Carbamidomethyl (C)` |  | fixed modifications |
| `-generation:fragment_mz` | double list | `200, 1800` |  | fragment mz |
| `-generation:fragments` | int list | `3, 12` |  | fragments |
| `-generation:free_cysteine_rt_correction` | string | `true` | true,false | free cysteine rt correction |
| `-generation:instrument` | string | `QE` |  | instrument |
| `-generation:irt_rescale` | bool | `false` |  | irt rescale |
| `-generation:max_fragment_charge` | int | `2` | 1:2 | max fragment charge |
| `-generation:max_variable_modifications` | int | `1` | 0: | max variable modifications |
| `-generation:min_relative_intensity` | double | `1.0e-04` | 0.0: | min relative intensity |
| `-generation:missed_cleavages` | int | `1` | 0: | missed cleavages |
| `-generation:ms2_model` | string | (empty) |  | ms2 model |
| `-generation:n_terminal_methionine_excision` | string | `true` | true,false | n terminal methionine excision |
| `-generation:nce` | double | `-1.0` | -1.0:100.0 | Normalized collision energy (>0 and <=100). -1 selects the instrument default. |
| `-generation:peptide_length` | int list | `7, 30` |  | peptide length |
| `-generation:precursor_charges` | int list | `1, 2, 3, 4` | 1:8 | precursor charges |
| `-generation:precursor_mz` | double list | `350, 1200` |  | precursor mz |
| `-generation:recompute_decoy_mz` | bool | `false` |  | recompute decoy mz |
| `-generation:reserved_doubly_charged` | int | `0` | 0: | reserved doubly charged |
| `-generation:rt_model` | string | (empty) |  | rt model |
| `-generation:variable_modifications` | string list | (empty) |  | variable modifications |
| `-filter:q_value` | double | `0.01` |  | Precursor Q.Value threshold |
| `-filter:min_charge` | int | `2` | 1:8 | CCS only: lowest precursor charge used (RT collapses all charges). z1 is censored at the mobility ramp top on timsTOF |
| `-filter:allow_z1` | bool | `false` |  | CCS: permit filter:min_charge 1 (censored observations enter training) |
| `-filter:rt_spread_max` | double | `0.2` |  | Drop a sequence whose charge states' RTs span more than this (minutes) |
| `-filter:rt_max_minutes` | double | `0.0` |  | rt_norm denominator; 0 = the run's maximum observed RT |
| `-cohort:train_size` | int | `0` | 0: | Training units (sequences for rt, sequence x charge for ccs); 0 = full pool |
| `-cohort:train_frac` | double | `0.0` |  | Alternative to train_size: fraction of the pool |
| `-cohort:full_fit` | bool | `false` |  | Train on EVERY unit of the run, the test and validation cohorts included: the closest fit the data allows. val and TEST are then in-sample; judge this mode by searching a DIFFERENT run with the result, never by its own numbers |
| `-cohort:no_inner_val` | bool | `false` |  | No inner validation cohort: its units rejoin the pool and checkpoints are selected on TEST, which is then no longer a held-out number |
| `-train:epochs` | int | `100` | 1: | Horizon of the cosine schedule (and the maximum epochs) |
| `-train:warmup` | int | `10` | 0: | Linear warmup epochs |
| `-train:lr` | double | `1.0e-04` |  | Peak learning rate (Adam) |
| `-train:batch_size` | int | `1024` | 1: | Batch size within a length group |
| `-stop:eval_every` | int | `1` | 1: | Validate every n epochs |
| `-stop:min_epochs` | int | `20` | 0: | Never stop before this epoch |
| `-stop:patience` | int | `10` | 0: | Stop after this many epochs without progress (never before max(min_epochs, warmup)) |
| `-stop:rel_tol` | double | `5.0e-03` |  | Progress = the selection metric beats the anchor by this fraction |
| `-stop:abs_tol` | double | `0.0` |  | Progress = beats the anchor by this absolute amount (0 = use rel_tol) |
| `-stop:max_seconds` | double | `0.0` |  | Wall-clock budget per head, checked at epoch boundaries; final evaluation/export may exceed it (0 = none) |
| `-stop:select` | string | `calibrated_sd` | calibrated_sd,rmse | Selection metric on the validation cohort |
| `-machine:device` | string | `cpu` |  | cpu or cuda[:N] |
| `-machine:threads` | int | `4` | 1: | CPU threads used for training |
| `-machine:no_cudnn` | bool | `false` |  | CUDA: do not use cuDNN (needed when only its loader shim is installed, as in pytorch.org's libtorch zips); slower |
| `-machine:seed` | int | `20260803` | 0: | Seed for the training subsample and batch order |
| `-search:candidates` | string | `evidence` | evidence,random | Candidate selection: evidence = target-decoy pairs of which the target or the decoy has fragment evidence in the run (search:prefilter_*; the same rule for both); random = a deterministic, label-blind random subset of pairs (search:subset) |
| `-search:subset` | int | `0` | 0: | search:candidates random: targets drawn from the library before the cap (0 = every eligible target) |
| `-search:max_pairs` | int | `200000` | 0: | Cap on the target-decoy pairs searched; time is linear in pairs (0 = no cap). With evidence candidates the cap is stratified by isolation window, library-RT decile and charge and keeps the pairs with the most evidence |
| `-search:prefilter_depth` | int | `5` | 1:6 | search:candidates evidence: keep a pair when its target or its decoy has this many of its 6 most intense predicted fragments among one MS2 spectrum's top peaks |
| `-search:prefilter_top_peaks` | int | `1000` | 1: | search:candidates evidence: the most intense peaks of each MS2 spectrum the prefilter matches |
| `-search:prefilter_ppm` | double | `10.0` | 0.1:100.0 | search:candidates evidence: fragment match tolerance, ppm either side |
| `-search:decoys` | string | `shuffle` | shuffle,pseudo_reverse | How the search's in-memory decoys are built from the selected targets; both methods keep the termini. Decoys in the library file are not searched and stay in the output |
| `-search:intensities` | string | `predicted` | predicted,library | Fragments and intensities of both members of each target-decoy pair: predicted = both predicted by the PeptDeep MS2 model (tune_models, $DIALIBGEN_MODEL_DIR or the bundled models) from their own sequences, each taking its own most intense fragments; library = the library's target assay, the decoy in its target's fragment slots with its target's intensities (makes decoys weaker than null targets; for comparison only) |
| `-search:instrument` | string | `auto` |  | search:intensities predicted: the MS2 model's instrument (QE, Lumos, timsTOF, SciexTOF, ThermoTOF or an alias); auto = the instrument DIALibGen generate recorded in the library, else timsTOF for an ion-mobility run and QE otherwise |
| `-search:nce` | double | `-1.0` | -1.0:100.0 | search:intensities predicted: the MS2 model's collision energy; -1 = the one recorded in the library with its instrument, else the instrument's default |
| `-search:seed` | int | `42` | 0: | Salt of the candidate draw: changes which pairs are searched, not how |
| `-search:passes` | int | `1` | 1:1 | Extraction passes (1 in this version) |
| `-search:rt_window` | double | `0.0` | 0.0: | Full RT extraction window, seconds (0 = from the calibration) |
| `-search:mz_ppm` | double | `0.0` | 0.0: | Full fragment m/z extraction window, ppm (0 = automatic) |
| `-search:im_window` | double | `0.0` | -1.0: | Full 1/K0 extraction window on ion-mobility (diaPASEF) runs: 0 = automatic, from the run's 1/K0 calibration; > 0 = this width; -1 = off (searched by m/z and RT only, no observed 1/K0: -write_im and the CCS head are refused) |
| `-search:ms1` | string | `true` | true,false | Extract MS1 traces and use the MS1 sub-scores |
| `-search:rt_im_scores` | string | `true` | true,false | Let the RT and 1/K0 deviation sub-scores into the classifier; false is an ablation for tuning, which exists to correct those deviations |
| `-search:calibration_min_rsq` | double | `0.7` | 0.0:1.0 | Least r^2 of the RT calibration on the seed assays |
| `-search:calibration_min_coverage` | double | `0.3` | 0.0:1.0 | Least fraction of the seed assays found in the run that the RT outlier removal must keep |
| `-search:allow_bootstrap` | bool | `false` |  | When the RT calibration fails, map the library RT range linearly onto the run instead of aborting. A test hook, recorded in the provenance |
| `-search:readoptions` | string | `auto` | auto,normal,cache | How the run is held: normal = in memory, cache = per-window cache files, auto = cache for runs above 3 GB |
| `-search:cache_dir` | string | (empty) |  | Directory for cache files (default: the directory of -out_ids); they are removed after the search |
| `-search:min_ids` | int | `200` | 0: | Abort when fewer target precursors pass q <= 0.01 |
| `-search:max_target_fraction` | double | `0.5` | 0.01:1.0 | Abort when more than this fraction of the scored target precursors passes q <= 0.01: no honest decoy set looks like that |
| `-search:report_max_q` | double | `0.1` | 0.01:1.0 | Precursors, targets and decoys, up to this precursor q-value go into -out_ids |
| `-search:entrapment_tag` | string | (empty) |  | Prefix of entrapment protein ids, at the start of an id or right after a '\|' (sp\|ENTRAP_P12345\|...): log and record the combined entrapment FDP estimate and the entrapment winner test. A validation aid |
| `-search:selftest` | string | `true` | true,false | Also score with swapped and with random pair labels, and abort unless both identify (almost) nothing. These catch a classifier that leaks labels, not decoys that are built weaker than null targets |
| `-search:batch_size` | int | `500` | 1: | Advanced: precursors per extraction batch within one isolation window |
| `-search:chunk` | int | `20000` | 2: | Advanced: precursors per extraction call; target-decoy pairs stay together. Each chunk takes batch-sized pieces from across the m/z range, so OpenSWATH, which parallelises over isolation windows, has work for its threads; it bounds extraction memory only, not the calibration's |
