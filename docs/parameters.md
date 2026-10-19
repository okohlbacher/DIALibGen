# Parameter reference

Generated from DIALibGen 0.12.1 with `scripts/update-parameters.py`. Do not edit this table by hand.

The schema lists all modes. Without `-mode`, a FASTA generates and a library converts; refinement and training options apply to their respective modes. `--help` shows the options marked in the column *--help* (`--help -mode <mode>` one mode's), `--helphelp` all of them. Generation booleans take `true`/`false`; list items are separate arguments. `-threads` defaults to all physical cores unless the command line or the INI sets it. See [usage](usage.md).

| Option | Type | Default | Constraints | --help | Description |
|---|---|---|---|---|---|
| `-mode` | string | `generate` | generate,refine,tune,convert,prepare | yes | Generate: FASTA to library; refine: write a run's observed RT (and 1/K0) into a library; tune: fine-tune RT (and CCS) on a run and re-predict the whole library; convert: change the format; prepare (macOS): precompile the GPU model. Without -mode: generate for a FASTA, convert for a library. |
| `-in` | input-file | (empty) |  | yes | Protein FASTA (generate) or spectral library (refine, tune, convert) |
| `-out` | output-file | (empty) |  | yes | Library to write: .tsv for DIA-NN; .pqp for OpenSWATH (add -generation:decoys); .parquet for DIALibGen (DIA-NN cannot read it). Other layouts: -out_format |
| `-out_format` | string | `auto` | auto,dialibgen_parquet,diann_tsv,diann2_parquet,openswath_pqp,openswath_tsv |  | Layout of -out. auto: by extension (.parquet dialibgen_parquet, .tsv diann_tsv, .pqp openswath_pqp). diann2_parquet (.parquet): DIA-NN 2.0's own Parquet library layout. openswath_tsv (.tsv): OpenSWATH's assay TSV, as TargetedFileConverter reads it |
| `-config` | input-file | (empty) |  | yes | JSON recipe for the mode (-write_config writes one); options on the command line or in the INI override it |
| `-write_config` | output-file | (empty) |  |  | Write the mode's effective recipe (JSON) and exit |
| `-models` | string | (empty) |  |  | Directory of the stock peptdeep_{rt,ms2,ccs}_dynamic.onnx models for generate, -tune and -run. Default: $DIALIBGEN_MODEL_DIR, else the bundled models. -generation:{rt,ms2,ccs}_model name single files |
| `-irt_standards` | input-file | (empty) |  |  | Generate: iRT standards for -generation:irt_rescale true; default: the bundled table |
| `-instrument_file` | input-file | (empty) |  |  | Generate: the instrument's defaults (NCE, ion mobility, m/z windows, charges) from this YAML file. Otherwise one file per instrument from $DIALIBGEN_INSTRUMENT_DIR, then the user's DIALibGen/instruments directory (~/.config, ~/Library/Application Support or %APPDATA%), then the shipped share/DIALibGen/instruments; default.yaml for the rest (docs/usage.md#instrument-defaults) |
| `-ids` | input-file | (empty) |  | yes | Refine, tune: the DIA-NN report.parquet of one run (or, with -empirical_library, a pre-filtered library) |
| `-no_filter` | bool | `false` |  | yes | Refine: keep the precursors the run did not identify (needs -tune or -no_write_rt, so that RT stays in one unit; JSON key filter) |
| `-no_write_rt` | bool | `false` |  | yes | Refine: keep the predicted RT (JSON key write_rt) |
| `-write_im` | bool | `false` |  | yes | Refine: also write the observed 1/K0 (ion-mobility runs), for charges >= -im_min_charge |
| `-write_intensity` | bool | `false` |  | yes | Refine: also write the observed fragment intensities (a report exported with --report-lib-info (DIA-NN 1.9) or --export-quant (2.x)); fragments are cross-checked on m/z, and a run that replaces nothing is an error |
| `-q_precursor` | double | `0.01` |  | yes | Refine: precursor q-value cut-off (>= 1 disables); -q_global and -q_protein add gates (--helphelp) |
| `-tune` | bool | `false` |  | yes | Refine: fine-tune RT (and CCS) on -ids and re-predict the whole library first, so precursors the run did not identify are corrected too (a training run plus whole-library prediction) |
| `-tune_heads` | string | `auto` | auto,rt,ccs,both | yes | Tune: models to fine-tune; auto = RT, plus CCS when the report has enough ion-mobility data |
| `-tune_out_models` | string | (empty) |  | yes | Tune: also keep the tuned models (ONNX, with their .tune.json and .trajectory.tsv) in this directory |
| `-run` | input-file | (empty) |  |  | Experimental: a centroided DIA run (mzML, or an mzPeak archive or its unpacked directory) that the built-in identification searches instead of reading -ids; not validated yet (-ids is). Settings: search: |
| `-out_ids` | output-file | (empty) |  |  | With -run: the identification report it writes (Parquet, DIA-NN column names), then used as -ids. Default: <out>.ids.parquet. Never overwritten |
| `-out_report` | output-file | (empty) |  |  | Refine, tune: per-axis residuals before the overwrite (TSV) |
| `-empirical_library` | bool | `false` |  |  | Refine: -ids is a pre-filtered empirical library, not a report: gates whose columns are missing are bypassed and each bypass is recorded (without it a missing gate column is an error; JSON key require_gates) |
| `-q_global` | double | `1.0` |  |  | Refine: global (peptide-level) q-value cut-off; 1 = off |
| `-q_protein` | double | `1.0` |  |  | Refine: protein q-value cut-off; 1 = off |
| `-min_fragments` | int | `0` | 0: |  | Refine: least distinct fragments the report identified per precursor (needs fragment identities in -ids); 0 = off |
| `-rt_unit` | string | `observed` | observed,minmax |  | Refine: observed = the run's own RT units; minmax = rescaled to 0..100 over the matched set (not with -no_filter or tuning) |
| `-dedup` | string | `lowest_q` | lowest_q,highest_evidence |  | Refine: which observation wins for a precursor seen more than once: lowest_q (ties: lower PEP, then first seen) or highest_evidence (needs an Evidence column) |
| `-min_match_fraction` | double | `0.0` |  |  | Refine: refuse when fewer than this fraction of the passing report precursors are in the library; 0 = only when none is |
| `-im_min_charge` | int | `2` |  |  | Refine, with -write_im: charges below this keep their predicted 1/K0 (z1 is censored at the ramp top on timsTOF diaPASEF) |
| `-im_ramp_top` | double | `0.0` |  |  | Refine, with -write_im: the mobility ramp's top, if known; observations within -im_ramp_margin of it are censored. 0 = unknown |
| `-im_ramp_margin` | double | `0.02` |  |  | Refine, with -write_im: see -im_ramp_top |
| `-intensity_min_correlation` | double | `0.0` |  |  | Refine, with -write_intensity: least fragment quality (correlation in DIA-NN 1.9, Score in 2.x); -1 disables the gate |
| `-intensity_no_restrict` | bool | `false` |  |  | Refine, with -write_intensity: replace a precursor only when every transition is trusted, else keep its predictions (neutral-loss transitions never match; JSON key intensity_restrict) |
| `-intensity_no_rerank` | bool | `false` |  |  | Refine, with -write_intensity: keep a replaced precursor's transitions in their order (JSON key intensity_rerank) |
| `-intensity_min_fragments` | int | `3` | 0: |  | Refine, with -write_intensity: a replaced precursor keeps at least this many transitions, or keeps its predictions |
| `-intensity_norm` | string | `library_max` | library_max,base_peak,sum,raw |  | Refine, with -write_intensity: how an observed area becomes an intensity; library_max scales to the precursor's previous maximum (the file's own convention) |
| `-intensity_min_relative` | double | `1.0e-04` |  |  | Refine, with -write_intensity: drop observed intensities below this fraction of the kept maximum |
| `-intensity_mz_tol_ppm` | double | `20.0` |  |  | Refine, with -write_intensity: a fragment matches only when the report's m/z agrees with the library's within this |
| `-intensity_max_mz_mismatch` | double | `0.01` |  |  | Refine, with -write_intensity: refuse when more than this fraction of matches fail the m/z check; 1 surveys a suspect pairing |
| `-allow_mixed_intensity` | bool | `false` |  |  | Refine: permit -write_intensity with -no_filter, which leaves observed and predicted intensities in one library (recorded) |
| `-tune_predict_sessions` | int | `0` | 0: |  | Tune: ONNX Runtime sessions for the re-prediction; 0 = -threads |
| `-tune_keep_free_cysteine_offset` | bool | `false` |  |  | Tune: keep the free-cysteine RT offset when re-predicting with the tuned RT model (off: the offset was fitted to the stock model, and a model tuned on the run has learned the effect) |
| `-log` | string | (empty) |  |  | Name of log file (created only when specified) |
| `-debug` | int | `0` |  |  | Sets the debug level |
| `-threads` | int | `0` |  | yes | Number of threads; 0 = all physical cores (the default when the command line and the INI do not set it) |
| `-no_progress` | bool | `false` |  |  | Disables progress logging to command line |
| `-force` | bool | `false` |  |  | Overrides tool-specific checks |
| `-test` | bool | `false` |  |  | Enables the test mode (needed for internal use only) |
| `-generation:instrument` | string | `QE` |  | yes | Instrument of the run: Astral, Exploris, Lumos, QE, timsTOF, or a vendor name such as 'Orbitrap Astral' or 'timsTOF Pro'. Its defaults (instrument files, see --helphelp -instrument_file) set the NCE, ion mobility (timsTOF) and m/z windows and charges (Astral, timsTOF). Default: default.yaml's (QE) |
| `-generation:nce` | double | `-1.0` | -1.0:100.0 | yes | Collision energy; -1 = the instrument's (Lumos and Astral 25, others 30) |
| `-generation:precursor_mz` | double list | (empty) |  | yes | Precursor m/z window; match it to the acquisition method. Default: the instrument's (Astral 380-980), else 350-1200 |
| `-generation:fragment_mz` | double list | (empty) |  | yes | Fragment m/z window. Default: the instrument's (Astral 150-2000), else 200-1800 |
| `-generation:precursor_charges` | int list | `1, 2, 3, 4` | 1:8 | yes | Precursor charges. Default: the instrument's (timsTOF 2-4), else 1-4 |
| `-generation:fixed_modifications` | string list | `Carbamidomethyl (C)` |  | yes | OpenMS names, e.g. 'Carbamidomethyl (C)'; give the option without a value for none |
| `-generation:variable_modifications` | string list | (empty) |  | yes | OpenMS names, e.g. 'Oxidation (M)'; at most -generation:max_variable_modifications (1) per peptide |
| `-generation:enzyme` | string | `Trypsin/P` |  | yes | Protease; Trypsin/P also cuts before P |
| `-generation:missed_cleavages` | int | `1` | 0: | yes | Missed cleavages per peptide |
| `-generation:peptide_length` | int list | `7, 30` |  |  | Peptide length window, in residues |
| `-generation:max_variable_modifications` | int | `1` | 0: |  | Most variable modifications on one peptide |
| `-generation:n_terminal_methionine_excision` | string | `true` | true,false |  | Also digest each protein without its initiator methionine |
| `-generation:fragments` | int list | `3, 12` |  |  | Fragments per precursor: fewer than min drops the precursor, the most intense up to max are kept |
| `-generation:max_fragment_charge` | int | `2` | 1:2 |  | Highest fragment charge (1 or 2) |
| `-generation:min_relative_intensity` | double | `1.0e-04` | 0.0: |  | Drop fragments below this fraction of the precursor's most intense one |
| `-generation:reserved_doubly_charged` | int | `0` | 0: |  | Fragment slots reserved for doubly charged ions; 0 ranks all fragments by predicted intensity |
| `-generation:predict_ion_mobility` | string | `auto` | auto,true,false |  | Predict CCS (and 1/K0, see derive_ion_mobility). auto: only for an instrument with ion mobility (timsTOF); true: always; false: never, the library then carries neither |
| `-generation:derive_ion_mobility` | string | `true` | true,false |  | Also write 1/K0, derived from the predicted CCS (when CCS is predicted) |
| `-generation:mobility_correction` | string | `auto` | auto,off |  | auto: correct the derived 1/K0 with the instrument file's mobility_correction (timsTOF; only for the CCS model it was fitted for); off: the plain Mason-Schamp 1/K0 of the predicted CCS |
| `-generation:irt_rescale` | bool | `false` |  |  | Rescale RT to iRT with the standards of -irt_standards; false: RT is the model's raw 0..1 scale (DIA-NN fits its own RT map) |
| `-generation:free_cysteine_rt_correction` | string | `true` | true,false |  | Add the measured RT offset of unmodified cysteines (nothing to do when cysteines are alkylated) |
| `-generation:decoys` | string | `none` | none,mutate,pseudo_reverse,reverse,shuffle |  | Decoys in the library: none (DIA-NN makes its own), mutate, pseudo_reverse, reverse, shuffle |
| `-generation:recompute_decoy_mz` | bool | `false` |  |  | Recompute a decoy's precursor m/z from its decoy sequence (with -generation:decoys) |
| `-generation:rt_model` | string | (empty) |  |  | RT model (ONNX); default: peptdeep_rt_dynamic.onnx in -models |
| `-generation:ms2_model` | string | (empty) |  |  | MS2 model (ONNX); default: peptdeep_ms2_dynamic.onnx in -models |
| `-generation:ccs_model` | string | (empty) |  |  | CCS model (ONNX); default: peptdeep_ccs_dynamic.onnx in -models |
| `-filter:q_value` | double | `0.01` |  |  | Precursor q-value cut-off of the training set |
| `-filter:min_charge` | int | `2` | 1:8 |  | CCS: lowest precursor charge (RT pools all charges; z1 is censored at the ramp top on timsTOF) |
| `-filter:allow_z1` | bool | `false` |  |  | CCS: permit filter:min_charge 1 (censored observations then train) |
| `-filter:rt_spread_max` | double | `0.2` |  |  | Drop a sequence whose charge states' RTs span more than this (minutes) |
| `-filter:rt_max_minutes` | double | `0.0` |  |  | RT normalisation; 0 = the run's latest observed RT |
| `-cohort:train_size` | int | `0` | 0: |  | Training units (sequences for RT, sequence x charge for CCS); 0 = the whole pool |
| `-cohort:train_frac` | double | `0.0` |  |  | Instead of train_size: the fraction of the pool |
| `-cohort:full_fit` | bool | `false` |  |  | Train on every unit of the run, test and validation included: validation and test are then in-sample, so judge the result on a different run |
| `-cohort:no_inner_val` | bool | `false` |  |  | No inner validation cohort: its units rejoin the pool and checkpoints are selected on the test cohort, which is then no longer held out |
| `-train:epochs` | int | `100` | 1: |  | Length of the cosine schedule (and the most epochs) |
| `-train:warmup` | int | `10` | 0: |  | Linear warm-up epochs |
| `-train:lr` | double | `1.0e-04` |  |  | Peak learning rate (Adam) |
| `-train:batch_size` | int | `1024` | 1: |  | Batch size within a length group |
| `-stop:eval_every` | int | `1` | 1: |  | Validate every n epochs |
| `-stop:min_epochs` | int | `20` | 0: |  | Never stop before this epoch |
| `-stop:patience` | int | `10` | 0: |  | Stop after this many epochs without progress (not before max(min_epochs, warmup)) |
| `-stop:rel_tol` | double | `5.0e-03` |  |  | Progress: the selection metric beats the best so far by this fraction |
| `-stop:abs_tol` | double | `0.0` |  |  | Progress: ... by this absolute amount (0 = use rel_tol) |
| `-stop:max_seconds` | double | `0.0` |  |  | Time budget per head, checked between epochs (the final evaluation may exceed it); 0 = none |
| `-stop:select` | string | `calibrated_sd` | calibrated_sd,rmse |  | Selection metric on the validation cohort |
| `-machine:device` | string | `auto` |  | yes | Auto: the GPU when one is usable and pays off (NVIDIA: CUDA, RT and CCS with cuDNN 9; Apple: CoreML for MS2 of large libraries, MPS for fine-tuning), else the CPU; cpu: CPU only (the reference); also cuda[:N], coreml, mps |
| `-machine:coreml_units` | string | `gpu` | gpu,all,ane,cpu |  | CoreML compute units: gpu (CPU and GPU, fastest), all (Core ML's choice, Neural Engine included), ane (CPU and Neural Engine), cpu (checks) |
| `-machine:train_threads` | int | `0` | 0: |  | CPU threads for fine-tuning (0: 16 on Linux and Windows, 8 on macOS). On Linux the tuned weights depend on the count (recorded); on macOS they do not |
| `-machine:cores` | string | `all` | all,performance |  | MacOS: which cores -threads 0 counts: all, or performance (saves energy, not time). Does not change the output |
| `-machine:no_cudnn` | bool | `false` |  |  | CUDA fine-tuning without cuDNN (when only its loader shim is installed, as in pytorch.org's libtorch zips); slower |
| `-machine:seed` | int | `20260803` | 0: |  | Seed of the training subsample and batch order |
| `-machine:tf32` | bool | `false` |  |  | TensorFloat-32 on NVIDIA GPUs, for prediction and fine-tuning: faster, but about 1 % of precursors select a different fragment set (off: GPU and CPU agree to float rounding) |
| `-machine:gpu_sessions` | int | `4` | 1:16 |  | CUDA: ONNX Runtime sessions per GPU and model, each with its own stream; CoreML: concurrent batches. Does not change the output |
| `-machine:gpu_devices` | int list | (empty) | 0: |  | CUDA devices to spread the sessions over (e.g. 0 1). Default: the first that works. Does not change the output |
| `-machine:require_gpu` | bool | `false` |  |  | Fail instead of predicting on the CPU when -machine:device asks for a GPU (or auto finds one) that a model cannot use |
| `-search:candidates` | string | `evidence` | evidence,random |  | Candidate pairs: evidence = target-decoy pairs of which either member has fragment evidence in the run (search:prefilter_*); random = a deterministic, label-blind random subset (search:subset) |
| `-search:subset` | int | `0` | 0: |  | Random candidates: targets drawn before the cap (0 = every eligible target) |
| `-search:max_pairs` | int | `200000` | 0: |  | Most target-decoy pairs searched; time is linear in pairs (0 = no cap). Evidence candidates are capped per isolation window, RT decile and charge, keeping the most evidence |
| `-search:prefilter_depth` | int | `5` | 1:6 |  | Evidence candidates: a pair is kept when either member has this many of its 6 most intense predicted fragments among one spectrum's top peaks |
| `-search:prefilter_top_peaks` | int | `1000` | 1: |  | Evidence candidates: the most intense peaks per MS2 spectrum that the prefilter matches |
| `-search:prefilter_ppm` | double | `10.0` | 0.1:100.0 |  | Evidence candidates: fragment tolerance, ppm either side |
| `-search:decoys` | string | `shuffle` | shuffle,pseudo_reverse |  | How the search builds its in-memory decoys (both keep the termini); decoys in the library are not searched and stay in the output |
| `-search:intensities` | string | `predicted` | predicted,library |  | Validation only: fragments of both pair members. predicted = predicted from each member's own sequence by the MS2 model (-models); library = the target's assay, the decoy in its target's fragment slots (weaker decoys) |
| `-search:instrument` | string | `auto` |  |  | Predicted intensities: the MS2 model's instrument; auto = the one the library records, else timsTOF for an ion-mobility run and QE otherwise |
| `-search:nce` | double | `-1.0` | -1.0:100.0 |  | Predicted intensities: the MS2 model's collision energy; -1 = the one the library records with its instrument, else the instrument's |
| `-search:seed` | int | `42` | 0: |  | Salt of the candidate draw: changes which pairs are searched, not how |
| `-search:rt_window` | double | `0.0` | 0.0: |  | Full RT extraction window, seconds (0 = from the calibration) |
| `-search:mz_ppm` | double | `0.0` | 0.0: |  | Full fragment m/z extraction window, ppm (0 = automatic) |
| `-search:im_window` | double | `0.0` | -1.0: |  | Full 1/K0 extraction window on ion-mobility runs: 0 = from the run's 1/K0 calibration; > 0 = this width; -1 = off (no observed 1/K0: -write_im and the CCS head are refused) |
| `-search:ms1` | string | `true` | true,false |  | Extract MS1 traces and use the MS1 sub-scores |
| `-search:rt_im_scores` | string | `true` | true,false |  | Validation only: false keeps the RT and 1/K0 deviation sub-scores out of the classifier (an ablation) |
| `-search:calibration_min_rsq` | double | `0.7` | 0.0:1.0 |  | Least r^2 of the RT calibration on the seed assays |
| `-search:calibration_min_coverage` | double | `0.3` | 0.0:1.0 |  | Least fraction of the seed assays found in the run that the RT outlier removal must keep |
| `-search:readoptions` | string | `auto` | auto,normal,cache |  | How the run is held: normal = in memory, cache = per-window cache files, auto = cache for runs above 3 GB |
| `-search:cache_dir` | string | (empty) |  |  | Directory of the cache files (default: that of -out_ids); removed after the search |
| `-search:ms1_in_memory` | string | `auto` | auto,on,off |  | Memory: hold the MS1 map in memory during extraction (auto = when the free memory is at least twice its size; on; off). Changes speed and memory, never the report |
| `-search:min_ids` | int | `200` | 0: |  | Abort when fewer target precursors pass q <= 0.01 |
| `-search:max_target_fraction` | double | `0.5` | 0.01:1.0 |  | Abort when more than this fraction of the scored targets passes q <= 0.01 (no honest decoy set looks like that) |
| `-search:report_max_q` | double | `0.1` | 0.01:1.0 |  | Precursors (targets and decoys) up to this q-value go into -out_ids |
| `-search:entrapment_tag` | string | (empty) |  |  | Validation only: prefix of entrapment protein ids (at the start or after a '\|'); logs and records the entrapment FDP estimate and winner test |
| `-search:selftest` | string | `true` | true,false |  | Also score with swapped and with random pair labels and abort unless both identify (almost) nothing (catches a classifier that leaks labels) |
| `-search:batch_size` | int | `100` | 1: |  | Memory: precursors per extraction batch within one isolation window (the report does not depend on it) |
| `-search:chunk` | int | `20000` | 2: |  | Memory: precursors per extraction call (pairs stay together); bounds extraction memory, not the calibration's |

## Deprecated options

Accepted with a WARNING in 0.12 and removed in 0.13; `-write_ini` does not write them.

| Option | Replacement |
|---|---|
| `-tune_models` | -models |
| `-tune_predict_gpu` | -machine:device (re-prediction follows it) |
| `-machine:threads` | -machine:train_threads |
| `-search:passes` | nothing (1 is its only value) |
| `-search:allow_bootstrap` | ODIA_SEARCH_ALLOW_BOOTSTRAP=1 |
