# Message review (0.12)

The user-facing messages were reviewed for 0.12: each says what went wrong,
the likely cause and the fix, in one or two lines, without development
history or measurements (those are in [usage](usage.md) and the
CHANGELOG). Details for diagnosis moved behind `-debug 1` (or
`ODIA_VERBOSE=1`). Every warning starts with `WARNING:`. Unchanged on
purpose: the provider lines (`[MS2] provider CUDA device 0 sessions 4`), the
`digest:`, `after fragment selection:`, `decoys:` and `wrote` lines, and the
phrases tests and scripts look for (`the GPU is required (-machine:require_gpu)`,
`unknown config key:`, `unknown instrument: '`, `fewer than 95 %`, ...).

The main changes (`...` abbreviates):

| where | before | after | verbose (`-debug 1`) adds |
|---|---|---|---|
| prediction device (new) | -- | `device: prediction on CUDA (cuDNN 9.1.0)`, `device: MS2 on the Apple GPU (CoreML), RT and CCS on the CPU: 3.4 M precursors, CoreML pays from 480 k (...)` | `device check:` providers, GPU visibility, cuDNN dlopen result and places looked in |
| no cuDNN | per model: `WARNING: RT: peptdeep_rt_dynamic.onnx runs on the CPU although CUDA was requested: CUDA device 0: cuDNN is unavailable ... Fix: install cuDNN 9 for CUDA 12, including the unversioned libcudnn.so (apt: ...; or pip ... plus a libcudnn.so symlink, on LD_LIBRARY_PATH); see https://github.com/...` | once: `WARNING: RT and CCS run on the CPU, not the GPU: no libcudnn.so on the library path. Fix: put cuDNN 9 for CUDA 12, with libcudnn.so, on LD_LIBRARY_PATH (one is in <dir>); docs/usage.md#gpu` | the loader error, the directories looked in |
| model fallback | `WARNING: RT: peptdeep_rt_dynamic.onnx runs on the CPU although CoreML was requested: CoreML has no operator for LSTM x2: ONNX Runtime keeps 5 of the model's nodes on the CPU (2 CoreML partitions). Fix: this model runs on the CPU path (...); see https://...` | `WARNING: RT: runs on the CPU, not CoreML: CoreML has no operator for LSTM x2. Fix: none needed: this model runs on the CPU as fast as without CoreML (-machine:device auto does not ask)` | model path, node and partition counts |
| fallback summary | `2 model(s) requested on the GPU ran on the CPU; recorded as odia.warnings (TSV: log only)` | `WARNING: part of the prediction ran on the CPU instead of the GPU (see the WARNINGs above); recorded in the library's metadata (odia.warnings)` | |
| iRT fit | `[IRT-FIT] provider ...` and a second copy of the RT fallback warning | nothing (the RT line says it) | the `[IRT-FIT]` line |
| CoreML shape | `WARNING: MS2: CoreML could not build the 2048x16 shape (<ORT error>); batches of that shape run on the CPU` | `WARNING: MS2: CoreML cannot compile batches of 2048 peptides of 14 residues; those run on the CPU (same output, slower)` | the ONNX Runtime error |
| CoreML summary | `[MS2] CoreML: 42 shapes built in 105.0 s (summed), ... padding rows (8.8 %)`, cache pruning line | -- | both lines |
| ONNX Runtime | `[W:onnxruntime:odia, file.cc:123 Func] ...` | errors only | its warnings |
| unknown instrument | `unknown instrument: 'X'. Known: QE, ... [about 700 characters on aliases and untrained slots]` | `unknown instrument: 'X'. Use QE, Lumos, timsTOF, SciexTOF or ThermoTOF, or a vendor name such as Astral, Exploris 480 or timsTOF Pro (docs/usage.md#generation-settings); Lumos means no instrument correction.` | |
| ThermoTOF, SciexTOF | `... carries no trained weights in the shipped MS2 checkpoint ... (median spectral cosine 0.993 ...)` | `WARNING: ThermoTOF has no trained MS2 weights; its spectra are effectively Lumos's. Use -generation:instrument Lumos to say so.` (SciexTOF alike) | |
| instrument alias, NCE | `instrument 'Astral' is 'Lumos' to this model (upstream's instrument_group)`; `nce not set; using 25.000000 for Lumos (recorded as instrument-default:Lumos)` | `instrument 'Astral' is predicted as Lumos`; `NCE not set; using Lumos's default, 25` | |
| ion mobility (new) | -- | `Astral has no ion mobility: CCS and 1/K0 not predicted (-generation:predict_ion_mobility true to predict them)` | |
| models not found | `no rt_model configured and peptdeep_rt_dynamic.onnx was not found. No tagged OpenMS release ships the PeptDeep models; ... Searched:` + list | `peptdeep_rt_dynamic.onnx not found; set DIALIBGEN_MODEL_DIR to the directory holding the three PeptDeep models, or give -generation:rt_model (-debug 1 lists the directories searched)` | the directories |
| model path | `rt_model does not exist: <path>` | `-generation:rt_model <path> does not exist` | |
| RT misses | `N precursors have no predicted RT. Refusing rather than shipping a library with NaN retention times, which extracts from nowhere and looks like a search bug.` | `N precursors have no predicted RT (a modification the RT model cannot encode?); no library written` | |
| iRT | `RT column is the model's RAW 0..1 output, NOT iRT (irt_rescale=false). The per-run map is scale-invariant, ...`; `rescaled to iRT: 1.0234 * raw + -2.1000 from 11 standards (worst off by 2.1000 iRT)`; `no iRT standards found; pass -irt_standards. Searched:` + list | `RT is the model's raw 0..1 scale, not iRT; -generation:irt_rescale true writes iRT for a library another tool reads as iRT`; `rescaled to iRT with 11 standards (largest error 2.1 iRT)`; `no iRT standards file found; give -irt_standards` | slope and intercept; the directories |
| overlap | `RT and CCS predicted on a thread of their own, 8 sessions, while the GPU predicts MS2 (ODIA_OVERLAP=0: one after the other)` | `RT and CCS predicted on the CPU (8 sessions) while the GPU predicts MS2` | |
| streaming | `no thread for the streaming writer (...); predicting in 0.11.1's order and writing after MS2`; `CCS of N precursors predicted again in their 0.11.1 batching` | `WARNING: cannot start the streaming writer (...); the library is written after prediction instead (same output, more memory)`; -- | the re-prediction count |
| TSV, OpenSWATH | warnings `TSV cannot carry the recipe: the config is NOT embedded and this file is not cache-eligible. ...` | info `the TSV does not embed the generation settings; .parquet keeps them with the library` | |
| -threads 0 | `-threads 0: 16 threads (hw.physicalcpu (12 performance and 4 efficiency cores))`; Linux `(physical cores of the 384 CPUs in the affinity mask, limited to 32 by the cgroup CPU quota)` | `-threads 0: 16 threads (physical cores: 12 performance, 4 efficiency)`; `(physical cores, limited by the CPU quota)` | |
| prepare | `prepared 42 CoreML shape(s) of peptdeep_ms2_dynamic.onnx for 7-30 residues in 105.0 s (summed over the builders); cache <dir>` | `prepared 42 CoreML shape(s) for peptide lengths 7-30; cache: <dir>` | model, seconds |
| library writers | `warning: N modification tokens have no UniMod accession and were written verbatim; ...`; `warning: dropped N transitions whose fragment ordinal or charge cannot describe a fragment of the stated peptide`; OpenSWATH decoys (4 lines) | `WARNING: N modifications have no UniMod accession and were written by name; DIA-NN may skip those precursors`; `WARNING: dropped N transitions with a fragment number or charge impossible for their peptide`; one line | |
| refine | `refine: search: ...`; `MIXED INTENSITY PROVENANCE: ...`; `NOTE: the RT column now holds the REFERENCE RUN's ...`; `no reference precursor matched the library. That is a join failure -- ...` | `search: ...`; `WARNING: mixed intensities: matched precursors carry observed, the others predicted intensities (-allow_mixed_intensity)`; `RT now holds the reference run's observed retention times (minutes), not iRT: use this library for that run's gradient`; `no reference precursor matched the library; check that -ids and -in use the same fixed modifications (e.g. carbamidomethyl C) and modification names` | |
| -run, no 1/K0 | `... is diaPASEF and is searched with its ion mobility, but only N of M ... (about 600 characters)` | `search: the run X is diaPASEF, but only N of M library targets (P) have a 1/K0 (IM or CCS), fewer than 95 %: the others would never be extracted. Generate the library with ion mobility (-generation:instrument timsTOF or -generation:predict_ion_mobility true), or search the run by m/z and RT only (search:im_window -1). Nothing was searched.` | |
| tune | `filter:rt_spread_max must be finite and in [0.000000, 1797693134862315708145274237317043567980...]`; `-out must not be the stock model itself` | `filter:rt_spread_max must be a finite number >= 0`; `-tune_out_models must not be the -tune_models directory (it would overwrite the stock model)` | |
| parameters | `machine:*` descriptions with measurements (`one H100, whole proteome: ~25 s less`, `+30-45 %`, `1.5x slower`, `-3 to +16 % wall`) | the behaviour only; numbers in usage.md | |

Not changed in this round (listed in the review, lower priority): the
built-in identification's (`-run`, experimental) long diagnostic messages
other than the one above, and the stock-model parse errors of the
fine-tuning ONNX reader. The trainer now logs every tenth epoch and the last
(every epoch with `-debug 1`; all of them are in the `.trajectory.tsv`).
