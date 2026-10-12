# The built-in identification step (-run): odia_search.
#
# The C++ tests build their own synthetic libraries and, for the pipeline test,
# a synthetic run behind the Identifier seams; none needs data, and only
# identify_predicted_assays_model needs a model (the PeptDeep MS2 model; the
# others use toy_fragment_model.h or search:intensities library).

add_executable(identify_subset_pairing identify_subset_pairing.cpp)
target_link_libraries(identify_subset_pairing PRIVATE odia_search)
add_test(NAME identify_subset_pairing COMMAND identify_subset_pairing)

# IonLadder: every fragment ion type bit-identical to OpenMS's getPrefix/getSuffix.
add_executable(identify_ion_ladder identify_ion_ladder.cpp)
target_link_libraries(identify_ion_ladder PRIVATE odia_search)
add_test(NAME identify_ion_ladder COMMAND identify_ion_ladder)

add_executable(identify_decoy_symmetry identify_decoy_symmetry.cpp)
target_link_libraries(identify_decoy_symmetry PRIVATE odia_search)
add_test(NAME identify_decoy_symmetry COMMAND identify_decoy_symmetry)

# Search decoys exchangeable with null targets: search:intensities predicted
# against library, on a toy fragment model (toy_fragment_model.h), and the
# entrapment tag rule.
add_executable(identify_decoy_exchangeability identify_decoy_exchangeability.cpp)
target_link_libraries(identify_decoy_exchangeability PRIVATE odia_search)
add_test(NAME identify_decoy_exchangeability COMMAND identify_decoy_exchangeability)
set_tests_properties(identify_decoy_exchangeability PROPERTIES TIMEOUT 600 ENVIRONMENT "OMP_NUM_THREADS=4")

# The same with the real PeptDeep MS2 model: the searched targets' assays are
# the generator's, the decoys' their own predictions.
add_executable(identify_predicted_assays_model identify_predicted_assays_model.cpp)
target_link_libraries(identify_predicted_assays_model PRIVATE odia_search)
if(EXISTS "${ODIA_MS2_MODEL_PATH}")
  add_test(NAME identify_predicted_assays_model COMMAND identify_predicted_assays_model "${ODIA_MS2_MODEL_PATH}")
  set_tests_properties(identify_predicted_assays_model PROPERTIES TIMEOUT 600 ENVIRONMENT "OMP_NUM_THREADS=2")
endif()

# search:ms1_in_memory's decision rule and the free-memory query.
add_executable(identify_ms1_memory identify_ms1_memory.cpp)
target_link_libraries(identify_ms1_memory PRIVATE odia_search)
add_test(NAME identify_ms1_memory COMMAND identify_ms1_memory)

# -threads 0 for -run: the physical cores, so more than one prediction session.
add_executable(identify_threads identify_threads.cpp)
target_link_libraries(identify_threads PRIVATE odia_search)
add_test(NAME identify_threads COMMAND identify_threads "${ODIA_MS2_MODEL_PATH}")

# The RT calibration's line-or-LOWESS choice on synthetic calibration points.
add_executable(identify_rt_model identify_rt_model.cpp)
target_link_libraries(identify_rt_model PRIVATE odia_search nlohmann_json::nlohmann_json)
add_test(NAME identify_rt_model COMMAND identify_rt_model)

add_executable(identify_report_roundtrip identify_report_roundtrip.cpp)
target_link_libraries(identify_report_roundtrip PRIVATE odia_search nlohmann_json::nlohmann_json)
if(TARGET odia_tune)
  # The trainer's report loader is part of the hand-off wherever tuning exists.
  target_link_libraries(identify_report_roundtrip PRIVATE odia_tune)
  target_compile_definitions(identify_report_roundtrip PRIVATE HAS_TUNE)
endif()
add_test(NAME identify_report_roundtrip
         COMMAND identify_report_roundtrip "${CMAKE_CURRENT_BINARY_DIR}/report-roundtrip")

add_executable(identify_pipeline_synthetic identify_pipeline_synthetic.cpp)
target_link_libraries(identify_pipeline_synthetic PRIVATE odia_search nlohmann_json::nlohmann_json)
add_test(NAME identify_pipeline_synthetic
         COMMAND identify_pipeline_synthetic "${CMAKE_CURRENT_BINARY_DIR}/pipeline-synthetic")
set_tests_properties(identify_pipeline_synthetic PROPERTIES TIMEOUT 600)

if(TARGET DIALibGen AND ODIA_PYTHON)
  add_test(NAME identify_cli_contract
           COMMAND ${CMAKE_COMMAND} -E env PYTHONDONTWRITEBYTECODE=1 ${ODIA_PYTHON}
                   "${CMAKE_CURRENT_SOURCE_DIR}/identify_cli_contract.py" $<TARGET_FILE:DIALibGen>)
  set_tests_properties(identify_cli_contract PROPERTIES TIMEOUT 300)
endif()

# The run-facing stages through stock OpenMS on a synthetic mzML with planted
# peaks (synthetic_run.h); the fixture is generated at test time, never stored.
add_executable(identify_extraction identify_extraction.cpp)
target_link_libraries(identify_extraction PRIVATE odia_search)
add_test(NAME identify_extraction
         COMMAND identify_extraction "${CMAKE_CURRENT_BINARY_DIR}/extraction")
# OpenMP would otherwise take every core of the machine (and spin on them).
set_tests_properties(identify_extraction PROPERTIES TIMEOUT 900 ENVIRONMENT "OMP_NUM_THREADS=4")

# Ion mobility (M2) on a synthetic diaPASEF mzML: two 1/K0 bands per window,
# a mis-calibrated library 1/K0, pair symmetry, and no change without it.
add_executable(identify_ion_mobility identify_ion_mobility.cpp)
target_link_libraries(identify_ion_mobility PRIVATE odia_search)
add_test(NAME identify_ion_mobility
         COMMAND identify_ion_mobility "${CMAKE_CURRENT_BINARY_DIR}/ion-mobility")
set_tests_properties(identify_ion_mobility PROPERTIES TIMEOUT 900 ENVIRONMENT "OMP_NUM_THREADS=4")

# The evidence prefilter and the chunking on the same kind of synthetic run.
add_executable(identify_prefilter identify_prefilter.cpp)
target_link_libraries(identify_prefilter PRIVATE odia_search)
add_test(NAME identify_prefilter
         COMMAND identify_prefilter "${CMAKE_CURRENT_BINARY_DIR}/prefilter")
set_tests_properties(identify_prefilter PROPERTIES TIMEOUT 900 ENVIRONMENT "OMP_NUM_THREADS=4")

# A measurement tool for the prefilter's settings on a real library and run
# (docs/design/built-in-identification.md, M3); built, not a test.
add_executable(identify_prefilter_probe identify_prefilter_probe.cpp)
target_link_libraries(identify_prefilter_probe PRIVATE odia_search nlohmann_json::nlohmann_json)

# The same fixture as files, for the command-line tests.
add_executable(identify_synth_run synth_run.cpp)
target_link_libraries(identify_synth_run PRIVATE odia_search)

if(TARGET DIALibGen AND ODIA_PYTHON)
  add_test(NAME identify_e2e_refine
           COMMAND ${CMAKE_COMMAND} -E env PYTHONDONTWRITEBYTECODE=1 ${ODIA_PYTHON}
                   "${CMAKE_CURRENT_SOURCE_DIR}/identify_e2e.py" refine $<TARGET_FILE:DIALibGen>
                   $<TARGET_FILE:identify_synth_run>)
  set_tests_properties(identify_e2e_refine PROPERTIES TIMEOUT 1800)
  if(TARGET odia_tune)
    add_test(NAME identify_e2e_tune
             COMMAND ${CMAKE_COMMAND} -E env PYTHONDONTWRITEBYTECODE=1 ${ODIA_PYTHON}
                     "${CMAKE_CURRENT_SOURCE_DIR}/identify_e2e.py" tune $<TARGET_FILE:DIALibGen>
                     $<TARGET_FILE:identify_synth_run> "${ODIA_MODEL_DIR}")
    set_tests_properties(identify_e2e_tune PROPERTIES TIMEOUT 2400 SKIP_RETURN_CODE 77)
  endif()
endif()

# mzPeak input (-run): an mzPeak archive and an mzML of the same data load into
# the same OpenSWATH maps (identify_mzpeak_run), and give the same report.
add_executable(identify_mzpeak_run identify_mzpeak_run.cpp)
target_link_libraries(identify_mzpeak_run PRIVATE odia_search)
# The fixture: the first 123 spectra (3 cycles: MS1 and 40 diaPASEF windows,
# 1/K0 bands and per-peak 1/K0) of identify_synth_run's "im" run (120 peptides,
# 0.3, seed 20261012), converted by mzpeak-convert 0.18.0 (MZPC_MAX_SPECTRA=123;
# chunked delta layout), and that archive's own mzML export (gzipped). No
# instrument data.
add_test(NAME identify_mzpeak_run
         COMMAND identify_mzpeak_run "${PROJECT_SOURCE_DIR}/test/data/mzpeak/synthetic-dia-pasef.mzpeak"
                                     "${PROJECT_SOURCE_DIR}/test/data/mzpeak/synthetic-dia-pasef.mzML.gz" 2 normal)
add_test(NAME identify_mzpeak_run_cache
         COMMAND identify_mzpeak_run "${PROJECT_SOURCE_DIR}/test/data/mzpeak/synthetic-dia-pasef.mzpeak"
                                     "${PROJECT_SOURCE_DIR}/test/data/mzpeak/synthetic-dia-pasef.mzML.gz" 2 cache)
set_tests_properties(identify_mzpeak_run identify_mzpeak_run_cache PROPERTIES TIMEOUT 300 SKIP_RETURN_CODE 77)
