# Tests for the library generator.
#
# Three tiers, by what they need:
#   1. nothing          -- always run
#   2. python3          -- the independent reference encoder
#   3. the ONNX models  -- prediction and calibration
#
# Tier 3 is gated rather than required because no tagged OpenMS release ships
# the PeptDeep models. A skipped tier-3 test is a missing model, not a failure.

find_program(ODIA_PYTHON python3)

set(ODIA_MODEL_DIR "$ENV{ODIA_OPENMS}/share/OpenMS/models"
    CACHE PATH "Directory holding peptdeep_{rt,ms2,ccs}_dynamic.onnx")
set(ODIA_RT_MODEL_PATH  "${ODIA_MODEL_DIR}/peptdeep_rt_dynamic.onnx"  CACHE FILEPATH "PeptDeep RT model")
set(ODIA_MS2_MODEL_PATH "${ODIA_MODEL_DIR}/peptdeep_ms2_dynamic.onnx" CACHE FILEPATH "PeptDeep MS2 model")
set(ODIA_CCS_MODEL_PATH "${ODIA_MODEL_DIR}/peptdeep_ccs_dynamic.onnx" CACHE FILEPATH "PeptDeep CCS model")

set(ODIA_FIXTURE_DIR "" CACHE PATH "Optional external fixtures (diann_library_small.tsv)")
set(ODIA_IRT_FIXTURE "$ENV{ODIA_OPENMS_SOURCE}/src/tests/class_tests/openms/data/peptdeep_irt_peptides_predicted.csv"
    CACHE FILEPATH "OpenMS's predicted iRT standards fixture")

foreach(helper odia_encode_dump odia_group_dump odia_mod_symbols odia_library_subset
               odia_predict_rt odia_predict_ms2 odia_predict_ccs
               odia_irt_calibration odia_ms2_sessions)
  add_executable(${helper} tools/${helper}.cpp)
  target_link_libraries(${helper} PRIVATE odia_library)
endforeach()

# --------------------------------------------------------------- tier 1

# The tool must answer --help and round-trip its own INI. This catches
# parameter-registration mistakes without needing any data.
if(TARGET DIALibraryGenerator)
  add_test(NAME DIALibraryGenerator_help COMMAND DIALibraryGenerator --help)
  add_test(NAME DIALibraryGenerator_write_ini
           COMMAND DIALibraryGenerator -write_ini "${CMAKE_CURRENT_BINARY_DIR}/DIALibraryGenerator.ini")

  # Everything that is right for a tool shipped INSIDE OpenMS and wrong for one
  # that merely links it: own version, update check off, bare-environment start,
  # input validation. See test/standalone_test.sh.
  add_test(NAME standalone_tool_contract
           COMMAND ${CMAKE_COMMAND} -E env bash
                   "${CMAKE_CURRENT_SOURCE_DIR}/standalone_test.sh"
                   $<TARGET_FILE:DIALibraryGenerator> "${PROJECT_VERSION}"
                   "${CMAKE_CURRENT_SOURCE_DIR}/../example/proteins.fasta")

  # -write_ctd and friends, which stock OpenMS refuses to a tool outside its own
  # tree. Skips (77) without python3.
  add_test(NAME write_ctd
           COMMAND ${CMAKE_COMMAND} -E env bash
                   "${CMAKE_CURRENT_SOURCE_DIR}/write_ctd_test.sh"
                   $<TARGET_FILE:DIALibraryGenerator>)
  set_tests_properties(write_ctd PROPERTIES SKIP_RETURN_CODE 77)

  # -write_config must emit the EFFECTIVE config: every default materialised, so
  # the embedded recipe is complete and a reader never has to know our defaults.
  # It takes no -in and no -out: seeing the defaults is useful on its own, and
  # requiring an input file nobody reads to get them was a needless refusal.
  add_test(NAME config_write_effective
           COMMAND DIALibraryGenerator
                   -write_config "${CMAKE_CURRENT_BINARY_DIR}/effective.json")
  set_tests_properties(config_write_effective PROPERTIES
    PASS_REGULAR_EXPRESSION "effective config")

  # An unknown config key is a typo that would otherwise silently become a
  # default, and a library built from a typo looks exactly like one built right.
  file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/typo.json" "{\"missed_cleavage\": 2}\n")
  add_test(NAME config_rejects_unknown_key
           COMMAND DIALibraryGenerator -in "${CMAKE_CURRENT_SOURCE_DIR}/../example/proteins.fasta"
                   -config "${CMAKE_CURRENT_BINARY_DIR}/typo.json"
                   -out "${CMAKE_CURRENT_BINARY_DIR}/typo.tsv")
  set_tests_properties(config_rejects_unknown_key PROPERTIES WILL_FAIL TRUE)
endif()

# groupByLength and the homogeneity guard. Collapsing every peptide into one
# group, or removing the guard, previously passed the whole suite -- and mixing
# lengths changes a prediction by more than most modifications do.
add_test(NAME peptdeep_length_grouping
         COMMAND odia_group_dump PEPTIDEK PEPTIDER PEPTIDEKK PEPTIDEKKK PEPTIDEM)
set_tests_properties(peptdeep_length_grouping PROPERTIES
  PASS_REGULAR_EXPRESSION "groups: \\[0,1,4\\] \\[2\\] \\[3\\]\nmixed_batch_refused: yes")

# A modification with no elemental composition encodes to all zeros, which is
# indistinguishable from no modification. OpenMS accepts bare mass shifts.
add_test(NAME peptdeep_rejects_massshift_modification
         COMMAND odia_encode_dump "PEPTIDEC[999]")
set_tests_properties(peptdeep_rejects_massshift_modification PROPERTIES WILL_FAIL TRUE)

# Every element symbol in every modification OpenMS knows must either map into
# AlphaPeptDeep's 109-element list or be a deliberate '?' bucket entry. The four
# symbols pinned here are halogen isotopes from PSI-MOD, which AlphaPeptDeep does
# not load at all; they are recorded so a new unmapped symbol fails instead of
# silently landing in '?'.
add_test(NAME peptdeep_modification_symbols COMMAND odia_mod_symbols)
set_tests_properties(peptdeep_modification_symbols PROPERTIES
  PASS_REGULAR_EXPRESSION "unmapped symbols: \\(35\\)Cl \\(37\\)Cl \\(79\\)Br \\(81\\)Br")

# A subset must own its strings. They were pointers into the parent arena's
# blocks, so a subset resolved into the PARENT's storage -- correct until the
# parent died. The test destroys the parent before reading a single string.
add_test(NAME library_subset_outlives_parent COMMAND odia_library_subset)

# --------------------------------------------------------------- tier 2

if(ODIA_PYTHON)
  # The PeptDeep reference encoder is the oracle the C++ encoder is checked
  # against, so it carries its own guard.
  add_test(NAME peptdeep_reference_encoding
           COMMAND ${CMAKE_COMMAND} -E env PYTHONDONTWRITEBYTECODE=1 ${ODIA_PYTHON}
                   ${CMAKE_CURRENT_SOURCE_DIR}/check_peptdeep_encoding.py)

  # The C++ header and the Python data file are the same table for two
  # languages, and nothing compiled or tested the header before.
  add_test(NAME peptdeep_element_table
           COMMAND ${CMAKE_COMMAND} -E env PYTHONDONTWRITEBYTECODE=1 ${ODIA_PYTHON}
                   ${CMAKE_CURRENT_SOURCE_DIR}/check_element_table.py)

  # The MS2 meta-input constants, pinned. A comparison between this code and the
  # oracle cannot see a constant that is wrong in both, and one of these -- the
  # instrument order -- is documented incorrectly in OpenMS's own headers.
  add_test(NAME peptdeep_meta_inputs
           COMMAND ${CMAKE_COMMAND} -E env PYTHONDONTWRITEBYTECODE=1 ${ODIA_PYTHON}
                   ${CMAKE_CURRENT_SOURCE_DIR}/check_meta_inputs.py)

  # The C++ encoder against the independent Python reference. Written from the
  # same spec but from different inputs -- OpenMS AASequence versus the raw
  # string -- so agreement is evidence the spec was implemented, not that one
  # copied the other.
  add_test(NAME peptdeep_encoder_matches_reference
           COMMAND ${CMAKE_COMMAND} -E env PYTHONDONTWRITEBYTECODE=1 ${ODIA_PYTHON}
                   ${CMAKE_CURRENT_SOURCE_DIR}/compare_encoders.py
                   $<TARGET_FILE:odia_encode_dump>)
endif()

# --------------------------------------------------------------- tier 3

# Predictions, not just tensors. Encoded values can match while dtype, layout or
# ordering assumptions differ, and the sequence list mixes lengths so a broken
# grouping or a lost reordering shows up as wrong values rather than an error.
if(EXISTS "${ODIA_RT_MODEL_PATH}" AND ODIA_PYTHON)
  add_test(NAME peptdeep_predictions_match_reference
           COMMAND ${CMAKE_COMMAND} -E env PYTHONDONTWRITEBYTECODE=1 ${ODIA_PYTHON}
                   ${CMAKE_CURRENT_SOURCE_DIR}/compare_predictions.py
                   $<TARGET_FILE:odia_predict_rt> "${ODIA_RT_MODEL_PATH}")
endif()

# MS2 intensities, covering the inputs the RT model does not have: charge, NCE
# and the rank-1 instrument index. Each is varied across the cases, because a
# defect that ignores one of them passes any fixture that holds it constant.
if(EXISTS "${ODIA_MS2_MODEL_PATH}" AND ODIA_PYTHON)
  add_test(NAME peptdeep_ms2_matches_reference
           COMMAND ${CMAKE_COMMAND} -E env PYTHONDONTWRITEBYTECODE=1 ${ODIA_PYTHON}
                   ${CMAKE_CURRENT_SOURCE_DIR}/compare_ms2.py
                   $<TARGET_FILE:odia_predict_ms2> "${ODIA_MS2_MODEL_PATH}")
endif()

# CCS, the third model kind. Its three inputs match neither the RT model's two
# nor the MS2 model's five, so the model guard rejected it outright until it
# grew a third arm.
if(EXISTS "${ODIA_CCS_MODEL_PATH}" AND ODIA_PYTHON)
  add_test(NAME peptdeep_ccs_matches_reference
           COMMAND ${CMAKE_COMMAND} -E env PYTHONDONTWRITEBYTECODE=1 ${ODIA_PYTHON}
                   ${CMAKE_CURRENT_SOURCE_DIR}/compare_ccs.py
                   $<TARGET_FILE:odia_predict_ccs> "${ODIA_CCS_MODEL_PATH}")
endif()

# The iRT calibration, against OpenMS's own fixture rather than against a
# constant this project wrote down. The fixture carries both the raw ONNX
# prediction and the calibrated iRT for the 11 Biognosys standards, so it tests
# the inference and the line separately.
if(EXISTS "${ODIA_IRT_FIXTURE}" AND EXISTS "${ODIA_RT_MODEL_PATH}" AND ODIA_PYTHON)
  add_test(NAME irt_calibration_matches_openms
           COMMAND ${CMAKE_COMMAND} -E env PYTHONDONTWRITEBYTECODE=1 ${ODIA_PYTHON}
                   ${CMAKE_CURRENT_SOURCE_DIR}/check_irt.py
                   $<TARGET_FILE:odia_irt_calibration> "${ODIA_RT_MODEL_PATH}"
                   "${CMAKE_CURRENT_SOURCE_DIR}/../data/irt_standards.tsv"
                   "${ODIA_IRT_FIXTURE}")
endif()

# A degenerate fit has infinite slope and would be applied silently to every
# precursor in the library.
if(EXISTS "${ODIA_RT_MODEL_PATH}")
  file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/degenerate_irt.tsv"
       "# every peptide identical, so every prediction is identical\nPEPTIDEK\t0.0\nPEPTIDEK\t50.0\nPEPTIDEK\t100.0\n")
  add_test(NAME irt_calibration_rejects_degenerate
           COMMAND odia_irt_calibration "${ODIA_RT_MODEL_PATH}"
                   "${CMAKE_CURRENT_BINARY_DIR}/degenerate_irt.tsv")
  set_tests_properties(irt_calibration_rejects_degenerate PROPERTIES WILL_FAIL TRUE)

  # Two standards is not a calibration, it is a line through two points with no
  # residual to check. Three is the documented floor.
  file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/too_few_irt.tsv"
       "LGGNEQVTR\t-24.92\nLFLQFGAQGSPFLK\t100.0\n")
  add_test(NAME irt_calibration_rejects_too_few
           COMMAND odia_irt_calibration "${ODIA_RT_MODEL_PATH}"
                   "${CMAKE_CURRENT_BINARY_DIR}/too_few_irt.tsv")
  set_tests_properties(irt_calibration_rejects_too_few PROPERTIES WILL_FAIL TRUE)

  # Standards that parse but cannot be encoded. One bad standard must cost
  # precision, not the calibration; three bad ones out of four must fail loudly
  # rather than fit a line through the survivor.
  file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/unencodable_irt.tsv"
       "LGGNEQVTR\t-24.92\nPEPTIDEC[999]K\t20.0\nSAMPLEC[999]R\t60.0\nAAAAC[999]K\t90.0\n")
  add_test(NAME irt_calibration_rejects_unencodable
           COMMAND odia_irt_calibration "${ODIA_RT_MODEL_PATH}"
                   "${CMAKE_CURRENT_BINARY_DIR}/unencodable_irt.tsv")
  set_tests_properties(irt_calibration_rejects_unencodable PROPERTIES WILL_FAIL TRUE)
endif()

# Spreading MS2 batches across sessions must be bit-for-bit identical to running
# them on one. A 4.7x speedup that quietly perturbs intensities would change
# which fragments rank into the top twelve.
if(EXISTS "${ODIA_MS2_MODEL_PATH}" AND ODIA_FIXTURE_DIR
   AND EXISTS "${ODIA_FIXTURE_DIR}/diann_library_small.tsv")
  add_test(NAME ms2_sessions_are_deterministic
           COMMAND odia_ms2_sessions "${ODIA_MS2_MODEL_PATH}"
                   "${ODIA_FIXTURE_DIR}/diann_library_small.tsv" 400)
endif()

# The whole pipeline, run the way a user runs it, with the result read back.
# Everything else checks a piece; this catches a library that is written
# successfully and is wrong. Skips (77) without the models or without pyarrow.
if(TARGET DIALibraryGenerator AND ODIA_PYTHON)
  add_test(NAME e2e
           COMMAND ${ODIA_PYTHON} "${CMAKE_CURRENT_SOURCE_DIR}/e2e_test.py"
                   $<TARGET_FILE:DIALibraryGenerator>
                   "${CMAKE_CURRENT_SOURCE_DIR}/../example/proteins.fasta"
                   "${ODIA_MODEL_DIR}"
                   "${CMAKE_CURRENT_BINARY_DIR}/e2e")
  set_tests_properties(e2e PROPERTIES SKIP_RETURN_CODE 77)
  file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/e2e")
endif()

# The free-cysteine RT correction is a table of measured constants, and a table
# of constants is exactly the kind of thing that drifts unnoticed. This also
# guards the self-gating property that lets it default to on: an alkylated
# library must be bit-identical whether the flag is set or not.
if(TARGET DIALibraryGenerator AND ODIA_PYTHON AND EXISTS "${ODIA_RT_MODEL_PATH}")
  add_test(NAME free_cysteine_correction
           COMMAND ${CMAKE_COMMAND}
                   -DTOOL=$<TARGET_FILE:DIALibraryGenerator>
                   -DPYTHON=${ODIA_PYTHON}
                   -DCHECKER=${CMAKE_CURRENT_SOURCE_DIR}/check_free_cysteine.py
                   -DFASTA=${CMAKE_CURRENT_SOURCE_DIR}/../example/proteins.fasta
                   -DWORKDIR=${CMAKE_CURRENT_BINARY_DIR}
                   -DRT_MODEL=${ODIA_RT_MODEL_PATH}
                   -DMS2_MODEL=${ODIA_MS2_MODEL_PATH}
                   -DCCS_MODEL=${ODIA_CCS_MODEL_PATH}
                   -P ${CMAKE_CURRENT_SOURCE_DIR}/free_cysteine.cmake)

  # The decoy methods must actually be different methods. They were not: the
  # generator hand-rolled its own name lookup and mapped every unknown name to
  # Mutate, so three of the four produced identical libraries in silence.
  add_test(NAME decoy_methods_distinct
           COMMAND ${CMAKE_COMMAND}
                   -DGEN=$<TARGET_FILE:DIALibraryGenerator>
                   -DPYTHON=${ODIA_PYTHON}
                   -DCHECKER=${CMAKE_CURRENT_SOURCE_DIR}/check_decoy_methods.py
                   -DFASTA=${CMAKE_CURRENT_SOURCE_DIR}/../example/proteins.fasta
                   -DWORKDIR=${CMAKE_CURRENT_BINARY_DIR}
                   -DRT_MODEL=${ODIA_RT_MODEL_PATH}
                   -DMS2_MODEL=${ODIA_MS2_MODEL_PATH}
                   -DCCS_MODEL=${ODIA_CCS_MODEL_PATH}
                   -P ${CMAKE_CURRENT_SOURCE_DIR}/decoy_methods.cmake)
endif()
