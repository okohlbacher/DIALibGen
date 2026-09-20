add_executable(dialibgen_canonical_modseq canonical_modseq.cpp)
target_link_libraries(dialibgen_canonical_modseq PRIVATE odia_refine)
add_test(NAME refinement_canonical_modseq COMMAND dialibgen_canonical_modseq)

add_executable(dialibgen_intensity_match intensity_match.cpp)
target_link_libraries(dialibgen_intensity_match PRIVATE odia_refine)
add_test(NAME refinement_intensity_match COMMAND dialibgen_intensity_match)

if(TARGET DIALibGen)
  if(ODIA_PYTHON)
    add_test(NAME refinement_report_contract COMMAND "${ODIA_PYTHON}"
      "${CMAKE_CURRENT_SOURCE_DIR}/report_contract.py" $<TARGET_FILE:DIALibGen>)
    set_tests_properties(refinement_report_contract PROPERTIES TIMEOUT 180)
  endif()
  add_test(NAME refinement_config_contract COMMAND ${CMAKE_COMMAND}
    "-DTOOL=$<TARGET_FILE:DIALibGen>" "-DWORK=${CMAKE_CURRENT_BINARY_DIR}/config-contract"
    -P "${CMAKE_CURRENT_SOURCE_DIR}/config_contract.cmake")
  add_test(NAME refinement_config_write_effective
    COMMAND DIALibGen -mode refine -write_config "${CMAKE_CURRENT_BINARY_DIR}/refine-effective.json")
  set_tests_properties(refinement_config_write_effective PROPERTIES PASS_REGULAR_EXPRESSION "effective config")
  file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/refine-typo.json" "{\"q_precursour\": 0.05}\n")
  add_test(NAME refinement_config_rejects_unknown_key
    COMMAND DIALibGen -mode refine -in x.parquet -ids y.parquet -out z.parquet
      -config "${CMAKE_CURRENT_BINARY_DIR}/refine-typo.json")
  set_tests_properties(refinement_config_rejects_unknown_key PROPERTIES WILL_FAIL TRUE)
endif()

if(TARGET odia_tune)
  add_executable(dialibgen_tune_params tune_params.cpp)
  target_link_libraries(dialibgen_tune_params PRIVATE odia_tune)
  add_test(NAME tune_parameter_contract COMMAND dialibgen_tune_params)
  add_executable(dialibgen_tune_parity tune_parity.cpp)
  target_link_libraries(dialibgen_tune_parity PRIVATE odia_tune)
  set(DIALIBGEN_PEPTDEEP_ONNX_DIR "${ODIA_MODEL_DIR}" CACHE PATH "Stock peptdeep RT/CCS models for refinement tests")
  if(EXISTS "${DIALIBGEN_PEPTDEEP_ONNX_DIR}/peptdeep_rt_dynamic.onnx" AND
     EXISTS "${DIALIBGEN_PEPTDEEP_ONNX_DIR}/peptdeep_ccs_dynamic.onnx")
    add_test(NAME tune_parity COMMAND dialibgen_tune_parity
      "${DIALIBGEN_PEPTDEEP_ONNX_DIR}/peptdeep_rt_dynamic.onnx"
      "${DIALIBGEN_PEPTDEEP_ONNX_DIR}/peptdeep_ccs_dynamic.onnx")
    if(TARGET DIALibGen)
      find_program(DIALIBGEN_REFINEMENT_PYTHON NAMES python3 python)
      set(DIALIBGEN_EXAMPLE_FASTA "-" CACHE FILEPATH "FASTA for refinement training test; - uses synthetic proteins")
      if(DIALIBGEN_REFINEMENT_PYTHON)
        add_test(NAME tune_training_contract COMMAND "${DIALIBGEN_REFINEMENT_PYTHON}"
          "${CMAKE_CURRENT_SOURCE_DIR}/training_contract.py" $<TARGET_FILE:DIALibGen> "${DIALIBGEN_PEPTDEEP_ONNX_DIR}")
        set_tests_properties(tune_training_contract PROPERTIES SKIP_RETURN_CODE 77 TIMEOUT 600)
        add_test(NAME tune_e2e COMMAND bash "${CMAKE_CURRENT_SOURCE_DIR}/tune_e2e.sh"
          $<TARGET_FILE:DIALibGen> "${DIALIBGEN_PEPTDEEP_ONNX_DIR}" "${DIALIBGEN_EXAMPLE_FASTA}" "${DIALIBGEN_REFINEMENT_PYTHON}")
        set_tests_properties(tune_e2e PROPERTIES SKIP_RETURN_CODE 77 TIMEOUT 2400)
      endif()
    endif()
  endif()
endif()
