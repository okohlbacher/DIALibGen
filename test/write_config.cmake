file(REMOVE "${OUT}")
execute_process(COMMAND "${TOOL}" -mode "${MODE}" -write_config "${OUT}"
  RESULT_VARIABLE rc OUTPUT_VARIABLE output ERROR_VARIABLE error)
if(NOT rc EQUAL 0 OR NOT EXISTS "${OUT}")
  message(FATAL_ERROR "write_config failed: ${output}\n${error}")
endif()
message(STATUS "${output}")
