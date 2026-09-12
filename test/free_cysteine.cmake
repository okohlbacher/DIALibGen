# Builds four small libraries and checks the free-cysteine RT correction.
# Expects: TOOL PYTHON CHECKER FASTA WORKDIR RT_MODEL MS2_MODEL CCS_MODEL

function(build name fixmod flag)
  set(cfg "${WORKDIR}/fcys_${name}.json")
  file(WRITE "${cfg}"
    "{ \"rt_model\":\"${RT_MODEL}\", \"ms2_model\":\"${MS2_MODEL}\","
    " \"ccs_model\":\"${CCS_MODEL}\", \"decoys\":\"none\","
    " \"fixed_modifications\":[${fixmod}],"
    " \"free_cysteine_rt_correction\":${flag} }")
  execute_process(COMMAND ${TOOL} -in ${FASTA} -config ${cfg}
                          -out ${WORKDIR}/fcys_${name}.parquet
                  RESULT_VARIABLE rc OUTPUT_VARIABLE out ERROR_VARIABLE err)
  if(NOT rc EQUAL 0)
    message(FATAL_ERROR "DIALibGen (${name}) failed: ${rc}\n${out}\n${err}")
  endif()
endfunction()

build(on  ""                          true)
build(off ""                          false)
build(camon  "\"Carbamidomethyl (C)\"" true)
build(camoff "\"Carbamidomethyl (C)\"" false)

execute_process(COMMAND ${PYTHON} ${CHECKER}
                        ${WORKDIR}/fcys_on.parquet ${WORKDIR}/fcys_off.parquet
                        ${WORKDIR}/fcys_camon.parquet ${WORKDIR}/fcys_camoff.parquet
                RESULT_VARIABLE rc OUTPUT_VARIABLE out ERROR_VARIABLE err)
message(STATUS "${out}${err}")
if(NOT rc EQUAL 0)
  message(FATAL_ERROR "free-cysteine correction check failed")
endif()
