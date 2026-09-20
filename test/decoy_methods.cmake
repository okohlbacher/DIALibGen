# Builds one small library per decoy method and asserts they differ -- from the
# target and from each other. Expects: GEN PYTHON CHECKER FASTA WORKDIR
# RT_MODEL MS2_MODEL CCS_MODEL
set(args "")
foreach(m mutate shuffle pseudo_reverse reverse)
  set(cfg "${WORKDIR}/decoy_${m}.json")
  file(WRITE "${cfg}"
    "{ \"rt_model\":\"${RT_MODEL}\", \"ms2_model\":\"${MS2_MODEL}\","
    " \"ccs_model\":\"${CCS_MODEL}\", \"decoys\":\"${m}\" }")
  file(REMOVE "${WORKDIR}/decoy_${m}.parquet")
  execute_process(COMMAND ${GEN} -in ${FASTA} -config ${cfg}
                          -out ${WORKDIR}/decoy_${m}.parquet
                  RESULT_VARIABLE rc OUTPUT_VARIABLE o ERROR_VARIABLE e)
  if(NOT rc EQUAL 0)
    message(FATAL_ERROR "generation failed for ${m}: ${rc}\n${o}\n${e}")
  endif()
  list(APPEND args "${m}=${WORKDIR}/decoy_${m}.parquet")
endforeach()
execute_process(COMMAND ${PYTHON} ${CHECKER}
                        ${WORKDIR}/decoy_mutate.parquet ${args}
                RESULT_VARIABLE rc OUTPUT_VARIABLE o ERROR_VARIABLE e)
message(STATUS "${o}${e}")
if(NOT rc EQUAL 0)
  message(FATAL_ERROR "decoy method check failed")
endif()
