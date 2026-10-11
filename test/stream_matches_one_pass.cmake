# generate predicts RT and CCS first and writes the library while MS2 runs
# (Ms2Stream, DIANNLibraryFile::StreamWriter). ODIA_STREAM=0 is the 0.11.1
# order (MS2, then CCS) with the one-pass writer. Both must write the same
# bytes, TSV and Parquet, at any -threads.
#   -DTOOL=DIALibGen -DCONFIG=... -DWORKDIR=...
#
# Small MS2 blocks, row groups and CCS calls make the stream cross block and
# row-group boundaries, and a strict intensity floor drops precursors after
# CCS has run on them, so CcsBatchGuard has rows to predict again: with
# 37-precursor CCS calls many rows are alone in their batch.
file(REMOVE_RECURSE "${WORKDIR}")
file(MAKE_DIRECTORY "${WORKDIR}")
set(FASTA "${WORKDIR}/synthetic.fasta")
set(text "")
foreach(i RANGE 1 70)
  string(RANDOM LENGTH 300 ALPHABET "ACDEFGHIKLMNPQRSTVWYKRKR" RANDOM_SEED ${i} sequence)
  string(APPEND text ">sp|SYN${i}|SYN${i}_TEST synthetic ${i}\nM${sequence}\n")
endforeach()
file(WRITE "${FASTA}" "${text}")
set(knobs ODIA_MS2_BLOCK=700 ODIA_ROW_GROUP_ROWS=500 ODIA_CCS_BLOCK=37)
set(resettled FALSE)
foreach(ext tsv parquet)
  foreach(threads 1 4)
    foreach(stream 0 1)
      set(out "${WORKDIR}/s${stream}_t${threads}.${ext}")
      execute_process(COMMAND ${CMAKE_COMMAND} -E env ${knobs} ODIA_STREAM=${stream}
                              "${TOOL}" -in "${FASTA}" -config "${CONFIG}" -out "${out}" -threads ${threads}
                      RESULT_VARIABLE rc OUTPUT_VARIABLE log ERROR_VARIABLE log)
      if(NOT rc EQUAL 0)
        message(FATAL_ERROR "ODIA_STREAM=${stream} -threads ${threads} -out ${ext} failed (${rc}):\n${log}")
      endif()
      if(log MATCHES "CCS of [1-9][0-9]* precursors predicted again")
        set(resettled TRUE)
      endif()
    endforeach()
    execute_process(COMMAND ${CMAKE_COMMAND} -E compare_files "${WORKDIR}/s0_t${threads}.${ext}"
                            "${WORKDIR}/s1_t${threads}.${ext}" RESULT_VARIABLE differ)
    if(differ)
      message(FATAL_ERROR "the streamed ${ext} at -threads ${threads} differs from the one-pass one")
    endif()
  endforeach()
  execute_process(COMMAND ${CMAKE_COMMAND} -E compare_files "${WORKDIR}/s1_t1.${ext}" "${WORKDIR}/s1_t4.${ext}"
                  RESULT_VARIABLE differ)
  if(differ)
    message(FATAL_ERROR "the streamed ${ext} differs between -threads 1 and 4")
  endif()
endforeach()
if(NOT resettled)
  message(FATAL_ERROR "no CCS was predicted again: the test did not exercise CcsBatchGuard")
endif()
