# Generation at different -threads must write the same bytes: digestion,
# enumeration, fragment selection and both writers run on host threads, and
# the ONNX session count follows -threads too.
#   -DTOOL=DIALibGen -DWORKDIR=... -DCONFIG=...
#
# The input is a synthetic FASTA, large enough for several parallel tasks per
# stage and, with ODIA_MS2_BLOCK=10000 (set by the test), for three MS2 blocks, so the pipelined
# fragment selection hands a block over at least once.
file(REMOVE_RECURSE "${WORKDIR}")
file(MAKE_DIRECTORY "${WORKDIR}")
set(fasta "${WORKDIR}/synthetic.fasta")
set(text "")
foreach(i RANGE 1 140)
  string(RANDOM LENGTH 360 ALPHABET "ACDEFGHIKLMNPQRSTVWYKRKR" RANDOM_SEED ${i} sequence)
  string(APPEND text ">sp|SYN${i}|SYN${i}_TEST synthetic ${i}\nM${sequence}\n")
endforeach()
file(WRITE "${fasta}" "${text}")

foreach(ext tsv parquet)
  foreach(threads 2 8)
    set(out "${WORKDIR}/t${threads}.${ext}")
    execute_process(COMMAND "${TOOL}" -in "${fasta}" -config "${CONFIG}" -out "${out}" -threads ${threads}
                    RESULT_VARIABLE rc OUTPUT_VARIABLE log ERROR_VARIABLE log)
    if(NOT rc EQUAL 0)
      message(FATAL_ERROR "-threads ${threads} -out ${ext} failed (${rc}):\n${log}")
    endif()
  endforeach()
  string(REGEX MATCH "digest: [^\n]*" digest "${log}")
  message(STATUS "${digest}")
  execute_process(COMMAND ${CMAKE_COMMAND} -E compare_files "${WORKDIR}/t2.${ext}" "${WORKDIR}/t8.${ext}"
                  RESULT_VARIABLE differ)
  if(differ)
    message(FATAL_ERROR "-threads 8 wrote a different ${ext} than -threads 2")
  endif()
endforeach()
