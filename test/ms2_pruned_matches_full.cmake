# The pruned MS2 graph must give the library the full graph gives.
#   -DTOOL=DIALibGen -DFULL=full.json -DPRUNED=pruned.json -DFASTA=... -DWORKDIR=...
#
# The pruned graph is the full one without its modloss branch: channels 0-3
# bit-identical, 4-7 zero, and generation reads only 0-3. The TSV holds no
# provenance, so it must be the same bytes. (The Parquet records which model
# ran -- its hash in the fingerprint, its path in the config -- so it differs in
# those metadata fields and nowhere else.)
file(REMOVE_RECURSE "${WORKDIR}")
file(MAKE_DIRECTORY "${WORKDIR}")
foreach(variant full pruned)
  if(variant STREQUAL "full")
    set(config "${FULL}")
  else()
    set(config "${PRUNED}")
  endif()
  execute_process(COMMAND "${TOOL}" -in "${FASTA}" -config "${config}" -out "${WORKDIR}/${variant}.tsv" -threads 4
                  RESULT_VARIABLE rc OUTPUT_VARIABLE log ERROR_VARIABLE log)
  if(NOT rc EQUAL 0)
    message(FATAL_ERROR "${variant} MS2 model failed (${rc}):\n${log}")
  endif()
  string(REGEX MATCH "MS2 model: [^\n]*" model "${log}")
  message(STATUS "${model}")
endforeach()
execute_process(COMMAND ${CMAKE_COMMAND} -E compare_files "${WORKDIR}/full.tsv" "${WORKDIR}/pruned.tsv"
                RESULT_VARIABLE differ)
if(differ)
  message(FATAL_ERROR "the pruned MS2 model wrote a different library than the full one")
endif()
