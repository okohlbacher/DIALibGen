# The library must not depend on how the prediction is cut into calls: small
# MS2 blocks (many pipelined blocks, the helper thread parsing, selecting and
# appending beside every call) and the 0.11.1 order (ODIA_STREAM=0) must write
# the bytes the defaults write, where the whole input is one MS2 block. MS2 is
# bit-identical for every batch composition, so any difference is the
# pipeline's. CCS is not -- a row alone in its batch gets a slightly different
# value (CcsBatchGuard) -- so ODIA_CCS_BLOCK is varied only together with the
# order, against the same CCS block.
#   -DTOOL=DIALibGen -DCONFIG=... -DWORKDIR=...
file(REMOVE_RECURSE "${WORKDIR}")
file(MAKE_DIRECTORY "${WORKDIR}")
set(fasta "${WORKDIR}/synthetic.fasta")
set(text "")
foreach(i RANGE 1 140)
  string(RANDOM LENGTH 360 ALPHABET "ACDEFGHIKLMNPQRSTVWYKRKR" RANDOM_SEED ${i} sequence)
  string(APPEND text ">sp|SYN${i}|SYN${i}_TEST synthetic ${i}\nM${sequence}\n")
endforeach()
file(WRITE "${fasta}" "${text}")

# name | environment (";"-free, space-separated) | reference it must equal
set(variants
  "default||"
  "ms2_997|ODIA_MS2_BLOCK=997|default"
  "one_pass|ODIA_STREAM=0|default"
  "ccs_4099|ODIA_CCS_BLOCK=4099|"
  "ccs_4099_one_pass_ms2_997|ODIA_CCS_BLOCK=4099 ODIA_STREAM=0 ODIA_MS2_BLOCK=997|ccs_4099")
foreach(ext tsv parquet)
  foreach(variant IN LISTS variants)
    string(REPLACE "|" ";" fields "${variant}")
    list(LENGTH fields nfields)
    list(GET fields 0 name)
    set(envs "")
    set(reference "")
    if(nfields GREATER 1)
      list(GET fields 1 envs)
      separate_arguments(envs UNIX_COMMAND "${envs}")
    endif()
    if(nfields GREATER 2)
      list(GET fields 2 reference)
    endif()
    set(out "${WORKDIR}/${name}.${ext}")
    execute_process(COMMAND ${CMAKE_COMMAND} -E env --unset=ODIA_MS2_BLOCK --unset=ODIA_CCS_BLOCK
                            --unset=ODIA_STREAM ${envs}
                            "${TOOL}" -in "${fasta}" -config "${CONFIG}" -out "${out}" -threads 4
                    RESULT_VARIABLE rc OUTPUT_VARIABLE log ERROR_VARIABLE log)
    if(NOT rc EQUAL 0)
      message(FATAL_ERROR "${name} -out ${ext} failed (${rc}):\n${log}")
    endif()
    if(name STREQUAL "default")
      string(REGEX MATCH "digest: [^\n]*" digest "${log}")
      message(STATUS "${digest}")
    endif()
    if(reference)
      execute_process(COMMAND ${CMAKE_COMMAND} -E compare_files "${WORKDIR}/${reference}.${ext}" "${out}"
                      RESULT_VARIABLE differ)
      if(differ)
        message(FATAL_ERROR "${name} (${envs}) wrote a different ${ext} than ${reference}")
      endif()
    endif()
  endforeach()
endforeach()
