# Generation with ODIA_PARSE_CHECK=1: every sequence parsePeptide parses is
# parsed by AASequence::fromString too, and a difference fails the run.
#   -DTOOL=DIALibGen -DFASTA=... -DCONFIG=... -DWORKDIR=...
file(REMOVE_RECURSE "${WORKDIR}")
file(MAKE_DIRECTORY "${WORKDIR}")
execute_process(COMMAND ${CMAKE_COMMAND} -E env ODIA_PARSE_CHECK=1 CUDA_VISIBLE_DEVICES=
                        "${TOOL}" -in "${FASTA}" -config "${CONFIG}" -out "${WORKDIR}/lib.tsv" -threads 4
                RESULT_VARIABLE rc OUTPUT_VARIABLE log ERROR_VARIABLE log)
if(NOT rc EQUAL 0)
  message(FATAL_ERROR "generation with ODIA_PARSE_CHECK=1 failed (${rc}):\n${log}")
endif()
string(REGEX MATCHALL "variable modification '[^\n]*" mods "${log}")
message(STATUS "${mods}")
