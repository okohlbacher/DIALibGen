add_executable(dialibgen_text_writer text_writer.cpp)
target_include_directories(dialibgen_text_writer PRIVATE "${PROJECT_SOURCE_DIR}/include")
add_test(NAME text_writer_contract COMMAND dialibgen_text_writer "${CMAKE_CURRENT_BINARY_DIR}/text-writer")

add_executable(dialibgen_atomic_file atomic_file.cpp)
target_include_directories(dialibgen_atomic_file PRIVATE "${PROJECT_SOURCE_DIR}/include")
add_test(NAME atomic_file_contract COMMAND dialibgen_atomic_file "${CMAKE_CURRENT_BINARY_DIR}/atomic-file")

add_executable(dialibgen_library_io library_io.cpp)
target_link_libraries(dialibgen_library_io PRIVATE odia_library)
foreach(mode compact-chunks compact-empty)
  add_test(NAME library_io_${mode} COMMAND dialibgen_library_io ${mode} "${CMAKE_CURRENT_BINARY_DIR}/${mode}.parquet")
endforeach()
# Multithreaded writers produce the single-threaded bytes; the parallel loop
# under them runs each index once and reports the lowest failing one.
foreach(mode threads-identical parallel-for-each)
  add_test(NAME library_io_${mode} COMMAND dialibgen_library_io ${mode} "${CMAKE_CURRENT_BINARY_DIR}/${mode}")
endforeach()
if(ODIA_PYTHON)
  add_test(NAME library_io_contract COMMAND ${ODIA_PYTHON} "${CMAKE_CURRENT_SOURCE_DIR}/library_io.py"
    $<TARGET_FILE:dialibgen_library_io> "${CMAKE_CURRENT_BINARY_DIR}/library-io")
  set_tests_properties(library_io_contract PROPERTIES SKIP_RETURN_CODE 77)
endif()
if(ODIA_PYTHON AND TARGET DIALibGen)
  add_test(NAME library_convert_contract COMMAND ${ODIA_PYTHON} "${CMAKE_CURRENT_SOURCE_DIR}/convert_contract.py"
    $<TARGET_FILE:dialibgen_library_io> $<TARGET_FILE:DIALibGen> "${CMAKE_CURRENT_BINARY_DIR}/convert")
endif()

# The opt-in layouts read back: DIA-NN 2.0's Parquet schema column by column.
add_executable(dialibgen_library_export library_export.cpp)
target_link_libraries(dialibgen_library_export PRIVATE odia_library)
foreach(mode diann2)
  add_test(NAME library_export_${mode} COMMAND dialibgen_library_export ${mode} "${CMAKE_CURRENT_BINARY_DIR}/export-${mode}")
endforeach()
