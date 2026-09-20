add_executable(dialibgen_text_writer text_writer.cpp)
target_include_directories(dialibgen_text_writer PRIVATE "${PROJECT_SOURCE_DIR}/include")
add_test(NAME text_writer_contract COMMAND dialibgen_text_writer "${CMAKE_CURRENT_BINARY_DIR}/text-writer")

add_executable(dialibgen_atomic_file atomic_file.cpp)
target_include_directories(dialibgen_atomic_file PRIVATE "${PROJECT_SOURCE_DIR}/include")
add_test(NAME atomic_file_contract COMMAND dialibgen_atomic_file "${CMAKE_CURRENT_BINARY_DIR}/atomic-file")

add_executable(dialibgen_library_io library_io.cpp)
target_link_libraries(dialibgen_library_io PRIVATE odia_library)
if(ODIA_PYTHON)
  add_test(NAME library_io_contract COMMAND ${ODIA_PYTHON} "${CMAKE_CURRENT_SOURCE_DIR}/library_io.py"
    $<TARGET_FILE:dialibgen_library_io> "${CMAKE_CURRENT_BINARY_DIR}/library-io")
  set_tests_properties(library_io_contract PROPERTIES SKIP_RETURN_CODE 77)
endif()
