# The classifier and FDR core vendored from OpenDIAlyzer (src/odia-core/).
#
# Plain C++ with main(), no OpenMS: a test returns 0 on success. They are built
# as C++17 although the rest of the tree is C++23, so the headers stay usable
# from OpenDIAlyzer and from any compiler the release matrix carries. With
# OpenMP the fold loop runs in parallel, which is what the thread-count
# determinism checks need to mean anything.

set(ODIA_CORE_TESTS
  odia_lda_test
  odia_lda_adversarial_test
  odia_fdr_test
  odia_pair_folds_test)

foreach(t IN LISTS ODIA_CORE_TESTS)
  add_executable(odia_core_${t} ${t}.cpp)
  target_link_libraries(odia_core_${t} PRIVATE odia_core)
  set_target_properties(odia_core_${t} PROPERTIES
    CXX_STANDARD 17 CXX_STANDARD_REQUIRED ON CXX_EXTENSIONS OFF)
  if(OpenMP_CXX_FOUND)
    target_link_libraries(odia_core_${t} PRIVATE OpenMP::OpenMP_CXX)
  endif()
  add_test(NAME odia_core_${t} COMMAND odia_core_${t})
endforeach()

# The vendored files match their manifest and carry their licence line.
if(ODIA_PYTHON)
  add_test(NAME odia_core_manifest
           COMMAND ${CMAKE_COMMAND} -E env PYTHONDONTWRITEBYTECODE=1 ${ODIA_PYTHON}
                   "${PROJECT_SOURCE_DIR}/scripts/sync-odia-core.py" --check)
endif()
