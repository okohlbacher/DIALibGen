# The built-in identification step (-run): odia_search.
#
# The C++ tests build their own synthetic libraries and, for the pipeline test,
# a synthetic run behind the Identifier seams; none needs data or models.

add_executable(identify_subset_pairing identify_subset_pairing.cpp)
target_link_libraries(identify_subset_pairing PRIVATE odia_search)
add_test(NAME identify_subset_pairing COMMAND identify_subset_pairing)

add_executable(identify_report_roundtrip identify_report_roundtrip.cpp)
target_link_libraries(identify_report_roundtrip PRIVATE odia_search nlohmann_json::nlohmann_json)
if(TARGET odia_tune)
  # The trainer's report loader is part of the hand-off wherever tuning exists.
  target_link_libraries(identify_report_roundtrip PRIVATE odia_tune)
  target_compile_definitions(identify_report_roundtrip PRIVATE HAS_TUNE)
endif()
add_test(NAME identify_report_roundtrip
         COMMAND identify_report_roundtrip "${CMAKE_CURRENT_BINARY_DIR}/report-roundtrip")

add_executable(identify_pipeline_synthetic identify_pipeline_synthetic.cpp)
target_link_libraries(identify_pipeline_synthetic PRIVATE odia_search nlohmann_json::nlohmann_json)
add_test(NAME identify_pipeline_synthetic
         COMMAND identify_pipeline_synthetic "${CMAKE_CURRENT_BINARY_DIR}/pipeline-synthetic")
set_tests_properties(identify_pipeline_synthetic PROPERTIES TIMEOUT 600)
