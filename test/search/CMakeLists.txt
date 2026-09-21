# The built-in identification step (-run): odia_search.
#
# The C++ tests build their own synthetic libraries and, for the pipeline test,
# a synthetic run behind the Identifier seams; none needs data or models.

add_executable(identify_subset_pairing identify_subset_pairing.cpp)
target_link_libraries(identify_subset_pairing PRIVATE odia_search)
add_test(NAME identify_subset_pairing COMMAND identify_subset_pairing)
