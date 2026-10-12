# The vendored mzpeak-openms reader (third_party/mzpeak-openms), as the static
# library odia_mzpeak. Its own CMakeLists.txt is not used: this one builds the
# same sources against the Arrow/Parquet targets DIALibGen already resolved,
# takes Boost.JSON header-only (no compiled Boost library to ship), and keeps
# the library's bundled MSNumpress apart from OpenMS's copy of it, which has
# the same namespace (ms::numpress::MSNumpress).
set(_mzp ${CMAKE_CURRENT_SOURCE_DIR}/third_party/mzpeak-openms)

if(NOT TARGET libzip::zip)
  find_package(libzip 1.11 CONFIG QUIET)
endif()
if(NOT TARGET libzip::zip)
  find_library(LIBZIP_LIB NAMES zip libzip)
  find_path(LIBZIP_INCLUDE_DIR zip.h)
  if(NOT LIBZIP_LIB OR NOT LIBZIP_INCLUDE_DIR)
    message(FATAL_ERROR "mzPeak input needs libzip (zip.h and the library); install libzip, set LIBZIP_LIB and "
                        "LIBZIP_INCLUDE_DIR, or configure with -DDIALIBGEN_WITH_MZPEAK=OFF")
  endif()
  add_library(libzip::zip UNKNOWN IMPORTED)
  set_target_properties(libzip::zip PROPERTIES IMPORTED_LOCATION "${LIBZIP_LIB}" INTERFACE_INCLUDE_DIRECTORIES "${LIBZIP_INCLUDE_DIR}")
endif()

if(NOT BOOST_INCLUDE_DIR OR NOT EXISTS "${BOOST_INCLUDE_DIR}/boost/json/src.hpp"
   OR NOT EXISTS "${BOOST_INCLUDE_DIR}/boost/compat/function_ref.hpp")
  message(FATAL_ERROR "mzPeak input needs Boost headers with Boost.JSON and Boost.Compat (Boost >= 1.86) in "
                      "'${BOOST_INCLUDE_DIR}'; set BOOST_INCLUDE_DIR, or configure with -DDIALIBGEN_WITH_MZPEAK=OFF")
endif()
file(STRINGS "${BOOST_INCLUDE_DIR}/boost/version.hpp" _mzp_boost_line REGEX "define BOOST_LIB_VERSION")
message(STATUS "mzPeak input: on (libzip, Boost ${_mzp_boost_line})")

set(_mzp_sources
  src/data/array_index.cpp src/data/signals.cpp src/index.cpp src/io/directory.cpp src/io/file.cpp src/io/zip.cpp
  src/metadata/scans.cpp src/metadata/spectrum.cpp src/metadata/table.cpp src/open.cpp src/run_metadata.cpp
  src/schema/buffer_format.cpp src/schema/cv.cpp src/schema/data_kind.cpp src/schema/entity_type.cpp src/schema/file.cpp
  src/schema/group.cpp src/schema/psi/array_type.cpp src/schema/psi/chunk_encoding.cpp src/schema/psi/data_type.cpp
  src/chromatogram.cpp src/chromatograms.cpp src/schema/psi/transform.cpp src/spectra.cpp src/spectrum.cpp
  src/wavelength_spectrum.cpp src/wavelength_spectra.cpp src/util/arrow.cpp src/util/decode_pool.cpp src/util/decoders.cpp
  src/util/executor.cpp src/util/manager.cpp src/util/parquet.cpp src/util/sha512.cpp src/util/json_writer.cpp
  src/util/key_runs.cpp src/util/metadata_model.cpp src/util/null_fill.cpp src/util/numpress.cpp src/util/parquet_writer.cpp
  src/util/planner.cpp src/util/row_group_cache.cpp src/util/projection.cpp src/util/query.cpp src/util/slice.cpp
  src/util/types.cpp src/writer.cpp
  subprojects/msnumpress/src/main/cpp/MSNumpress.cpp)
list(TRANSFORM _mzp_sources PREPEND "${_mzp}/")

add_library(odia_mzpeak STATIC ${_mzp_sources} ${CMAKE_CURRENT_SOURCE_DIR}/src/search/MzPeakBoostJson.cpp)
target_compile_features(odia_mzpeak PUBLIC cxx_std_23)
target_include_directories(odia_mzpeak SYSTEM PUBLIC ${_mzp}/include)
target_include_directories(odia_mzpeak PRIVATE ${_mzp}/subprojects/msnumpress/src/main/cpp)
# Arrow's headers before any other copy (as for odia_library), then Boost.
target_include_directories(odia_mzpeak SYSTEM BEFORE PUBLIC
  "$<TARGET_PROPERTY:${DIALIBGEN_ARROW_TARGET},INTERFACE_INCLUDE_DIRECTORIES>"
  "$<TARGET_PROPERTY:${DIALIBGEN_PARQUET_TARGET},INTERFACE_INCLUDE_DIRECTORIES>")
target_include_directories(odia_mzpeak SYSTEM PUBLIC ${BOOST_INCLUDE_DIR})
target_link_libraries(odia_mzpeak PUBLIC ${DIALIBGEN_ARROW_TARGET} ${DIALIBGEN_PARQUET_TARGET} libzip::zip)
# Boost.JSON is compiled into MzPeakBoostJson.cpp; never auto-link a Boost library (MSVC).
target_compile_definitions(odia_mzpeak PRIVATE BOOST_JSON_NO_LIB BOOST_CONTAINER_NO_LIB BOOST_ALL_NO_LIB)
# The bundled MSNumpress in a namespace of its own (OpenMS links the same one).
set_source_files_properties(${_mzp}/src/util/numpress.cpp ${_mzp}/subprojects/msnumpress/src/main/cpp/MSNumpress.cpp
  PROPERTIES COMPILE_DEFINITIONS "numpress=odia_mzpeak_numpress")
set_target_properties(odia_mzpeak PROPERTIES POSITION_INDEPENDENT_CODE ON)
if(MSVC)
  target_compile_options(odia_mzpeak PRIVATE /W0 /bigobj /utf-8)
  target_compile_definitions(odia_mzpeak PRIVATE NOMINMAX _CRT_SECURE_NO_WARNINGS)
else()
  target_compile_options(odia_mzpeak PRIVATE -w)
endif()
