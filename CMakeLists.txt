cmake_minimum_required(VERSION 3.21 FATAL_ERROR)

project(DIALibraryGenerator
  VERSION 0.1.0
  DESCRIPTION "In-silico DIA spectral library generation from a FASTA, built on OpenMS"
  HOMEPAGE_URL "https://github.com/okohlbacher/DIALibraryGenerator"
  LANGUAGES CXX)

# OpenMS's imported target declares INTERFACE_COMPILE_FEATURES cxx_std_23, so
# consumers are raised to C++23 whatever is set here. The sources themselves
# compile at C++20 (and at C++2b with a stock GCC 13); the binding constraint is
# the OpenMS install, not this code.
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Built standalone, or pulled in by a consumer via add_subdirectory/FetchContent.
# In the latter case the tests and the install rules are the consumer's business,
# not ours.
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  set(DIALIBGEN_IS_TOP_LEVEL ON)
else()
  set(DIALIBGEN_IS_TOP_LEVEL OFF)
endif()
option(DIALIBGEN_BUILD_TOOL  "Build the DIALibraryGenerator executable" ${DIALIBGEN_IS_TOP_LEVEL})
option(DIALIBGEN_BUILD_TESTS "Build the test suite" ${DIALIBGEN_IS_TOP_LEVEL})
option(DIALIBGEN_INSTALL     "Generate install rules" ${DIALIBGEN_IS_TOP_LEVEL})

include(GNUInstallDirs)

# ---------------------------------------------------------------------------
# Dependencies.
#
# OpenMS is consumed as an installed, read-only package and is never modified.
# Note the version reality: no tagged OpenMS release ships the PeptDeep ONNX
# models (WITH_ONNX exists only on develop and defaults OFF), so against a
# released OpenMS the three .onnx files must be supplied separately -- see
# README, "Models".
# ---------------------------------------------------------------------------
find_package(OpenMS REQUIRED)
message(STATUS "OpenMS: ${OpenMS_DIR}")

find_package(Arrow CONFIG REQUIRED)
find_package(Parquet CONFIG REQUIRED)
if(Arrow_VERSION VERSION_LESS 23)
  message(FATAL_ERROR "Apache Arrow >= 23 required (found ${Arrow_VERSION}).")
endif()
message(STATUS "Arrow: ${Arrow_VERSION}")

if(TARGET Arrow::arrow_shared)
  set(DIALIBGEN_ARROW_TARGET Arrow::arrow_shared)
else()
  set(DIALIBGEN_ARROW_TARGET Arrow::arrow_static)
endif()
if(TARGET Parquet::parquet_shared)
  set(DIALIBGEN_PARQUET_TARGET Parquet::parquet_shared)
else()
  set(DIALIBGEN_PARQUET_TARGET Parquet::parquet_static)
endif()

# ONNX Runtime drives the PeptDeep models. OpenMS's own binding is not usable
# here: it rejects modified peptides -- PeptDeepInputBuilder's only entry points
# are buildUnmodified*Batch() and it zero-fills the mod_x tensor -- and its
# headers are not installed.
find_path(ONNXRuntime_INCLUDE_DIR NAMES onnxruntime_cxx_api.h
          PATH_SUFFIXES include onnxruntime/core/session)
find_library(ONNXRuntime_LIBRARY NAMES onnxruntime PATH_SUFFIXES lib)
if(NOT ONNXRuntime_INCLUDE_DIR OR NOT ONNXRuntime_LIBRARY)
  message(FATAL_ERROR "ONNX Runtime not found. Set ONNXRuntime_INCLUDE_DIR and "
                      "ONNXRuntime_LIBRARY, or install onnxruntime-cpp.")
endif()
message(STATUS "ONNX Runtime: ${ONNXRuntime_LIBRARY}")

# nlohmann/json parses the library config. Declared explicitly: it was
# previously satisfied only because one conda prefix happened to hold it, which
# is invisible until someone builds against a system OpenMS.
find_package(nlohmann_json CONFIG QUIET)
if(NOT nlohmann_json_FOUND)
  find_path(NLOHMANN_JSON_INCLUDE_DIR NAMES nlohmann/json.hpp)
  if(NOT NLOHMANN_JSON_INCLUDE_DIR)
    message(FATAL_ERROR "nlohmann/json not found. Install nlohmann_json, or set "
                        "NLOHMANN_JSON_INCLUDE_DIR.")
  endif()
  add_library(nlohmann_json::nlohmann_json INTERFACE IMPORTED)
  set_target_properties(nlohmann_json::nlohmann_json PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${NLOHMANN_JSON_INCLUDE_DIR}")
endif()

# Boost reaches this code only through OpenMS's public headers, and OpenMS's
# imported target does not export Boost's include path. Same latent defect as
# nlohmann above: it resolves by accident in a single-prefix environment.
# Located by header rather than by find_package: only the include path is
# needed, and FindBoost is removed in CMake 4 (policy CMP0167).
find_path(BOOST_INCLUDE_DIR NAMES boost/version.hpp)
if(BOOST_INCLUDE_DIR)
  message(STATUS "Boost: ${BOOST_INCLUDE_DIR} (needed transitively by OpenMS headers)")
endif()

# Optional: the ONNX sessions parallelise over peptide batches. Without OpenMP
# the pragmas are ignored and the build is slower, not broken.
find_package(OpenMP)

# ---------------------------------------------------------------------------
# The library. Consumers link this; the executable is a thin driver over it.
# ---------------------------------------------------------------------------
add_library(odia_library STATIC
  src/library/Library.cpp
  src/library/LibraryGenerator.cpp
  src/library/DIANNLibraryFile.cpp
  src/predict/PeptDeepEncoder.cpp
  src/predict/PeptDeepPredictor.cpp
)
add_library(DIALibraryGenerator::odia_library ALIAS odia_library)

target_include_directories(odia_library PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_include_directories(odia_library SYSTEM PUBLIC ${ONNXRuntime_INCLUDE_DIR})
if(BOOST_INCLUDE_DIR)
  target_include_directories(odia_library SYSTEM PUBLIC ${BOOST_INCLUDE_DIR})
endif()

target_link_libraries(odia_library PUBLIC
  OpenMS
  ${DIALIBGEN_ARROW_TARGET} ${DIALIBGEN_PARQUET_TARGET}
  ${ONNXRuntime_LIBRARY})
if(OpenMP_CXX_FOUND)
  target_link_libraries(odia_library PUBLIC OpenMP::OpenMP_CXX)
endif()
set_target_properties(odia_library PROPERTIES POSITION_INDEPENDENT_CODE ON)

# Where the shipped data tables live. Kept identical in name and meaning to the
# definition this code was extracted with, so no source change was needed.
# KNOWN LIMITATION: this is an absolute path into the source tree, so an
# installed binary still resolves its iRT standards relative to the build.
# -irt_standards overrides it, and the default path (irt_rescale=false) reads no
# data file at all. See README, "Known limitations".
set(DIALIBGEN_DATA_DIR "${CMAKE_CURRENT_SOURCE_DIR}/data")

# ---------------------------------------------------------------------------
# The TOPP-compatible tool.
# ---------------------------------------------------------------------------
if(DIALIBGEN_BUILD_TOOL)
  add_executable(DIALibraryGenerator src/DIALibraryGenerator.cpp)
  target_link_libraries(DIALibraryGenerator PRIVATE
    odia_library nlohmann_json::nlohmann_json)
  target_compile_definitions(DIALibraryGenerator PRIVATE
    ODIA_DATA_DIR="${DIALIBGEN_DATA_DIR}")
endif()

# ---------------------------------------------------------------------------
# Install and export, so a consumer can find_package(DIALibraryGenerator).
# ---------------------------------------------------------------------------
if(DIALIBGEN_INSTALL)
  include(CMakePackageConfigHelpers)

  install(TARGETS odia_library EXPORT DIALibraryGeneratorTargets
          ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
  if(DIALIBGEN_BUILD_TOOL)
    install(TARGETS DIALibraryGenerator RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
  endif()
  install(DIRECTORY include/odia DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
  install(DIRECTORY data/ DESTINATION ${CMAKE_INSTALL_DATADIR}/dialibrarygenerator)

  install(EXPORT DIALibraryGeneratorTargets
          FILE DIALibraryGeneratorTargets.cmake
          NAMESPACE DIALibraryGenerator::
          DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/DIALibraryGenerator)

  configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/DIALibraryGeneratorConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/DIALibraryGeneratorConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/DIALibraryGenerator)
  write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/DIALibraryGeneratorConfigVersion.cmake
    COMPATIBILITY SameMinorVersion)
  install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/DIALibraryGeneratorConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/DIALibraryGeneratorConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/DIALibraryGenerator)
endif()

if(DIALIBGEN_BUILD_TESTS)
  enable_testing()
  add_subdirectory(test)
endif()
