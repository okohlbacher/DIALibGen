cmake_minimum_required(VERSION 3.21 FATAL_ERROR)

project(DIALibraryGenerator
  VERSION 0.2.2
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

# Without this the macOS floor is whatever SDK the build machine happens to
# have, and it MOVES when a CI runner image is retired. dyld refuses to load a
# Mach-O on an older macOS than its minos, so the failure lands on the user as
# an app that installs, launches and then cannot run -- and a Homebrew cask
# cannot state a floor it does not know.
#
# 13.3 is MEASURED, not chosen: at 13.0 the build fails with
#   TextWriter.h:105: 'to_chars' is unavailable: introduced in macOS 13.3
# because libc++ shipped the floating-point std::to_chars in 13.3. That is the
# whole value of setting a target at all -- without one this was not an error,
# it was a binary that linked fine on the runner and refused to start for a
# user on 13.2.
if(APPLE AND NOT CMAKE_OSX_DEPLOYMENT_TARGET)
  set(CMAKE_OSX_DEPLOYMENT_TARGET "13.3" CACHE STRING "Minimum macOS version" FORCE)
endif()

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

# Every file carrying this project's version must agree with project(VERSION).
# Runs at configure time, so the first cmake after a bump fails when one was
# missed -- which is how a release ends up announcing the wrong number.
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/version-guard.cmake)

# ---------------------------------------------------------------------------
# Dependencies.
#
# OpenMS is consumed as an installed, read-only package and is never modified.
# Note the version reality: no tagged OpenMS release ships the PeptDeep ONNX
# models (WITH_ONNX exists only on develop and defaults OFF), so against a
# released OpenMS the three .onnx files must be supplied separately -- see
# README, "Models".
# ---------------------------------------------------------------------------
# OpenMSConfig.cmake does find_dependency(Eigen3 3.4.0...<6), and Eigen 3.4's
# own config-version file rejects any range whose maximum crosses 4 -- so the
# find_package below fails on a perfectly good Eigen 3.4.0, on every platform,
# before OpenMS's fallback can run. Probed rather than assumed: the shim is
# selected only when a range-respecting lookup has actually failed, so an Eigen
# that behaves (or one the user pointed at themselves) is left alone.
# The probe runs unconditionally, and is NOT skipped when the caller set
# Eigen3_DIR: a caller who points at a perfectly good Eigen 3.4.0 hits the same
# rejection, and skipping the probe for them means the shim exists but never
# engages. find_package honours their Eigen3_DIR here, so an Eigen that
# satisfies the range is left exactly as they chose it.
find_package(Eigen3 3.4.0...<6 QUIET CONFIG)
if(NOT Eigen3_FOUND)
  set(Eigen3_DIR "${CMAKE_CURRENT_SOURCE_DIR}/cmake/eigen3-range-shim"
      CACHE PATH "Eigen3 package directory" FORCE)
  message(STATUS "Eigen3: using cmake/eigen3-range-shim "
                 "(the installed Eigen rejects OpenMS's version range)")
endif()

find_package(OpenMS REQUIRED)
message(STATUS "OpenMS: ${OpenMS_DIR}")

find_package(Arrow CONFIG REQUIRED)
find_package(Parquet CONFIG REQUIRED)
# 19 is the real floor, not the 23 this once declared: the newest API used is
# parquet::arrow::OpenFile's Result-returning overload, added in Arrow 19.0.0.
# The old number excluded bioconda's libarrow 21 -- which is what OpenMS itself
# pins, so a consumer linking both could not satisfy it.
if(Arrow_VERSION VERSION_LESS 19)
  message(FATAL_ERROR "Apache Arrow >= 19 required (found ${Arrow_VERSION}).")
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
# The header lands in a different place in each distribution: bare include/ in
# Microsoft's release tarballs, include/onnxruntime/ in conda-forge's
# onnxruntime-cpp, include/onnxruntime/core/session/ in a source build.
find_path(ONNXRuntime_INCLUDE_DIR NAMES onnxruntime_cxx_api.h
          PATH_SUFFIXES include onnxruntime include/onnxruntime
                        onnxruntime/core/session include/onnxruntime/core/session)
# onnxruntime_conda is conda-forge's name for it on Windows (the import library
# is onnxruntime_conda.lib, the DLL onnxruntime_conda.dll); every other
# distribution calls it onnxruntime. Searching for the plain name alone made a
# perfectly good conda install look absent.
find_library(ONNXRuntime_LIBRARY NAMES onnxruntime onnxruntime_conda
             PATH_SUFFIXES lib)
if(NOT ONNXRuntime_INCLUDE_DIR OR NOT ONNXRuntime_LIBRARY)
  # Name the half that failed: the header and the library are packaged in
  # different places by different distributions, and "not found" for both is a
  # different problem from "not found" for one.
  set(_ort_missing "")
  if(NOT ONNXRuntime_INCLUDE_DIR)
    string(APPEND _ort_missing "\n  onnxruntime_cxx_api.h (set ONNXRuntime_INCLUDE_DIR)")
  endif()
  if(NOT ONNXRuntime_LIBRARY)
    string(APPEND _ort_missing "\n  the onnxruntime library (set ONNXRuntime_LIBRARY)")
  endif()
  message(FATAL_ERROR
    "ONNX Runtime not found. Missing:${_ort_missing}\n"
    "Searched CMAKE_PREFIX_PATH: ${CMAKE_PREFIX_PATH}\n"
    "Install onnxruntime-cpp, or set the variables above.")
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
# BUILD_INTERFACE, not a bare PUBLIC path. These are where THIS machine happens
# to keep ONNX Runtime and Boost; baking them into the exported target would
# hand a consumer an absolute path from the build machine. CMake refuses it
# outright when the path lies inside the source tree -- which is how this
# surfaced, on a Windows runner where contrib is built under the checkout --
# but exporting it was wrong everywhere, just silently. The installed config
# re-finds both for the consumer.
target_include_directories(odia_library SYSTEM PUBLIC
  $<BUILD_INTERFACE:${ONNXRuntime_INCLUDE_DIR}>)
if(BOOST_INCLUDE_DIR)
  target_include_directories(odia_library SYSTEM PUBLIC
    $<BUILD_INTERFACE:${BOOST_INCLUDE_DIR}>)
endif()

target_link_libraries(odia_library PUBLIC
  OpenMS
  ${DIALIBGEN_ARROW_TARGET} ${DIALIBGEN_PARQUET_TARGET}
  ${ONNXRuntime_LIBRARY})
if(OpenMP_CXX_FOUND)
  target_link_libraries(odia_library PUBLIC OpenMP::OpenMP_CXX)
endif()
set_target_properties(odia_library PROPERTIES POSITION_INDEPENDENT_CODE ON)

# The shipped data tables (irt_standards.tsv). The tool looks for them relative
# to its own executable first -- $DIALIBGEN_DATA_DIR, then ../share/DIALibraryGenerator,
# then ../data -- so an installed, packaged or relocated binary finds them. This
# absolute source path is the LAST candidate and exists only so an uninstalled
# build tree works; it is deliberately not the only answer any more.
set(DIALIBGEN_DATA_DIR "${CMAKE_CURRENT_SOURCE_DIR}/data")

# ---------------------------------------------------------------------------
# The TOPP-compatible tool.
# ---------------------------------------------------------------------------
if(DIALIBGEN_BUILD_TOOL)
  add_executable(DIALibraryGenerator src/DIALibraryGenerator.cpp)
  target_link_libraries(DIALibraryGenerator PRIVATE
    odia_library nlohmann_json::nlohmann_json)
  # PRIVATE, and therefore per-target: a second tool in this repo would need
  # its own copy of both defines, and grep -L would find the one that missed.
  target_compile_definitions(DIALibraryGenerator PRIVATE
    ODIA_DATA_DIR="${DIALIBGEN_DATA_DIR}"
    # Without this the tool announces the version of whatever OpenMS tree it was
    # built against -- a version the user cannot file a bug against.
    DIALIBGEN_VERSION="${PROJECT_VERSION}")
  # An installed layout the tool can actually resolve from argv[0]:
  #   <prefix>/bin/DIALibraryGenerator
  #   <prefix>/share/DIALibraryGenerator/irt_standards.tsv
  set_target_properties(DIALibraryGenerator PROPERTIES
    INSTALL_RPATH_USE_LINK_PATH ON)
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
  # share/DIALibraryGenerator, matching what the tool searches for. The old
  # lower-case directory name was reachable by nothing.
  install(DIRECTORY data/ DESTINATION ${CMAKE_INSTALL_DATADIR}/DIALibraryGenerator)

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
