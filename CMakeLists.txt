cmake_minimum_required(VERSION 3.21 FATAL_ERROR)

project(DIALibGen
  VERSION 0.11.0
  DESCRIPTION "In-silico DIA spectral library generation from a FASTA, built on OpenMS"
  HOMEPAGE_URL "https://github.com/okohlbacher/DIALibGen"
  LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(APPLE AND NOT CMAKE_OSX_DEPLOYMENT_TARGET)
  set(CMAKE_OSX_DEPLOYMENT_TARGET "13.3" CACHE STRING "Minimum macOS version" FORCE)
endif()

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  set(DIALIBGEN_IS_TOP_LEVEL ON)
else()
  set(DIALIBGEN_IS_TOP_LEVEL OFF)
endif()
option(DIALIBGEN_BUILD_FINETUNE "Build RT/CCS model training with LibTorch" ON)
option(DIALIBGEN_BUILD_TOOL  "Build the DIALibGen executable" ${DIALIBGEN_IS_TOP_LEVEL})

set(ODIA_MODEL_DIR "$ENV{ODIA_OPENMS}/share/OpenMS/models"
    CACHE PATH "Directory holding peptdeep_{rt,ms2,ccs}_dynamic.onnx")
option(DIALIBGEN_BUILD_TESTS "Build the test suite" ${DIALIBGEN_IS_TOP_LEVEL})
option(DIALIBGEN_INSTALL     "Generate install rules" ${DIALIBGEN_IS_TOP_LEVEL})

include(GNUInstallDirs)

include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/version-guard.cmake)

find_package(Eigen3 3.4.0...<6 QUIET CONFIG)
if(NOT Eigen3_FOUND)
  set(Eigen3_DIR "${CMAKE_CURRENT_SOURCE_DIR}/cmake/eigen3-range-shim")
  message(STATUS "Eigen3: using cmake/eigen3-range-shim "
                 "(the installed Eigen rejects OpenMS's version range)")
endif()

# OpenMS's package compatibility check can reject a newer minor version.
# Enforce our source minimum after discovery; the installed SDK pins its ABI.
find_package(OpenMS CONFIG REQUIRED)
if(NOT OpenMS_VERSION OR OpenMS_VERSION VERSION_LESS 3.5)
  message(FATAL_ERROR "OpenMS >= 3.5 required (found ${OpenMS_VERSION}).")
endif()
message(STATUS "OpenMS: ${OpenMS_VERSION} (${OpenMS_DIR})")

find_package(Arrow CONFIG REQUIRED)
find_package(Parquet CONFIG REQUIRED)
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

find_path(ONNXRuntime_INCLUDE_DIR NAMES onnxruntime_cxx_api.h
          PATH_SUFFIXES include onnxruntime include/onnxruntime
                        onnxruntime/core/session include/onnxruntime/core/session)
find_library(ONNXRuntime_LIBRARY NAMES onnxruntime onnxruntime_conda
             PATH_SUFFIXES lib)
if(NOT ONNXRuntime_INCLUDE_DIR OR NOT ONNXRuntime_LIBRARY)
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
if(NOT TARGET ONNXRuntime::ONNXRuntime)
  add_library(ONNXRuntime::ONNXRuntime UNKNOWN IMPORTED)
  set_target_properties(ONNXRuntime::ONNXRuntime PROPERTIES
    IMPORTED_LOCATION "${ONNXRuntime_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${ONNXRuntime_INCLUDE_DIR}")
endif()

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

find_path(BOOST_INCLUDE_DIR NAMES boost/version.hpp)
if(BOOST_INCLUDE_DIR)
  message(STATUS "Boost: ${BOOST_INCLUDE_DIR} (needed transitively by OpenMS headers)")
endif()

find_package(OpenMP)

add_library(odia_library STATIC
  src/library/Library.cpp
  src/library/LibraryGenerator.cpp
  src/library/DIANNLibraryFile.cpp
  src/predict/PeptDeepEncoder.cpp
  src/predict/PeptDeepPredictor.cpp
)
add_library(DIALibGen::odia_library ALIAS odia_library)
target_compile_features(odia_library PUBLIC cxx_std_23)

target_include_directories(odia_library PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
if(BOOST_INCLUDE_DIR)
  target_include_directories(odia_library SYSTEM PUBLIC
    $<BUILD_INTERFACE:${BOOST_INCLUDE_DIR}>)
endif()
# OpenMS contrib may also contain older Arrow headers alongside Boost. Keep the
# headers of the linked Arrow/Parquet targets first, including for SDK consumers.
target_include_directories(odia_library SYSTEM BEFORE PUBLIC
  "$<TARGET_PROPERTY:${DIALIBGEN_ARROW_TARGET},INTERFACE_INCLUDE_DIRECTORIES>"
  "$<TARGET_PROPERTY:${DIALIBGEN_PARQUET_TARGET},INTERFACE_INCLUDE_DIRECTORIES>")

target_link_libraries(odia_library PUBLIC
  OpenMS
  ${DIALIBGEN_ARROW_TARGET} ${DIALIBGEN_PARQUET_TARGET}
  ONNXRuntime::ONNXRuntime)
if(OpenMP_CXX_FOUND)
  target_link_libraries(odia_library PUBLIC OpenMP::OpenMP_CXX)
endif()
set_target_properties(odia_library PROPERTIES POSITION_INDEPENDENT_CODE ON)

add_library(odia_refine STATIC src/refine/LibraryRefiner.cpp)
target_link_libraries(odia_refine PUBLIC odia_library)
add_library(DIALibGen::odia_refine ALIAS odia_refine)
set_target_properties(odia_refine PROPERTIES POSITION_INDEPENDENT_CODE ON)
set(DIALIBGEN_TARGETS odia_library odia_refine)
if(DIALIBGEN_BUILD_FINETUNE)
  find_package(Torch REQUIRED CONFIG)
  find_package(Qt6 REQUIRED COMPONENTS Core)
  if(TORCH_CXX_FLAGS MATCHES "_GLIBCXX_USE_CXX11_ABI=0")
    message(FATAL_ERROR "LibTorch must use the CXX11 ABI to link with OpenMS and Arrow.")
  endif()
  add_library(odia_tune STATIC src/finetune/OnnxWeights.cpp src/finetune/Trainer.cpp)
  target_link_libraries(odia_tune PUBLIC odia_library torch PRIVATE nlohmann_json::nlohmann_json Qt6::Core)
  separate_arguments(DIALIBGEN_TORCH_FLAGS NATIVE_COMMAND "${TORCH_CXX_FLAGS}")
  target_compile_options(odia_tune PUBLIC ${DIALIBGEN_TORCH_FLAGS})
  add_library(DIALibGen::odia_tune ALIAS odia_tune)
  set_target_properties(odia_tune PROPERTIES POSITION_INDEPENDENT_CODE ON)
  list(APPEND DIALIBGEN_TARGETS odia_tune)
endif()

set(DIALIBGEN_DATA_DIR "${CMAKE_CURRENT_SOURCE_DIR}/data")

if(DIALIBGEN_BUILD_TOOL)
  add_executable(DIALibGen src/DIALibGen.cpp src/Generate.cpp src/Refine.cpp)
  target_link_libraries(DIALibGen PRIVATE
    odia_refine nlohmann_json::nlohmann_json)
  if(DIALIBGEN_BUILD_FINETUNE)
    target_link_libraries(DIALibGen PRIVATE odia_tune)
    target_compile_definitions(DIALibGen PRIVATE DIALIBGEN_WITH_FINETUNE)
    # CMake omits in-tree SDKs from INSTALL_RPATH_USE_LINK_PATH. A source
    # install still needs Torch; portable packages later replace this RPATH.
    if(UNIX)
      set_property(TARGET DIALibGen APPEND PROPERTY INSTALL_RPATH "$<TARGET_FILE_DIR:torch>")
    endif()
  endif()
  target_compile_definitions(DIALibGen PRIVATE
    ODIA_DATA_DIR="${DIALIBGEN_DATA_DIR}"
    DIALIBGEN_VERSION="${PROJECT_VERSION}")
  set_target_properties(DIALibGen PROPERTIES
    INSTALL_RPATH_USE_LINK_PATH ON)
endif()

if(DIALIBGEN_INSTALL)
  include(CMakePackageConfigHelpers)
  install(FILES LICENSE THIRD-PARTY-NOTICES.md DESTINATION ${CMAKE_INSTALL_DATADIR}/DIALibGen)
  install(DIRECTORY licenses/ DESTINATION ${CMAKE_INSTALL_DATADIR}/DIALibGen/licenses)

  install(TARGETS ${DIALIBGEN_TARGETS} EXPORT DIALibGenTargets
          ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
  if(DIALIBGEN_BUILD_TOOL)
    install(TARGETS DIALibGen RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

    set(_dlg_models "")
    foreach(_m IN ITEMS peptdeep_rt_dynamic.onnx
                        peptdeep_ms2_dynamic.onnx
                        peptdeep_ccs_dynamic.onnx)
      if(EXISTS "${ODIA_MODEL_DIR}/${_m}")
        list(APPEND _dlg_models "${ODIA_MODEL_DIR}/${_m}")
      endif()
    endforeach()
    list(LENGTH _dlg_models _dlg_model_count)
    if(_dlg_model_count EQUAL 3)
      install(FILES ${_dlg_models}
              DESTINATION ${CMAKE_INSTALL_DATADIR}/DIALibGen/models)
      message(STATUS "DIALibGen: bundling the PeptDeep models from ${ODIA_MODEL_DIR}")
    else()
      message(STATUS "DIALibGen: ${_dlg_model_count}/3 PeptDeep models in "
                     "'${ODIA_MODEL_DIR}' -- the install will not carry them, "
                     "run dialibgen-fetch-models after installing")
    endif()

    install(PROGRAMS scripts/fetch-models.sh
            DESTINATION ${CMAKE_INSTALL_BINDIR}
            RENAME dialibgen-fetch-models)
  endif()
  install(DIRECTORY include/odia DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
  install(DIRECTORY data/ DESTINATION ${CMAKE_INSTALL_DATADIR}/DIALibGen)

  install(DIRECTORY cmake/eigen3-range-shim DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/DIALibGen)

  install(EXPORT DIALibGenTargets
          FILE DIALibGenTargets.cmake
          NAMESPACE DIALibGen::
          DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/DIALibGen)

  configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/DIALibGenConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/DIALibGenConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/DIALibGen)
  write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/DIALibGenConfigVersion.cmake
    COMPATIBILITY SameMinorVersion)
  install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/DIALibGenConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/DIALibGenConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/DIALibGen)
endif()

if(DIALIBGEN_BUILD_TESTS)
  enable_testing()
  add_subdirectory(test)
  add_subdirectory(test/refinement)
endif()
