# Shim for Eigen 3.4's config-version file, which rejects CMake version RANGES.
#
# OpenMSConfig.cmake does `find_dependency(Eigen3 3.4.0...<6)`. Eigen 3.4's own
# Eigen3ConfigVersion.cmake enforces a same-major-version policy and marks any
# range whose maximum crosses 4 as incompatible -- so a perfectly good Eigen
# 3.4.0 is refused and `find_package(OpenMS)` fails before OpenMS's own fallback
# can run. It is the first thing anyone building against a conda OpenMS hits,
# on every platform, and it is not fixable in OpenMS or in Eigen from here.
#
# This shim declares the imported target against the real headers and nothing
# else. It is selected automatically by the top-level CMakeLists, but only when
# a plain range-respecting find_package(Eigen3) has already failed -- so where
# Eigen behaves, this directory is never consulted.
#
# Use directly with: -DEigen3_DIR=<this directory> [-DEIGEN3_INCLUDE_DIR=<headers>]
if(NOT DEFINED EIGEN3_INCLUDE_DIR OR NOT EXISTS "${EIGEN3_INCLUDE_DIR}/signature_of_eigen3_matrix_library")
  find_path(EIGEN3_INCLUDE_DIR signature_of_eigen3_matrix_library
            PATH_SUFFIXES eigen3 include/eigen3)
endif()

if(NOT EIGEN3_INCLUDE_DIR)
  set(Eigen3_FOUND FALSE)
  return()
endif()

if(NOT TARGET Eigen3::Eigen)
  add_library(Eigen3::Eigen INTERFACE IMPORTED)
  set_target_properties(Eigen3::Eigen PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${EIGEN3_INCLUDE_DIR}")
endif()

set(EIGEN3_INCLUDE_DIRS "${EIGEN3_INCLUDE_DIR}")
set(Eigen3_FOUND TRUE)
