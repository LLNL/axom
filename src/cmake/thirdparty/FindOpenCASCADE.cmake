# Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Setup OpenCASCADE
#------------------------------------------------------------------------------
# This file defines:
#  OPENCASCADE_FOUND - If OpenCASCADE was found
#  OPENCASCADE_INCLUDE_DIRS - The OpenCASCADE include directories
#  OPENCASCADE_LIBRARIES - The OpenCASCADE libraries
#------------------------------------------------------------------------------

# Note: Axom is using OpenCASCADE for its data exchange component
axom_assert_is_directory(DIR_VARIABLE OPENCASCADE_DIR)
find_package(OpenCASCADE CONFIG QUIET NO_DEFAULT_PATH
    HINTS ${OPENCASCADE_DIR}/lib/cmake
          ${OPENCASCADE_DIR}
)

#Update the include dir to not include 'opencascade'
string(REPLACE "include/opencascade" "include" OpenCASCADE_INCLUDE_DIR ${OpenCASCADE_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set OpenCASCADE_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(OpenCASCADE DEFAULT_MSG
                                  OpenCASCADE_INCLUDE_DIR
                                  OpenCASCADE_LIBRARIES )

if(NOT OpenCASCADE_FOUND)
    message(FATAL_ERROR "OPENCASCADE_DIR is not a path to a valid OpenCASCADE install: ${OPENCASCADE_DIR}")
endif()

message(STATUS "OpenCASCADE includes: ${OpenCASCADE_INCLUDE_DIR}")
message(STATUS "OpenCASCADE libraries: ${OpenCASCADE_LIBRARIES}")
