# Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Setup Axom
#------------------------------------------------------------------------------
# Expects AXOM_DIR to point to an Axom installation.
#
# This file defines the following CMake variables:
#  AXOM_FOUND - If Axom was found
#  AXOM_INCLUDE_DIRS - The Axom include directories
#  AXOM_LIBRARIES - The Axom libraries
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Check for AXOM_DIR
#------------------------------------------------------------------------------
if(NOT AXOM_DIR)
    MESSAGE(FATAL_ERROR "Could not find Axom. Axom requires explicit AXOM_DIR.")
endif()

set(AXOM_INCLUDE_DIRS ${AXOM_DIR}/include)

# NOTE: fmt, CLI11 and sparsehash are useful open-source projects that Axom
# uses internally and we export for other codes use as well
set(AXOM_LIBRARIES sparsehash cli11 fmt axom )

foreach(_library ${AXOM_LIBRARIES})
    
    set(_target_file ${AXOM_DIR}/lib/cmake/${_library}-targets.cmake)

    if(NOT EXISTS ${_target_file})
        MESSAGE(FATAL_ERROR "Could not find Axom CMake exported target file (${_target_file})")
    endif()

    include(${_target_file})

endforeach()

#------------------------------------------------------------------------------
# Sets AXOM_FOUND if AXOM_INCLUDE_DIRS and AXOM_LIBRARIES are not empty
#------------------------------------------------------------------------------
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(AXOM  DEFAULT_MSG
                                  AXOM_INCLUDE_DIRS
                                  AXOM_LIBRARIES )
