#------------------------------------------------------------------------------
# Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
#
# Produced at the Lawrence Livermore National Laboratory
#
# LLNL-CODE-741217
#
# All rights reserved.
#
# This file is part of Axom.
#
# For details about use and distribution, please read axom/LICENSE.
#------------------------------------------------------------------------------
###############################################################################
#
# Setup Axom
#
###############################################################################
#
# Expects AXOM_DIR to point to an Axom installation.
#
# This file defines the following CMake variables:
#  AXOM_FOUND - If Axom was found
#  AXOM_INCLUDE_DIRS - The Axom include directories
#  AXOM_LIBRARIES - The Axom libraries
#
###############################################################################

###############################################################################
# Check for AXOM_DIR
###############################################################################
if(NOT AXOM_DIR)
    MESSAGE(FATAL_ERROR "Could not find Axom. Axom requires explicit AXOM_DIR.")
endif()

set(AXOM_INCLUDE_DIRS ${AXOM_DIR}/include)

# NOTE: This turns on all components of Axom which may not be true for your code
set(AXOM_LIBRARIES sparsehash fmt core )
if (ENABLE_MPI)
    list(APPEND AXOM_LIBRARIES lumberjack)
endif()
list(APPEND AXOM_LIBRARIES slic primal mint slam quest sidre)

foreach(_component ${AXOM_LIBRARIES})
    
    set(_target_file ${AXOM_DIR}/lib/cmake/${_component}-targets.cmake)

    if(NOT EXISTS ${_target_file})
        MESSAGE(FATAL_ERROR "Could not find Axom CMake include file (${_target_file})")
    endif()

    include(${_target_file})

endforeach()

###############################################################################
# Sets AXOM_FOUND if AXOM_INCLUDE_DIRS and AXOM_LIBRARIES are not empty
###############################################################################
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(AXOM  DEFAULT_MSG
                                  AXOM_INCLUDE_DIRS
                                  AXOM_LIBRARIES )
