# Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)


#------------------------------------------------------------------------------
# Setup Conduit
#
# This file defines:
#  CONDUIT_FOUND - If Conduit was found
#  CONDUIT_INCLUDE_DIRS - The Conduit include directories
#  
#  If found, the conduit CMake targets will also be imported
#------------------------------------------------------------------------------

# first Check for CONDUIT_DIR

if(NOT CONDUIT_DIR)
    MESSAGE(FATAL_ERROR "Could not find Conduit. Conduit requires explicit CONDUIT_DIR.")
endif()

set(_conduit_config "${CONDUIT_DIR}/lib/cmake/ConduitConfig.cmake")
if(NOT EXISTS ${_conduit_config})
    MESSAGE(FATAL_ERROR "Could not find Conduit cmake include file ${_conduit_config}")
endif()

find_package(Conduit REQUIRED
             NO_DEFAULT_PATH
             PATHS ${CONDUIT_DIR}/lib/cmake)

# handle the QUIETLY and REQUIRED arguments and set CONDUIT_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(CONDUIT  DEFAULT_MSG
                                  CONDUIT_INCLUDE_DIRS
                                  )





