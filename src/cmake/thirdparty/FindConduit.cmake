# Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)


#------------------------------------------------------------------------------
# Setup Conduit
#
# This file defines:
#  CONDUIT_FOUND - If Conduit was found
#  
#  If found, the Conduit CMake targets will also be imported
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
