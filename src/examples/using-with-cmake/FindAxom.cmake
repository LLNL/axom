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
#
###############################################################################

###############################################################################
# Check for AXOM_DIR
###############################################################################
if(NOT AXOM_DIR)
    MESSAGE(FATAL_ERROR "Could not find Axom. Axom requires explicit AXOM_DIR.")
endif()

###############################################################################
#TODO Axom needs to export a main "axom.cmake" to make this easier and more 
# useful
###############################################################################

if(NOT EXISTS ${AXOM_DIR}/lib/cmake/axom_utils-targets.cmake)
    MESSAGE(FATAL_ERROR "Could not find Axom CMake include file (${AXOM_DIR}/lib/cmake/axom_utils-targets.cmake)")
endif()

###############################################################################
# Import Axom's CMake targets
###############################################################################
include(${AXOM_DIR}/lib/cmake/axom_utils-targets.cmake)

###############################################################################
# Set remaning CMake variables 
###############################################################################
# we found Axom
set(AXOM_FOUND TRUE)
# provide location of the headers in AXOM_INCLUDE_DIRS
set(AXOM_INCLUDE_DIRS ${AXOM_DIR}/include/)




