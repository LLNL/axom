# Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Example that shows how to use an installed instance of Axom 
# in a BLT-based build system.
#
# Required configuration variables are stored in a CMake cache file 
# called 'host-config.cmake' in the same directory
#
# This file defines paths to BLT, axom and optionally conduit.
# It also contains information about the compiler used to build axom.
#------------------------------------------------------------------------------
#
# To build:
#  mkdir build
#  cd build
#  cmake -C ../host-config.cmake ..
#  make
#  ./bin/example
#
#------------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.8)

project(using_with_blt)

# Define some config options
option(EXAMPLE_VERBOSE_OUTPUT "Verbose output about imported targets" OFF)

# Simple function to check validity of a path variable, e.g. `FOO_DIR`
# and a required file under it of the form `${FOO_DIR}/some/subpath/file`
function(check_var_and_path)
    set(_variable ${ARGV0})
    set(_file     ${ARGV1})

    if(NOT DEFINED ${_variable})
        message(FATAL_ERROR "Required variable '${_variable}' not provided."
                            " Please supply '${_variable}' to configuration." )
    endif()

    if(NOT EXISTS ${_file})
        message (FATAL_ERROR "Path '${_variable}' does not contain required"
                             " file '${_file}'.")
    endif()
endfunction()



#------------------------------------------------------------------------------
# Set up BLT with validity checks
#------------------------------------------------------------------------------

# Check that path to BLT is provided and valid
check_var_and_path(BLT_SOURCE_DIR ${BLT_SOURCE_DIR}/SetupBLT.cmake)

include(${BLT_SOURCE_DIR}/SetupBLT.cmake)


#------------------------------------------------------------------------------
# Check for AXOM_DIR and use CMake's find_package to import axom's targets
#------------------------------------------------------------------------------
check_var_and_path(AXOM_DIR ${AXOM_DIR}/lib/cmake/axom-config.cmake)

find_package(axom REQUIRED
             NO_DEFAULT_PATH 
             PATHS ${AXOM_DIR}/lib/cmake)

#------------------------------------------------------------------------------
# Check for CONDUIT_DIR -- required if Axom was built with Conduit
# Note: 'AXOM_USE_CONDUIT' variable is provided by Axom's config file
#------------------------------------------------------------------------------
if(AXOM_USE_CONDUIT)
    check_var_and_path(CONDUIT_DIR ${CONDUIT_DIR}/lib/cmake/ConduitConfig.cmake)

    find_package(Conduit REQUIRED
                 NO_DEFAULT_PATH 
                 PATHS ${CONDUIT_DIR}/lib/cmake)
endif()

#------------------------------------------------------------------------------
# Set up example target that depends on axom
#------------------------------------------------------------------------------

blt_add_executable(NAME       example 
                   SOURCES    example.cpp
                   DEPENDS_ON axom fmt)


#------------------------------------------------------------------------------
# Optionally, print out information about imported targets
#------------------------------------------------------------------------------
if(EXAMPLE_VERBOSE_OUTPUT)
    foreach(_target axom cli11 fmt sparsehash conduit::conduit example)
        blt_print_target_properties(TARGET ${_target})  
    endforeach()
endif()

