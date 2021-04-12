# Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Setup SCR
#------------------------------------------------------------------------------
# This file defines:
#  SCR_FOUND - If SCR was found
#  SCR_INCLUDE_DIRS - The SCR include directories
#  SCR_LIBRARY - The SCR library
#------------------------------------------------------------------------------

# first Check for SCR_DIR

if(NOT SCR_DIR)
    message(FATAL_ERROR "Could not find SCR. SCR support needs explicit SCR_DIR")
endif()

#find includes
find_path( SCR_INCLUDE_DIRS scr.h
           PATHS  ${SCR_DIR}/include/
           NO_DEFAULT_PATH
           NO_CMAKE_ENVIRONMENT_PATH
           NO_CMAKE_PATH
           NO_SYSTEM_ENVIRONMENT_PATH
           NO_CMAKE_SYSTEM_PATH)

find_library( SCR_LIBRARY NAMES scr libscr
              PATHS ${SCR_DIR}/lib
              NO_DEFAULT_PATH
              NO_CMAKE_ENVIRONMENT_PATH
              NO_CMAKE_PATH
              NO_SYSTEM_ENVIRONMENT_PATH
              NO_CMAKE_SYSTEM_PATH)


include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set SCR_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(SCR  DEFAULT_MSG
                                  SCR_INCLUDE_DIRS
                                  SCR_LIBRARY )
