# Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Setup MFEM
#------------------------------------------------------------------------------
# This file defines:
#  MFEM_FOUND        - If mfem was found
#  MFEM_INCLUDE_DIRS - The mfem include directories
#  MFEM_LIBRARY      - The mfem library
#------------------------------------------------------------------------------

if (NOT MFEM_DIR)
  message(FATAL_ERROR "Cannot find MFEM. MFEM_DIR is not defined. ")
endif()

# Using the CMake build system on Windows, and make-based system otherwise
set(_use_mfem_config FALSE)
if(WIN32)
    set(_use_mfem_config TRUE)
endif()

if(_use_mfem_config)
    # Allow for several different configurations of MFEM
    find_package(mfem CONFIG 
        REQUIRED
        HINTS ${MFEM_DIR}/cmake/mfem 
              ${MFEM_DIR}/lib/cmake/mfem
              ${MFEM_DIR}/share/cmake/mfem
              ${MFEM_DIR}/share/mfem
              ${MFEM_DIR}/mfem)

    if(NOT MEFM_LIBRARY)
       set(MFEM_LIBRARY ${MFEM_LIBRARIES})
    endif()

else()

    find_path( MFEM_INCLUDE_DIRS mfem.hpp
               PATHS 
                ${MFEM_DIR}/include/
                ${MFEM_DIR}
               NO_DEFAULT_PATH
               NO_CMAKE_ENVIRONMENT_PATH
               NO_CMAKE_PATH
               NO_SYSTEM_ENVIRONMENT_PATH
               NO_CMAKE_SYSTEM_PATH
               )

    find_library( MFEM_LIBRARY NAMES mfem
                  PATHS 
                    ${MFEM_DIR}/lib
                    ${MFEM_DIR}
                  NO_DEFAULT_PATH
                  NO_CMAKE_ENVIRONMENT_PATH
                  NO_CMAKE_PATH
                  NO_SYSTEM_ENVIRONMENT_PATH
                  NO_CMAKE_SYSTEM_PATH
                  )
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MFEM  DEFAULT_MSG
                                  MFEM_INCLUDE_DIRS
                                  MFEM_LIBRARY )
