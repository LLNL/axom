#------------------------------------------------------------------------------
# Copyright (c) 2017, Lawrence Livermore National Security, LLC.
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


#------------------------------------------------------------------------------
# Setup MFEM
#
# This file defines:
#  MFEM_FOUND        - If mfem was found
#  MFEM_INCLUDE_DIRS - The mfem include directories
#  MFEM_LIBRARY      - The mfem library
#------------------------------------------------------------------------------

if (NOT MFEM_DIR)
  message(FATAL_ERROR "Cannot find MFEM. MFEM_DIR is not defined. ")
endif()

find_path( MFEM_INCLUDE_DIRS mfem.hpp
           PATHS ${MFEM_DIR}/include/
           NO_DEFAULT_PATH
           NO_CMAKE_ENVIRONMENT_PATH
           NO_CMAKE_PATH
           NO_SYSTEM_ENVIRONMENT_PATH
           NO_CMAKE_SYSTEM_PATH
           )

find_library( MFEM_LIBRARY NAMES mfem
              PATHS ${MFEM_DIR}/lib
              NO_DEFAULT_PATH
              NO_CMAKE_ENVIRONMENT_PATH
              NO_CMAKE_PATH
              NO_SYSTEM_ENVIRONMENT_PATH
              NO_CMAKE_SYSTEM_PATH
              )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MFEM  DEFAULT_MSG
                                  MFEM_INCLUDE_DIRS
                                  MFEM_LIBRARY )
