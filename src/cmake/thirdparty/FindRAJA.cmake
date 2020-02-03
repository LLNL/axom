# Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Setup RAJA
#------------------------------------------------------------------------------

if (NOT RAJA_DIR)
  message(FATAL_ERROR "Could not find RAJA. RAJA_DIR must be explicitly specified when configuring CMake" )
endif()

find_package(RAJA REQUIRED)
if(TARGET RAJA)
   message(STATUS "Found RAJA: ${RAJA_DIR}")
endif()
