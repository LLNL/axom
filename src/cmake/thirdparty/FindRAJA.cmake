# Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

if (NOT RAJA_DIR)
  message(FATAL_ERROR "Could not find RAJA. RAJA_DIR must be explicitly specified when configuring CMake" )
endif()

find_package(RAJA REQUIRED)
find_package_handle_standard_args( RAJA DEFAULT_MSG
                                   RAJA_INCLUDE_DIR
                                   RAJA_LIB_DIR
                                   )
