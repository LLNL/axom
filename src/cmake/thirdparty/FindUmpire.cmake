# Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

if (NOT UMPIRE_DIR)
  message(FATAL_ERROR "Could not find Umpire. UMPIRE_DIR must be explicitly specified when configuring CMake" )
endif()

find_package(umpire REQUIRED PATHS ${UMPIRE_DIR} )
find_package_handle_standard_args( UMPIRE DEFAULT_MSG
                                   UMPIRE_INCLUDE_DIRS
                                   )