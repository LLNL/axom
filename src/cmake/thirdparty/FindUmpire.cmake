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

if (NOT UMPIRE_DIR)
  message(FATAL_ERROR "Could not find Umpire. UMPIRE_DIR must be explicitly specified when configuring CMake" )
endif()

find_package(umpire REQUIRED PATHS ${UMPIRE_DIR} )
find_package_handle_standard_args( UMPIRE DEFAULT_MSG
                                   UMPIRE_INCLUDE_DIRS
                                   )