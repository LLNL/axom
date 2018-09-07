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

if (NOT RAJA_DIR)
  message(FATAL_ERROR "Could not find RAJA. RAJA_DIR must be explicitly specified when configuring CMake" )
endif()

find_package(RAJA REQUIRED)
find_package_handle_standard_args( RAJA DEFAULT_MSG
                                   RAJA_INCLUDE_DIR
                                   RAJA_LIB_DIR
                                   )
