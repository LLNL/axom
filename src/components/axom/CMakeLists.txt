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

## list of axom headers
set( axom_installed_headers
     ${CMAKE_BINARY_DIR}/include/axom/config.hpp
     cstdint_wrapper.hpp
     Macros.hpp
     Types.hpp
     )

################################
# Add code checks
################################
axom_add_code_checks(PREFIX axom_includes)

## add rules for the make install target
install( FILES ${axom_installed_headers} DESTINATION {CMAKE_INSTALL_PREFIX}/include/axom )
