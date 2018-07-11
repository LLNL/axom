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

################################
# Check necessary dependencies
################################
axom_ensure_dependencies(NAME       Mint
                         COMPONENTS Axom_utils SLIC)

################################
# Add mint options
################################
option( MINT_USE_64BIT_INDEXTYPE "Use 64-bit indices for Mint" OFF )
cmake_dependent_option( MINT_USE_SIDRE "Enables Sidre support in Mint" ON
                       "ENABLE_SIDRE" OFF)

################################
# Add mint sources
################################
add_subdirectory(src)

################################
# Add examples
################################
if (ENABLE_EXAMPLES)
  add_subdirectory(examples)
endif()

################################
# Add tests
################################
if (ENABLE_TESTS)
  add_subdirectory(tests)
endif()

################################
# Add docs
################################
if (ENABLE_DOCS)
  add_subdirectory(docs)
endif()

################################
# Add code checks
################################
axom_add_code_checks(PREFIX mint)

################################
# Create CMake importable
# exports for all of our targets
################################
install(EXPORT mint-targets DESTINATION lib/cmake)
