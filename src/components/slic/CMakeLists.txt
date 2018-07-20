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
axom_component_requires(NAME       SLIC
                        COMPONENTS Axom_utils)

################################
# Add the SLIC sources
################################
add_subdirectory(src)

################################
# Add examples
################################
if (AXOM_ENABLE_EXAMPLES)
  add_subdirectory(examples)
endif()

################################
# Add tests
################################
if (AXOM_ENABLE_TESTS)
  add_subdirectory(tests)
endif()

################################
# Add docs
################################
if (AXOM_ENABLE_DOCS)
  add_subdirectory(docs)
endif()

################################
# Add code checks
################################
axom_add_code_checks(PREFIX slic)

################################
# Create CMake importable
# exports for all of our targets
################################
install(EXPORT slic-targets DESTINATION lib/cmake)
