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
# add sources
################################
add_subdirectory(src)

################################
# add examples
################################
if (ENABLE_EXAMPLES)
  add_subdirectory(examples)
endif()

################################
# add tests
################################
if (ENABLE_TESTS)
  add_subdirectory(tests)
endif()

################################
# add docs
################################
if (ENABLE_DOCS)
  add_subdirectory(docs)
endif()

################################
# Add code checks
################################
axom_add_code_checks(PREFIX primal)

################################
# Create CMake importable
# exports for all of our targets
################################
install(EXPORT primal-targets DESTINATION lib/cmake)
