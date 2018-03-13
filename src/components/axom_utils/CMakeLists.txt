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
# Add the axom_utils sources
################################
add_subdirectory(src)

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
axom_add_code_checks(PREFIX axom_utils)

################################
# Create CMake importable
# exports for all of our targets
################################
install(EXPORT axom_utils-targets DESTINATION lib/cmake)
