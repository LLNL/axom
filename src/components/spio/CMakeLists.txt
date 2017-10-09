###############################################################################
# Copyright (c) 2017, Lawrence Livermore National Security, LLC.
#
# Produced at the Lawrence Livermore National Laboratory
#
# LLNL-CODE-xxxxxx
#
# All rights reserved.
#
# This file is part of Axom.
#
# For details about use and distribution, please read axom/LICENSE.
#
##############################################################################

################################
# IOParallel
################################
project(spio)

################################
# Add the sources
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
# Create CMake importable
# exports for all of our targets
################################
install(EXPORT ${PROJECT_NAME}-targets DESTINATION lib/cmake)
