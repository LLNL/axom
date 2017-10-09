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
# primal
################################
project(primal)

################################
# add sources
################################
add_subdirectory(src)

################################
# add examples
################################
# if (ENABLE_EXAMPLES)
#   add_subdirectory(examples)
# endif()

################################
# add tests
################################
if (ENABLE_TESTS)
  add_subdirectory(tests)
endif()

################################
# add uncrustify
################################
blt_add_code_check_targets(uncrustify.cfg)

################################
# add docs
################################
if (ENABLE_DOCS)
  add_subdirectory(docs)
endif()

################################
# Create CMake importable
# exports for all of our targets
################################
install(EXPORT ${PROJECT_NAME}-targets DESTINATION lib/cmake)
