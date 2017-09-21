#-------------------------------------------------------------------------------
# Copyright (c) 2015, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
#
# All rights reserved.
#
# This source code cannot be distributed without permission and further
# review from Lawrence Livermore National Laboratory.
#-------------------------------------------------------------------------------

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
