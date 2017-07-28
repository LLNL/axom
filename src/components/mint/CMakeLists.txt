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
# mint
################################
project(mint)

################################
# Add mint sources
################################
add_subdirectory(src)

################################
# Add examples
################################
#if (ENABLE_EXAMPLES)
#  add_subdirectory(examples)
#endif()

################################
# Add tests
################################
if (ENABLE_TESTS)
  add_subdirectory(tests)
endif()

################################
# add uncrustify
################################
add_code_check_targets(uncrustify.cfg)

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
