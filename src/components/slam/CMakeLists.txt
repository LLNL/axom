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

##################################################################################
# SLAM component -- Set-theoretic Lightweight API for Meshes (low level mesh API) 
##################################################################################

### Some additional flags that might be useful
# To avoid warnings: -Wno-missing-braces  
# Other: -fno-strict-aliasing -Wextra

blt_append_custom_compiler_flag(
            FLAGS_VAR SLAM_VECTORIZER_FLAGS
            GNU "-ftree-vectorize -fdump-tree-vect -ftree-vectorizer-verbose=7 -fdump-tree-optimized-details " #-march=native "
            )

################################
# Add the slam (meshapi) sources
################################
add_subdirectory(src)

################################
# Add examples
################################
if (ENABLE_EXAMPLES)
  add_subdirectory(examples)
endif()

################################
# Add tests and benchmarks
################################
if (ENABLE_TESTS)
  add_subdirectory(tests)
  if (ENABLE_BENCHMARKS)
    add_subdirectory(benchmarks)
  endif()
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
axom_add_code_checks(
        PREFIX slam
        EXCLUDES 
            slam/examples/lulesh2.0.3
            slam/examples/tinyHydro )

################################
# Create CMake importable
# exports for all of our targets
################################
install(EXPORT slam-targets DESTINATION lib/cmake)
