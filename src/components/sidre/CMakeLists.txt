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
axom_component_requires(NAME       Sidre
                        COMPONENTS Axom_utils SLIC
                        TPLS       Conduit)

# Note: Sidre's unordered_map is either sparsehash OR C++11's unordered_map
#       Default to std::map and warn if neither is available
if(NOT SPARSEHASH_FOUND)
  if(NOT CMAKE_CXX_STANDARD EQUAL 11)
    if(NOT CMAKE_CXX_STANDARD EQUAL 14)
      message(WARNING "No hashmap implementation available. 
        Defaulting to using an ordered map for Sidre groups.
        For improved performance, set SPARSEHASH_DIR to location of a sparsehash install
        or enable C++11 (e.g. set BLT_CXX_STD to 'c++11').")
    endif()
  endif()
endif()

################################
# Add the sidre sources
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
axom_add_code_checks(PREFIX   sidre
                     EXCLUDES sidre/examples/lulesh2)

################################
# Create CMake importable
# exports for all of our targets
################################
install(EXPORT sidre-targets DESTINATION lib/cmake)
