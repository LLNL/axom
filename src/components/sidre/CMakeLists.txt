################################
# Sidre
################################
project(sidre)

################################
# Check necessary dependencies
################################
if(NOT CONDUIT_FOUND)
  message(FATAL_ERROR "Sidre requires Conduit. Set CONDUIT_DIR to location of a Conduit install.")
endif()

if(NOT HDF5_FOUND)
  message(FATAL_ERROR "Sidre requires HDF5. Set HDF5_DIR to location of a HDF5 install.")
endif()

# Note: sidre depends on boost OR C++11 b/c USE_UNORDERED_MAP is defined in DataGroup.hpp

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
if (ENABLE_EXAMPLES)
  add_subdirectory(examples)
endif()

################################
# Add tests
################################
if (ENABLE_TESTS)
  add_subdirectory(tests)
endif()


blt_add_code_check_targets(uncrustify.cfg)

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
