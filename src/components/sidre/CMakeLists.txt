################################
# Datastore
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

if(NOT BOOST_FOUND AND NOT ATK_USE_CXX11)
  message(FATAL_ERROR "Sidre requires Boost, or C++11 support. 
           Set BOOST_DIR to location of a Boost install if specifying Boost.")
endif()

# Enable when sparse hash is required ....   USE_DENSE_HASH_MAP 
# 
#if(NOT SPARSEHASH_FOUND)
#  message(FATAL_ERROR "Sidre requires SparseHash. Set SPARSEHASH_DIR to location of a SparseHash install.")
#endif()


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
