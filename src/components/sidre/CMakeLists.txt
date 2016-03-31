################################
# Datastore
################################
project(sidre)

################################
# Add the conduit include path
################################
include_directories(${CONDUIT_INCLUDE_DIRS})
include_directories(${HDF5_INCLUDE_DIRS})

################################
# Add the sparsehash headers include
################################
include_directories( ${SPARSEHASH_INCLUDE_DIRS} )

################################
# Add the boost headers include
################################
include_directories( ${Boost_INCLUDE_DIRS} )

################################
# Add the sidre sources
################################
add_subdirectory(src)

################################
# Add examples
################################
if (BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()

################################
# Add tests
################################
if (BUILD_TESTING)
  add_subdirectory(tests)
endif()


add_code_check_targets(uncrustify.cfg)

################################
# Add docs
################################
if (BUILD_DOCS)
  add_subdirectory(docs)
endif()

################################
# Create CMake importable
# exports for all of our targets
################################
install(EXPORT ${PROJECT_NAME}-targets DESTINATION lib/cmake)
