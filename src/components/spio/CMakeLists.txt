################################
# IOParallel
################################
project(spio)

################################
# Add the conduit include path
################################
include_directories(${CONDUIT_INCLUDE_DIRS})

################################
# Add the HDF5 include path
################################
include_directories(${HDF5_INCLUDE_DIRS})

################################
# Add the boost headers include
################################
include_directories( ${Boost_INCLUDE_DIRS} )

################################
# Add the sources
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
if (BUILD_TESTS)
  add_subdirectory(tests)
endif()


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
