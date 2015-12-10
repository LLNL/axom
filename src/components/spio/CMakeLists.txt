################################
# IOParallel
################################
project(spio)

################################
# Add the conduit include path
################################
include_directories(${CONDUIT_INCLUDE_DIRS})

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
add_subdirectory(examples)

################################
# Add tests
################################
#if (BUILD_TESTING)
#  add_subdirectory(tests)
#endif()


################################
# Add docs
################################
#add_subdirectory(docs)

################################
# Create CMake importable
# exports for all of our targets
################################
install(EXPORT ${PROJECT_NAME}-targets DESTINATION lib/cmake)
