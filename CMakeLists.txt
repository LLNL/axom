cmake_minimum_required(VERSION 3.1)

################################
# Datastore
################################

project(Datastore)

################################
# Prevent in-source builds
################################
include(cmake/PreventInSourceBuilds.cmake)

################################
# Set some standard cmake opts
################################
include(cmake/CMakeBasics.cmake)

################################
# Setup 3rd Party Libs
################################
include(cmake/Setup3rdParty.cmake)

################################
# Add the datastore lib
################################
add_subdirectory(datastore)

################################
# Add docs
################################
add_subdirectory(docs)

################################
# Add tests
################################
add_subdirectory(tests)

################################
# Create CMake importable
# exports for all of our targets
################################
#install(EXPORT datastore DESTINATION lib/cmake)