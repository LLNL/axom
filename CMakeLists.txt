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
# Add the datastore lib
################################
# add datastore include dirs so units tests have access to the headers
add_subdirectory(datastore)

################################
# Add docs
################################
add_subdirectory(docs)

################################
# Create CMake importable
# exports for all of our targets
################################
#install(EXPORT datastore DESTINATION lib/cmake)