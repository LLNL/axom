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
# Add examples
################################
add_subdirectory(examples)

################################
# Add tests
################################
#add_subdirectory(tests)

################################
# install targets from build
################################s
install(DIRECTORY ${PROJECT_BINARY_DIR}/bin DESTINATION . )
install(DIRECTORY ${PROJECT_BINARY_DIR}/include DESTINATION .)

################################
# Create CMake importable
# exports for all of our targets
################################
install(EXPORT ${PROJECT_NAME}-targets DESTINATION lib/cmake)