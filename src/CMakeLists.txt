cmake_minimum_required(VERSION 3.1)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

################################
# ASCTK ToolKit
################################
project(asctk LANGUAGES C CXX)

################################
# Prevent in-source builds
################################
include(PreventInSourceBuilds)

################################
# Set some standard cmake opts
################################
include(CMakeBasics)

################################
# Setup 3rd Party Libs
################################
include(Setup3rdParty)

################################
# Setup toolkit docs targets
################################
include(SetupDocs)

# Setup toolkit source checks
################################
include(SetupCodeChecks)

################################
# Add components
################################
add_subdirectory(components)

################################
# Add docs
################################
add_subdirectory(docs)

################################
# install targets from build
################################s
install(DIRECTORY ${PROJECT_BINARY_DIR}/bin DESTINATION . )
install(DIRECTORY ${PROJECT_BINARY_DIR}/include DESTINATION .)
