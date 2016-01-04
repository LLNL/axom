cmake_minimum_required(VERSION 3.1)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

################################
# ASCTK ToolKit
################################
project(asctk LANGUAGES C CXX)

if (ENABLE_FORTRAN)
   enable_language(Fortran)
endif()

################################
# Include standard build system logic and options
################################
include(CMakeBasics)

################################
# Add Thirdparty Builtin Libs
################################
add_subdirectory(thirdparty)

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
