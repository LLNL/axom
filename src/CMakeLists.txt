cmake_minimum_required(VERSION 3.1)

################################
# ASCTK ToolKit
################################
project(asctk LANGUAGES C CXX)

if (ENABLE_FORTRAN)
   enable_language(Fortran)
endif()

################################
# Include blt
################################
include(blt/SetupBLT.cmake)

################################
# Include standard build system logic and options
################################
include(cmake/CMakeBasics.cmake)

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
if (ENABLE_DOCS)
  add_subdirectory(docs)
endif()

################################
# install targets from build
################################s
# NOTE: Commented out installation of bin directory for now since contents (shroud) 
#       don't currently work in other directories.
#install(DIRECTORY ${PROJECT_BINARY_DIR}/bin DESTINATION . OPTIONAL)
install(DIRECTORY ${PROJECT_BINARY_DIR}/include DESTINATION . OPTIONAL)
install(DIRECTORY ${PROJECT_BINARY_DIR}/lib/fortran DESTINATION lib OPTIONAL)

