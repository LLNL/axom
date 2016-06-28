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
if (NOT EXISTS ${PROJECT_SOURCE_DIR}/blt/SetupBLT.cmake)
    message(FATAL_ERROR "\
The BLT submodule is not present. \
If in git repository run the following two commands:\n \
git submodule init\n \
git submodule update")

endif()
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

