cmake_minimum_required(VERSION 3.1)

################################
# axom project
################################
project(axom LANGUAGES C CXX)

if (ENABLE_FORTRAN)
   enable_language(Fortran)
endif()

################################
# BLT
################################
if (DEFINED BLT_SOURCE_DIR)
    # Support having a shared BLT outside of the repository if given a BLT_SOURCE_DIR

    if (NOT EXISTS ${BLT_SOURCE_DIR}/SetupBLT.cmake)
        message(FATAL_ERROR "Given BLT_SOURCE_DIR does not contain SetupBLT.cmake")
    endif()
else()
    # Use internal BLT if no BLT_SOURCE_DIR is given

    set(BLT_SOURCE_DIR "${PROJECT_SOURCE_DIR}/blt" CACHE PATH "")
    if (NOT EXISTS ${BLT_SOURCE_DIR}/SetupBLT.cmake)
        message(FATAL_ERROR
            "The BLT submodule is not present. "
            "Run the following two commands in your git repository: \n"
            "    git submodule init \n "
            "    git submodule update" )
    endif()
endif()

set(ENABLE_COPY_HEADERS ON CACHE BOOL "")

include(${BLT_SOURCE_DIR}/SetupBLT.cmake)

################################
# Include standard build system logic and options
################################
include(cmake/CMakeBasics.cmake)

################################
# Add Thirdparty Builtin Libs
################################
add_subdirectory(thirdparty)

#################################
# Add top-level includes
#################################
add_subdirectory(include)

################################
# Add components
################################
add_subdirectory(components)

################################
# Adds directory for utility code
# that uses or supports axom
################################
add_subdirectory(util)

################################
# Add docs
################################
if (ENABLE_DOCS)
  add_subdirectory(docs)
endif()

################################
# Create header of configuration options
################################s
include(cmake/AxomConfig.cmake)
