# Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Axom project
#------------------------------------------------------------------------------

if (ENABLE_CUDA)
    cmake_minimum_required(VERSION 3.9)
else()
    cmake_minimum_required(VERSION 3.8.2)
endif()

project(axom LANGUAGES C CXX)

if (ENABLE_FORTRAN)
   enable_language(Fortran)
endif()

#------------------------------------------------------------------------------
# Initialize BLT build system
#------------------------------------------------------------------------------
if (DEFINED BLT_SOURCE_DIR)
    # Support having a shared BLT outside of the repository if given a BLT_SOURCE_DIR

    if (NOT EXISTS ${BLT_SOURCE_DIR}/SetupBLT.cmake)
        message(FATAL_ERROR "Given BLT_SOURCE_DIR does not contain SetupBLT.cmake")
    endif()
else()
    # Use internal 'blt' submodule path if BLT_SOURCE_DIR not provided
    set(BLT_SOURCE_DIR "${PROJECT_SOURCE_DIR}/cmake/blt" CACHE PATH "")
    if (NOT EXISTS ${BLT_SOURCE_DIR}/SetupBLT.cmake)
        message(FATAL_ERROR
            "Cannot locate BLT. "
            "Either run the following two commands in your git repository: \n"
            "    git submodule init\n"
            "    git submodule update\n"
            "Or add -DBLT_SOURCE_DIR=/path/to/blt to your CMake command." )
    endif()
endif()

set(ENABLE_GMOCK ON CACHE STRING "")

if (“${PROJECT_SOURCE_DIR}” STREQUAL “${CMAKE_SOURCE_DIR}”)
    # Set some default BLT options before loading BLT only if not included in
    # another project
    if (NOT BLT_CXX_STD)
        set(BLT_CXX_STD "c++11" CACHE STRING "")
    endif()

    # These are not used in Axom, turn them off
    set(_unused_blt_tools
        CLANGQUERY
        CLANGTIDY
        VALGRIND
        ASTYLE
        CMAKEFORMAT
        UNCRUSTIFY
        YAPF)
    foreach(_tool ${_unused_blt_tools})
        set(ENABLE_${_tool} OFF CACHE BOOL "")
    endforeach()

    # These are only used by Axom developers, so turn them off
    # unless an explicit executable path is given
    set(_used_blt_tools
        CLANGFORMAT
        CPPCHECK
        DOXYGEN
        SPHINX)
    foreach(_tool ${_used_blt_tools})
        if(NOT ${_tool}_EXECUTABLE)
            set(ENABLE_${_tool} OFF CACHE BOOL "")
        else()
            set(ENABLE_${_tool} ON CACHE BOOL "")
        endif()
    endforeach()

    set(BLT_REQUIRED_CLANGFORMAT_VERSION  "10" CACHE STRING "")

    # Mark BLT's built-in third-party targets as EXPORTABLE so they can be added
    # to the axom-targets export set
    set(BLT_EXPORT_THIRDPARTY ON CACHE BOOL "")
endif()

include(${BLT_SOURCE_DIR}/SetupBLT.cmake)

#------------------------------------------------------------------------------
# Attempt to set Axom's data directory (used for Axom's tests and examples)
#------------------------------------------------------------------------------
if(NOT AXOM_DATA_DIR)
    # Use internal 'axom_dir' submodule path if AXOM_DATA_DIR not provided
    get_filename_component(_data_dir "${PROJECT_SOURCE_DIR}/../data" ABSOLUTE)
    if(EXISTS ${_data_dir}/README.md)
        set(AXOM_DATA_DIR ${_data_dir} CACHE PATH "")
    endif()
endif()

#------------------------------------------------------------------------------
# Include build system logic and options
#------------------------------------------------------------------------------
include(CMakeDependentOption)
include(cmake/CMakeBasics.cmake)

#------------------------------------------------------------------------------
# Global includes for the Axom projects (restrict these as much as possible)
#------------------------------------------------------------------------------
include_directories(${CMAKE_BINARY_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR})

#------------------------------------------------------------------------------
# Add source directories
#------------------------------------------------------------------------------
add_subdirectory(thirdparty)

include(cmake/AxomVersion.cmake)
add_subdirectory(axom)

if(AXOM_ENABLE_TOOLS)
    add_subdirectory(tools)
endif()

# Using Axom install examples
add_subdirectory(examples)

if(AXOM_ENABLE_DOCS)
    if(SPHINX_FOUND)
        blt_add_sphinx_target( axom_docs )
    endif()
    add_subdirectory(docs)
endif()

#------------------------------------------------------------------------------
# Generate header file with configuration options
#------------------------------------------------------------------------------
include(cmake/AxomConfig.cmake)
