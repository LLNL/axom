# Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Example that shows how to use an installed instance of Axom 
# in a CMake-based build system.
#
# Configuration variables are stored in a CMake cache file 'host-config.cmake'
# which defines paths to axom and possibly other TPLs.
# It also contains information about the compiler used to build axom.
#------------------------------------------------------------------------------
#
# To build:
#  mkdir build
#  cd build
#  cmake -C ../host-config.cmake ..
#  make
#  ./example
#
#------------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.14)

project(using_with_cmake)

message(STATUS "CMake Version: ${CMAKE_VERSION}")

#------------------------------------------------------------------------------
# Check for AXOM_DIR and use CMake's find_package to import axom's targets
#------------------------------------------------------------------------------
if(NOT DEFINED AXOM_DIR OR NOT EXISTS ${AXOM_DIR}/lib/cmake/axom-config.cmake)
    message(FATAL_ERROR "Missing required 'AXOM_DIR' variable pointing to an installed axom")
endif()

if (ENABLE_CUDA)
    enable_language(CUDA)
endif()

if (ENABLE_HIP)
    if (NOT ROCM_PATH)
        find_path(ROCM_PATH
            hip
            ENV{ROCM_DIR}
            ENV{ROCM_PATH}
            ENV{HIP_PATH}
            ${HIP_PATH}/..
            ${HIP_ROOT_DIR}/../
            ${ROCM_ROOT_DIR}
            /opt/rocm)
    endif()
    set(CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH};${ROCM_PATH}")
    find_package(hip REQUIRED CONFIG PATHS ${ROCM_PATH})
endif()

# Required to create Conduit's MPI dependencies
if (ENABLE_MPI)
    find_package(MPI REQUIRED)
    
    if(NOT TARGET mpi AND TARGET MPI::MPI_CXX)
        add_library(mpi INTERFACE)
        target_link_libraries(mpi INTERFACE MPI::MPI_CXX)
    endif()
endif()

find_package(axom REQUIRED
             NO_DEFAULT_PATH 
             PATHS ${AXOM_DIR}/lib/cmake)

# BEGIN FIXME:
# Create fake empty target, this stops CMake from adding -lcuda_runtime to the link line
add_library(cuda_runtime INTERFACE)

# Create fake empty target, this stops CMake from adding -lblt_hip_runtime to the link line
add_library(blt_hip_runtime INTERFACE)
# END FIXME

# Remove implicitly added link directories added by CMake that are problematic when
# the default system libraries are older than the ones used by the compiler
if(BLT_CMAKE_IMPLICIT_LINK_DIRECTORIES_EXCLUDE)
    list(REMOVE_ITEM CMAKE_C_IMPLICIT_LINK_DIRECTORIES
                     ${BLT_CMAKE_IMPLICIT_LINK_DIRECTORIES_EXCLUDE})
    list(REMOVE_ITEM CMAKE_CXX_IMPLICIT_LINK_DIRECTORIES
                     ${BLT_CMAKE_IMPLICIT_LINK_DIRECTORIES_EXCLUDE})
    list(REMOVE_ITEM CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES
                     ${BLT_CMAKE_IMPLICIT_LINK_DIRECTORIES_EXCLUDE})
endif()

#------------------------------------------------------------------------------
# Set up example target that depends on axom
#------------------------------------------------------------------------------
add_executable(example example.cpp)

# setup the axom include path
target_include_directories(example PRIVATE ${AXOM_INCLUDE_DIRS})

# link to axom targets
target_link_libraries(example axom)
target_link_libraries(example axom::fmt)
