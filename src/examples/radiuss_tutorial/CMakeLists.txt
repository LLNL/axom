# Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Tutorial that demonstrates the usage of several Axom components
# in a BLT-based build system.
#
# Configuration variables are stored in a CMake cache file 'host-config.cmake'
# which defines paths to blt, axom and possibly other TPLs.
# It also contains information about the compiler used to build axom.
#------------------------------------------------------------------------------
#
# To build:
#  mkdir build
#  cd build
#  cmake -C ../host-config.cmake ..
#  make
#  ./bin/example
#
#------------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.14)

project(radiuss_tutorial)

#------------------------------------------------------------------------------
# Set up BLT with validity checks
#------------------------------------------------------------------------------

# Check that path to BLT is provided and valid
if(NOT DEFINED BLT_SOURCE_DIR OR NOT EXISTS ${BLT_SOURCE_DIR}/SetupBLT.cmake)
    message(FATAL_ERROR "Missing required 'BLT_SOURCE_DIR' variable pointing to a valid blt")
endif()
include(${BLT_SOURCE_DIR}/SetupBLT.cmake)


#------------------------------------------------------------------------------
# Check for AXOM_DIR and use CMake's find_package to import axom's targets
#------------------------------------------------------------------------------
if(NOT DEFINED AXOM_DIR OR NOT EXISTS ${AXOM_DIR}/lib/cmake/axom-config.cmake)
    message(FATAL_ERROR "Missing required 'AXOM_DIR' variable pointing to an installed axom")
endif()

include(CMakeFindDependencyMacro)

find_dependency(axom REQUIRED
                NO_DEFAULT_PATH 
                PATHS ${AXOM_DIR}/lib/cmake)

#------------------------------------------------------------------------------
# Set up target executables for the tutorial
#------------------------------------------------------------------------------

set(tutorial_deps axom axom::fmt axom::cli11)
blt_list_append(TO tutorial_deps ELEMENTS cuda IF CUDA_FOUND)

blt_add_executable(NAME       lesson_00_check_axom_configuration 
                   SOURCES    lesson_00/check_axom_configuration.cpp
                   DEPENDS_ON ${tutorial_deps})

blt_add_executable(NAME       lesson_01_load_stl_mesh 
                   SOURCES    lesson_01/load_stl_mesh.cpp
                   DEPENDS_ON ${tutorial_deps})

blt_add_executable(NAME       lesson_02_naive_self_intersections
                   SOURCES    lesson_02/naive_self_intersections.cpp
                   DEPENDS_ON ${tutorial_deps})

if(umpire_FOUND AND RAJA_FOUND)
    blt_list_append(TO tutorial_deps ELEMENTS RAJA IF RAJA_FOUND)
    blt_list_append(TO tutorial_deps ELEMENTS umpire IF umpire_FOUND)

    blt_add_executable(NAME       lesson_03_device_self_intersections
                       SOURCES    lesson_03/device_self_intersections.cpp
                       DEPENDS_ON ${tutorial_deps})

    blt_add_executable(NAME       lesson_04_device_spatial_indexes
                       SOURCES    lesson_04/device_spatial_indexes.cpp
                       DEPENDS_ON ${tutorial_deps})
endif()

#------------------------------------------------------------------------------
# Optionally, print out information about imported targets
#------------------------------------------------------------------------------
if(EXAMPLE_VERBOSE_OUTPUT)
    foreach(_target axom axom::cli11 axom::fmt axom::sparsehash conduit::conduit umpire RAJA)
        blt_print_target_properties(TARGET ${_target})  
    endforeach()
endif()

