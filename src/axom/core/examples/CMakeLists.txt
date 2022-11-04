# Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Core containers, numerics, and utilities examples
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# List of single source file examples
#------------------------------------------------------------------------------
set(example_sources
    core_numerics.cpp
    core_containers.cpp 
    core_acceleration.cpp)

set( core_examples_depends axom )
blt_list_append( TO core_examples_depends ELEMENTS blt::hip IF ${ENABLE_HIP} )
blt_list_append( TO core_examples_depends ELEMENTS cuda IF ${ENABLE_CUDA} )


#------------------------------------------------------------------------------
# Add targets and tests for serial examples
#------------------------------------------------------------------------------
foreach(example_source ${example_sources})
    get_filename_component(exe_name ${example_source} NAME_WE)
    blt_add_executable(
        NAME        ${exe_name}_ex
        SOURCES     ${example_source}
        OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
        DEPENDS_ON  ${core_examples_depends}
        FOLDER      axom/core/examples )

    if(AXOM_ENABLE_TESTS)
        axom_add_test(
            NAME    ${exe_name}
            COMMAND ${exe_name}_ex )
    endif()
endforeach()

if(ENABLE_CUDA)
    blt_add_executable(
          NAME        core_acceleration_cuda_on_ex
          SOURCES     core_acceleration.cpp
          OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
          DEPENDS_ON  ${core_examples_depends} cuda
          FOLDER      axom/core/examples )

    blt_add_executable(
          NAME        core_containers_cuda_on_ex
          SOURCES     core_containers.cpp
          OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
          DEPENDS_ON  axom cuda
          FOLDER      axom/core/examples )
endif()

blt_add_executable(
    NAME        core_utilities_ex
    SOURCES     core_utilities.cpp
    OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
    DEPENDS_ON  ${core_examples_depends}
    FOLDER      axom/core/examples )

if(AXOM_ENABLE_TESTS)
    axom_add_test(
        NAME    core_utilities
        COMMAND core_utilities_ex ${CMAKE_CURRENT_SOURCE_DIR} )
endif()
