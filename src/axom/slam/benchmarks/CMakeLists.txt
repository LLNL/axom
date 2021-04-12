# Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# C++ Benchmarks for Slam component
#------------------------------------------------------------------------------

set(slam_benchmark_files
    slam_array.cpp
    slam_sets.cpp
    )

if (ENABLE_BENCHMARKS)
    foreach(test ${slam_benchmark_files})
        get_filename_component( test_name ${test} NAME_WE )
        set(test_name "${test_name}_benchmark")
        
        blt_add_executable(
            NAME        ${test_name}
            SOURCES     ${test}
            OUTPUT_DIR  ${TEST_OUTPUT_DIRECTORY}
            DEPENDS_ON  slic slam gbenchmark
            FOLDER      axom/slam/benchmarks
            )

        blt_add_benchmark( 
            NAME        ${test_name} 
            COMMAND     ${test_name}
            )
    endforeach()
endif()
