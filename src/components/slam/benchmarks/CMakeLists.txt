#------------------------------------------------------------------------------
# Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
#
# Produced at the Lawrence Livermore National Laboratory
#
# LLNL-CODE-741217
#
# All rights reserved.
#
# This file is part of Axom.
#
# For details about use and distribution, please read axom/LICENSE.
#------------------------------------------------------------------------------

################
# C++ Benchmarks
################
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
            DEPENDS_ON  slic slam gbenchmark boost
            FOLDER      axom/slam/benchmarks
            )

        blt_add_benchmark( 
            NAME        ${test_name} 
            COMMAND     ${test_name}
            )
    endforeach()
endif()
