# Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Unit tests for Slic component
#------------------------------------------------------------------------------

#
# Specify list of tests
#

set(gtest_slic_tests
    slic_asserts.cpp
    slic_interface.cpp
    slic_macros.cpp
    )

set(bench_slic_tests
    slic_benchmark_asserts.cpp
   )


#
# Add gtest based tests
#
foreach(test ${gtest_slic_tests})
    get_filename_component( test_name ${test} NAME_WE )
    blt_add_executable( NAME ${test_name}_test
                        SOURCES ${test}
                        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                        DEPENDS_ON axom gtest
                        FOLDER     axom/slic/tests )
    axom_add_test( NAME ${test_name}
                   COMMAND ${test_name}_test )
endforeach()

blt_add_executable( NAME slic_fmt_test
                    SOURCES slic_fmt.cpp
                    OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                    DEPENDS_ON axom gtest fmt
                    FOLDER     axom/slic/tests )
axom_add_test( NAME slic_fmt
               COMMAND slic_fmt_test )



if (ENABLE_BENCHMARKS)
    foreach(test ${bench_slic_tests})
        get_filename_component( test_name ${test} NAME_WE )

        blt_add_executable(
            NAME ${test_name}
            SOURCES ${test}
            OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
            DEPENDS_ON axom gbenchmark
            FOLDER     axom/slic/benchmarks )

        blt_add_benchmark(
            NAME ${test_name}
            COMMAND ${test_name} --benchmark_min_time=0.0001
            )
    endforeach()
endif()
