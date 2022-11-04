# Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Unit tests for Slic component
#------------------------------------------------------------------------------

#
# Serial Unit Tests
#
set(serial_slic_tests
    slic_asserts.cpp
    slic_fmt.cpp
    slic_interface.cpp
    slic_macros.cpp
    slic_uninit.cpp )

set(slic_device_depends "")
blt_list_append( TO slic_device_depends ELEMENTS cuda IF ${ENABLE_CUDA} )

foreach(test ${serial_slic_tests})
    get_filename_component( test_name ${test} NAME_WE )
    blt_add_executable(NAME       ${test_name}_test
                       SOURCES    ${test}
                       OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                       DEPENDS_ON axom gtest fmt ${slic_device_depends}
                       FOLDER     axom/slic/tests )
    axom_add_test(NAME    ${test_name}
                  COMMAND ${test_name}_test )
endforeach()


#
# Parallel Unit Tests
#
if (ENABLE_MPI AND AXOM_ENABLE_LUMBERJACK)
    set(parallel_slic_tests
        slic_macros_parallel.cpp )

    foreach(test ${parallel_slic_tests})
        get_filename_component( test_name ${test} NAME_WE )
        blt_add_executable(NAME       ${test_name}_test
                           SOURCES    ${test}
                           OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                           DEPENDS_ON axom gtest fmt ${slic_device_depends}
                           FOLDER     axom/slic/tests )
        axom_add_test(NAME          ${test_name}
                      COMMAND       ${test_name}_test 
                      NUM_MPI_TASKS 4)
    endforeach()
endif()


#
# Benchmarks
#
if (ENABLE_BENCHMARKS)
    set(slic_benchmarks
        slic_benchmark_asserts.cpp )

    foreach(test ${slic_benchmarks})
        get_filename_component( test_name ${test} NAME_WE )

        blt_add_executable(NAME       ${test_name}
                           SOURCES    ${test}
                           OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                           DEPENDS_ON axom gbenchmark ${slic_device_depends}
                           FOLDER     axom/slic/benchmarks )

        blt_add_benchmark(
            NAME    ${test_name}
            COMMAND ${test_name} --benchmark_min_time=0.0001 )
    endforeach()
endif()


#
# Fortran tests
#
if (ENABLE_FORTRAN)
    set(fruit_slic_tests
        slic_interface_F.f )

    foreach(test ${fruit_slic_tests})
        get_filename_component( test_name ${test} NAME_WE )
        blt_add_executable(NAME ${test_name}_test
                           SOURCES ${test}
                           OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                           DEPENDS_ON slic fruit ${slic_device_depends}
                           FOLDER axom/slic/tests )
        axom_add_test(NAME    ${test_name}
                      COMMAND ${test_name}_test )
    endforeach()
endif()
