# Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

#
# Specify list of tests
#

# Add Serial GTests based tests
set(gtest_lumberjack_tests
    lumberjack_Lumberjack.cpp
    lumberjack_Message.cpp
    lumberjack_TextEqualityCombiner.cpp
    )

foreach(test ${gtest_lumberjack_tests})
    get_filename_component( test_name ${test} NAME_WE )
    blt_add_executable( NAME ${test_name}_test
                        SOURCES ${test}
                        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                        DEPENDS_ON lumberjack gtest 
                        FOLDER axom/lumberjack/tests )
    blt_add_test( NAME ${test_name} 
                  COMMAND ${test_name}_test )
endforeach()


# Add MPI GTest based tests
set(mpi_gtest_lumberjack_tests
    lumberjack_BinaryCommunicator.cpp
    lumberjack_RootCommunicator.cpp
    )

foreach(test ${mpi_gtest_lumberjack_tests})
    get_filename_component( test_name ${test} NAME_WE )
    blt_add_executable( NAME ${test_name}_test
                        SOURCES ${test}
                        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                        DEPENDS_ON lumberjack gtest 
                        FOLDER axom/lumberjack/tests )
    foreach(task_count RANGE 1 8)
        blt_add_test( NAME ${test_name}_${task_count}
                      COMMAND ${test_name}_test
                      NUM_MPI_TASKS ${task_count})
    endforeach()
endforeach()

#Unique Tests
blt_add_executable(NAME lumberjack_speed_test
                   SOURCES lumberjack_speedTest.cpp
                   OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                   DEPENDS_ON lumberjack
                   FOLDER axom/lumberjack/tests )

set(lumberjack_sample_input_dir ${PROJECT_SOURCE_DIR}/axom/lumberjack/tests/sampleInput)
blt_add_test( NAME lumberjack_speedTest_binary
              COMMAND lumberjack_speed_test b 10 ${lumberjack_sample_input_dir}/loremIpsum02
              NUM_MPI_TASKS 4)
              
blt_add_test( NAME lumberjack_speedTest_root
              COMMAND lumberjack_speed_test r 10 ${lumberjack_sample_input_dir}/loremIpsum02
              NUM_MPI_TASKS 4)
              
