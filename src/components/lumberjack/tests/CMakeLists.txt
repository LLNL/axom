#-------------------------------------------------------------------------------
# Copyright (c) 2015, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
#
# All rights reserved.
#
# This source code cannot be distributed without permission and further
# review from Lawrence Livermore National Laboratory.
#-------------------------------------------------------------------------------

#
# Specify list of tests
#

set(gtest_lumberjack_tests
    lumberjack_Lumberjack.cpp
    lumberjack_Message.cpp
    lumberjack_TextEqualityCombiner.cpp
    lumberjack_RootCommunicator.cpp
    )

#
# Add gtest based tests
#
foreach(test ${gtest_lumberjack_tests})
    get_filename_component( test_name ${test} NAME_WE )
    blt_add_executable( NAME ${test_name}
                        SOURCES ${test}
                        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                        DEPENDS_ON lumberjack gtest )
    blt_add_test( NAME ${test_name} 
                  COMMAND ${test_name} )
endforeach()

blt_add_executable(NAME lumberjack_speedTest
                   SOURCES lumberjack_speedTest.cpp
                   OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                   DEPENDS_ON lumberjack )
blt_add_test( NAME lumberjack_speedTest_binary
              COMMAND lumberjack_speedTest b 10 ${lumberjack_SOURCE_DIR}/tests/sampleInput/loremIpsum02)
blt_add_test( NAME lumberjack_speedTest_root
              COMMAND lumberjack_speedTest r 10 ${lumberjack_SOURCE_DIR}/tests/sampleInput/loremIpsum02)
