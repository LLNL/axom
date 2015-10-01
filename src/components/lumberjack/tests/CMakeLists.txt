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
    add_gtest(TEST_SOURCE ${test} DEPENDS_ON lumberjack)
endforeach()
