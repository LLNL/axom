#
# Specify list of tests
#

set(gtest_lumberjack_tests
    lumberjack_Logger.cpp
    lumberjack_MessageInfo.cpp
    lumberjack_MessageEqualityCombiner.cpp
    lumberjack_RootCommunicator.cpp
    )

#
# Add gtest based tests
#
foreach(test ${gtest_lumberjack_tests})
    add_gtest(TEST_SOURCE ${test} DEPENDS_ON lumberjack)
endforeach()
