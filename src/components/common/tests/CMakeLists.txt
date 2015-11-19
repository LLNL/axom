#
# Specify list of tests
#

set(gtest_common_tests
    gtest_smoke.cpp)

#
# Add gtest based tests
#
foreach(test ${test_common_tests})
    add_gtest( TEST_SOURCE ${test})
endforeach()

if(ENABLE_OMP)
    make_executable( EXECUTABLE_NAME openmpTester EXECUTABLE_SOURCE testOpenMP.cpp WITH_OPENMP ADD_CTEST)
endif()
