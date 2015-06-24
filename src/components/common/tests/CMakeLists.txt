#
# Specify list of tests
#

set(gtest_common_tests
   common_utilities.cpp
   )

#
# Add gtest based tests
#
foreach(test ${gtest_common_tests})
    add_gtest( TEST_SOURCE ${test} )
    get_filename_component(test_name ${arg_TEST_SOURCE} NAME_WE)
    add_dependencies(${test_name} copy_headers_common)
endforeach()

if(ENABLE_OMP)
    make_executable( EXECUTABLE_NAME openmpTester EXECUTABLE_SOURCE testOpenMP.cpp WITH_OPENMP)
endif()
