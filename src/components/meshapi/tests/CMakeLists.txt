#
# Specify list of tests
#

set(gtest_meshapi_tests
   testHello.cxx
   )


#
# Add gtest based tests
#
foreach(test ${gtest_meshapi_tests})
    add_gtest( TEST_SOURCE ${test} DEPENDS_ON slic meshapi )
endforeach()


