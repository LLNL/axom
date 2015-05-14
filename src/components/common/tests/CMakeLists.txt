#
# Specify list of tests
#

set(gtest_common_tests
   common_utilities.cxx
   )

#
# Add gtest based tests
#
foreach(test ${gtest_datastore_tests})
    add_gtest( TEST_SOURCE ${test} )
endforeach()
