#
# Specify list of tests
#

set(gtest_common_tests
   common_fileUtilities.cpp
   common_types.cpp
   )


#
# Add gtest based tests
#
foreach(test ${gtest_common_tests})
    add_gtest( TEST_SOURCE ${test} DEPENDS_ON common slic)
endforeach()



