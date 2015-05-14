#
# Specify list of tests
#

set(gtest_meshapi_tests
   testOrderedSet.cxx
   testDynamicVariableRelation.cxx
   testStaticVariableRelation.cxx
   testStaticConstantRelation.cxx
   testSet.cxx
   )


#
# Add gtest based tests
#
foreach(test ${gtest_meshapi_tests})
    add_gtest( TEST_SOURCE ${test} DEPENDS_ON meshapi)
endforeach()



