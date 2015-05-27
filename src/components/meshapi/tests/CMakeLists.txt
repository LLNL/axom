#
# Specify list of tests
#

set(gtest_meshapi_tests
   # test sets
   testSet.cxx
   testNullSet.cpp
   testOrderedSet.cxx

   # test relations
   testDynamicVariableRelation.cxx
   testStaticVariableRelation.cxx
   testStaticConstantRelation.cxx
   
   # test maps
   testMap.cxx
   )


#
# Add gtest based tests
#
foreach(test ${gtest_meshapi_tests})
    add_gtest( TEST_SOURCE ${test} DEPENDS_ON meshapi)
endforeach()



