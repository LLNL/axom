#
# Specify list of tests
#

set(gtest_meshapi_tests
   # test sets
   testSet.cpp
   testNullSet.cpp
   testOrderedSet.cpp

   # test relations
   testDynamicVariableRelation.cpp
   testStaticVariableRelation.cpp
   testStaticConstantRelation.cpp
   
   # test maps
   testMap.cpp
   )


#
# Add gtest based tests
#
foreach(test ${gtest_meshapi_tests})
    add_gtest( TEST_SOURCE ${test} DEPENDS_ON meshapi)
endforeach()



