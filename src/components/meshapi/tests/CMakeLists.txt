#
# Specify list of tests
#

set(gtest_meshapi_tests
   # test sets
   testSet.cpp
   testNullSet.cpp
   testRangeSet.cpp
   testWindowedRangeSet.cpp
   testIndirectionSet.cpp

   # test relations
   testDynamicVariableRelation.cpp
   testStaticVariableRelation.cpp
   testStaticConstantRelation.cpp
   
   # test maps
   testMap.cpp
   
   # test integration of sets, relations and maps
   testAccessingRelationDataInMap.cpp
   )


#
# Add gtest based tests
#
foreach(test ${gtest_meshapi_tests})
    add_gtest( TEST_SOURCE ${test} DEPENDS_ON meshapi slic)
endforeach()



