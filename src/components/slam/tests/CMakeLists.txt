#
# Specify list of tests
#

set(gtest_slam_tests
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

   # aux tests
   testModularInt.cpp
   )


#
# Add gtest based tests
#
foreach(test ${gtest_slam_tests})
    add_gtest( TEST_SOURCE ${test} DEPENDS_ON slam slic)
endforeach()



