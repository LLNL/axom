#
# Specify list of tests
#

set(gtest_slam_tests
   # test sets
   slam_Set.cpp
   slam_NullSet.cpp
   slam_RangeSet.cpp
   slam_WindowedRangeSet.cpp
   slam_IndirectionSet.cpp

   # test relations
   slam_DynamicVariableRelation.cpp
   slam_StaticVariableRelation.cpp
   slam_StaticConstantRelation.cpp
   
   # test maps
   slam_Map.cpp
   
   # test integration of sets, relations and maps
   slam_AccessingRelationDataInMap.cpp

   # aux tests
   slam_ModularInt.cpp
   )


#
# Add gtest based tests
#
foreach(test ${gtest_slam_tests})
    add_gtest( TEST_SOURCE ${test} DEPENDS_ON slam slic)
endforeach()



