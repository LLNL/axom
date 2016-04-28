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

    # testing utility functions
    slam_utilities.cpp
    )


#
# Add gtest based tests
#
foreach(test ${gtest_slam_tests})
    get_filename_component( test_name ${test} NAME_WE )
    blt_add_executable( NAME ${test_name}
                        SOURCES ${test}
                        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                        DEPENDS_ON slam slic gtest )
    blt_add_test( NAME ${test_name} 
                  COMMAND ${test_name} )
endforeach()
