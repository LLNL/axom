#
# Specify list of tests
#

set(gtest_slam_tests
    # test sets
    slam_set_Set.cpp
    slam_set_NullSet.cpp
    slam_set_PositionSet.cpp
    slam_set_RangeSet.cpp
    slam_set_IndirectionSet.cpp

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


# blt_add_executable depends_on handles our include dir settings
# these teste are using boost so boost needs to be a dep
set(slam_tests_depends_on axom_utils slic slam gtest)
if (BOOST_FOUND)
    list(APPEND slam_tests_depends_on boost)
endif()


#
# Add gtest based tests
#
foreach(test ${gtest_slam_tests})
    get_filename_component( test_name ${test} NAME_WE )
    blt_add_executable( NAME ${test_name}
                        SOURCES ${test}
                        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                        DEPENDS_ON ${slam_tests_depends_on} )
    blt_add_test( NAME ${test_name}
                  COMMAND ${test_name} )
endforeach()

# The slam utilities test requires a directory
blt_add_executable(
    NAME    slam_utilities
    SOURCES slam_utilities.cpp
    OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
    DEPENDS_ON ${slam_tests_depends_on}
    )
blt_add_test(
    NAME    slam_utilities
    COMMAND slam_utilities ${CMAKE_CURRENT_SOURCE_DIR}
    )
