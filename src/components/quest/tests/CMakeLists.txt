#-------------------------------------------------------------------------------
# Copyright (c) 2015, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
#
# All rights reserved.
#
# This source code cannot be distributed without permission and further
# review from Lawrence Livermore National Laboratory.
#-------------------------------------------------------------------------------

set(quest_tests
    quest_boundingBox.cpp
    quest_bucket_tree.cpp
    quest_inout_octree.cpp
    quest_intersection.cpp
    quest_mesh_extent.cpp
    quest_morton.cpp
    quest_numeric_array.cpp
    quest_octree.cpp
    quest_orient.cpp
    quest_point.cpp
    quest_signed_distance.cpp
    quest_spatial_octree.cpp
    quest_squared_distance.cpp
    quest_vector.cpp
   )

set(quest_tests_depends
    quest
    slic
    gtest
   )

if (ENABLE_BOOST)
    list(APPEND quest_tests_depends boost)
endif()

foreach(test ${quest_tests})
    get_filename_component( test_name ${test} NAME_WE )
    blt_add_executable(
        NAME ${test_name}
        SOURCES ${test}
        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
        DEPENDS_ON ${quest_tests_depends}
    )
    blt_add_test(
        NAME ${test_name}
        COMMAND ${test_name}
    )
endforeach()

