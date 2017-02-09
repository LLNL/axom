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
    quest_intersection_impl.cpp
    quest_morton.cpp
    quest_numeric_array.cpp
    quest_octree.cpp
    quest_orient.cpp
    quest_point.cpp
    quest_signed_distance.cpp
    quest_spatial_octree.cpp
    quest_squared_distance.cpp
    quest_triangle.cpp
    quest_vector.cpp
   )

set(quest_tests_depends
    quest
    slic
    fmt
    gtest
   )

if (BOOST_FOUND)
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


# Regression tests for quest signed distance and InOutOctree point containment queries
if (ENABLE_MPI AND ENABLE_SIDRE)

    blt_add_executable(
        NAME quest_regression
        SOURCES quest_regression.cpp
        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
        DEPENDS_ON ${quest_tests_depends} conduit sidre spio hdf5 sparsehash
    )
    set(quest_data_dir  ${CMAKE_CURRENT_SOURCE_DIR}/../data)
    
    # Note: The following are all valid regression tests.
    #       Most are commented out to avoid excessive test runtimes, especially in debug mode.
    set(quest_regression_datasets
        sphere
        #naca0012
        #plane_simp
        #plane_binary
    )
    set(quest_regression_resolutions
        16
        #32
    )
    
    # Add resolution tests for each dataset and resolution
    foreach(dataset ${quest_regression_datasets})
        foreach(res ${quest_regression_resolutions})
            blt_add_test(
                NAME quest_regression_${dataset}_${res}
                COMMAND quest_regression 
                        --mesh ${quest_data_dir}/${dataset}.stl 
                        --baseline ${quest_data_dir}/${dataset}_${res}_baseline.root
            )
        endforeach()
    endforeach()
endif()
