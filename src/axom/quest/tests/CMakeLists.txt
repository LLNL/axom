#------------------------------------------------------------------------------
# Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
#
# Produced at the Lawrence Livermore National Laboratory
#
# LLNL-CODE-741217
#
# All rights reserved.
#
# This file is part of Axom.
#
# For details about use and distribution, please read axom/LICENSE.
#------------------------------------------------------------------------------

set(quest_tests
    quest_all_nearest_neighbors.cpp
    quest_implicit_grid.cpp
    quest_inout_octree.cpp
    quest_meshtester.cpp
    quest_octree.cpp
    quest_signed_distance.cpp
    quest_signed_distance_interface.cpp
    quest_spatial_octree.cpp
    quest_stl_reader.cpp
    quest_vertex_weld.cpp
   )

set(quest_tests_depends
    quest
    slic
    fmt
    gtest
   )

if (ENABLE_OPENMP)
    list(APPEND quest_tests_depends openmp)
endif()

set(quest_vertex_weld_extra_depends mint)


foreach(test ${quest_tests})
    get_filename_component( test_name ${test} NAME_WE )
    blt_add_executable(
        NAME ${test_name}_test
        SOURCES ${test}
        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
        DEPENDS_ON ${quest_tests_depends} ${${test_name}_extra_depends}
        FOLDER axom/quest/tests
        )

    blt_add_test(
        NAME ${test_name}
        COMMAND ${test_name}_test
        )

endforeach()


if(MFEM_FOUND)
    set(test_name quest_point_in_cell_mfem)
    blt_add_executable(
        NAME ${test_name}
        SOURCES ${test_name}.cpp
        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
        DEPENDS_ON ${quest_tests_depends} mfem core
        )

    set_property(
        TARGET ${test_name}
        APPEND_STRING PROPERTY COMPILE_FLAGS "${AXOM_DISABLE_UNUSED_PARAMETER_WARNINGS}"
        )

    blt_add_test(
        NAME ${test_name}
        COMMAND ${test_name}
        )
endif()

# Test quest::initialize() passing a Mesh * instead of a file name
set(test_name quest_initialize)
set(quest_initialize_depends ${quest_tests_depends})
if (ENABLE_MPI)
    list(APPEND quest_initialize_depends mpi)
endif()
blt_add_executable(
    NAME ${test_name}_test
    SOURCES ${test_name}.cpp
    OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
    DEPENDS_ON ${quest_initialize_depends}
    FOLDER axom/quest/tests
)
if (ENABLE_MPI)
    blt_add_test(
        NAME ${test_name}
        COMMAND ${test_name}_test
        NUM_MPI_TASKS 4
    )
else()
    blt_add_test(
        NAME ${test_name}
        COMMAND ${test_name}_test
    )
endif()

# Regression tests for quest signed distance and InOutOctree point containment queries
if (ENABLE_MPI AND AXOM_ENABLE_SIDRE)

    set(quest_regression_depends
        ${quest_tests_depends}
        conduit
        sidre
        )

    if (HDF5_FOUND)
       list(APPEND quest_regression_depends hdf5)
    endif()

    if (SCR_FOUND)
       list(APPEND quest_regression_depends scr)
    endif()

    set(test_name quest_regression_test)

    blt_add_executable(
        NAME ${test_name}
        SOURCES quest_regression.cpp
        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
        DEPENDS_ON ${quest_regression_depends}
        FOLDER axom/quest/tests
        )

    set(quest_data_dir  ${CMAKE_CURRENT_SOURCE_DIR}/../data)

    if ( NOT ENABLE_OPENMP )
        set_property(
            TARGET ${test_name}
            APPEND_STRING PROPERTY COMPILE_FLAGS "${AXOM_DISABLE_OMP_PRAGMA_WARNINGS}" )
    endif()

    if ( "$ENV{SYS_TYPE}" STREQUAL "bgqos_0" )
       set_property(
           TARGET ${test_name}
           APPEND_STRING PROPERTY LINK_FLAGS ${AXOM_ALLOW_MULTIPLE_DEFINITIONS}
       )

    endif()


    # Note: By default we only run a small subset of the quest regression tests.
    #       To run additional regression tests, configure the code
    #       with ENABLE_QUEST_EXTRA_REGRESSION_TESTS set to True.

    # set(ENABLE_QUEST_EXTRA_REGRESSION_TESTS True)

    set(quest_regression_datasets
        sphere
    )
    set(quest_regression_resolutions
        16
    )

    if( ${ENABLE_QUEST_EXTRA_REGRESSION_TESTS} )
       list(APPEND quest_regression_datasets
           naca0012
           plane_simp
           plane_binary
           aatbase_3_binary
       )
       list(APPEND quest_regression_resolutions
           32
       )
    endif()

    # Add resolution tests for each dataset and resolution
    foreach(dataset ${quest_regression_datasets})
        foreach(res ${quest_regression_resolutions})
            blt_add_test(
                NAME quest_regression_${dataset}_${res}
                COMMAND ${test_name}
                        --mesh ${quest_data_dir}/${dataset}.stl
                        --baseline ${quest_data_dir}/${dataset}_${res}_baseline.root
                        --resolution ${res} ${res} ${res}
                NUM_MPI_TASKS 1
            )
        endforeach()
    endforeach()
endif()
