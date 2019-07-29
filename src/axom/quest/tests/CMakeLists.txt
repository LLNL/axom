# Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Quest tests
#------------------------------------------------------------------------------

set(quest_tests
    quest_all_nearest_neighbors.cpp
    quest_inout_octree.cpp
    quest_inout_quadtree.cpp
    quest_signed_distance.cpp
    quest_stl_reader.cpp
    quest_vertex_weld.cpp
   )

# Optionally, add tests that require AXOM_DATA_DIR
blt_list_append(TO       quest_tests
                IF       AXOM_DATA_DIR
                ELEMENTS quest_meshtester.cpp)

set(quest_tests_depends
    axom
    fmt
    gtest
   )

blt_list_append( TO quest_tests_depends ELEMENTS mpi IF ENABLE_MPI )
blt_list_append( TO quest_tests_depends ELEMENTS cuda IF ENABLE_CUDA )

foreach(test ${quest_tests})
    get_filename_component( test_name ${test} NAME_WE )
    blt_add_executable(
        NAME ${test_name}_test
        SOURCES ${test}
        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
        DEPENDS_ON ${quest_tests_depends}
        FOLDER axom/quest/tests
        )

    axom_add_test(
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
        DEPENDS_ON ${quest_tests_depends} mfem
        FOLDER axom/quest/tests
        )

    string(STRIP "${AXOM_DISABLE_UNUSED_PARAMETER_WARNINGS}" MFEM_COMPILE_FLAGS)

    if (ENABLE_CUDA)
        set(MFEM_COMPILE_FLAGS "-Xcompiler=${MFEM_COMPILE_FLAGS}")
    endif()

    blt_add_target_compile_flags( TO ${test_name} FLAGS "${MFEM_COMPILE_FLAGS}" )

    axom_add_test(
        NAME ${test_name}
        COMMAND ${test_name}
        )
endif()


#------------------------------------------------------------------------------
# Tests that use MPI when available
#------------------------------------------------------------------------------

set(quest_mpi_tests
    quest_initialize.cpp
    quest_signed_distance_interface.cpp
    )

# Optionally, add tests that require AXOM_DATA_DIR
blt_list_append(TO       quest_mpi_tests 
                IF       AXOM_DATA_DIR
                ELEMENTS quest_inout_interface.cpp )

foreach(test ${quest_mpi_tests})
    get_filename_component( test_name ${test} NAME_WE )
    blt_add_executable(
        NAME ${test_name}_test
        SOURCES ${test_name}.cpp
        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
        DEPENDS_ON ${quest_tests_depends}
        FOLDER axom/quest/tests
    )

    set(_numMPITasks 4)
    if (${test_name} STREQUAL "quest_signed_distance_interface")
        set(_numMPITasks 1)
    endif()
    if (${test_name} STREQUAL "quest_inout_interface")
        set(_numMPITasks 1)
    endif()

    if (ENABLE_MPI)
        axom_add_test(
            NAME ${test_name}
            COMMAND ${test_name}_test
            NUM_MPI_TASKS ${_numMPITasks}
        )
    else()
        axom_add_test(
            NAME ${test_name}
            COMMAND ${test_name}_test
        )
    endif()
endforeach()


#------------------------------------------------------------------------------
# Regression tests for quest signed distance and inout queries
#
# Note: By default we only run a small subset of the quest regression tests.
#       To run additional regression tests, configure the code
#       with AXOM_QUEST_ENABLE_EXTRA_REGRESSION_TESTS set to ON.
#------------------------------------------------------------------------------

cmake_dependent_option(
    AXOM_QUEST_ENABLE_EXTRA_REGRESSION_TESTS
    "Enable extra regression tests for quest" OFF
    "ENABLE_MPI;AXOM_ENABLE_SIDRE;AXOM_DATA_DIR" OFF )
mark_as_advanced(AXOM_QUEST_ENABLE_EXTRA_REGRESSION_TESTS)

if (ENABLE_MPI AND AXOM_ENABLE_SIDRE)

    set(quest_regression_depends
        ${quest_tests_depends}
        cli11
        )

    blt_list_append(TO quest_regression_depends ELEMENTS hdf5 IF HDF5_FOUND)
    blt_list_append(TO quest_regression_depends ELEMENTS scr IF SCR_FOUND)

    set(test_name quest_regression_test)

    blt_add_executable(
        NAME ${test_name}
        SOURCES quest_regression.cpp
        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
        DEPENDS_ON ${quest_regression_depends}
        FOLDER axom/quest/tests
        )

    if ( NOT ENABLE_OPENMP )
        set(extra_compile_flags "${AXOM_DISABLE_OMP_PRAGMA_WARNINGS}")

        if ( ENABLE_CUDA )
            string(STRIP "${extra_compile_flags}" extra_compile_flags)
            set(extra_compile_flags "-Xcompiler=${extra_compile_flags}")
        endif()

        blt_add_target_compile_flags( TO ${test_name} FLAGS "${extra_compile_flags}" )

    endif()

    # Add resolution tests for each dataset and resolution
    if(AXOM_DATA_DIR)
        set(quest_data_dir  ${AXOM_DATA_DIR}/quest)
        
        set(quest_regression_datasets
            sphere
        )
        set(quest_regression_resolutions
            16
        )

        blt_list_append(
            TO          quest_regression_datasets
            ELEMENTS    aatbase_3_binary
                        cap
                        naca0012
                        plane_binary
                        plane_simp
            IF          AXOM_QUEST_ENABLE_EXTRA_REGRESSION_TESTS
           )
        blt_list_append(
            TO          quest_regression_resolutions
            ELEMENTS    32
            IF          AXOM_QUEST_ENABLE_EXTRA_REGRESSION_TESTS
           )

        foreach(dataset ${quest_regression_datasets})
            foreach(res ${quest_regression_resolutions})
                axom_add_test(
                    NAME quest_regression_${dataset}_${res}
                    COMMAND ${test_name}
                            --mesh ${quest_data_dir}/${dataset}.stl
                            --baseline ${quest_data_dir}/regression/${dataset}_${res}_baseline.root
                    NUM_MPI_TASKS 1
                )
            endforeach()
        endforeach()
    endif()
endif()
