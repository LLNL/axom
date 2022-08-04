# Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Quest examples
#------------------------------------------------------------------------------

set(quest_example_depends axom fmt cli11)

blt_list_append(TO quest_example_depends ELEMENTS mpi IF ENABLE_MPI)
blt_list_append(TO quest_example_depends ELEMENTS openmp IF ENABLE_OPENMP)
blt_list_append(TO quest_example_depends ELEMENTS cuda IF ENABLE_CUDA)
blt_list_append(TO quest_example_depends ELEMENTS blt::hip IF ENABLE_HIP)
blt_list_append(TO quest_example_depends ELEMENTS RAJA IF RAJA_FOUND)

set(quest_data_dir  ${AXOM_DATA_DIR}/quest)

# In/out octree containment example -------------------------------------------
blt_add_executable(
    NAME        quest_containment_driver_ex
    SOURCES     containment_driver.cpp
    OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
    DEPENDS_ON  ${quest_example_depends}
    FOLDER      axom/quest/examples
    )

# Shaping example -------------------------------------------------------------
if(ENABLE_MPI AND MFEM_FOUND AND MFEM_USE_MPI
              AND AXOM_ENABLE_SIDRE AND AXOM_ENABLE_MFEM_SIDRE_DATACOLLECTION
              AND AXOM_ENABLE_KLEE)
    blt_add_executable(
        NAME        quest_shaping_driver_ex
        SOURCES     shaping_driver.cpp
        OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
        DEPENDS_ON  ${quest_example_depends} mfem mpi
        FOLDER      axom/quest/examples
        )
endif()

# Distributed closest point example -------------------------------------------
if(ENABLE_MPI AND MFEM_FOUND AND MFEM_USE_MPI
              AND AXOM_ENABLE_SIDRE AND AXOM_ENABLE_MFEM_SIDRE_DATACOLLECTION)
    blt_add_executable(
            NAME        quest_distributed_distance_query_ex
            SOURCES     quest_distributed_distance_query_example.cpp
            OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
            DEPENDS_ON  ${quest_example_depends} mfem mpi
            FOLDER      axom/quest/examples
            )

    if(AXOM_ENABLE_TESTS AND AXOM_DATA_DIR)
        set(_nranks 3)

        # Bugfix: Serialize ordering of DCP tests as workaround for MPI issues 
        # encountered when these tests ran in parallel
        set(_dcp_previous_test)

        # Run the distributed closest point example on N ranks for each enabled policy
        set(_policies "seq")
        blt_list_append(TO _policies ELEMENTS "omp" IF ENABLE_OPENMP)
        blt_list_append(TO _policies ELEMENTS "cuda" IF ENABLE_CUDA)
        blt_list_append(TO _policies ELEMENTS "hip" IF ENABLE_HIP)
        foreach(_pol ${_policies})
            # Add test with data on all ranks
            set(_test "quest_distributed_closest_point_run_2D_${_pol}")
            axom_add_test(
                NAME    ${_test}
                COMMAND quest_distributed_distance_query_ex
                            --mesh-file ${quest_data_dir}/box_2D_r${_nranks}.root
                            --num-samples 500
                            --dist-threshold .45
                            --policy ${_pol}
                            --object-file dcp_object_mesh_2d_${_pol}
                            --distance-file dcp_closest_point_2d_${_pol}
                NUM_MPI_TASKS ${_nranks})
            if(${_pol} STREQUAL "omp")
                set_property(TEST ${_test} APPEND PROPERTY ENVIRONMENT OMP_NUM_THREADS=4)
            endif()

            if(_dcp_previous_test)
                set_tests_properties(${_test} PROPERTIES DEPENDS ${_dcp_previous_test})
            endif()
            set(_dcp_previous_test ${_test})

            # Add test with data on a subset of the ranks
            set(_test "quest_distributed_closest_point_run_2D_empty_ranks_${_pol}")
            axom_add_test(
                NAME    ${_test}
                COMMAND quest_distributed_distance_query_ex
                            --mesh-file ${quest_data_dir}/box_2D_r${_nranks}.root
                            --num-samples 500
                            --dist-threshold .45
                            --policy ${_pol}
                            --object-file dcp_object_mesh_2d_${_pol}
                            --distance-file dcp_closest_point_2d_${_pol}
                            --empty-rank-probability .5
                NUM_MPI_TASKS ${_nranks})
            if(${_pol} STREQUAL "omp")
                set_property(TEST ${_test} APPEND PROPERTY ENVIRONMENT OMP_NUM_THREADS=4)
            endif()

            set_tests_properties(${_test} PROPERTIES DEPENDS ${_dcp_previous_test})
            set(_dcp_previous_test ${_test})
    
        endforeach()

        unset(_nranks)
        unset(_policies)
        unset(_test)
        unset(_dcp_previous_test)
    endif()
endif()

# Point in cell example -------------------------------------------------------
if(MFEM_FOUND)
    blt_add_executable(
        NAME        quest_point_in_cell_benchmark_ex
        SOURCES     point_in_cell_benchmark.cpp
        OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
        DEPENDS_ON  ${quest_example_depends} mfem
        FOLDER      axom/quest/examples
        )
endif()

# Delaunay triangulation example ----------------------------------------------
blt_add_executable(
    NAME        quest_delaunay_triangulation_ex
    SOURCES     delaunay_triangulation.cpp
    OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
    DEPENDS_ON  ${quest_example_depends}
    FOLDER      axom/quest/examples
    )

if(AXOM_ENABLE_TESTS)
    foreach(d 2 3)
        axom_add_test(
            NAME quest_delaunay_${d}_test
            COMMAND quest_delaunay_triangulation_ex -d ${d} -o delaunay_mesh_${d}d -n 100
            )
    endforeach()
endif()

# Scattered interpolation example --------------------------------------------
if(AXOM_ENABLE_SIDRE)
    blt_add_executable(
        NAME        quest_scattered_interpolation_ex
        SOURCES     scattered_interpolation.cpp
        OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
        DEPENDS_ON  ${quest_example_depends}
        FOLDER      axom/quest/examples
        )

    if(AXOM_ENABLE_TESTS)
        foreach(d 2 3)
            axom_add_test(
                NAME quest_scattered_interpolation_${d}_test
                COMMAND quest_scattered_interpolation_ex -d ${d} -n 10000 -q 20000
                )
        endforeach()
    endif()
endif()


# Quest signed distance and inout interface examples (C++ and Fortran) --------

blt_add_executable(
    NAME       quest_signed_distance_interface_ex
    SOURCES    quest_signed_distance_interface.cpp
    OUTPUT_DIR ${EXAMPLE_OUTPUT_DIRECTORY}
    DEPENDS_ON ${quest_example_depends}
    FOLDER      axom/quest/examples
    )

blt_add_executable(
    NAME        quest_inout_interface_ex
    SOURCES     quest_inout_interface.cpp
    OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
    DEPENDS_ON  ${quest_example_depends}
    FOLDER      axom/quest/examples
    )

if(AXOM_ENABLE_TESTS AND AXOM_DATA_DIR)
    # Add a test for the quest interface; Set up for MPI, when available
    if (ENABLE_MPI)
        axom_add_test(
            NAME quest_inout_interface_3D_mpi_test
            COMMAND quest_inout_interface_ex -i ${quest_data_dir}/sphere_binary.stl
            NUM_MPI_TASKS 2
            )
        if(C2C_FOUND)
            axom_add_test(
                NAME quest_inout_interface_2D_mpi_test
                COMMAND quest_inout_interface_ex -i ${AXOM_DATA_DIR}/contours/unit_circle.contour
                NUM_MPI_TASKS 2
                )
        endif()
    else()
        axom_add_test(
            NAME quest_inout_interface_3D_test
            COMMAND quest_inout_interface_ex -i ${quest_data_dir}/sphere_binary.stl
            )
        if(C2C_FOUND)
            axom_add_test(
                NAME quest_inout_interface_2D_test
                COMMAND quest_inout_interface_ex -i ${AXOM_DATA_DIR}/contours/unit_circle.contour
                )
        endif()
    endif()
endif()

if (ENABLE_FORTRAN)

    blt_list_append(TO quest_example_fortran_depends ELEMENTS mpi IF ENABLE_MPI)

    blt_add_executable(
        NAME       quest_signed_distance_interface_F_ex
        SOURCES    quest_signed_distance_interface.F
        OUTPUT_DIR ${EXAMPLE_OUTPUT_DIRECTORY}
        DEPENDS_ON ${quest_example_depends} ${quest_example_fortran_depends}
        FOLDER     axom/quest/examples
        )

    # This target fails to compile with hipcc/amdflang on debugs
    if(NOT "${CMAKE_Fortran_COMPILER}" MATCHES "amdflang$" OR NOT "${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
        blt_add_executable(
            NAME        quest_inout_interface_F_ex
            SOURCES     quest_inout_interface.F
            OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
            DEPENDS_ON  ${quest_example_depends} ${quest_example_fortran_depends}
            FOLDER      axom/quest/examples
            )

        if (NOT ENABLE_CUDA)
            if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
                set_target_properties(quest_inout_interface_F_ex
                                      PROPERTIES LINKER_LANGUAGE Fortran)
            else()
                set_target_properties(quest_inout_interface_F_ex
                                      PROPERTIES LINKER_LANGUAGE CXX)
            endif()
        endif()
    endif()

    # When CUDA is enabled, BLT will determine the correct linker, so don't override it here
    if (NOT ENABLE_CUDA)
        # When using the Intel compiler we need to link with the Fortran compiler to get openmp to work correctly.
        if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
            set_target_properties( quest_signed_distance_interface_F_ex
                                   PROPERTIES LINKER_LANGUAGE Fortran)
        else()
            set_target_properties( quest_signed_distance_interface_F_ex
                                   PROPERTIES LINKER_LANGUAGE CXX)
        endif()
    endif()

    if(AXOM_ENABLE_TESTS AND AXOM_DATA_DIR)
        if(NOT "${CMAKE_Fortran_COMPILER}" MATCHES "amdflang$" OR NOT "${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
            if (ENABLE_MPI)
                axom_add_test(
                    NAME quest_inout_interface_mpi_F_test
                    COMMAND quest_inout_interface_F_ex ${quest_data_dir}/sphere_binary.stl
                    NUM_MPI_TASKS 2
                    )
            else()
                axom_add_test(
                    NAME quest_inout_interface_F_test
                    COMMAND quest_inout_interface_F_ex ${quest_data_dir}/sphere_binary.stl
                    )
            endif()
        endif()
    endif()
endif()
