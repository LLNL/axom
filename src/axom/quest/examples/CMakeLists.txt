# Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Quest examples
#------------------------------------------------------------------------------

set(quest_example_depends quest fmt cli11)

blt_list_append(TO quest_example_depends ELEMENTS RAJA IF RAJA_FOUND)

set(quest_data_dir  ${AXOM_DATA_DIR}/quest)
set(shaping_data_dir  ${AXOM_DATA_DIR}/shaping)

# In/out octree containment example -------------------------------------------
axom_add_executable(
    NAME        quest_containment_driver_ex
    SOURCES     containment_driver.cpp
    OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
    DEPENDS_ON  ${quest_example_depends}
    FOLDER      axom/quest/examples
    )

# BVH two pass example --------------------------------------------------------
if (RAJA_FOUND AND UMPIRE_FOUND)
    axom_add_executable(
        NAME        quest_bvh_two_pass_ex
        SOURCES     quest_bvh_two_pass.cpp
        OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
        DEPENDS_ON  ${quest_example_depends}
        FOLDER      axom/quest/examples
        )
endif()

# Shaping example -------------------------------------------------------------
if(AXOM_ENABLE_MPI AND MFEM_FOUND AND MFEM_USE_MPI
                   AND AXOM_ENABLE_SIDRE AND AXOM_ENABLE_MFEM_SIDRE_DATACOLLECTION
                   AND AXOM_ENABLE_KLEE)
    axom_add_executable(
        NAME        quest_shaping_driver_ex
        SOURCES     shaping_driver.cpp
        OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
        DEPENDS_ON  ${quest_example_depends} mfem
        FOLDER      axom/quest/examples
        )

    if(AXOM_ENABLE_TESTS AND AXOM_DATA_DIR)
        # 2D shaping tests
        set(_nranks 1)
        if(C2C_FOUND)
            set(_testname quest_shaping_driver_ex_sampling_circles)
            axom_add_test(
                NAME ${_testname}
                COMMAND  quest_shaping_driver_ex 
                            -i ${shaping_data_dir}/circles.yaml
                            --method sampling
                            inline_mesh --min -6 -6 --max 6 6 --res 25 25 -d 2
                NUM_MPI_TASKS ${_nranks})
            # Analytic area for annulus w/ outer/inner radii 5 and 2.5 is ~58.905
            set_tests_properties(${_testname} PROPERTIES 
                PASS_REGULAR_EXPRESSION  "Volume of material 'steel' is 58.")

            set(_testname quest_shaping_driver_ex_sampling_circles_opposite)
            axom_add_test(
                NAME ${_testname}
                COMMAND  quest_shaping_driver_ex 
                        -i ${shaping_data_dir}/circles_opposite.yaml
                        --method sampling
                        inline_mesh --min -6 -6 --max 6 6 --res 25 25 -d 2
                NUM_MPI_TASKS ${_nranks})
            # Analytic area for annulus w/ outer/inner radii 5 and 2.5 is ~65.97
            set_tests_properties(${_testname} PROPERTIES 
                PASS_REGULAR_EXPRESSION  "Volume of material 'steel' is 65.")

            set(_testname quest_shaping_driver_ex_sampling_balls_and_jacks)
            axom_add_test(
                NAME ${_testname}
                COMMAND  quest_shaping_driver_ex 
                        -i ${shaping_data_dir}/balls_and_jacks.yaml
                        --method sampling
                        --background-material air
                        inline_mesh --min -10 0 --max 60 50 --res 25 25 -d 2
                NUM_MPI_TASKS ${_nranks})
            # background: 3500; balls: ~373; jacks: ~230; air ~2897
            set_tests_properties(${_testname} PROPERTIES 
                PASS_REGULAR_EXPRESSION  "Volume of material 'air' is 2895.")

            set(_testname quest_shaping_driver_ex_sampling_ball_impact)
            axom_add_test(
                NAME ${_testname}
                COMMAND  quest_shaping_driver_ex 
                        -i ${shaping_data_dir}/ball_impact.yaml
                        --method sampling
                        --background-material air
                        inline_mesh --min -10 0 --max 40 30 --res 25 25 -d 2
                NUM_MPI_TASKS ${_nranks})
            # steel ball: semi-cirle of radius 5.5 has area ~47.516
            set_tests_properties(${_testname} PROPERTIES 
                PASS_REGULAR_EXPRESSION  "Volume of material 'steel' is 47.")

            # Tests 3D mesh against 2D contours
            set(_testname quest_shaping_driver_ex_sampling_ball_impact_3d)
            axom_add_test(
                NAME ${_testname}
                COMMAND  quest_shaping_driver_ex 
                        -i ${shaping_data_dir}/ball_impact.yaml
                        --method sampling
                        --background-material air
                        inline_mesh --min -30 -30 -10 --max 30 30 40 --res 16 16 16 -d3
                NUM_MPI_TASKS ${_nranks})
            # steel ball: cirle of radius 5.5 has analytic volume ~696.91
            set_tests_properties(${_testname} PROPERTIES 
                PASS_REGULAR_EXPRESSION  "Volume of material 'steel' is 6[89][0-9].")          
        endif()

        # 3D shaping tests
        set(_nranks 3)
        set(_testname quest_shaping_driver_ex_sampling_sphere)
        axom_add_test(
            NAME ${_testname}
            COMMAND  quest_shaping_driver_ex 
                    -i ${shaping_data_dir}/sphere.yaml
                    --method sampling
                    inline_mesh --min -6 -6 -6 --max 6 6 6 --res 16 16 16 -d 3
            NUM_MPI_TASKS ${_nranks})
        # sphere of radius 5 has volume: 4/3 pi * r^3 ~523.6
        # input sphere is discretized and mesh is coarse; accuracy is within ~1%
        set_tests_properties(${_testname} PROPERTIES 
            PASS_REGULAR_EXPRESSION  "Volume of material 'steel' is 519.")

        set(_testname quest_shaping_driver_ex_sampling_spheres)
        axom_add_test(
            NAME ${_testname}
            COMMAND  quest_shaping_driver_ex 
                    -i ${shaping_data_dir}/spheres.yaml
                    --method sampling
                    --background-material void
                    inline_mesh --min -6 -6 -6 --max 6 6 6 --res 16 16 16 -d 3
            NUM_MPI_TASKS ${_nranks})
        # bbox volume: 12^3 = 1728; sphere(r=5): ~523.6; sphere(r=2): 33.5
        # expected analytic volume when fully resolved: ~1237.9
        set_tests_properties(${_testname} PROPERTIES 
            PASS_REGULAR_EXPRESSION  "Volume of material 'void' is 1239.")

        set(_testname quest_shaping_driver_ex_sampling_plane)
        axom_add_test(
            NAME ${_testname}
            COMMAND  quest_shaping_driver_ex 
                    -i ${shaping_data_dir}/plane.yaml
                    --method sampling
                    --background-material air
                    inline_mesh --min -430 -45 -35 --max 430 855 210 --res 16 16 16 -d 3
            NUM_MPI_TASKS ${_nranks})
        # background: 860*900*245=189,630,000; airplane volume from meshlab: 7,420,578.5
        # expected volume when fully resolved: 182,209,421.5
        # actual value is accurate to about 1% at coarse resolution
        set_tests_properties(${_testname} PROPERTIES 
            PASS_REGULAR_EXPRESSION  "Volume of material 'air' is 182227963")
    endif()
endif()

# Distributed closest point example -------------------------------------------
if(AXOM_ENABLE_MPI AND AXOM_ENABLE_SIDRE AND HDF5_FOUND)
    axom_add_executable(
            NAME        quest_distributed_distance_query_ex
            SOURCES     quest_distributed_distance_query_example.cpp
            OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
            DEPENDS_ON  ${quest_example_depends}
            FOLDER      axom/quest/examples
            )

    if(AXOM_ENABLE_TESTS AND AXOM_DATA_DIR)
        set(_nranks 3)

        # Run the distributed closest point example on N ranks for each enabled policy
        set(_policies "seq")
        if(RAJA_FOUND)
            blt_list_append(TO _policies ELEMENTS "cuda" IF AXOM_ENABLE_CUDA)
            blt_list_append(TO _policies ELEMENTS "hip"  IF AXOM_ENABLE_HIP)
            blt_list_append(TO _policies ELEMENTS "omp"  IF AXOM_ENABLE_OPENMP)
        endif()

        # Non-zero empty-rank probability tests domain underloading case
        set(_meshes "mdmesh.2x1" "mdmesh.2x3")
        # The mdmesh.* files were generated by these commands:
        # src/tools/gen-multidom-structured-mesh.py -ml=0,0 -mu=2,2 -ms=100,100 -dc=2,1 -o mdmesh.2x1
        # src/tools/gen-multidom-structured-mesh.py -ml=0,0 -mu=2,2 -ms=100,100 -dc=2,3 -o mdmesh.2x3
        foreach(_pol ${_policies})
            foreach(_mesh ${_meshes})
                # Add test with data on all ranks
                set(_test "quest_distributed_closest_point_run_2D_${_pol}_${_mesh}")
                axom_add_test(
                    NAME    ${_test}
                    COMMAND quest_distributed_distance_query_ex
                                --mesh-file ${quest_data_dir}/${_mesh}.root
                                --num-samples 500
                                --center 1.2 1.5
                                --radius 0.75
                                --obj-domain-count-range 0 2
                                --dist-threshold .3
                                --no-random-spacing
                                --check-results
                                --policy ${_pol}
                                --object-file dcp_object_mesh_2d_${_pol}_${_mesh}
                                --distance-file dcp_closest_point_2d_${_pol}_${_mesh}
                    NUM_MPI_TASKS ${_nranks})
                if(${_pol} STREQUAL "omp")
                    set_property(TEST ${_test} APPEND PROPERTY ENVIRONMENT OMP_NUM_THREADS=4)
                endif()
            endforeach()
        endforeach()

        unset(_nranks)
        unset(_policies)
        unset(_test)
    endif()
endif()

# Marching cubes example -------------------------------------------
if(CONDUIT_FOUND)
    axom_add_executable(
            NAME        quest_marching_cubes_ex
            SOURCES     quest_marching_cubes_example.cpp
            OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
            DEPENDS_ON  ${quest_example_depends}
            FOLDER      axom/quest/examples
            )

    # These examples currently segfault on windows without MPI when loading the test data into conduit
    if(AXOM_ENABLE_TESTS AND AXOM_DATA_DIR AND NOT WIN32)
        if (ENABLE_MPI)
            set(_nranks 3)
        endif()

        # Run the marching cubes example on N ranks for each enabled policy
        set(_policies "seq")
        if(RAJA_FOUND)
            blt_list_append(TO _policies ELEMENTS "omp" IF AXOM_ENABLE_OPENMP)
            blt_list_append(TO _policies ELEMENTS "cuda" IF AXOM_ENABLE_CUDA)
            blt_list_append(TO _policies ELEMENTS "hip" IF AXOM_ENABLE_HIP)
        endif()

        # Non-zero empty-rank probability tests domain underloading case
        set(_meshes "mdmesh.2x1" "mdmesh.2x3" "mdmesh.2x2x1")
        # The amc.* files were generated by these commands:
        # src/tools/gen-multidom-structured-mesh.py -ml=0,0 -mu=2,2 -ms=20,20 -dc=2,2 -o amc_mesh.2x2
        # src/tools/gen-multidom-structured-mesh.py -ml=0,0,0 -mu=2,2,2 -ms=20,20,15 -dc=2,2,1 -o mdmesh.2x2x1
        foreach(_pol ${_policies})
            foreach(_mesh ${_meshes})
                # Determine problem dimension by mesh filename.
                # and set dimension-dependent arguments.
                string(REGEX MATCH "\\.[0-9]+(x[0-9]+)+$" _sizes "${_mesh}")
                string(REGEX MATCHALL "[0-9]+" _sizes ${_sizes})
                list(LENGTH _sizes _ndim)

                if(_ndim EQUAL 2)
                  set(_dir 1.0 0.4)
                  set(_center 1.0 0.4)
                elseif(_ndim EQUAL 3)
                  set(_dir 1.0 0.4 1.2)
                  set(_center 1.0 0.4 1.2)
                endif()

                set(_test "quest_marching_cubes_run_${_ndim}D_${_pol}_${_mesh}")
                axom_add_test(
                    NAME    ${_test}
                    COMMAND quest_marching_cubes_ex
                                --policy ${_pol}
                                --mesh-file ${quest_data_dir}/${_mesh}.root
                                --fields-file ${_test}.field
                                --dir ${_dir}
                                --center ${_center}
                                --contourVal 1.25
                                --check-results
                                --verbose
                    NUM_MPI_TASKS ${_nranks})
            endforeach()
        endforeach()

        unset(_nranks)
        unset(_policies)
        unset(_test)
    endif()

endif()

# Point in cell example -------------------------------------------------------
if(MFEM_FOUND)
    axom_add_executable(
        NAME        quest_point_in_cell_benchmark_ex
        SOURCES     point_in_cell_benchmark.cpp
        OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
        DEPENDS_ON  ${quest_example_depends} mfem
        FOLDER      axom/quest/examples
        )
endif()

# Delaunay triangulation example ----------------------------------------------
axom_add_executable(
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
    axom_add_executable(
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

axom_add_executable(
    NAME       quest_signed_distance_interface_ex
    SOURCES    quest_signed_distance_interface.cpp
    OUTPUT_DIR ${EXAMPLE_OUTPUT_DIRECTORY}
    DEPENDS_ON ${quest_example_depends}
    FOLDER      axom/quest/examples
    )

axom_add_executable(
    NAME        quest_inout_interface_ex
    SOURCES     quest_inout_interface.cpp
    OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
    DEPENDS_ON  ${quest_example_depends}
    FOLDER      axom/quest/examples
    )

if(AXOM_ENABLE_TESTS AND AXOM_DATA_DIR)
    # Add a test for the quest interface; Set up for MPI, when available
    if (AXOM_ENABLE_MPI)
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

    set(quest_fortran_examples quest_signed_distance_interface)
    
    # The inout Fortran example fails to compile with hipcc/amdflang in debug configurations
    if("${CMAKE_Fortran_COMPILER}" MATCHES "amdflang$" AND "${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
        set(_has_inout_fortran_example FALSE)
    else()
        set(_has_inout_fortran_example TRUE)
        blt_list_append(TO quest_fortran_examples ELEMENTS quest_inout_interface)
    endif()

    foreach(_t ${quest_fortran_examples})
        set(_example_name "${_t}_F_ex")
        axom_add_executable(
            NAME       ${_example_name}
            SOURCES    ${_t}.F
            OUTPUT_DIR ${EXAMPLE_OUTPUT_DIRECTORY}
            DEPENDS_ON ${quest_example_depends}
            FOLDER     axom/quest/examples)

        # When CUDA is enabled, BLT will determine the correct linker, so don't override it here
        if(NOT AXOM_ENABLE_CUDA)
            if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "IntelLLVM")
                set_target_properties(${_example_name}
                                      PROPERTIES LINKER_LANGUAGE Fortran)
            else()
                set_target_properties(${_example_name}
                                      PROPERTIES LINKER_LANGUAGE CXX)
            endif()
        endif()
    endforeach()

    if(AXOM_ENABLE_TESTS AND AXOM_DATA_DIR AND ${_has_inout_fortran_example})
        if (AXOM_ENABLE_MPI)
            axom_add_test(
                NAME quest_inout_interface_mpi_F_test
                COMMAND quest_inout_interface_F_ex ${quest_data_dir}/sphere_binary.stl
                NUM_MPI_TASKS 2)
        else()
            axom_add_test(
                NAME quest_inout_interface_F_test
                COMMAND quest_inout_interface_F_ex ${quest_data_dir}/sphere_binary.stl)
        endif()
    endif()
endif()
