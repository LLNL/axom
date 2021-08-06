# Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Quest examples
#------------------------------------------------------------------------------

set(quest_example_depends axom fmt cli11)

blt_list_append(TO quest_example_depends ELEMENTS cuda IF ENABLE_CUDA)
blt_list_append(TO quest_example_depends ELEMENTS RAJA IF RAJA_FOUND)

blt_add_executable(
    NAME        quest_containment_driver_ex
    SOURCES     containment_driver.cpp
    OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
    DEPENDS_ON  ${quest_example_depends}
    FOLDER      axom/quest/examples
    )

if(MFEM_FOUND AND AXOM_ENABLE_MFEM_SIDRE_DATACOLLECTION AND ENABLE_MPI)
    blt_add_executable(
        NAME        quest_shaping_driver_ex
        SOURCES     shaping_driver.cpp
        OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
        DEPENDS_ON  ${quest_example_depends} mfem mpi
        FOLDER      axom/quest/examples
        )
endif()

## Quest interface examples

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

# Add a test for the quest interface; Set up for MPI, when available
set(quest_data_dir  ${AXOM_DATA_DIR}/quest)

if(AXOM_ENABLE_TESTS AND AXOM_DATA_DIR)
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

    blt_add_executable(
        NAME        quest_inout_interface_F_ex
        SOURCES     quest_inout_interface.F
        OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
        DEPENDS_ON  ${quest_example_depends} ${quest_example_fortran_depends}
        FOLDER      axom/quest/examples
        )

    # When CUDA is enabled, BLT will determine the correct linker, so don't override it here
    if (NOT ENABLE_CUDA)
        # When using the Intel compiler we need to link with the Fortran compiler to get openmp to work correctly.
        if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
            set_target_properties( quest_signed_distance_interface_F_ex
                                   PROPERTIES LINKER_LANGUAGE Fortran)

            set_target_properties(quest_inout_interface_F_ex
                                  PROPERTIES LINKER_LANGUAGE Fortran)
        else()
            set_target_properties( quest_signed_distance_interface_F_ex
                                   PROPERTIES LINKER_LANGUAGE CXX)

            set_target_properties(quest_inout_interface_F_ex
                                  PROPERTIES LINKER_LANGUAGE CXX)
        endif()
    endif()

    if(AXOM_ENABLE_TESTS AND AXOM_DATA_DIR)
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
