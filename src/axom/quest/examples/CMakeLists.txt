# Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Quest examples
#------------------------------------------------------------------------------

set(quest_example_depends
        core
        slic
        mint
        slam
        quest
        )

blt_list_append(TO quest_example_depends ELEMENTS openmp IF ENABLE_OPENMP)

blt_list_append(TO quest_example_depends ELEMENTS cuda IF ENABLE_CUDA)

blt_add_executable(
    NAME        quest_containment_driver_ex
    SOURCES     containment_driver.cpp
    OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
    DEPENDS_ON  ${quest_example_depends} fmt cli11
    FOLDER      axom/quest/examples
    )

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

# Add a test for the quest interface.  Set up for MPI, when available
set(quest_data_dir  ${AXOM_DATA_DIR}/quest)

if(AXOM_ENABLE_TESTS AND AXOM_DATA_DIR)
    if (ENABLE_MPI)
        blt_add_test(
            NAME quest_inout_interface_mpi_test
            COMMAND quest_inout_interface_ex ${quest_data_dir}/sphere_binary.stl
            NUM_MPI_TASKS 2
            )
    else()
        blt_add_test(
            NAME quest_inout_interface_test
            COMMAND quest_inout_interface_ex ${quest_data_dir}/sphere_binary.stl
            )
    endif()
endif()

if (ENABLE_FORTRAN)

    blt_add_executable(
        NAME       quest_signed_distance_interface_F_ex
        SOURCES    quest_signed_distance_interface.F
        OUTPUT_DIR ${EXAMPLE_OUTPUT_DIRECTORY}
        DEPENDS_ON ${quest_example_depends}
        FOLDER     axom/quest/examples
        )

    blt_add_executable(
        NAME        quest_inout_interface_F_ex
        SOURCES     quest_inout_interface.F
        OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
        DEPENDS_ON  ${quest_example_depends}
        FOLDER      axom/quest/examples
        )

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

    if(AXOM_ENABLE_TESTS AND AXOM_DATA_DIR)
        if (ENABLE_MPI)
            blt_add_test(
                NAME quest_inout_interface_mpi_F_test
                COMMAND quest_inout_interface_F_ex ${quest_data_dir}/sphere_binary.stl
                NUM_MPI_TASKS 2
                )
        else()
            blt_add_test(
                NAME quest_inout_interface_F_test
                COMMAND quest_inout_interface_F_ex ${quest_data_dir}/sphere_binary.stl
                )
        endif()
    endif()
endif()
