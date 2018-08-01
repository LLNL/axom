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

set(quest_example_depends
        core
        slic
        mint
        slam
        quest
        )

if (ENABLE_OPENMP)
    list(APPEND quest_example_depends openmp)
endif()

blt_add_executable(
    NAME        quest_containment_driver_ex
    SOURCES     containment_driver.cpp
    OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
    DEPENDS_ON  ${quest_example_depends} fmt
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
    NAME        quest_interface_ex
    SOURCES     quest_interface.cpp
    OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
    DEPENDS_ON  ${quest_example_depends}
    FOLDER      axom/quest/examples
    )

# Add a test for the quest interface.  Set up for MPI, when available
set(quest_data_dir  ${CMAKE_CURRENT_SOURCE_DIR}/../data)
if (ENABLE_MPI)
    blt_add_test(
        NAME quest_interface_mpi_test
        COMMAND quest_interface_ex ${quest_data_dir}/sphere_binary.stl
        NUM_MPI_TASKS 2
        )
else()
    blt_add_test(
        NAME quest_interface_test
        COMMAND quest_interface_ex ${quest_data_dir}/sphere_binary.stl
        )
endif()


if (ENABLE_FORTRAN)

    blt_add_executable(
        NAME       quest_signed_distance_interface_F_ex
        SOURCES    quest_signed_distance_interface.F
        OUTPUT_DIR ${EXAMPLE_OUTPUT_DIRECTORY}
        DEPENDS_ON ${quest_example_depends}
        FOLDER     axom/quest/examples
        )

    set_target_properties( quest_signed_distance_interface_F_ex
                           PROPERTIES LINKER_LANGUAGE Fortran)

    blt_add_executable(
        NAME        quest_interface_F_ex
        SOURCES     quest_interface_fortran.F
        OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
        DEPENDS_ON  ${quest_example_depends}
        FOLDER      axom/quest/examples
        )

    set_target_properties(quest_interface_F_ex
                          PROPERTIES LINKER_LANGUAGE Fortran)

    if (ENABLE_MPI)
        blt_add_test(
            NAME quest_interface_mpi_F_test
            COMMAND quest_interface_F_ex ${quest_data_dir}/sphere_binary.stl
            NUM_MPI_TASKS 2
            )
    else()
        blt_add_test(
            NAME quest_interface_F_test
            COMMAND quest_interface_F_ex ${quest_data_dir}/sphere_binary.stl
            )
    endif()
endif()
