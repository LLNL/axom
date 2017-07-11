#-------------------------------------------------------------------------------
# Copyright (c) 2015, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
#
# All rights reserved.
#
# This source code cannot be distributed without permission and further
# review from Lawrence Livermore National Laboratory.
#-------------------------------------------------------------------------------

set(quest_example_depends
        axom_utils
        slic
        mint
        slam
        quest
        )

if (BOOST_FOUND)
    list(APPEND quest_example_depends boost)
endif()

if (ENABLE_OPENMP)
    list(APPEND quest_example_depends openmp)
endif()

blt_add_executable(
    NAME        quest_distance_driver_ex
    SOURCES     distance_driver.cpp
    OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
    DEPENDS_ON  ${quest_example_depends}
    )

blt_add_executable(
    NAME        quest_sphere_distance_ex
    SOURCES     sphere_distance.cpp
    OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
    DEPENDS_ON  ${quest_example_depends}
    )

blt_add_executable(
    NAME        quest_containment_driver_ex
    SOURCES     containment_driver.cpp
    OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
    DEPENDS_ON  ${quest_example_depends} fmt
    )

blt_add_executable(
    NAME        quest_mesh_tester_ex
    SOURCES     mesh_tester.cpp
    OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
    DEPENDS_ON  ${quest_example_depends}
    )

if (ENABLE_MPI)

  set(quest_interface_extra_depends)
  if (ENABLE_LUMBERJACK)
    list(APPEND quest_interface_extra_depends lumberjack)
  endif()

  blt_add_executable(
    NAME        quest_interface_ex
    SOURCES     quest_interface.cpp
    OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
    DEPENDS_ON  ${quest_example_depends} ${quest_interface_extra_depends}
    )

  if (ENABLE_FORTRAN)
    blt_add_executable(
        NAME        quest_interface_F_ex
        SOURCES     quest_interface_fortran.F
        OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
        DEPENDS_ON  ${quest_example_depends} ${quest_interface_extra_depends}
        )

    set_target_properties(quest_interface_F_ex
                          PROPERTIES LINKER_LANGUAGE Fortran)

  endif()

endif()
