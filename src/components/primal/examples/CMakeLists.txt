#-------------------------------------------------------------------------------
# Copyright (c) 2015, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
#
# All rights reserved.
#
# This source code cannot be distributed without permission and further
# review from Lawrence Livermore National Laboratory.
#-------------------------------------------------------------------------------

set(primal_example_depends
        axom_utils
        slic
        mint
        slam
        quest
        )

if (BOOST_FOUND)
    list(APPEND primal_example_depends boost)
endif()

blt_add_executable(
    NAME        mesh_tester
    SOURCES     mesh_tester.cpp
    OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
    DEPENDS_ON  ${primal_example_depends}
    USE_OPENMP  ${ENABLE_OPENMP}
    )

