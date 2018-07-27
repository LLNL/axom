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

######################################
# Slam (MeshAPI) version of lulesh 2.0
######################################


set(lulesh_headers
    lulesh.hpp
    # lulesh_tuple.hpp
    )

set(lulesh_sources
    lulesh.cpp
    lulesh-comm.cpp
    lulesh-init.cpp
    lulesh-util.cpp
    lulesh-viz.cpp
    )


# Ignore omp pragmas when openMP is not defined
if ( NOT ENABLE_OPENMP )
    set_source_files_properties(${lulesh_sources}
        PROPERTIES COMPILE_FLAGS "${AXOM_DISABLE_OMP_PRAGMA_WARNINGS}" )
endif()

set(slam_tests_depends_on axom_utils slic slam)

if (ENABLE_OPENMP)
    list(APPEND slam_tests_depends_on openmp)
endif()

if (ENABLE_MPI)
    list(APPEND slam_tests_depends_on mpi)
endif()

blt_add_executable(
    NAME        slam_lulesh_ex
    SOURCES     ${lulesh_sources} ${lulesh_headers}
    OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
    DEPENDS_ON  ${slam_tests_depends_on}
    FOLDER      axom/slam/examples
    )

if(ENABLE_MPI)
    blt_add_test(NAME        slam_lulesh
                 COMMAND     slam_lulesh_ex
                 NUM_MPI_TASKS   8
                )
else()
    blt_add_test(NAME        slam_lulesh
                 COMMAND     slam_lulesh_ex
                )
endif()

