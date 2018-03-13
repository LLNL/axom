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

#
# Build lulesh example
#

set(lulesh_headers
    lulesh.h
    lulesh_tuple.h
    )

set(lulesh_sources
    lulesh.cc
    lulesh-comm.cc
    lulesh-init.cc
    lulesh-util.cc
    lulesh-viz.cc
    )

blt_add_executable(
    NAME sidre_lulesh2_ex
    SOURCES ${lulesh_sources} ${lulesh_headers}
    OUTPUT_DIR ${EXAMPLE_OUTPUT_DIRECTORY}
    DEPENDS_ON axom_utils sidre conduit conduit_relay hdf5 slic
    FOLDER axom/sidre/examples
    )

if(ENABLE_MPI)
    blt_add_test(NAME sidre_lulesh2
                 COMMAND sidre_lulesh2_ex -s 5
                 NUM_MPI_TASKS 8
                )
else()
    blt_add_test(NAME sidre_lulesh2
                 COMMAND sidre_lulesh2_ex -s 5
                )
endif()
