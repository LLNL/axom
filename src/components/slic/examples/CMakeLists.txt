#-------------------------------------------------------------------------------
# Copyright (c) 2015, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
#
# All rights reserved.
#
# This source code cannot be distributed without permission and further
# review from Lawrence Livermore National Laboratory.
#-------------------------------------------------------------------------------

## Build examples using logapi
make_executable( EXECUTABLE_SOURCE logging_example.cc DEPENDS_ON slic )

## Build MPI examples using logapi
if (ENABLE_MPI)

   make_executable( LINK_MPI
                    EXECUTABLE_SOURCE parallel_logging_example.cc
                    DEPENDS_ON slic )

   make_executable( LINK_MPI
                    EXECUTABLE_SOURCE seperate_file_per_rank.cc
                    DEPENDS_ON slic )

endif()