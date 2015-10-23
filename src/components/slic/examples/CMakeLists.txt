#-------------------------------------------------------------------------------
# Copyright (c) 2015, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
#
# All rights reserved.
#
# This source code cannot be distributed without permission and further
# review from Lawrence Livermore National Laboratory.
#-------------------------------------------------------------------------------

## Build examples using slic
make_executable( EXECUTABLE_SOURCE basic/logging_example.cc DEPENDS_ON slic )
make_executable( EXECUTABLE_SOURCE basic/cpplogger_example.cc DEPENDS_ON slic )
make_executable( EXECUTABLE_SOURCE multicode/driver.cpp DEPENDS_ON slic )

## Build MPI examples using slic
if (ENABLE_MPI)

   make_executable( WITH_MPI
                    EXECUTABLE_SOURCE basic/parallel_logging_example.cc
                    DEPENDS_ON slic )

   make_executable( WITH_MPI
                    EXECUTABLE_SOURCE basic/seperate_file_per_rank.cc
                    DEPENDS_ON slic )

endif()