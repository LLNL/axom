#-------------------------------------------------------------------------------
# Copyright (c) 2015, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
#
# All rights reserved.
#
# This source code cannot be distributed without permission and further
# review from Lawrence Livermore National Laboratory.
#-------------------------------------------------------------------------------

if (ENABLE_MPI)
   make_executable( ADD_CTEST
                    EXECUTABLE_SOURCE basic/parallel_logging_example.cc
                    DEPENDS_ON slic )

   make_executable( ADD_CTEST
                    EXECUTABLE_SOURCE basic/seperate_file_per_rank.cc
                    DEPENDS_ON slic )
endif()

make_executable( ADD_CTEST
                 EXECUTABLE_SOURCE basic/cpplogger_example.cc
                 DEPENDS_ON slic )

make_executable( ADD_CTEST
                 EXECUTABLE_SOURCE basic/logging_example.cc
                 DEPENDS_ON slic )

make_executable( ADD_CTEST
                 EXECUTABLE_SOURCE multicode/driver.cpp
                 DEPENDS_ON slic )
