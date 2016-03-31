#-------------------------------------------------------------------------------
# Copyright (c) 2015, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
#
# All rights reserved.
#
# This source code cannot be distributed without permission and further
# review from Lawrence Livermore National Laboratory.
#-------------------------------------------------------------------------------

make_executable(EXECUTABLE_NAME lumberjack_basicExample
                EXECUTABLE_SOURCE basicExample.cpp
                DEPENDS_ON lumberjack
                IS_EXAMPLE )
