#-------------------------------------------------------------------------------
# Copyright (c) 2015, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
#
# All rights reserved.
#
# This source code cannot be distributed without permission and further
# review from Lawrence Livermore National Laboratory.
#-------------------------------------------------------------------------------


make_executable( EXECUTABLE_SOURCE quest_driver.cpp DEPENDS_ON quest slic)
make_executable( EXECUTABLE_NAME quest_sphere EXECUTABLE_SOURCE sphere.cpp DEPENDS_ON quest slic )
