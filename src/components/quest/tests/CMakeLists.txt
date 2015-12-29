#-------------------------------------------------------------------------------
# Copyright (c) 2015, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
#
# All rights reserved.
#
# This source code cannot be distributed without permission and further
# review from Lawrence Livermore National Laboratory.
#-------------------------------------------------------------------------------

set(quest_tests
    quest_mesh_extent.cpp
    quest_boundingBox.cpp
    quest_intersection.cpp
    quest_numeric_array.cpp
    quest_orient.cpp
    quest_point.cpp
    quest_squared_distance.cpp
    quest_vector.cpp
   )

foreach(test ${quest_tests})
  add_gtest( TEST_SOURCE ${test} DEPENDS_ON quest slic )
endforeach()