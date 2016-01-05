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
    test_extent.cpp
    test_quest_boundingBox.cpp
    test_quest_orient.cpp
    test_quest_point.cpp
    test_quest_squared_distance.cpp
    test_quest_vector.cpp
    test_quest_numeric_array.cpp
   )

foreach(test ${quest_tests})
  add_gtest( TEST_SOURCE ${test} DEPENDS_ON quest slic )
endforeach()