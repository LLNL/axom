#-------------------------------------------------------------------------------
# Copyright (c) 2015, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
#
# All rights reserved.
#
# This source code cannot be distributed without permission and further
# review from Lawrence Livermore National Laboratory.
#-------------------------------------------------------------------------------

set( primal_tests
    primal_boundingbox.cpp
    primal_bvhtree.cpp
    primal_clip.cpp
    primal_intersect_impl.cpp
    primal_intersect.cpp
    primal_numeric_array.cpp
    primal_orientation.cpp
    primal_point.cpp
    primal_squared_distance.cpp
    primal_triangle.cpp
    primal_vector.cpp
    primal_uniform_grid.cpp
   )

foreach ( test ${primal_tests} )

   get_filename_component( test_name ${test} NAME_WE )

   blt_add_executable(
        NAME       ${test_name}_test
        SOURCES    ${test}
        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
        DEPENDS_ON primal slic gtest
       )

   blt_add_test(
        NAME    ${test_name}
        COMMAND ${test_name}_test
       )

endforeach()
