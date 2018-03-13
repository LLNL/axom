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

set( primal_tests
    primal_boundingbox.cpp
    primal_orientedboundingbox.cpp
    primal_bvhtree.cpp
    primal_clip.cpp
    primal_closest_point.cpp
    primal_intersect_impl.cpp
    primal_intersect.cpp
    primal_morton.cpp
    primal_numeric_array.cpp
    primal_compute_bounding_box.cpp
    primal_orientation.cpp
    primal_point.cpp
    primal_rectangular_lattice.cpp
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
        DEPENDS_ON fmt primal slic gtest
        FOLDER     axom/primal/tests
       )

   blt_add_test(
        NAME    ${test_name}
        COMMAND ${test_name}_test
       )

endforeach()
