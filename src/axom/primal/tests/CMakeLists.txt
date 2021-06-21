# Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Primal unit tests
#------------------------------------------------------------------------------

set( primal_tests
    primal_bezier_curve.cpp
    primal_bezier_intersect.cpp
    primal_boundingbox.cpp
    primal_clip.cpp
    primal_closest_point.cpp
    primal_compute_bounding_box.cpp
    primal_in_sphere.cpp
    primal_intersect.cpp
    primal_intersect_impl.cpp
    primal_numeric_array.cpp
    primal_orientation.cpp
    primal_orientedboundingbox.cpp
    primal_plane.cpp
    primal_point.cpp
    primal_ray_intersect.cpp
    primal_bounding_box_intersect.cpp
    primal_segment.cpp
    primal_sphere.cpp
    primal_squared_distance.cpp
    primal_tetrahedron.cpp
    primal_triangle.cpp
    primal_vector.cpp
    primal_curved_polygon.cpp
    primal_curved_polygon_intersect.cpp
   )

set(primal_test_depends axom fmt gtest)

blt_list_append( TO primal_test_depends ELEMENTS cuda IF ${ENABLE_CUDA} )

foreach ( test ${primal_tests} )

   get_filename_component( test_name ${test} NAME_WE )

   blt_add_executable(
        NAME       ${test_name}_test
        SOURCES    ${test}
        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
        DEPENDS_ON ${primal_test_depends}
        FOLDER     axom/primal/tests
       )

   axom_add_test(
        NAME    ${test_name}
        COMMAND ${test_name}_test
       )

endforeach()
