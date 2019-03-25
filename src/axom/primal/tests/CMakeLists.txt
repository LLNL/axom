# Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

set( primal_tests
    primal_boundingbox.cpp
    primal_orientedboundingbox.cpp
    primal_bvh.cpp
    primal_clip.cpp
    primal_closest_point.cpp
    primal_in_sphere.cpp
    primal_intersect_impl.cpp
    primal_intersect.cpp
    primal_numeric_array.cpp
    primal_compute_bounding_box.cpp
    primal_orientation.cpp
    primal_plane.cpp
    primal_point.cpp
    primal_sphere.cpp
    primal_squared_distance.cpp
    primal_triangle.cpp
    primal_tetrahedron.cpp
    primal_vector.cpp
   )

foreach ( test ${primal_tests} )

   get_filename_component( test_name ${test} NAME_WE )

   set(primal_test_depends fmt primal slic gtest)
  
   blt_list_append( TO primal_test_depends ELEMENTS raja IF RAJA_FOUND )
   blt_list_append( TO primal_test_depends ELEMENTS umpire IF UMPIRE_FOUND )
   blt_list_append( TO primal_test_depends ELEMENTS openmp  IF ENABLE_OPENMP )

   blt_add_executable(
        NAME       ${test_name}_test
        SOURCES    ${test}
        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
        DEPENDS_ON ${primal_test_depends}
        FOLDER     axom/primal/tests
       )

   blt_add_test(
        NAME    ${test_name}
        COMMAND ${test_name}_test
       )

endforeach()
