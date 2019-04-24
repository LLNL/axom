# Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

set( spin_tests
    spin_bvhtree.cpp
    spin_implicit_grid.cpp
    spin_morton.cpp
    spin_octree.cpp
    spin_rectangular_lattice.cpp
    spin_spatial_octree.cpp
    spin_uniform_grid.cpp
   )

foreach ( test ${spin_tests} )

   get_filename_component( test_name ${test} NAME_WE )

   blt_add_executable(
        NAME       ${test_name}_test
        SOURCES    ${test}
        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
        DEPENDS_ON core slam fmt spin primal slic gtest
        FOLDER     axom/spin/tests
       )

   blt_add_test(
        NAME    ${test_name}
        COMMAND ${test_name}_test
       )

endforeach()
