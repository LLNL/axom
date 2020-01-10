# Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
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

set( spin_depends
     core
     slic
     slam
     primal
     spin
     gtest
     fmt
    )

blt_list_append( TO spin_depends ELEMENTS raja IF ${RAJA_FOUND} )
blt_list_append( TO spin_depends ELEMENTS umpire IF ${UMPIRE_FOUND} )
blt_list_append( TO spin_depends ELEMENTS openmp  IF ${ENABLE_OPENMP} )
blt_list_append( TO spin_depends ELEMENTS cuda IF ${ENABLE_CUDA} )

foreach ( test ${spin_tests} )

   get_filename_component( test_name ${test} NAME_WE )

   blt_add_executable(
        NAME       ${test_name}_test
        SOURCES    ${test}
        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
        DEPENDS_ON ${spin_depends}
        FOLDER     axom/spin/tests
       )

   blt_add_test(
        NAME    ${test_name}
        COMMAND ${test_name}_test
       )

endforeach()

if ( RAJA_FOUND AND UMPIRE_FOUND )

   blt_add_executable(
      NAME       spin_bvh_test
      SOURCES    spin_bvh.cpp
      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
      DEPENDS_ON mint ${spin_depends}
      FOLDER     axom/spin/tests
     )

   blt_add_test(
      NAME    spin_bvh
      COMMAND spin_bvh_test
      NUM_OMP_THREADS 4
     )

endif()

