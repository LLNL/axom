# Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Unit tests for Spin component
#------------------------------------------------------------------------------

set(spin_tests
    spin_bvhtree.cpp
    spin_implicit_grid.cpp
    spin_morton.cpp
    spin_octree.cpp
    spin_rectangular_lattice.cpp
    spin_spatial_octree.cpp
    spin_uniform_grid.cpp
   )

set(spin_depends
    axom
    gtest
    fmt
    )

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

    axom_add_test(
      NAME    ${test_name}
      COMMAND ${test_name}_test
      )

endforeach()

if ( RAJA_FOUND AND UMPIRE_FOUND )

    blt_add_executable(
      NAME       spin_bvh_test
      SOURCES    spin_bvh.cpp
      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
      DEPENDS_ON ${spin_depends}
      FOLDER     axom/spin/tests
      )

    axom_add_test(
      NAME            spin_bvh
      COMMAND         spin_bvh_test
      NUM_OMP_THREADS 4
      )

    # fix for rzmanta of following runtime error:
    #   CUDAassert: peer access is not supported between these two devices /unimportant/path/raja-0.12.1/include/RAJA/policy/cuda/sort.hpp 341
    #   unknown file: Failure
    #   C++ exception with description "CUDAassert" thrown in the test body.
    set_property(TEST spin_bvh
                 APPEND
                 PROPERTY ENVIRONMENT  "CUDA_VISIBLE_DEVICES=0") # Zero based indexed, so only first GPU is visible

endif()
