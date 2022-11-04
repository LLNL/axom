# Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Unit tests for Klee component
#------------------------------------------------------------------------------

#
# Specify list of tests
#

set(gtest_klee_tests
    klee_config.cpp
    klee_geometry.cpp
    klee_geometry_operators.cpp
    klee_geometry_operators_io.cpp
    klee_io.cpp
    klee_io_util.cpp
    klee_shape.cpp
    klee_shape_set.cpp
    klee_units.cpp
   )

blt_add_library(
        NAME        klee_test_utils
        SOURCES     KleeTestUtils.cpp
        HEADERS     KleeTestUtils.hpp KleeMatchers.hpp
        DEPENDS_ON  axom gtest gmock
        FOLDER      axom/klee/tests)

set(klee_device_test_depends)
blt_list_append( TO klee_device_test_depends ELEMENTS cuda IF ${ENABLE_CUDA} )

#
# Add gtest based tests
#
foreach(test ${gtest_klee_tests})
    get_filename_component( test_name ${test} NAME_WE )
    blt_add_executable( NAME ${test_name}_test
                        SOURCES ${test}
                        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                        DEPENDS_ON gtest klee_test_utils ${klee_device_test_depends}
                        FOLDER      axom/klee/tests )
    blt_add_test( NAME ${test_name}
                  COMMAND ${test_name}_test )
endforeach()

