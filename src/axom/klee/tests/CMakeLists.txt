# Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
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
        klee_io.cpp
        klee_shape.cpp
   )

#
# Add gtest based tests
#
foreach(test ${gtest_klee_tests})
    get_filename_component( test_name ${test} NAME_WE )
    blt_add_executable( NAME ${test_name}_test
                        SOURCES ${test}
                        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                        DEPENDS_ON klee gtest
                        FOLDER      axom/klee/tests )
    blt_add_test( NAME ${test_name}
                  COMMAND ${test_name}_test )
endforeach()

