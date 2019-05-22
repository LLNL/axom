# Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Mir unit tests
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# Specify list of tests
#------------------------------------------------------------------------------

set(gtest_mir_tests
    mir_smoke.cpp
    )


set(mir_tests_depends_on core slic mir gtest)

#------------------------------------------------------------------------------
# Add gtest based tests
#------------------------------------------------------------------------------
foreach(test ${gtest_mir_tests})
    get_filename_component( test_name ${test} NAME_WE )
    blt_add_executable( NAME ${test_name}_test
                        SOURCES ${test}
                        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                        DEPENDS_ON ${mir_tests_depends_on}
                        FOLDER axom/mir/tests )
    blt_add_test( NAME ${test_name}
                  COMMAND ${test_name}_test )
endforeach()

