# Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Multimat unit tests
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Specify list of tests
#------------------------------------------------------------------------------

set(gtest_multimat_tests
    multimat_test.cpp
    )


set(multimat_tests_depends_on core slic slam multimat gtest)

#------------------------------------------------------------------------------
# Add gtest based tests
#------------------------------------------------------------------------------
foreach(test ${gtest_multimat_tests})
    get_filename_component( test_name ${test} NAME_WE )
    blt_add_executable( NAME ${test_name}_test
                        SOURCES ${test}
                        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                        DEPENDS_ON ${multimat_tests_depends_on} 
                        FOLDER axom/multimat/tests )
    blt_add_test( NAME ${test_name}
                  COMMAND ${test_name}_test )
endforeach()
