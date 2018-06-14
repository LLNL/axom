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

#
# Specify list of tests
#

set(gtest_multimat_tests
    # test sets
    multimat_test.cpp
    )


set(multimat_tests_depends_on axom_utils slic slam multimat gtest)

#
# Add gtest based tests
#
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
