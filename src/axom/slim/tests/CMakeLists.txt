# Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Tests for SLIM component
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Specify list of tests
#------------------------------------------------------------------------------

# Add Serial GTests based tests
set(gtest_slim_tests
    slim_LuaMap.cpp
    )

foreach(test ${gtest_slim_tests})
    get_filename_component( test_name ${test} NAME_WE )
    blt_add_executable( NAME ${test_name}_test
                        SOURCES ${test}
                        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                        DEPENDS_ON slim gtest 
                        FOLDER axom/slim/tests )
    blt_add_test( NAME    ${test_name} 
                  COMMAND ${test_name}_test )
endforeach()


