# Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Tests for Inlet component
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Specify list of tests
#------------------------------------------------------------------------------

# Add Serial GTests based tests
set(gtest_inlet_tests
    inlet_LuaMap.cpp
    inlet_Structure.cpp )

foreach(test ${gtest_inlet_tests})
    get_filename_component( test_name ${test} NAME_WE )
    blt_add_executable( NAME       ${test_name}_test
                        SOURCES    ${test}
                        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                        DEPENDS_ON inlet gtest 
                        FOLDER axom/inlet/tests )
    blt_add_test( NAME    ${test_name} 
                  COMMAND ${test_name}_test )
endforeach()
