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
if (SOL_FOUND)
        set(gtest_inlet_tests inlet_Reader.cpp inlet_object.cpp inlet_Inlet.cpp)
        blt_add_library(
                NAME        inlet_test_utils
                SOURCES     inlet_test_utils.cpp
                HEADERS     inlet_test_utils.hpp
                DEPENDS_ON  inlet gtest
                )
endif()

foreach(test ${gtest_inlet_tests})
    get_filename_component( test_name ${test} NAME_WE )
    blt_add_executable( NAME       ${test_name}_test
                        SOURCES    ${test}
                        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                        DEPENDS_ON inlet gtest inlet_test_utils
                        FOLDER     axom/inlet/tests )
    blt_add_test( NAME    ${test_name} 
                    COMMAND ${test_name}_test )
endforeach()
