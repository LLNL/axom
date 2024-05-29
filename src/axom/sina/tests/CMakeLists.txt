# Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Sina unit tests
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Specify list of tests
#------------------------------------------------------------------------------

set(gtest_sina_tests
    src/sina_ConduitUtil.cpp
    src/sina_CppBridge.cpp
    src/sina_Curve.cpp
    src/sina_CurveSet.cpp
    src/sina_DataHolder.cpp
    src/sina_Datum.cpp
    src/sina_Document.cpp
    src/sina_File.cpp
    src/sina_ID.cpp
    src/sina_Record.cpp
    src/sina_Relationship.cpp
    src/sina_Run.cpp
   )

# Define Sina test utility sources and headers
set(sina_test_utils_sources
    src/ConduitTestUtils.cpp
    src/TestRecord.cpp
)
set(sina_test_utils_headers
    include/ConduitTestUtils.hpp
    include/TestRecord.hpp
)

# Define Sina test dependencies
set(sina_gtests_depends_on sina gtest gmock conduit::conduit)

# Add tests using Adiak if necessary and Adiak dependency
blt_list_append( TO gtest_sina_tests ELEMENTS src/sina_AdiakWriter.cpp IF AXOM_SINA_USE_ADIAK)
blt_list_append( TO sina_gtests_depends_on ELEMENTS adiak::adiak IF AXOM_SINA_USE_ADIAK )

# Create a library for the test utilities so they can be used
axom_add_library(
    NAME       sina_test_utils
    SOURCES    ${sina_test_utils_sources}
    HEADERS    ${sina_test_utils_headers}
    DEPENDS_ON ${sina_gtests_depends_on}
    FOLDER     axom/sina/tests)

#------------------------------------------------------------------------------
# Add gtest C++ tests
#------------------------------------------------------------------------------
foreach(test ${gtest_sina_tests})
    get_filename_component( test_name ${test} NAME_WE )
    axom_add_executable(NAME ${test_name}_test
                        SOURCES ${test}
                        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                        DEPENDS_ON sina_test_utils
                        FOLDER axom/sina/tests
                        )

    axom_add_test( NAME    ${test_name}
                   COMMAND ${test_name}_test
                   )
endforeach()