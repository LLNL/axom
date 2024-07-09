# Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Sina unit tests
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# Add gtest C++ tests
#------------------------------------------------------------------------------
set(gtest_sina_tests
    src/sina_CppBridge.cpp
    src/sina_File.cpp
   )

set(sina_gtests_depends_on gtest sina conduit::conduit)

foreach(test ${gtest_sina_tests})
    get_filename_component( test_name ${test} NAME_WE )
    axom_add_executable(NAME ${test_name}_test
                        SOURCES ${test}
                        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                        DEPENDS_ON ${sina_gtests_depends_on}
                        FOLDER axom/sina/tests
                        )

    axom_add_test( NAME    ${test_name}
                   COMMAND ${test_name}_test
                   )
endforeach()


if (ENABLE_GMOCK)

    #------------------------------------------------------------------------------
    # Create test utilities library for Sina
    #------------------------------------------------------------------------------

    # Define Sina test utility sources and headers
    set(sina_test_utils_sources
        src/ConduitTestUtils.cpp
        src/TestRecord.cpp
        )
    set(sina_test_utils_headers
        include/ConduitTestUtils.hpp
        include/TestRecord.hpp
        )

    set(sina_test_utils_depends_on gmock ${sina_gtests_depends_on})

    # Create a library for the test utilities so they can be used
    axom_add_library(
        NAME       sina_test_utils
        SOURCES    ${sina_test_utils_sources}
        HEADERS    ${sina_test_utils_headers}
        DEPENDS_ON ${sina_test_utils_depends_on}
        FOLDER     axom/sina/tests)

    #------------------------------------------------------------------------------
    # Add gmock C++ tests
    #------------------------------------------------------------------------------

    set(gmock_sina_tests
        src/sina_ConduitUtil.cpp
        src/sina_Curve.cpp
        src/sina_CurveSet.cpp
        src/sina_DataHolder.cpp
        src/sina_Datum.cpp
        src/sina_Document.cpp
        src/sina_ID.cpp
        src/sina_Record.cpp
        src/sina_Relationship.cpp
        src/sina_Run.cpp
        )

    # set(sina_gmock_depends_on ${sina_gtests_depends_on} gmock)

    # Add tests using Adiak if necessary and Adiak dependency
    blt_list_append( TO gmock_sina_tests ELEMENTS src/sina_AdiakWriter.cpp IF AXOM_SINA_USE_ADIAK)
    blt_list_append( TO sina_gmock_depends_on ELEMENTS adiak::adiak IF AXOM_SINA_USE_ADIAK )

    foreach(test ${gmock_sina_tests})
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
endif()

#------------------------------------------------------------------------------
# Add fortran integration test
#------------------------------------------------------------------------------
if (ENABLE_FORTRAN)
    find_package(PythonInterp REQUIRED)
    configure_file(${CMAKE_SOURCE_DIR}/axom/sina/tests/test_fortran_integration.py ${TEST_OUTPUT_DIRECTORY}/test_fortran_integration.py COPYONLY)
    configure_file(${CMAKE_SOURCE_DIR}/axom/sina/interface/src/sina_schema.json ${TEST_OUTPUT_DIRECTORY}/sina_schema.json COPYONLY)

    axom_add_test( NAME     sina_fortran_integration_test
                   COMMAND  ${PYTHON_EXECUTABLE} ${TEST_OUTPUT_DIRECTORY}/test_fortran_integration.py -bd ${PROJECT_BINARY_DIR})
endif()
