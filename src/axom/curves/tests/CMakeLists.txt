# Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

## list curves component tests
set( curves_tests

     curves_dummy_test.cpp
   )

## lists curces component test dependencies
set( curves_test_dependencies
     core
     slic
     curves
     gtest
     )


foreach( test ${curves_tests} )

   get_filename_component( test_name ${test} NAME_WE )

   blt_add_executable(
        NAME       ${test_name}_test
        SOURCES    ${test}
        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
        DEPENDS_ON ${curves_test_dependencies} gtest
        FOLDER     axom/curves/tests
        )

   blt_add_test( NAME    ${test_name}
                 COMMAND ${test_name}_test
               )

endforeach()