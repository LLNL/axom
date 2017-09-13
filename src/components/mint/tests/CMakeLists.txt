#-------------------------------------------------------------------------------
# Copyright (c) 2015, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
#
# All rights reserved.
#
# This source code cannot be distributed without permission and further
# review from Lawrence Livermore National Laboratory.
#-------------------------------------------------------------------------------

set( mint_tests
     mint_mesh_extent.cpp
     mint_shape_functions.cpp
     mint_single_fe.cpp
     mint_write_vtk.cpp
   )

foreach( test ${mint_tests} )

   get_filename_component( test_name ${test} NAME_WE )

   blt_add_executable(
        NAME       ${test_name}_test
        SOURCES    ${test}
        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
        DEPENDS_ON axom_utils slic mint gtest
        )

   blt_add_test( NAME    ${test_name}
                 COMMAND ${test_name}_test
               )

endforeach()
