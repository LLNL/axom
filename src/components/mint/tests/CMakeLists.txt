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

set( mint_tests
     mint_array.cpp
     mint_mesh_extent.cpp
     mint_shape_functions.cpp
     mint_single_fe.cpp
     mint_write_vtk.cpp
# turn off for now
#     mint_connectivity.cpp
   )

foreach( test ${mint_tests} )

   get_filename_component( test_name ${test} NAME_WE )

   blt_add_executable(
        NAME       ${test_name}_test
        SOURCES    ${test}
        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
        DEPENDS_ON axom_utils slic mint gtest
        FOLDER     axom/mint/tests
        )

   blt_add_test( NAME    ${test_name}
                 COMMAND ${test_name}_test
               )

endforeach()
