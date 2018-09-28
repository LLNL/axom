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

     ## core tests
     mint_core_array.cpp

     ## execution model tests
     mint_execution_cell_traversals.cpp
     mint_execution_forall_traversals.cpp
     mint_execution_node_traversals.cpp
     mint_execution_policy.cpp
     mint_execution_xargs.cpp

     ## mesh tests
     mint_mesh.cpp
     mint_mesh_blueprint.cpp
     mint_mesh_cell_types.cpp
     mint_mesh_connectivity_array.cpp
     mint_mesh_coordinates.cpp
     mint_mesh_curvilinear_mesh.cpp
     mint_mesh_field.cpp
     mint_mesh_field_data.cpp
     mint_mesh_field_types.cpp
     mint_mesh_field_variable.cpp
     mint_mesh_particle_mesh.cpp
     mint_mesh_rectilinear_mesh.cpp
     mint_mesh_uniform_mesh.cpp
     mint_mesh_unstructured_mesh.cpp
     mint_mesh_face_relation.cpp

     ## fem tests
     mint_fem_shape_functions.cpp
     mint_fem_single_fe.cpp

     ## util tests
     mint_util_write_vtk.cpp
   )

## lists all axom depencies used for testing
set( mint_test_dependencies
     core
     slic
     mint
     )

if ( ${AXOM_MINT_USE_SIDRE} )
  list(APPEND mint_test_dependencies sidre conduit)
endif()

if ( ${RAJA_FOUND} )
  list(APPEND mint_test_dependencies raja )
endif()

if ( ${ENABLE_OPENMP} )
  list(APPEND mint_test_dependencies openmp)
endif()


foreach( test ${mint_tests} )

   get_filename_component( test_name ${test} NAME_WE )

   blt_add_executable(
        NAME       ${test_name}_test
        SOURCES    ${test}
        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
        DEPENDS_ON ${mint_test_dependencies} gtest
        FOLDER     axom/mint/tests
        )

   blt_add_test( NAME    ${test_name}
                 COMMAND ${test_name}_test
               )

endforeach()
