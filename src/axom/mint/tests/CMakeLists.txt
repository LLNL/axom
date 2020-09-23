# Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Tests for Mint component
#------------------------------------------------------------------------------

set( mint_tests

     ## execution model tests
     mint_execution_cell_traversals.cpp
     mint_execution_node_traversals.cpp
     mint_execution_face_traversals.cpp
     mint_execution_xargs.cpp

     ## mesh tests
     mint_mesh.cpp
     mint_mesh_cell_types.cpp
     mint_mesh_blueprint.cpp
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
     mint_su2_io.cpp
     mint_util_write_vtk.cpp
   )

## lists all axom depencies used for testing
set( mint_test_dependencies
     core
     slic
     mint
     )

blt_list_append( TO mint_test_dependencies
                 ELEMENTS sidre conduit::conduit
                 IF AXOM_MINT_USE_SIDRE )
blt_list_append( TO mint_test_dependencies ELEMENTS RAJA IF ${RAJA_FOUND} )
blt_list_append( TO mint_test_dependencies ELEMENTS openmp IF ${ENABLE_OPENMP} )
blt_list_append( TO mint_test_dependencies ELEMENTS cuda IF ${ENABLE_CUDA} )

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
                 NUM_OMP_THREADS 4
               )

endforeach()
