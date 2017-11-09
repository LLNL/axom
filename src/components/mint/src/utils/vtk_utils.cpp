/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-741217
 *
 * All rights reserved.
 *
 * This file is part of Axom.
 *
 * For details about use and distribution, please read axom/LICENSE.
 *
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

#include "axom/Macros.hpp"            // for AXOM_DEBUG_VAR

#include "mint/vtk_utils.hpp"         // file header

#include "axom_utils/Utilities.hpp"   // for utilities::max
#include "mint/CellType.hpp"          // for cell::vtk_types
#include "mint/Field.hpp"             // for Field
#include "mint/FieldData.hpp"         // for FieldData
#include "mint/FieldTypes.hpp"        // for *_FIELD_TYPE
#include "mint/Mesh.hpp"              // for Mesh
#include "mint/MeshType.hpp"          // for MINT_*_*_MESH
#include "mint/RectilinearMesh.hpp"   // for RectilinearMesh
#include "mint/StructuredMesh.hpp"    // for StructuredMesh
#include "mint/UniformMesh.hpp"       // for UniformMesh
#include "slic/slic.hpp"              // for slic macros

#include <fstream>                     // for std::ofstream
#include <iomanip>                     // for std::setfill, std::setw
#include <limits>                      // for std::numeric_limits
#include <string>                      // for std::string

namespace axom
{
namespace mint
{

//------------------------------------------------------------------------------
//  INTERNAL HELPER METHODS
//------------------------------------------------------------------------------
namespace internal
{

/*!
 * \brief Writes mesh node locations to a VTK file using the legacy
 *  ASCII format.
 * \param [in] mesh the mesh whose nodes will be written.
 * \param [in] file the stream to write to.
 * \pre mesh != AXOM_NULLPTR
 */
void write_points( const Mesh * mesh, std::ofstream& file )
{
  SLIC_ASSERT( mesh != AXOM_NULLPTR );
  const localIndex num_nodes = mesh->getMeshNumberOfNodes();
  const int mesh_dim = mesh->getDimension();

  file << "POINTS " << num_nodes << " double\n";
  for ( localIndex nodeIdx = 0 ; nodeIdx < num_nodes ; ++nodeIdx )
  {
    file << mesh->getMeshNodeCoordinate( nodeIdx, 0 );
    for ( int dim = 1 ; dim < mesh_dim ; ++dim )
    {
      file << " " << mesh->getMeshNodeCoordinate( nodeIdx, dim );
    }
    for ( int dim = 0 ; dim < 3 - mesh_dim ; ++dim )
    {
      file << " 0.0";
    }
    file << std::endl;
  }
}

/*!
 * \brief Writes mesh cell connectivity and type to a VTK file
 *  using the legacy ASCII format.
 * \param [in] mesh the mesh whose cells will be written.
 * \param [in] file the stream to write to.
 * \pre mesh != AXOM_NULLPTR
 */
void write_cells( const Mesh * mesh, std::ofstream& file )
{
  SLIC_ASSERT( mesh != AXOM_NULLPTR );
  const localIndex num_cells = mesh->getMeshNumberOfCells();

  /* First need to get total size of the connectivity array. */
  /* If the mesh only has one cell type we can calculate this directly. */
  int max_cell_nodes = mesh->getMeshNumberOfCellNodes( 0 );
  localIndex total_size = ( max_cell_nodes + 1 ) * num_cells;

  /* If the mesh has mixed cells then we need to loop over the elements. */
  if ( mesh->getMeshType() == MINT_UNSTRUCTURED_MIXED_ELEMENT_MESH )
  {
    total_size = num_cells;
    for ( localIndex cellIdx = 0 ; cellIdx < num_cells ; ++cellIdx )
    {
      const int num_cell_nodes = mesh->getMeshNumberOfCellNodes( cellIdx );
      max_cell_nodes = utilities::max(num_cell_nodes, max_cell_nodes);
      total_size += num_cell_nodes;
    }
  }
  file << "CELLS " << num_cells << " " << total_size << std::endl;

  /* Write out the mesh cell connectivity. */
  localIndex cell_nodes[ max_cell_nodes ];
  for ( localIndex cellIdx = 0 ; cellIdx < num_cells ; ++cellIdx )
  {
    const int num_cell_nodes = mesh->getMeshNumberOfCellNodes( cellIdx );
    mesh->getMeshCell( cellIdx, cell_nodes );

    file << num_cell_nodes;
    for ( int i = 0 ; i < num_cell_nodes ; ++i )
    {
      file << " " << cell_nodes[ i ];
    }
    file << std::endl;
  }

  /* Write out the mesh cell types. */
  file << "CELL_TYPES " << num_cells << std::endl;
  for ( localIndex cellIdx = 0 ; cellIdx < num_cells ; ++cellIdx )
  {
    int cell_type = mesh->getMeshCellType( cellIdx );
    file << cell::vtk_types[ cell_type ] << std::endl;
  }
}

/*!
 * \brief Writes the dimensions of a structured mesh to a VTK file using the
 *  legacy ASCII format.
 * \param [in] mesh the structured mesh to write out.
 * \param [in] file the stream to write to.
 * \pre mesh != AXOM_NULLPTR
 */
void write_dimensions( const StructuredMesh * mesh, std::ofstream& file )
{
  SLIC_ASSERT( mesh != AXOM_NULLPTR );

  localIndex ext[ 3 ];
  mesh->getExtentSize( ext );
  file << "DIMENSIONS ";
  file << ext[ 0 ] << " " << ext[ 1 ] << " " << ext[ 2 ] << std::endl;
}

/*!
 * \brief Writes a rectilinear mesh to a VTK file using the legacy
 *  ASCII format.
 * \param [in] mesh the rectilinear mesh to write out.
 * \param [in] file the stream to write to.
 * \pre mesh != AXOM_NULLPTR
 */
void write_rectilinear_mesh( const RectilinearMesh * mesh, std::ofstream& file )
{
  SLIC_ASSERT( mesh != AXOM_NULLPTR );

  write_dimensions( mesh, file );

  localIndex ext[ 3 ];
  mesh->getExtentSize( ext );
  std::string coord_names[3] = { "X_COORDINATES ", "Y_COORDINATES ",
                                 "Z_COORDINATES " };

  for ( int dim = 0 ; dim < mesh->getDimension() ; ++dim )
  {
    file << coord_names[ dim ] << ext[ dim ] << " double\n";
    const double * coords = mesh->getCoordinateArray( dim );
    file << coords[0];
    for (globalIndex i = 1 ; i < ext[ dim ] ; ++i )
    {
      file << " " << coords[i];
    }
    file << std::endl;
  }
  for ( int dim = mesh->getDimension() ; dim < 3 ; ++dim )
  {
    file << coord_names[ dim ] << ext[ dim ] << " double\n";
    file << 0.0 << std::endl;
  }
}

/*!
 * \brief Writes a uniform mesh to a VTK file using the legacy
 *  ASCII format.
 * \param [in] mesh the uniform mesh to write out.
 * \param [in] file the stream to write to.
 * \pre mesh != AXOM_NULLPTR
 */
void write_uniform_mesh( const UniformMesh * mesh, std::ofstream& file )
{
  SLIC_ASSERT( mesh != AXOM_NULLPTR );

  write_dimensions( mesh, file );

  double temp[3];
  mesh->getOrigin(temp);
  file << "ORIGIN ";
  file << temp[ 0 ] << " " << temp[ 1 ] << " " << temp[ 2 ] << std::endl;

  mesh->getSpacing(temp);
  file << "SPACING ";
  file << temp[ 0 ] << " " << temp[ 1 ] << " " << temp[ 2 ] << std::endl;
}

/*!
 * \brief Writes a scalar field to a VTK file using the legacy
 *  ASCII format.
 * \param [in] field the scalar field to write out.
 * \param [in] file the stream to write to.
 * \pre field != AXOM_NULLPTR
 * \pre field->getNumComponents() == 1
 */
void write_scalar_data( const Field * field, std::ofstream& file )
{
  SLIC_ASSERT(  field != AXOM_NULLPTR );
  SLIC_ASSERT(  field->getNumComponents() == 1 );
  const localIndex num_values = field->getNumTuples();

  file << "SCALARS " << field->getName() << " ";
  if ( field->getType() == DOUBLE_FIELD_TYPE )
  {
    file << "double\n";
    file << "LOOKUP_TABLE default\n";

    const double * data_ptr = field->getDoublePtr();
    SLIC_ASSERT( data_ptr != AXOM_NULLPTR );

    for ( localIndex i = 0 ; i < num_values ; ++i )
    {
      file << data_ptr[ i ] << std::endl;
    }
  }
  else if ( field->getType() == INTEGER_FIELD_TYPE )
  {
    file << "int\n";
    file << "LOOKUP_TABLE default\n";

    const int * data_ptr = field->getIntPtr();
    SLIC_ASSERT( data_ptr != AXOM_NULLPTR );

    for ( localIndex i = 0 ; i < num_values ; ++i )
    {
      file << data_ptr[ i ] << std::endl;
    }
  }
}

/*!
 * \brief Writes a vector field to a VTK file using the legacy
 *  ASCII format.
 * \param [in] field the vector field to write out.
 * \param [in] file the stream to write to.
 * \pre field != AXOM_NULLPTR
 * \pre field->getNumComponents() == 2 || field->getNumComponents() == 3
 */
void write_vector_data( const Field * field, std::ofstream& file )
{
  SLIC_ASSERT(  field != AXOM_NULLPTR );
  const int num_components = field->getNumComponents();
  const localIndex num_values = field->getNumTuples();
  SLIC_ASSERT(  num_components == 2 || num_components == 3 );

  file << "VECTORS " << field->getName() << " ";
  if ( field->getType() == DOUBLE_FIELD_TYPE )
  {
    file << "double\n";

    const double * data_ptr = field->getDoublePtr();
    SLIC_ASSERT( data_ptr != AXOM_NULLPTR );

    for ( localIndex i = 0 ; i < num_values ; ++i )
    {
      file << data_ptr[ num_components * i + 0 ] << " ";
      file << data_ptr[ num_components * i + 1 ] << " ";
      if ( num_components == 2 )
      {
        file << 0.0 << std::endl;
      }
      else
      {
        file << data_ptr[ num_components * i + 2 ] << std::endl;
      }
    }
  }
  else if ( field->getType() == INTEGER_FIELD_TYPE )
  {
    file << "int\n";

    const int * data_ptr = field->getIntPtr();
    SLIC_ASSERT( data_ptr != AXOM_NULLPTR );

    for ( localIndex i = 0 ; i < num_values ; ++i )
    {
      file << data_ptr[ num_components * i + 0 ] << " ";
      file << data_ptr[ num_components * i + 1 ] << " ";
      if ( num_components == 2 )
      {
        file << 0 << std::endl;
      }
      else
      {
        file << data_ptr[ num_components * i + 2 ] << std::endl;
      }
    }
  }
}

/*!
 * \brief Writes a multidimensional field to a VTK file using the legacy
 *  ASCII format.
 * \param [in] field the multidimensional field to write out.
 * \param [in] file the stream to write to.
 * \pre field != AXOM_NULLPTR
 * \pre field->getNumComponents > 3
 */
void write_multidim_data( const Field * field, std::ofstream& file )
{
  SLIC_ASSERT( field != AXOM_NULLPTR );
  const int field_type = field->getType();
  const int num_components = field->getNumComponents();
  const localIndex num_values = field->getNumTuples();
  SLIC_ASSERT( num_components > 3 );

  if ( field_type == DOUBLE_FIELD_TYPE )
  {
    for ( int cur_comp = 0 ; cur_comp < num_components ; ++cur_comp )
    {
      file << "SCALARS " << field->getName() << "_";
      file << std::setfill('0') << std::setw(3) << cur_comp;
      file << " double\n";
      file << "LOOKUP_TABLE default\n";

      const double * data_ptr = field->getDoublePtr();
      SLIC_ASSERT( data_ptr != AXOM_NULLPTR );

      for ( localIndex i = 0 ; i < num_values ; ++i )
      {
        file << data_ptr[ num_components * i + cur_comp ] << std::endl;
      }
    }
  }
  else if ( field_type == INTEGER_FIELD_TYPE )
  {
    for ( int cur_comp = 0 ; cur_comp < num_components ; ++cur_comp )
    {
      file << "SCALARS " << field->getName() << "_";
      file << std::setfill('0') << std::setw(3) << cur_comp;
      file << " int\n";
      file << "LOOKUP_TABLE default\n";

      const int * data_ptr = field->getIntPtr();
      SLIC_ASSERT( data_ptr != AXOM_NULLPTR );

      for ( localIndex i = 0 ; i < num_values ; ++i )
      {
        file << data_ptr[ num_components * i + cur_comp ] << std::endl;
      }
    }
  }
}

/*!
 * \brief Writes mesh FieldData to a VTK file using the legacy ASCII format.
 * \param [in] field_data the data to write out.
 * \param [in] num_values the number of tuples each field is expected to have.
 * \param [in] file the stream to write to.
 * \pre field_data != AXOM_NULLPTR
 */
void write_data( const FieldData& field_data, localIndex num_values,
                 std::ofstream& file )
{
  for ( int i = 0 ; i < field_data.getNumberOfFields() ; ++i )
  {
    const Field * field = field_data.getField( i );
    SLIC_ASSERT( field != AXOM_NULLPTR );
    const int num_components = field->getNumComponents();
    SLIC_ASSERT( field->getNumTuples() == num_values );

    if ( field->getType() != DOUBLE_FIELD_TYPE &&
         field->getType() != INTEGER_FIELD_TYPE )
    {
      SLIC_WARNING( "Field " << field->getName() <<
                    " type not double or integer." );
      continue;
    }

    if ( num_components == 1 )
    {
      write_scalar_data( field, file );
    }
    else if ( num_components == 2 || num_components == 3 )
    {
      write_vector_data( field, file );
    }
    else if ( num_components > 3 )
    {
      write_multidim_data( field, file );
    }
    else
    {
      SLIC_WARNING( "Field has an improper number of components.");
    }
  }
}

} /* namespace internal */

//------------------------------------------------------------------------------
int write_vtk( const Mesh * mesh, const std::string& file_path )
{
  SLIC_ASSERT( mesh != AXOM_NULLPTR );
  int mesh_type = mesh->getMeshType();

  std::ofstream file( file_path.c_str() );
  if ( !file.good() )
  {
    SLIC_WARNING( "Could not open file at path " << file_path );
    return -1;
  }

  file.setf(file.scientific);
#if __cplusplus >= 201103L
  file.precision(  std::numeric_limits< double >::max_digits10 );
#else
  file.precision(  std::numeric_limits< double >::digits10 + 2 );
#endif

  /* Write the VTK header */
  file << "# vtk DataFile Version 3.0\n";
  file << "Mesh generated by axom::mint::write_vtk\n";
  file << "ASCII\n";

  /* Write out the mesh node and cell coordinates. */
  if ( mesh_type == MINT_UNSTRUCTURED_SEGMENT_MESH ||
       mesh_type == MINT_UNSTRUCTURED_TRIANGLE_MESH ||
       mesh_type == MINT_UNSTRUCTURED_QUAD_MESH ||
       mesh_type == MINT_UNSTRUCTURED_TET_MESH ||
       mesh_type == MINT_UNSTRUCTURED_HEX_MESH ||
       mesh_type == MINT_UNSTRUCTURED_MIXED_ELEMENT_MESH ||
       mesh_type == MINT_PARTICLE_MESH )
  {
    file << "DATASET UNSTRUCTURED_GRID\n";
    internal::write_points( mesh, file );
    internal::write_cells( mesh, file );
  }
  else if ( mesh_type == MINT_STRUCTURED_CURVILINEAR_MESH )
  {
    file << "DATASET STRUCTURED_GRID\n";
    const StructuredMesh * struc_mesh =
      dynamic_cast< const StructuredMesh * >( mesh );
    internal::write_dimensions( struc_mesh, file );
    internal::write_points( struc_mesh, file );
  }
  else if ( mesh_type == MINT_STRUCTURED_RECTILINEAR_MESH )
  {
    file << "DATASET RECTILINEAR_GRID\n";
    const RectilinearMesh * rect_mesh =
      dynamic_cast< const RectilinearMesh * >( mesh );
    internal::write_rectilinear_mesh( rect_mesh, file );
  }
  else if ( mesh_type == MINT_STRUCTURED_UNIFORM_MESH )
  {
    file << "DATASET STRUCTURED_POINTS\n";
    const UniformMesh * uniform_mesh =
      dynamic_cast< const UniformMesh * >( mesh );
    internal::write_uniform_mesh( uniform_mesh, file );
  }
  else
  {
    SLIC_WARNING( "Mesh does not have a proper type (" << mesh_type << ") " <<
                  "write aborted." );
    file.close();
    remove(file_path.c_str());
    return -1;
  }

  /* Write out the node data if any. */
  const localIndex num_nodes = mesh->getMeshNumberOfNodes();
  const FieldData& node_data = mesh->getNodeFieldData();
  if ( node_data.getNumberOfFields() > 0 )
  {
    file << "POINT_DATA " << num_nodes << std::endl;
    internal::write_data( node_data, num_nodes, file );
  }

  /* Write out the cell data if any. */
  const localIndex num_cells = mesh->getMeshNumberOfCells();
  const FieldData& cell_data = mesh->getCellFieldData();
  if ( cell_data.getNumberOfFields() > 0 )
  {
    file << "CELL_DATA " << num_cells << std::endl;
    internal::write_data( cell_data, num_cells, file );
  }

  file.close();
  return 0;
}

} /* namespace mint */
} /* namespace axom */
