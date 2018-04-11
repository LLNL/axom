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
#include "axom_utils/Matrix.hpp"      // for numerics::Matrix

#include "mint/CellTypes.hpp"          // for cell::vtk_types
#include "mint/Field.hpp"             // for Field
#include "mint/FieldData.hpp"         // for FieldData
#include "mint/FieldTypes.hpp"        // for *_FIELD_TYPE
#include "mint/FiniteElement.hpp"     // for mint::FiniteElement
#include "mint/Mesh.hpp"              // for Mesh
#include "mint/MeshTypes.hpp"          // for MINT_*_*_MESH
#include "mint/ParticleMesh.hpp"      // for ParticleMesh
#include "mint/RectilinearMesh.hpp"   // for RectilinearMesh
#include "mint/StructuredMesh.hpp"    // for StructuredMesh
#include "mint/UniformMesh.hpp"       // for UniformMesh

#include "slic/slic.hpp"              // for slic macros

// C/C++ includes
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
 * \brief Computes the maximum number of nodes of a cell on the given mesh.
 *
 * \param [in] mesh the mesh object
 * \return max_cell_nodes the max number of nodes of a cell
 *
 * \pre mesh != AXOM_NULLPTR
 * \post max_cell_nodes >= 1
 */
int get_max_cell_nodes( const Mesh* mesh )
{
  SLIC_ASSERT( mesh != AXOM_NULLPTR );

  int max_cell_nodes = 0;

  if( !mesh->hasMixedCellTypes( ) )
  {
    // short-circuit
    max_cell_nodes = mint::cell_info[ mesh->getMeshCellType() ].num_nodes;
    return max_cell_nodes;
  }

  const mint::IndexType numCells = mesh->getNumberOfCells();
  for ( mint::IndexType icell=0; icell < numCells; ++icell )
  {
    CellType cell_type  = mesh->getMeshCellType( icell );
    const int num_nodes = mint::cell_info[ cell_type ].num_nodes;
    if ( num_nodes > max_cell_nodes )
    {
      max_cell_nodes = num_nodes;
    }

  }

  return ( max_cell_nodes );
}

/*!
 * \brief Writes mesh node locations to a VTK file in legacy ASCII format.
 * \param [in] mesh the mesh whose nodes will be written.
 * \param [in] file the stream to write to.
 * \pre mesh != AXOM_NULLPTR
 */
void write_points( const Mesh* mesh, std::ofstream& file )
{
  SLIC_ASSERT( mesh != AXOM_NULLPTR );
  const IndexType num_nodes = mesh->getNumberOfNodes();
  const int mesh_dim = mesh->getDimension();

  const double* x = mesh->getCoordinateArray( X_COORDINATE );
  SLIC_ASSERT( x != AXOM_NULLPTR );

  const double* y = ( mesh_dim > 1 ) ?
      mesh->getCoordinateArray( Y_COORDINATE ) : AXOM_NULLPTR;
  const double* z = ( mesh_dim > 2 ) ?
      mesh->getCoordinateArray( Z_COORDINATE ) : AXOM_NULLPTR;

  file << "POINTS " << num_nodes << " double\n";
  for ( IndexType nodeIdx = 0 ; nodeIdx < num_nodes ; ++nodeIdx )
  {
    double xx = x[ nodeIdx ];
    double yy = ( y != AXOM_NULLPTR ) ? y[ nodeIdx ] : 0.0;
    double zz = ( z != AXOM_NULLPTR ) ? z[ nodeIdx ] : 0.0;
    file << xx << " " << yy << " " << zz << std::endl;
  }

}

/*!
 * \brief Writes mesh cell connectivity and type to a VTK file
 *  using the legacy ASCII format.
 * \param [in] mesh the mesh whose cells will be written.
 * \param [in] file the stream to write to.
 * \pre mesh != AXOM_NULLPTR
 */
void write_cells( const Mesh* mesh, std::ofstream& file )
{
  SLIC_ASSERT( mesh != AXOM_NULLPTR );
  const IndexType num_cells = mesh->getNumberOfCells();

  /* First need to get total size of the connectivity array. */
  /* If the mesh only has one cell type we can calculate this directly. */
  int max_cell_nodes = get_max_cell_nodes( mesh );
  IndexType total_size = ( max_cell_nodes + 1 ) * num_cells;

  file << "CELLS " << num_cells << " " << total_size << std::endl;

  /* Write out the mesh cell connectivity. */
  IndexType cell_nodes[ max_cell_nodes ];
  for ( IndexType cellIdx = 0 ; cellIdx < num_cells ; ++cellIdx )
  {
    int num_nodes = mint::cell_info[mesh->getMeshCellType(cellIdx)].num_nodes;
    mesh->getMeshCell( cellIdx, cell_nodes );

    file << num_nodes;
    for ( int i = 0 ; i < num_nodes ; ++i )
    {
      file << " " << cell_nodes[ i ];
    }
    file << std::endl;
  }

  /* Write out the mesh cell types. */
  file << "CELL_TYPES " << num_cells << std::endl;
  for ( IndexType cellIdx = 0 ; cellIdx < num_cells ; ++cellIdx )
  {
    CellType cell_type = mesh->getMeshCellType( cellIdx );
    file << mint::cell_info[ cell_type ].vtk_type << std::endl;
  }

}

/*!
 * \brief Writes the dimensions of a structured mesh to a VTK file using the
 *  legacy ASCII format.
 * \param [in] mesh the structured mesh to write out.
 * \param [in] file the stream to write to.
 * \pre mesh != AXOM_NULLPTR
 */
void write_dimensions( const StructuredMesh* mesh, std::ofstream& file )
{
  SLIC_ASSERT( mesh != AXOM_NULLPTR );

  IndexType ext[ 3 ];
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
void write_rectilinear_mesh( const RectilinearMesh* mesh, std::ofstream& file )
{
  SLIC_ASSERT( mesh != AXOM_NULLPTR );

  write_dimensions( mesh, file );

  IndexType ext[ 3 ];
  mesh->getExtentSize( ext );
  std::string coord_names[3] = { "X_COORDINATES ", "Y_COORDINATES ",
                                 "Z_COORDINATES " };

  for ( int dim = 0 ; dim < mesh->getDimension() ; ++dim )
  {
    file << coord_names[ dim ] << ext[ dim ] << " double\n";
    const double* coords = mesh->getCoordinateArray( dim );
    file << coords[0];
    for (int64 i = 1 ; i < ext[ dim ] ; ++i )
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
void write_uniform_mesh( const UniformMesh* mesh, std::ofstream& file )
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
void write_scalar_data( const Field* field, std::ofstream& file )
{
  SLIC_ASSERT(  field != AXOM_NULLPTR );
  SLIC_ASSERT(  field->getNumComponents() == 1 );
  const IndexType num_values = field->getNumTuples();

  file << "SCALARS " << field->getName() << " ";
  if ( field->getType() == DOUBLE_FIELD_TYPE )
  {
    file << "double\n";
    file << "LOOKUP_TABLE default\n";

    const double* data_ptr = Field::getDataPtr< double >( field );
    SLIC_ASSERT( data_ptr != AXOM_NULLPTR );

    for ( IndexType i = 0 ; i < num_values ; ++i )
    {
      file << data_ptr[ i ] << std::endl;
    }
  }
  else if ( field->getType() == INT32_FIELD_TYPE )
  {
    file << "int\n";
    file << "LOOKUP_TABLE default\n";

    const int* data_ptr = Field::getDataPtr< int >( field );
    SLIC_ASSERT( data_ptr != AXOM_NULLPTR );

    for ( IndexType i = 0 ; i < num_values ; ++i )
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
void write_vector_data( const Field* field, std::ofstream& file )
{
  SLIC_ASSERT(  field != AXOM_NULLPTR );
  const int num_components = field->getNumComponents();
  const IndexType num_values = field->getNumTuples();
  SLIC_ASSERT(  num_components == 2 || num_components == 3 );

  file << "VECTORS " << field->getName() << " ";
  if ( field->getType() == DOUBLE_FIELD_TYPE )
  {
    file << "double\n";

//    const double* data_ptr = field->getDoublePtr();
    const double* data_ptr = Field::getDataPtr< double >( field );
    SLIC_ASSERT( data_ptr != AXOM_NULLPTR );

    for ( IndexType i = 0 ; i < num_values ; ++i )
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
  else if ( field->getType() == INT32_FIELD_TYPE )
  {
    file << "int\n";

    const int* data_ptr = Field::getDataPtr< int >( field );
    SLIC_ASSERT( data_ptr != AXOM_NULLPTR );

    for ( IndexType i = 0 ; i < num_values ; ++i )
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
void write_multidim_data( const Field* field, std::ofstream& file )
{
  SLIC_ASSERT( field != AXOM_NULLPTR );
  const int field_type = field->getType();
  const int num_components = field->getNumComponents();
  const IndexType num_values = field->getNumTuples();
  SLIC_ASSERT( num_components > 3 );

  if ( field_type == DOUBLE_FIELD_TYPE )
  {
    for ( int cur_comp = 0 ; cur_comp < num_components ; ++cur_comp )
    {
      file << "SCALARS " << field->getName() << "_";
      file << std::setfill('0') << std::setw(3) << cur_comp;
      file << " double\n";
      file << "LOOKUP_TABLE default\n";

//      const double* data_ptr = field->getDoublePtr();
      const double* data_ptr = Field::getDataPtr< double >( field );
      SLIC_ASSERT( data_ptr != AXOM_NULLPTR );

      for ( IndexType i = 0 ; i < num_values ; ++i )
      {
        file << data_ptr[ num_components * i + cur_comp ] << std::endl;
      }
    }
  }
  else if ( field_type == INT32_FIELD_TYPE )
  {
    for ( int cur_comp = 0 ; cur_comp < num_components ; ++cur_comp )
    {
      file << "SCALARS " << field->getName() << "_";
      file << std::setfill('0') << std::setw(3) << cur_comp;
      file << " int\n";
      file << "LOOKUP_TABLE default\n";

      const int* data_ptr = Field::getDataPtr< int >( field );
      SLIC_ASSERT( data_ptr != AXOM_NULLPTR );

      for ( IndexType i = 0 ; i < num_values ; ++i )
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
void write_data( const FieldData* field_data, IndexType num_values,
                 std::ofstream& file )
{
  const int numFields = field_data->getNumFields( );
  for ( int i = 0 ; i < numFields ; ++i )
  {
    const Field* field = field_data->getField( i );
    SLIC_ASSERT( field != AXOM_NULLPTR );
    const int num_components = field->getNumComponents();
    SLIC_ASSERT( field->getNumTuples() == num_values );

    const bool invalidType = ( field->getType() >= NUMBER_OF_FIELD_TYPES );
    SLIC_ERROR_IF( invalidType,
                   "Field [" << field->getName() << "] has invalid type" );

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
int write_vtk( const Mesh* mesh, const std::string& file_path )
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
  if ( mesh_type == mint::UNSTRUCTURED_MESH ||
       mesh_type == mint::PARTICLE_MESH )
  {
    file << "DATASET UNSTRUCTURED_GRID\n";
    internal::write_points( mesh, file );
    internal::write_cells( mesh, file );
  }
  else if ( mesh_type == mint::STRUCTURED_MESH )
  {
    file << "DATASET STRUCTURED_GRID\n";
    const StructuredMesh* struc_mesh =
      dynamic_cast< const StructuredMesh* >( mesh );
    internal::write_dimensions( struc_mesh, file );
    internal::write_points( struc_mesh, file );
  }
  else if ( mesh_type == mint::RECTILINEAR_MESH )
  {
    file << "DATASET RECTILINEAR_GRID\n";
    const RectilinearMesh* rect_mesh =
      dynamic_cast< const RectilinearMesh* >( mesh );
    internal::write_rectilinear_mesh( rect_mesh, file );
  }
  else if ( mesh_type == mint::UNIFORM_MESH )
  {
    file << "DATASET STRUCTURED_POINTS\n";
    const UniformMesh* uniform_mesh =
      dynamic_cast< const UniformMesh* >( mesh );
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
  const IndexType num_nodes = mesh->getNumberOfNodes();
  const FieldData* node_data = mesh->getFieldData( mint::NODE_CENTERED );
  if ( node_data->getNumFields() > 0 )
  {
    file << "POINT_DATA " << num_nodes << std::endl;
    internal::write_data( node_data, num_nodes, file );
  }

  /* Write out the cell data if any. */
  if ( mesh->getMeshType() != mint::PARTICLE_MESH )
  {
    const IndexType num_cells = mesh->getNumberOfCells();
    const FieldData* cell_data = mesh->getFieldData( mint::CELL_CENTERED );
    if ( cell_data->getNumFields() > 0 )
    {
      file << "CELL_DATA " << num_cells << std::endl;
      internal::write_data( cell_data, num_cells, file );
    }
  }

  file.close();
  return 0;
}

//------------------------------------------------------------------------------
int write_vtk( mint::FiniteElement& fe, const std::string& file_path )
{
  if ( file_path.empty() )
  {
    return -1;
  }

  std::ofstream ofs ( file_path.c_str( ) );
  if ( !ofs.is_open() )
  {
    SLIC_WARNING( "Could not open file at path " << file_path );
    return -1;
  }

  ofs.setf( ofs.scientific );
#if __cplusplus >= 201103L
  ofs.precision(  std::numeric_limits< double >::max_digits10 );
#else
  ofs.precision(  std::numeric_limits< double >::digits10 + 2 );
#endif

  const bool zero_copy = true;
  const int cell_type  = fe.getCellType( );
  const int ndims      = fe.getPhysicalDimension( );
  const int nnodes     = fe.getNumNodes( );
  numerics::Matrix< double > nodes( ndims, nnodes, fe.getPhysicalNodes(),
                                    zero_copy );

  // write VTK header
  ofs << "# vtk DataFile Version 3.0\n";
  ofs << " FiniteElement\n";
  ofs << "ASCII\n";
  ofs << "DATASET UNSTRUCTURED_GRID\n";
  ofs << "POINTS " << nnodes << " double\n";

  // write the cell coordinates
  for ( int i=0; i < nnodes; ++i )
  {
    const double* pt = nodes.getColumn( i );
    const double x   = pt[ 0 ];
    const double y   = ( ndims > 1 ) ? pt[ 1 ] : 0.0;
    const double z   = ( ndims > 2 ) ? pt[ 2 ] : 0.0;

    ofs << x << " " << y << " " << z << std::endl;

  } // END for all nodes

  // write cell connectivity
  ofs << "CELLS 1 " << nnodes+1 << std::endl;
  ofs << nnodes << " ";
  for ( int i=0; i < nnodes; ++i ) {
    ofs << i << " ";
  }
  ofs << std::endl;

  // write cell type information
  ofs << "CELL_TYPES 1\n";
  ofs << mint::cell_info[ cell_type ].vtk_type << std::endl;

  // close the file
  ofs.close();

  // return success
  return 0;
}

} /* namespace mint */
} /* namespace axom */
