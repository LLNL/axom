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

// Axom utils
#include "axom_utils/Utilities.hpp"     /* for utilities::max */

// Mint includes
#include "mint/config.hpp"              /* for IndexType, int64 */
#include "mint/CellTypes.hpp"           /* for cell::vtk_types */
#include "mint/CurvilinearMesh.hpp"     /* for CurvilinearMesh */
#include "mint/Field.hpp"               /* for Field */
#include "mint/FieldData.hpp"           /* for FieldData */
#include "mint/FieldTypes.hpp"          /* for *_FIELD_TYPE */
#include "mint/FieldVariable.hpp"       /* for FieldVariable */
#include "mint/Mesh.hpp"                /* for Mesh */
#include "mint/ParticleMesh.hpp"        /* for ParticleMesh */
#include "mint/RectilinearMesh.hpp"     /* for RectilinearMesh */
#include "mint/UniformMesh.hpp"         /* for UniformMesh */
#include "mint/UnstructuredMesh.hpp"    /* for UnstructuredMesh */
#include "mint/vtk_utils.hpp"           /* for write_vtk */

// Slic includes
#include "slic/slic.hpp"                /* for slic macros */
#include "slic/UnitTestLogger.hpp"      /* for UnitTestLogger */

// C/C++ includes
#include <cmath>                        /* for std::exp */
#include <cstdio>                       /* for std::remove */
#include <fstream>                      /* for std::ifstream */
#include <iomanip>                      /* for std::setfill, std::setw */
#include <string>                       /* for std::string */
#include <sstream>                      /* for std::stringstream */
#include <set>                          /* for std::set */

// gtest includes
#include "gtest/gtest.h"                /* for TEST and EXPECT_* macros */

#ifndef DELETE_VTK_FILES
  #define DELETE_VTK_FILES 1
#endif

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
 * \brief Creates artificial 1 dimensional mesh data.
 * \param [in] mesh the mesh to populate.
 * \pre mesh != AXOM_NULLPTR
 */
void create_scalar_data( Mesh* mesh )
{
  const IndexType mesh_num_nodes = mesh->getNumberOfNodes();
  const IndexType mesh_num_cells = mesh->getNumberOfCells();

  double* double_ptr = AXOM_NULLPTR;
  int* int_ptr    = AXOM_NULLPTR;

  double_ptr = mesh->createField< double >( "node_scalars_double",
                                            mint::NODE_CENTERED );
  int_ptr = mesh->createField< int >( "node_scalars_int",
                                      mint::NODE_CENTERED );

  for ( int idx = 0 ; idx < mesh_num_nodes ; ++idx )
  {
    double coords[] = { 0, 0, 0 };
    mesh->getNode( idx, coords );

    const double x = coords[0];
    const double y = coords[1];
    const double z = coords[2];

    double r2 = (x * x) + (y * y) + (z * z);
    r2 = std::exp(-r2 / 100.0);
    double temp = (x * x) - 4 * y + (2 * y * z) - (5 * x * y * z + 1) * r2;
    double_ptr[ idx ] = temp;
    int_ptr[ idx ] = static_cast<int>(temp);
  }

  /* Particle meshes can only have node centered fields. */
  if ( mesh->getMeshType() == PARTICLE_MESH )
  {
    return;
  }

  double_ptr = mesh->createField< double >( "cell_scalars_double",
                                            mint::CELL_CENTERED  );

  int_ptr    = mesh->createField< int >( "cell_scalars_int",
                                         mint::CELL_CENTERED );

  for ( int idx = 0 ; idx < mesh_num_cells ; ++idx )
  {
    double_ptr[ idx ] = idx;
    int_ptr[ idx ] = idx;
  }

}

/*!
 * \brief Creates artificial 3 dimensional mesh data.
 * \param [in] mesh the mesh to populate.
 * \pre mesh != AXOM_NULLPTR
 */
void create_vector_data( Mesh* mesh )
{
  const IndexType num_nodes = mesh->getNumberOfNodes();
  const IndexType num_cells = mesh->getNumberOfCells();

  double* double_ptr3 = AXOM_NULLPTR;
  double* double_ptr2 = AXOM_NULLPTR;
  int* int_ptr3    = AXOM_NULLPTR;
  int* int_ptr2    = AXOM_NULLPTR;

  double_ptr3 = mesh->createField< double >( "node_vectors_3double",
                                             mint::NODE_CENTERED, 3 );
  int_ptr3    = mesh->createField< int >( "node_vectors_3int",
                                          mint::NODE_CENTERED, 3 );
  double_ptr2 = mesh->createField< double >( "node_vectors_2double",
                                             mint::NODE_CENTERED, 2);
  int_ptr2    =  mesh->createField< int >( "node_vectors_2int",
                                           mint::NODE_CENTERED, 2 );

  for ( int idx = 0 ; idx < num_nodes ; ++idx )
  {
    double coords[] = { 0, 0, 0 };
    mesh->getNode( idx, coords );

    const double x = coords[0];
    const double y = coords[1];
    const double z = coords[2];

    double r2 = (x * x) + (y * y) + (z * z);
    r2 = std::exp(-r2 / 100.0);
    double v1 = (x * x) - 4 * y + (2 * y * z) - (5 * x * y * z + 1) * r2;
    double v2 = (y * y) - 4 * z + (2 * x * z) - (5 * x * y * z + 1) * r2;
    double v3 = (z * z) - 4 * x + (2 * x * y) - (5 * x * y * z + 1) * r2;

    double_ptr3[ 3 * idx ] = v1;
    double_ptr3[ 3 * idx + 1 ] = v2;
    double_ptr3[ 3 * idx + 2 ] = v3;

    int_ptr3[ 3 * idx ] = static_cast<int>(v1);
    int_ptr3[ 3 * idx + 1 ] = static_cast<int>(v2);
    int_ptr3[ 3 * idx + 2 ] = static_cast<int>(v3);

    double_ptr2[ 2 * idx ] = v1;
    double_ptr2[ 2 * idx + 1 ] = v2;

    int_ptr2[ 2 * idx ] = static_cast<int>(v1);
    int_ptr2[ 2 * idx + 1 ] = static_cast<int>(v2);
  }

  /* Particle meshes can only have node centered fields. */
  if ( mesh->getMeshType() == PARTICLE_MESH )
  {
    return;
  }

  double_ptr3 = mesh->createField< double >( "cell_vectors_3double",
                                             mint::CELL_CENTERED,3 );
  int_ptr3    = mesh->createField< int >( "cell_vectors_3int",
                                          mint::CELL_CENTERED, 3 );
  double_ptr2 = mesh->createField< double >( "cell_vectors_2double",
                                             mint::CELL_CENTERED,2);
  int_ptr2    = mesh->createField< int >( "cell_vectors_2int",
                                          mint::CELL_CENTERED, 2 );

  for ( int idx = 0 ; idx < num_cells ; ++idx )
  {
    double_ptr3[ 3 * idx ] = idx;
    double_ptr3[ 3 * idx + 1 ] = idx + 1;
    double_ptr3[ 3 * idx + 2 ] = idx + 2;

    int_ptr3[ 3 * idx ] = idx;
    int_ptr3[ 3 * idx + 1 ] = idx + 1;
    int_ptr3[ 3 * idx + 2 ] = idx + 2;

    double_ptr2[ 2 * idx ] = idx;
    double_ptr2[ 2 * idx + 1 ] = idx + 1;

    int_ptr2[ 2 * idx ] = idx;
    int_ptr2[ 2 * idx + 1 ] = idx + 1;
  }

}

/*!
 * \brief Creates artificial 4 dimensional mesh data.
 * \param [in] mesh the mesh to populate.
 * \pre mesh != AXOM_NULLPTR
 */
void create_multidim_data( Mesh* mesh )
{
  const IndexType num_nodes = mesh->getNumberOfNodes();
  const IndexType num_cells = mesh->getNumberOfCells();

  double* double_ptr = AXOM_NULLPTR;
  int* int_ptr    = AXOM_NULLPTR;

  double_ptr = mesh->createField< double >( "node_multidim_double",
                                            mint::NODE_CENTERED, 4 );
  int_ptr    = mesh->createField< int >( "node_multidim_int",
                                         mint::NODE_CENTERED, 4 );

  for ( int idx = 0 ; idx < num_nodes ; ++idx )
  {
    double coords[] = { 0, 0, 0 };
    mesh->getNode( idx, coords );

    const double x = coords[0];
    const double y = coords[1];
    const double z = coords[2];

    double r2 = (x * x) + (y * y) + (z * z);
    r2 = std::exp(-r2 / 100.0);
    double v1 = (x * x) - 4 * y + (2 * y * z) - (5 * x * y * z + 1) * r2;
    double v2 = (y * y) - 4 * z + (2 * x * z) - (5 * x * y * z + 1) * r2;
    double v3 = (z * z) - 4 * x + (2 * x * y) - (5 * x * y * z + 1) * r2;
    double v4 = (x * x + y * y + z * z) * r2;

    double_ptr[4 * idx + 0] = v1;
    double_ptr[4 * idx + 1] = v2;
    double_ptr[4 * idx + 2] = v3;
    double_ptr[4 * idx + 3] = v4;

    int_ptr[4 * idx + 0] = static_cast<int>(v1);
    int_ptr[4 * idx + 1] = static_cast<int>(v2);
    int_ptr[4 * idx + 2] = static_cast<int>(v3);
    int_ptr[4 * idx + 3] = static_cast<int>(v4);
  }

  /* Particle meshes can only have node centered fields. */
  if ( mesh->getMeshType() == PARTICLE_MESH )
  {
    return;
  }

  double_ptr = mesh->createField< double >( "cell_multidim_double",
                                            mint::CELL_CENTERED, 4 );
  int_ptr    = mesh->createField< int >( "cell_multidim_int",
                                         mint::CELL_CENTERED, 4 );

  for ( int idx = 0 ; idx < num_cells ; ++idx )
  {
    double_ptr[ 4 * idx  + 0] = idx;
    double_ptr[ 4 * idx  + 1] = idx + 1;
    double_ptr[ 4 * idx  + 2] = idx + 2;
    double_ptr[ 4 * idx  + 3] = idx + 3;

    int_ptr[ 4 * idx  + 0] = idx;
    int_ptr[ 4 * idx  + 1] = idx + 1;
    int_ptr[ 4 * idx  + 2] = idx + 2;
    int_ptr[ 4 * idx  + 3] = idx + 3;
  }

}

/*!
 * \brief Creates artificial mesh data and then writes the mesh out to disk.
 * \param [in] mesh the mesh to write out.
 * \param [in] path the path of the file to be written.
 * \pre mesh != AXOM_NULLPTR
 */
void populate_and_write( Mesh* mesh, const std::string& path )
{
  create_scalar_data( mesh );
  create_vector_data( mesh );
  create_multidim_data( mesh );
  write_vtk( mesh, path );
}

/*!
 * \brief Checks that the VTK header was written correctly.
 * \param [in] file the file to parse.
 */
void check_header( std::ifstream& file )
{
  std::string buffer;
  std::getline(file, buffer);
  EXPECT_EQ(buffer, "# vtk DataFile Version 3.0");
  std::getline(file, buffer);
  EXPECT_EQ(buffer, "Mesh generated by axom::mint::write_vtk");
  std::getline(file, buffer);
  EXPECT_EQ(buffer, "ASCII");
}

/*!
 * \brief Checks that a scalar field was written correctly.
 * \param [in] field the field to check against.
 * \param [in] file the file to parse.
 * \param [in] offset the offset into the field to start at.
 * \pre field != AXOM_NULLPTR
 */
void check_scalar( const Field* field, std::ifstream& file,
                   uint offset = 0 )
{
  const int num_components   = field->getNumComponents( );
  const IndexType num_values = field->getNumTuples( );

  std::string buffer;
  file >> buffer;
  EXPECT_EQ(  buffer, "LOOKUP_TABLE" );
  file >> buffer;
  EXPECT_EQ(  buffer, "default" );

  if ( field->getType() == DOUBLE_FIELD_TYPE )
  {
    const double* field_data = Field::getDataPtr< double >( field );
    for ( IndexType idx = 0 ; idx < num_values ; ++idx )
    {
      double temp;
      file >> temp;
      EXPECT_DOUBLE_EQ( temp, field_data[ num_components * idx + offset ] );
    }
  }
  else if ( field->getType() == INT32_FIELD_TYPE )
  {
    const int* field_data = Field::getDataPtr< int >( field );
    for ( IndexType idx = 0 ; idx < num_values ; ++idx )
    {
      int temp;
      file >> temp;
      EXPECT_EQ( field_data[ num_components * idx + offset ], temp );
    }
  }
}

/*!
 * \brief Checks that a vector field was written correctly.
 * \param [in] field the vector field to check against.
 * \param [in] file the file to parse.
 * \pre mesh != AXOM_NULLPTR
 * \pre field->getNumComponents() == 2 || field->getNumComponents() == 3
 */
void check_vector_data( const Field* field, std::ifstream& file )
{
  const int num_components   = field->getNumComponents();
  const IndexType num_values = field->getNumTuples();

  if ( field->getType() == DOUBLE_FIELD_TYPE )
  {
    const double* field_data = Field::getDataPtr< double >( field );
    double temp;
    for ( IndexType idx = 0 ; idx < num_values ; ++idx )
    {
      for ( int dim = 0 ; dim < num_components ; ++dim )
      {
        file >> temp;
        EXPECT_DOUBLE_EQ( temp, field_data[ num_components * idx + dim ] );
      }

      if ( num_components == 2 )
      {
        file >> temp;
        EXPECT_EQ( temp, 0.0 );
      }
    }
  }
  else if ( field->getType() == INT32_FIELD_TYPE )
  {
    const int* field_data = Field::getDataPtr< int >( field );
    int temp;
    for ( IndexType idx = 0 ; idx < num_values ; ++idx )
    {
      for ( int dim = 0 ; dim < num_components ; ++dim )
      {
        file >> temp;
        EXPECT_EQ( temp, field_data[ num_components * idx + dim ] );
      }

      if ( num_components == 2 )
      {
        file >> temp;
        EXPECT_EQ( temp, 0 );
      }
    }
  }
}

/*!
 * \brief Checks that a multidimensional field was written correctly.
 * \param [in] field the multidimensional field to check against.
 * \param [in] file the file to parse.
 * \pre mesh != AXOM_NULLPTR
 * \pre field->getNumComponents() > 3
 */
void check_multidim_data( const Field* field, std::ifstream& file )
{
  const int num_components = field->getNumComponents();
  const int field_type = field->getType();
  check_scalar( field, file, 0 );

  std::string type, name, d_type;
  for ( int comp = 1 ; comp < num_components ; ++comp )
  {
    file >> type >> name >> d_type;
    EXPECT_EQ( type, "SCALARS" );

    std::stringstream temp;
    temp << field->getName() << "_";
    temp << std::setfill('0') << std::setw(3) << comp;
    EXPECT_EQ( name, temp.str() );

    if ( d_type == "double" )
    {
      EXPECT_EQ( field_type, DOUBLE_FIELD_TYPE );
    }
    else if ( d_type == "int" )
    {
      EXPECT_EQ( field_type, INT32_FIELD_TYPE );
    }
    else
    {
      EXPECT_TRUE( false ) << "Unknown field data type: " << d_type;
    }

    check_scalar( field, file, comp );
  }
}

/*!
 * \brief Checks that field data was written correctly.
 * \param [in] field_data the field data to check against.
 * \param [in] file the file to parse.
 * \pre field_data != AXOM_NULLPTR
 */
void check_fieldData( const FieldData* field_data, std::ifstream& file )
{
  std::set< std::string > fields_read;
  std::string type, name, d_type;
  size_t cur_pos;
  while ( file.good() )
  {
    file >> type >> name >> d_type;

    if ( field_data->hasField( name ) )
    {
      fields_read.insert( name );
      const Field* field = field_data->getField( name );

      if ( d_type == "double" )
      {
        EXPECT_EQ( field->getType(), DOUBLE_FIELD_TYPE );
      }
      else if ( d_type == "int" )
      {
        EXPECT_EQ( field->getType(), INT32_FIELD_TYPE );
      }
      else
      {
        EXPECT_TRUE( false ) << "Unknown field data type: " << d_type;
      }

      if ( type == "SCALARS" )
      {
        EXPECT_TRUE( field->getNumComponents() == 1 );
        check_scalar( field, file );
      }
      else if ( type == "VECTORS" )
      {
        EXPECT_TRUE( field->getNumComponents() == 2 ||
                     field->getNumComponents() == 3 );
        check_vector_data( field, file );
      }
    }
    else
    {
      size_t underscore_pos = name.size() - 4;
      std::string true_name = name.substr( 0, underscore_pos );
      EXPECT_EQ( name.substr(underscore_pos), "_000" );
      ASSERT_TRUE( field_data->hasField( true_name ) ) << true_name;
      fields_read.insert( true_name );

      const Field* field = field_data->getField( true_name );

      if ( d_type == "double" )
      {
        EXPECT_EQ( field->getType(), DOUBLE_FIELD_TYPE );
      }
      else if ( d_type == "int" )
      {
        EXPECT_EQ( field->getType(), INT32_FIELD_TYPE );
      }
      else
      {
        EXPECT_TRUE(false) << "Unknown field data type: " << d_type;
      }

      EXPECT_TRUE( field->getNumComponents() > 3 );
      check_multidim_data( field, file );
    }

    cur_pos = file.tellg();
    file >> type;
    file.seekg( cur_pos );
    if ( type == "CELL_DATA" || type == "POINT_DATA" )
    {
      break;
    }
  }

  EXPECT_EQ( static_cast< int >( fields_read.size() ),
             field_data->getNumFields() );
  for ( std::set< std::string >::iterator it = fields_read.begin() ;
        it != fields_read.end() ; ++it )
  {
    EXPECT_TRUE( field_data->hasField( *it ) );
  }
}

/*!
 * \brief Checks that the mesh cell and node data were written correctly.
 * \param [in] mesh the mesh to check against.
 * \param [in] file the file to parse.
 * \pre mesh != AXOM_NULLPTR
 */
void check_data( const Mesh* mesh, std::ifstream& file )
{
  std::string data_type;
  IndexType data_size;
  while ( file.good() )
  {
    file >> data_type >> data_size;
    if ( data_type == "POINT_DATA" )
    {
      EXPECT_EQ( data_size, mesh->getNumberOfNodes() );
      const FieldData* node_data = mesh->getFieldData( mint::NODE_CENTERED );
      check_fieldData( node_data, file );
    }
    else if ( data_type == "CELL_DATA" )
    {
      EXPECT_EQ( data_size, mesh->getNumberOfCells() );
      const FieldData* cell_data = mesh->getFieldData( mint::CELL_CENTERED );
      check_fieldData( cell_data, file );
    }
  }
}

/*!
 * \brief Checks that the uniform mesh header was written correctly.
 * \param [in] mesh the mesh to check against.
 * \param [in] file the file to parse.
 * \pre mesh != AXOM_NULLPTR
 */
void check_uniform_mesh( const UniformMesh* u_mesh, std::ifstream& file )
{
  std::string buffer;
  file >> buffer;
  EXPECT_EQ( buffer, "DATASET" );
  file >> buffer;
  EXPECT_EQ( buffer, "STRUCTURED_POINTS" );

  IndexType ext_size[3];
  u_mesh->getExtentSize( ext_size );
  file >> buffer;
  EXPECT_EQ( buffer, "DIMENSIONS" );
  for ( int i = 0 ; i < 3 ; ++i )
  {
    IndexType temp;
    file >> temp;
    EXPECT_EQ( temp, ext_size[ i ] );
  }

  const double* origin = u_mesh->getOrigin( );
  file >> buffer;
  EXPECT_EQ( buffer, "ORIGIN" );
  for ( int i = 0 ; i < 3 ; ++i )
  {
    double temp;
    file >> temp;
    EXPECT_DOUBLE_EQ( temp, origin[ i ] );
  }

  const double* spacing = u_mesh->getSpacing( );
  file >> buffer;
  EXPECT_EQ( buffer, "SPACING" );
  for ( int i = 0 ; i < 3 ; ++i )
  {
    double temp;
    file >> temp;
    EXPECT_DOUBLE_EQ( temp, spacing[ i ]);
  }
}

/*!
 * \brief Checks that the rectilinear mesh header was written correctly.
 * \param [in] mesh the mesh to check against.
 * \param [in] file the file to parse.
 * \pre mesh != AXOM_NULLPTR
 */
void check_rectilinear_mesh( const RectilinearMesh* r_mesh,
                             std::ifstream& file )
{
  std::string buffer;
  file >> buffer;
  EXPECT_EQ(  buffer, "DATASET" );
  file >> buffer;
  EXPECT_EQ(  buffer, "RECTILINEAR_GRID" );

  IndexType ext_size[3];
  r_mesh->getExtentSize( ext_size );
  file >> buffer;
  EXPECT_EQ( buffer, "DIMENSIONS" );
  for ( int i = 0 ; i < 3 ; ++i )
  {
    IndexType temp;
    file >> temp;
    EXPECT_EQ( temp, ext_size[ i ] );
  }

  std::string coord_names[3] = { "X_COORDINATES", "Y_COORDINATES",
                                 "Z_COORDINATES" };
  std::string extracted_name, extracted_type;
  IndexType extracted_size;
  double extracted_coord;
  for ( int dim = 0 ; dim < r_mesh->getDimension() ; ++dim )
  {
    file >> extracted_name >> extracted_size >> extracted_type;
    EXPECT_EQ(  extracted_name, coord_names[ dim ] );
    EXPECT_EQ(  extracted_size, ext_size[ dim ] );
    EXPECT_EQ(  extracted_type, "double" );

    const double* coord_array = r_mesh->getCoordinateArray( dim );
    for (IndexType i = 0 ; i < ext_size[ dim ] ; ++i )
    {
      file >> extracted_coord;
      EXPECT_EQ( extracted_coord, coord_array[ i ] );
    }
  }
  for (int dim = r_mesh->getDimension() ; dim < 3 ; ++dim )
  {
    file >> extracted_name >> extracted_size >> extracted_type;
    file >> extracted_coord;

    EXPECT_EQ(  extracted_name,   coord_names[ dim ] );
    EXPECT_EQ(  extracted_size,   ext_size[ dim ] );
    EXPECT_EQ(  extracted_type,   "double" );
    EXPECT_EQ(  extracted_coord,  0.0 );
  }
}

/*!
 * \brief Checks that the mesh node coordinates were written correctly.
 * \param [in] mesh the mesh to check against.
 * \param [in] file the file to parse.
 * \pre mesh != AXOM_NULLPTR
 */
void check_points( const Mesh* mesh, std::ifstream& file )
{
  const IndexType num_nodes = mesh->getNumberOfNodes();

  std::string extracted_name, extracted_type;
  IndexType extracted_size;
  file >> extracted_name >> extracted_size >> extracted_type;
  EXPECT_EQ(  extracted_name, "POINTS" );
  EXPECT_EQ(  extracted_size, num_nodes );
  EXPECT_EQ(  extracted_type, "double" );

  double extracted_coord = 0.0;
  for ( IndexType idx = 0 ; idx < num_nodes ; ++idx )
  {
    double coords[] = { 0, 0, 0 };
    mesh->getNode( idx, coords );

    const double x = coords[0];
    const double y = coords[1];
    const double z = coords[2];

    file >> extracted_coord;
    EXPECT_EQ( extracted_coord, x );

    file >> extracted_coord;
    EXPECT_EQ( extracted_coord, y );

    file >> extracted_coord;
    EXPECT_EQ( extracted_coord, z );

  }
}

/*!
 * \brief Checks that the mesh cells were written correctly.
 * \param [in] mesh the mesh to check against.
 * \param [in] file the file to parse.
 * \pre mesh != AXOM_NULLPTR
 */
void check_cells( const Mesh* mesh, std::ifstream& file )
{
  const IndexType num_cells = mesh->getNumberOfCells();

  /* First need to get total size of the connectivity array. */
  /* If the mesh only has one cell type we can calculate this directly. */
  int max_cell_nodes = mesh->getNumberOfCellNodes( 0 );
  IndexType total_size = ( max_cell_nodes + 1 ) * num_cells;

  /* If the mesh has mixed cells then we need to loop over the elements. */
  if ( mesh->hasMixedCellTypes() )
  {
    total_size = num_cells;
    for ( IndexType cellIdx = 0 ; cellIdx < num_cells ; ++cellIdx )
    {
      const int num_cell_nodes = mesh->getNumberOfCellNodes( cellIdx );
      max_cell_nodes = utilities::max(num_cell_nodes, max_cell_nodes);
      total_size += num_cell_nodes;
    }
  }

  std::string type;
  IndexType extracted_cells, extracted_size;
  file >> type >> extracted_cells >> extracted_size;
  EXPECT_EQ( type, "CELLS" );
  EXPECT_EQ( extracted_cells, num_cells );
  EXPECT_EQ( extracted_size, total_size );

  /* Write out the mesh cell connectivity. */
  IndexType temp;
  IndexType cell_nodes[ max_cell_nodes ];
  for ( IndexType cellIdx = 0 ; cellIdx < num_cells ; ++cellIdx )
  {
    const int num_cell_nodes = mesh->getNumberOfCellNodes( cellIdx );
    mesh->getCell( cellIdx, cell_nodes );

    file >> temp;
    EXPECT_EQ( temp, num_cell_nodes );
    for ( int i = 0 ; i < num_cell_nodes ; ++i )
    {
      file >> temp;
      EXPECT_EQ( temp, cell_nodes[ i ] );
    }
  }

  /* Write out the mesh cell types. */
  file >> type >> extracted_cells;
  EXPECT_EQ( type, "CELL_TYPES" );
  EXPECT_EQ( extracted_cells,  num_cells );
  for ( IndexType cellIdx = 0 ; cellIdx < num_cells ; ++cellIdx )
  {
    int cell_type = mesh->getCellType( cellIdx );
    file >> temp;
    EXPECT_EQ( temp, cell_info[ cell_type ].vtk_type );
  }
}

/*!
 * \brief Checks that the curvilinear mesh header was written correctly.
 * \param [in] mesh the mesh to check against.
 * \param [in] file the file to parse.
 * \pre mesh != AXOM_NULLPTR
 */
void check_curvilinear_mesh( const CurvilinearMesh* c_mesh,
                             std::ifstream& file )
{
  std::string buffer;
  file >> buffer;
  EXPECT_EQ(  buffer, "DATASET" );
  file >> buffer;
  EXPECT_EQ(  buffer, "STRUCTURED_GRID" );

  IndexType ext_size[3];
  c_mesh->getExtentSize( ext_size );
  file >> buffer;
  EXPECT_EQ( buffer, "DIMENSIONS" );
  for ( int i = 0 ; i < 3 ; ++i )
  {
    IndexType temp;
    file >> temp;
    EXPECT_EQ( temp, ext_size[ i ] );
  }

  check_points( c_mesh, file );
}

/*!
 * \brief Checks that the unstructured mesh header was written correctly.
 * \param [in] mesh the mesh to check against.
 * \param [in] file the file to parse.
 * \pre mesh != AXOM_NULLPTR
 */
void check_unstructured_mesh( const Mesh* mesh, std::ifstream& file )
{
  std::string buffer;
  file >> buffer;
  EXPECT_EQ(  buffer, "DATASET" );
  file >> buffer;
  EXPECT_EQ(  buffer, "UNSTRUCTURED_GRID" );

  check_points( mesh, file );
  check_cells( mesh, file );
}

} /* end namespace internal */

//------------------------------------------------------------------------------
// UNIT TESTS
//------------------------------------------------------------------------------

/*!
 * \brief Creates a 3D UniformMesh and writes it out to disk using
 *  mint::write_vtk and then reads the file back in to check for correctness.
 */
TEST( mint_util_write_vtk, UniformMesh3D )
{
  const std::string path = "uniformMesh3D.vtk";
  const int64 ext[6] = { 0, 10, 0, 10, 0, 10 };
  const double origin[3] = { -5.0, -10.0, -15.0 };
  const double corner[3] = { 5.0, 10.0, 15.0 };
  UniformMesh* u_mesh = new UniformMesh( 3, ext, origin, corner );

  internal::populate_and_write( u_mesh, path );
  std::ifstream file( path.c_str() );
  ASSERT_TRUE( file );
  internal::check_header( file );
  internal::check_uniform_mesh( u_mesh, file );
  internal::check_data( u_mesh, file );

  file.close();
  delete u_mesh;
#if DELETE_VTK_FILES
  std::remove( path.c_str() );
#endif
}

/*!
 * \brief Creates a 2D UniformMesh and writes it out to disk using
 *  mint::write_vtk and then reads the file back in to check for correctness.
 */
TEST( mint_util_write_vtk, UniformMesh2D )
{
  const std::string path = "uniformMesh2D.vtk";
  const int64 ext[4] = { 0, 10, 0, 10 };
  const double origin[2] = { -5.0, -10.0 };
  const double corner[2] = { 5.0, 10.0 };
  UniformMesh* u_mesh = new UniformMesh( 2, ext, origin, corner );

  internal::populate_and_write( u_mesh, path );
  std::ifstream file( path.c_str() );
  ASSERT_TRUE( file );
  internal::check_header( file );
  internal::check_uniform_mesh( u_mesh, file );
  internal::check_data( u_mesh, file );

  file.close();
  delete u_mesh;
#if DELETE_VTK_FILES
  std::remove( path.c_str() );
#endif
}

/*!
 * \brief Creates a 1D UniformMesh and writes it out to disk using
 *  mint::write_vtk and then reads the file back in to check for correctness.
 */
TEST( mint_util_write_vtk, UniformMesh1D )
{
  const std::string path = "uniformMesh1D.vtk";
  const int64 ext[2] = { 0, 10 };
  const double origin[1] = { -5.0 };
  const double corner[1] = { 5.0 };
  UniformMesh* u_mesh = new UniformMesh( 1, ext, origin, corner );

  internal::populate_and_write( u_mesh, path );
  std::ifstream file( path.c_str() );
  ASSERT_TRUE( file );
  internal::check_header( file );
  internal::check_uniform_mesh( u_mesh, file );
  internal::check_data( u_mesh, file );

  file.close();
  delete u_mesh;
#if DELETE_VTK_FILES
  std::remove( path.c_str() );
#endif
}

/*!
 * \brief Creates a 3D RectilinearMesh and writes it out to disk using
 *  mint::write_vtk and then reads the file back in to check for correctness.
 */
TEST( mint_util_write_vtk, RectilinearMesh3D )
{
  const std::string path = "rectilinearMesh3D.vtk";
  int64 ext[6] = { 0, 10, 0, 11, 0, 12 };
  RectilinearMesh* r_mesh = new RectilinearMesh( 3, ext );

  IndexType ext_size[3];
  r_mesh->getExtentSize( ext_size );
  for ( int dim = 0 ; dim < 3 ; ++dim )
  {
    double* x = r_mesh->getCoordinateArray( dim );
    for ( IndexType i = 0 ; i < ext_size[ dim ] ; ++i )
    {
      x[ i ] = i * i / 10.0;
    }
  }

  internal::populate_and_write( r_mesh, path );
  std::ifstream file( path.c_str() );
  ASSERT_TRUE( file );
  internal::check_header( file );
  internal::check_rectilinear_mesh( r_mesh, file );
  internal::check_data( r_mesh, file );

  file.close();
  delete r_mesh;
#if DELETE_VTK_FILES
  std::remove( path.c_str() );
#endif
}

/*!
 * \brief Creates a 2D RectilinearMesh and writes it out to disk using
 *  mint::write_vtk and then reads the file back in to check for correctness.
 */
TEST( mint_util_write_vtk, RectilinearMesh2D )
{
  const std::string path = "rectilinearMesh2D.vtk";
  int64 ext[4] = { 0, 10, 0, 11 };
  RectilinearMesh* r_mesh = new RectilinearMesh( 2, ext );

  IndexType ext_size[3];
  r_mesh->getExtentSize( ext_size );
  for ( int dim = 0 ; dim < 2 ; ++dim )
  {
    double* x = r_mesh->getCoordinateArray( dim );
    for ( IndexType i = 0 ; i < ext_size[ dim ] ; ++i )
    {
      x[ i ] = i * i / 10.0;
    }
  }

  internal::populate_and_write( r_mesh, path );
  std::ifstream file( path.c_str() );
  ASSERT_TRUE( file );
  internal::check_header( file );
  internal::check_rectilinear_mesh( r_mesh, file );
  internal::check_data( r_mesh, file );

  file.close();
  delete r_mesh;
#if DELETE_VTK_FILES
  std::remove( path.c_str() );
#endif
}

/*!
 * \brief Creates a 1D RectilinearMesh and writes it out to disk using
 *  mint::write_vtk and then reads the file back in to check for correctness.
 */
TEST( mint_util_write_vtk, RectilinearMesh1D )
{
  const std::string path = "rectilinearMesh1D.vtk";
  int64 ext[2] = { 0, 10 };
  RectilinearMesh* r_mesh = new RectilinearMesh( 1, ext );

  IndexType ext_size[3];
  r_mesh->getExtentSize( ext_size );
  for ( int dim = 0 ; dim < 1 ; ++dim )
  {
    double* x = r_mesh->getCoordinateArray( dim );
    for ( IndexType i = 0 ; i < ext_size[ dim ] ; ++i )
    {
      x[ i ] = i * i / 10.0;
    }
  }

  internal::populate_and_write( r_mesh, path );
  std::ifstream file( path.c_str() );
  ASSERT_TRUE( file );
  internal::check_header( file );
  internal::check_rectilinear_mesh( r_mesh, file );
  internal::check_data( r_mesh, file );

  file.close();
  delete r_mesh;
#if DELETE_VTK_FILES
  std::remove( path.c_str() );
#endif
}

/*!
 * \brief Creates a 3D CurvilinearMesh and writes it out to disk using
 *  mint::write_vtk and then reads the file back in to check for correctness.
 */
TEST( mint_util_write_vtk, CurvilinearMesh3D )
{
  const std::string path = "curvilinearMesh3D.vtk";
  int64 ext[6] = { 0, 10, 0, 11, 0, 12 };
  CurvilinearMesh* c_mesh = new CurvilinearMesh( 3, ext );

  IndexType ext_size[3];
  c_mesh->getExtentSize( ext_size );
  double* x_coords = c_mesh->getCoordinateArray( X_COORDINATE );
  double* y_coords = c_mesh->getCoordinateArray( Y_COORDINATE );
  double* z_coords = c_mesh->getCoordinateArray( Z_COORDINATE );
  for ( IndexType i = 0 ; i < ext_size[0] ; ++i )
  {
    for ( IndexType j = 0 ; j < ext_size[1] ; ++j )
    {
      for ( IndexType k = 0 ; k < ext_size[2] ; ++k )
      {
        IndexType idx = c_mesh->getLinearIndex( i, j, k );
        double x = i + utilities::random_real( -0.45, 0.45 );
        double y = j + utilities::random_real( -0.45, 0.45 );
        double z = k + utilities::random_real( -0.45, 0.45 );
        x_coords[ idx ] = x;
        y_coords[ idx ] = y;
        z_coords[ idx ] = z;
      }
    }
  }

  internal::populate_and_write( c_mesh, path );
  std::ifstream file( path.c_str() );
  ASSERT_TRUE( file );
  internal::check_header( file );
  internal::check_curvilinear_mesh( c_mesh, file );
  internal::check_data( c_mesh, file );

  file.close();
  delete c_mesh;
#if DELETE_VTK_FILES
  std::remove( path.c_str() );
#endif
}

/*!
 * \brief Creates a 2D CurvilinearMesh and writes it out to disk using
 *  mint::write_vtk and then reads the file back in to check for correctness.
 */
TEST( mint_util_write_vtk, CurvilinearMesh2D )
{
  const std::string path = "curvilinearMesh2D.vtk";
  int64 ext[6] = { 0, 2, 0, 2 };
  CurvilinearMesh* c_mesh = new CurvilinearMesh( 2, ext );

  IndexType ext_size[3];
  c_mesh->getExtentSize( ext_size );
  double* x_coords = c_mesh->getCoordinateArray( X_COORDINATE );
  double* y_coords = c_mesh->getCoordinateArray( Y_COORDINATE );
  for ( IndexType i = 0 ; i < ext_size[0] ; ++i )
  {
    for ( IndexType j = 0 ; j < ext_size[1] ; ++j )
    {
      IndexType idx = c_mesh->getLinearIndex( i, j );
      double x = i + utilities::random_real( -0.45, 0.45 );
      double y = j + utilities::random_real( -0.45, 0.45 );
      x_coords[ idx ] = x;
      y_coords[ idx ] = y;
    }
  }

  internal::populate_and_write( c_mesh, path );
  std::ifstream file( path.c_str() );
  ASSERT_TRUE( file );
  internal::check_header( file );
  internal::check_curvilinear_mesh( c_mesh, file );
  internal::check_data( c_mesh, file );

  file.close();
  delete c_mesh;
#if DELETE_VTK_FILES
  std::remove( path.c_str() );
#endif
}

/*!
 * \brief Creates a 1D CurvilinearMesh and writes it out to disk using
 *  mint::write_vtk and then reads the file back in to check for correctness.
 * \note Mesh, but not mesh data will display in Visit.
 */
TEST( mint_util_write_vtk, CurvilinearMesh1D )
{
  const std::string path = "curvilinearMesh1D.vtk";
  int64 ext[6] = { 0, 1 };
  CurvilinearMesh* c_mesh = new CurvilinearMesh( 1, ext );

  IndexType ext_size[3];
  c_mesh->getExtentSize( ext_size );
  double* x_coords = c_mesh->getCoordinateArray( X_COORDINATE );
  for ( IndexType idx = 0 ; idx < ext_size[0] ; ++idx )
  {
    double x = idx + utilities::random_real( -0.45, 0.45 );
    x_coords[ idx ] = x;
  }

  internal::populate_and_write( c_mesh, path );
  std::ifstream file( path.c_str() );
  ASSERT_TRUE( file );
  internal::check_header( file );
  internal::check_curvilinear_mesh( c_mesh, file );
  internal::check_data( c_mesh, file );

  file.close();
  delete c_mesh;
#if DELETE_VTK_FILES
  std::remove( path.c_str() );
#endif
}

/*!
 * \brief Creates a 3D UnstructuredMesh and writes it out to disk using
 *  mint::write_vtk and then reads the file back in to check for correctness.
 */
TEST( mint_util_write_vtk, UnstructuredMesh3D )
{
  const std::string path = "unstructuredMesh3D.vtk";
  const IndexType nx = 11;
  const IndexType ny = 12;
  const IndexType nz = 13;
  const IndexType nNodes = nx * ny * nz;
  const IndexType nCells = (nx-1) * (ny-1) * (nz-1);
  UnstructuredMesh< Topology::SINGLE >* u_mesh =
    new UnstructuredMesh< Topology::SINGLE >( 3, HEX, nNodes, nCells );

  for ( IndexType idx = 0 ; idx < nx ; ++idx )
  {
    for ( IndexType idy = 0 ; idy < ny ; ++idy )
    {
      for ( IndexType idz = 0 ; idz < nz ; ++idz )
      {
        double x = idx + utilities::random_real( -0.45, 0.45 );
        double y = idy + utilities::random_real( -0.45, 0.45 );
        double z = idz + utilities::random_real( -0.45, 0.45 );
        u_mesh->appendNode( x, y, z );
      }
    }
  }

  IndexType cell[8];
  for ( IndexType idx = 0 ; idx < nx - 1 ; ++idx )
  {
    for ( IndexType idy = 0 ; idy < ny - 1 ; ++idy )
    {
      for ( IndexType idz = 0 ; idz < nz - 1 ; ++idz )
      {
        IndexType index = idx * ny * nz + idy * nz + idz;
        cell[0] = index;
        cell[1] = index + ny * nz;
        cell[2] = index + nz + ny * nz;
        cell[3] = index + nz;
        cell[4] = cell[0] + 1;
        cell[5] = cell[1] + 1;
        cell[6] = cell[2] + 1;
        cell[7] = cell[3] + 1;

        u_mesh->appendCell( cell );
      }
    }
  }

  internal::populate_and_write( u_mesh, path );
  std::ifstream file( path.c_str() );
  ASSERT_TRUE( file );
  internal::check_header( file );
  internal::check_unstructured_mesh( u_mesh, file );
  internal::check_data( u_mesh, file );

  file.close();
  delete u_mesh;
#if DELETE_VTK_FILES
  std::remove( path.c_str() );
#endif
}

/*!
 * \brief Creates a 2D UnstructuredMesh and writes it out to disk using
 *  mint::write_vtk and then reads the file back in to check for correctness.
 */
TEST( mint_util_write_vtk, UnstructuredMesh2D )
{
  const std::string path = "unstructuredMesh2D.vtk";
  const IndexType nx = 11;
  const IndexType ny = 12;
  const IndexType nNodes = nx * ny;
  const IndexType nCells = (nx-1) * (ny-1);
  UnstructuredMesh< Topology::SINGLE >* u_mesh =
    new UnstructuredMesh< Topology::SINGLE >( 2, QUAD, nNodes, nCells );

  for ( IndexType idx = 0 ; idx < nx ; ++idx )
  {
    for ( IndexType idy = 0 ; idy < ny ; ++idy )
    {
      double x = idx + utilities::random_real( -0.45, 0.45 );
      double y = idy + utilities::random_real( -0.45, 0.45 );
      u_mesh->appendNode( x, y );
    }
  }

  IndexType cell[4];
  for ( IndexType idx = 0 ; idx < nx - 1 ; ++idx )
  {
    for ( IndexType idy = 0 ; idy < ny - 1 ; ++idy )
    {
      IndexType index = idx * ny + idy;
      cell[0] = index;
      cell[1] = index + ny;
      cell[2] = index + ny + 1;
      cell[3] = index + 1;

      u_mesh->appendCell( cell );
    }
  }

  internal::populate_and_write( u_mesh, path );
  std::ifstream file( path.c_str() );
  ASSERT_TRUE( file );
  internal::check_header( file );
  internal::check_unstructured_mesh( u_mesh, file );
  internal::check_data( u_mesh, file );

  file.close();
  delete u_mesh;
#if DELETE_VTK_FILES
  std::remove( path.c_str() );
#endif
}

/*!
 * \brief Creates a 1D UnstructuredMesh and writes it out to disk using
 *  mint::write_vtk and then reads the file back in to check for correctness.
 */
TEST( mint_util_write_vtk, UnstructuredMesh1D )
{
  const std::string path = "unstructuredMesh1D.vtk";
  const IndexType nx = 11;
  UnstructuredMesh< Topology::SINGLE >* u_mesh =
    new UnstructuredMesh< Topology::SINGLE >( 1, SEGMENT, nx, nx - 1 );

  for ( IndexType idx = 0 ; idx < nx ; ++idx )
  {
    double x = idx + utilities::random_real( -0.45, 0.45 );
    u_mesh->appendNode( x );
  }

  IndexType cell[2];
  for ( IndexType idx = 0 ; idx < nx - 1 ; ++idx )
  {
    cell[0] = idx;
    cell[1] = idx + 1;
    u_mesh->appendCell( cell, mint::SEGMENT );
  }

  internal::populate_and_write( u_mesh, path );
  std::ifstream file( path.c_str() );
  ASSERT_TRUE( file );
  internal::check_header( file );
  internal::check_unstructured_mesh( u_mesh, file );
  internal::check_data( u_mesh, file );

  file.close();
  delete u_mesh;
#if DELETE_VTK_FILES
  std::remove( path.c_str() );
#endif
}

/*!
 * \brief Creates a 3D UnstructuredMesh with mixed elements and writes it out to
 *  disk using mint::write_vtk and then reads the file back in to check for
 *  correctness. The UnstructuredMesh consists of one hexahedron with pyramids
 *  on the faces.
 */
TEST( mint_util_write_vtk, UnstructuredMixedMesh3D )
{
  const std::string path = "unstructuredMixedMesh3D.vtk";
  const IndexType nx = 2;
  const IndexType ny = 2;
  const IndexType nz = 2;
  const IndexType nNodes = nx * ny * nz;
  const IndexType nCells = (nx-1) * (ny-1) * (nz-1);

  UnstructuredMesh< Topology::MIXED >* u_mesh =
    new UnstructuredMesh< Topology::MIXED >( 3, nNodes, nCells );

  /* Create the nodes for the hexahedron. */
  for ( IndexType idx = 0 ; idx < nx ; ++idx )
  {
    for ( IndexType idy = 0 ; idy < ny ; ++idy )
    {
      for ( IndexType idz = 0 ; idz < nz ; ++idz )
      {
        double x = idx + utilities::random_real( -0.45, 0.45 );
        double y = idy + utilities::random_real( -0.45, 0.45 );
        double z = idz + utilities::random_real( -0.45, 0.45 );
        u_mesh->appendNode( x, y, z );
      }
    }
  }

  /* Create the nodes for the pyramids. */
  u_mesh->appendNode( -1.0, 0.5, 0.5 );
  u_mesh->appendNode( 2.0, 0.5, 0.5 );
  u_mesh->appendNode( 0.5, -1.0, 0.5 );
  u_mesh->appendNode( 0.5, 2.0, 0.5 );
  u_mesh->appendNode( 0.5, 0.5, -1.0 );
  u_mesh->appendNode( 0.5, 0.5, 2.0 );

  /* Create the hexahedron. */
  IndexType hex[8];
  hex[0] = 0;
  hex[1] = 4;
  hex[2] = 6;
  hex[3] = 2;
  hex[4] = 1;
  hex[5] = 5;
  hex[6] = 7;
  hex[7] = 3;
  u_mesh->appendCell( hex, HEX );

  /* Create the pyramid for the -x face. */
  IndexType pyramid[5];
  pyramid[0] = hex[0];
  pyramid[1] = hex[4];
  pyramid[2] = hex[7];
  pyramid[3] = hex[3];
  pyramid[4] = 8;
  u_mesh->appendCell( pyramid, PYRAMID );

  /* Create the pyramid for the +x face. */
  pyramid[0] = hex[1];
  pyramid[1] = hex[2];
  pyramid[2] = hex[6];
  pyramid[3] = hex[5];
  pyramid[4] = 9;
  u_mesh->appendCell( pyramid, PYRAMID );

  /* Create the pyramid for the -y face. */
  pyramid[0] = hex[0];
  pyramid[1] = hex[1];
  pyramid[2] = hex[5];
  pyramid[3] = hex[4];
  pyramid[4] = 10;
  u_mesh->appendCell( pyramid, PYRAMID );

  /* Create the pyramid for the +x face. */
  pyramid[0] = hex[2];
  pyramid[1] = hex[3];
  pyramid[2] = hex[7];
  pyramid[3] = hex[6];
  pyramid[4] = 11;
  u_mesh->appendCell( pyramid, PYRAMID );

  /* Create the pyramid for the -z face. */
  pyramid[0] = hex[3];
  pyramid[1] = hex[2];
  pyramid[2] = hex[1];
  pyramid[3] = hex[0];
  pyramid[4] = 12;
  u_mesh->appendCell( pyramid, PYRAMID );

  /* Create the pyramid for the +z face. */
  pyramid[0] = hex[4];
  pyramid[1] = hex[5];
  pyramid[2] = hex[6];
  pyramid[3] = hex[7];
  pyramid[4] = 13;
  u_mesh->appendCell( pyramid, PYRAMID );

  internal::populate_and_write( u_mesh, path );
  std::ifstream file( path.c_str() );
  ASSERT_TRUE( file );
  internal::check_header( file );
  internal::check_unstructured_mesh( u_mesh, file );
  internal::check_data( u_mesh, file );

  file.close();
  delete u_mesh;
#if DELETE_VTK_FILES
  std::remove( path.c_str() );
#endif
}

/*!
 * \brief Creates a 2D UnstructuredMesh with mixed elements and writes it out to
 *  disk using mint::write_vtk and then reads the file back in to check for
 *  correctness. The UnstructuredMesh consists of one quadrilateral with
 *  triangles on the faces.
 */
TEST( mint_util_write_vtk, UnstructuredMixedMesh2D )
{
  const std::string path = "unstructuredMixedMesh2D.vtk";
  const IndexType nx = 2;
  const IndexType ny = 2;
  const IndexType nNodes = nx * ny;
  const IndexType nCells = (nx-1) * (ny-1);
  UnstructuredMesh< Topology::MIXED >* u_mesh =
    new UnstructuredMesh< Topology::MIXED >( 2, nNodes, nCells );

  /* Create the nodes for the quad. */
  for ( IndexType idx = 0 ; idx < nx ; ++idx )
  {
    for ( IndexType idy = 0 ; idy < ny ; ++idy )
    {
      double x = idx + utilities::random_real( -0.45, 0.45 );
      double y = idy + utilities::random_real( -0.45, 0.45 );
      u_mesh->appendNode( x, y );
    }
  }

  /* Create the nodes for the triangles. */
  u_mesh->appendNode( -1.0, 0.5 );
  u_mesh->appendNode( 2.0, 0.5 );
  u_mesh->appendNode( 0.5, -1.0 );
  u_mesh->appendNode( 0.5, 2.0 );

  /* Create the quad. */
  IndexType quad[4];
  quad[0] = 0;
  quad[1] = 2;
  quad[2] = 3;
  quad[3] = 1;
  u_mesh->appendCell( quad, QUAD );

  /* Create the triangle for the -x face. */
  IndexType triangle[3];
  triangle[0] = quad[0];
  triangle[1] = quad[3];
  triangle[2] = 4;
  u_mesh->appendCell( triangle, TRIANGLE );

  /* Create the triangle for the +x face. */
  triangle[0] = quad[2];
  triangle[1] = quad[1];
  triangle[2] = 5;
  u_mesh->appendCell( triangle, TRIANGLE );

  /* Create the triangle for the -y face. */
  triangle[0] = quad[1];
  triangle[1] = quad[0];
  triangle[2] = 6;
  u_mesh->appendCell( triangle, TRIANGLE );

  /* Create the triangle for the +y face. */
  triangle[0] = quad[3];
  triangle[1] = quad[2];
  triangle[2] = 7;
  u_mesh->appendCell( triangle, TRIANGLE );

  internal::populate_and_write( u_mesh, path );
  std::ifstream file( path.c_str() );
  ASSERT_TRUE( file );
  internal::check_header( file );
  internal::check_unstructured_mesh( u_mesh, file );
  internal::check_data( u_mesh, file );

  file.close();
  delete u_mesh;
#if DELETE_VTK_FILES
  std::remove( path.c_str() );
#endif
}

/*!
 * \brief Creates a 3D ParticleMesh and writes it out to disk using
 *  mint::write_vtk and then reads the file back in to check for correctness.
 */
TEST( mint_util_write_vtk, ParticleMesh3D )
{
  const std::string path = "particleMesh3D.vtk";
  const IndexType nParticles = 1000;
  ParticleMesh* p_mesh = new ParticleMesh( 3, nParticles );

  double* x = p_mesh->getCoordinateArray( X_COORDINATE );
  double* y = p_mesh->getCoordinateArray( Y_COORDINATE );
  double* z = p_mesh->getCoordinateArray( Z_COORDINATE );

  for ( IndexType i = 0 ; i < nParticles ; ++i )
  {
    x[ i ] = 10.0 * utilities::random_real( 0.0, 1.0 );
    y[ i ] = 10.0 * utilities::random_real( 0.0, 1.0 );
    z[ i ] = 10.0 * utilities::random_real( 0.0, 1.0 );
  }

  internal::populate_and_write( p_mesh, path );
  std::ifstream file( path.c_str() );
  ASSERT_TRUE( file );
  internal::check_header( file );
  internal::check_unstructured_mesh( p_mesh, file );
  internal::check_data( p_mesh, file );

  file.close();
  delete p_mesh;
#if DELETE_VTK_FILES
  std::remove( path.c_str() );
#endif
}

/*!
 * \brief Creates a 2D ParticleMesh and writes it out to disk using
 *  mint::write_vtk and then reads the file back in to check for correctness.
 */
TEST( mint_util_write_vtk, ParticleMesh2D )
{
  const std::string path = "particleMesh2D.vtk";
  const IndexType nParticles = 1000;
  ParticleMesh* p_mesh = new ParticleMesh( 2, nParticles );

  double* x = p_mesh->getCoordinateArray( X_COORDINATE );
  double* y = p_mesh->getCoordinateArray( Y_COORDINATE );

  for ( IndexType i = 0 ; i < nParticles ; ++i )
  {
    x[ i ] = 10.0 * utilities::random_real( 0.0, 1.0 );
    y[ i ] = 10.0 * utilities::random_real( 0.0, 1.0 );
  }

  internal::populate_and_write( p_mesh, path );
  std::ifstream file( path.c_str() );
  ASSERT_TRUE( file );
  internal::check_header( file );
  internal::check_unstructured_mesh( p_mesh, file );
  internal::check_data( p_mesh, file );

  file.close();
  delete p_mesh;
#if DELETE_VTK_FILES
  std::remove( path.c_str() );
#endif
}

/*!
 * \brief Creates a 1D ParticleMesh and writes it out to disk using
 *  mint::write_vtk and then reads the file back in to check for correctness.
 */
TEST( mint_util_write_vtk, ParticleMesh1D )
{
  const std::string path = "particleMesh1D.vtk";
  const IndexType nParticles = 1000;
  ParticleMesh* p_mesh = new ParticleMesh( 1, nParticles );

  double* x = p_mesh->getCoordinateArray( X_COORDINATE );

  for ( IndexType i = 0 ; i < nParticles ; ++i )
  {
    x[ i ] = 10.0 * utilities::random_real( 0.0, 1.0 );
  }

  internal::populate_and_write( p_mesh, path );
  std::ifstream file( path.c_str() );
  ASSERT_TRUE( file );
  internal::check_header( file );
  internal::check_unstructured_mesh( p_mesh, file );
  internal::check_data( p_mesh, file );

  file.close();
  delete p_mesh;
#if DELETE_VTK_FILES
  std::remove( path.c_str() );
#endif
}

} /* end namespace mint */
} /* end namespace axom */

//------------------------------------------------------------------------------
using axom::slic::UnitTestLogger;

int main( int argc, char* argv[] )
{
  int result = 0;
  ::testing::InitGoogleTest( &argc, argv );
  UnitTestLogger logger;
  result = RUN_ALL_TESTS();
  return result;
}
