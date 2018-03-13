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

#include "axom_utils/Utilities.hpp"     // for utilities::max
#include "gtest/gtest.h"                // for TEST and EXPECT_* macros
#include "mint/CellType.hpp"            // for cell::vtk_types
#include "mint/CurvilinearMesh.hpp"     // for CurvilinearMesh
#include "mint/Field.hpp"               // for Field
#include "mint/FieldData.hpp"           // for FieldData
#include "mint/FieldTypes.hpp"          // for *_FIELD_TYPE
#include "mint/FieldVariable.hpp"       // for FieldVariable
#include "mint/Mesh.hpp"                // for Mesh
#include "mint/ParticleMesh.hpp"        // for ParticleMesh
#include "mint/RectilinearMesh.hpp"     // for RectilinearMesh
#include "mint/UniformMesh.hpp"         // for UniformMesh
#include "mint/UnstructuredMesh.hpp"    // for UnstructuredMesh
#include "mint/vtk_utils.hpp"           // for write_vtk
#include "slic/slic.hpp"                // for slic macros
#include "slic/UnitTestLogger.hpp"      // for UnitTestLogger

#include <cmath>                        // for std::exp
#include <cstdio>                       // for std::remove
#include <fstream>                      // for std::ifstream
#include <iomanip>                      // for std::setfill, std::setw
#include <string>                       // for std::string
#include <sstream>                      // for std::stringstream
#include <set>                          // for std::set

#if __cplusplus >= 201103L
#include <random>                       // for random number generator
#else
#include <cstdlib>                      // for rand
#endif

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
 * \brief Return a random double between min and max.
 * \param [in] min the lower bound.
 * \param [in] max the upper bound.
 */
double randomD(double min, double max)
{
#if __cplusplus >= 201103L
  static std::random_device rd;
  static std::mt19937_64 mt( rd() );
  static std::uniform_real_distribution< double > dist(0.0, 1.0);
  double temp = dist(mt);
#else
  double temp = ( double( rand() ) / RAND_MAX );
#endif
  return temp * (max - min) + min;
}

/*!
 * \brief Creates artificial 1 dimensional mesh data.
 * \param [in] mesh the mesh to populate.
 * \pre mesh != AXOM_NULLPTR
 */
void create_scalar_data( Mesh* mesh )
{
  const int mesh_dim = mesh->getDimension();
  const int mesh_num_nodes = mesh->getMeshNumberOfNodes();
  const int mesh_num_cells = mesh->getMeshNumberOfCells();

  Field* node_field_d = new FieldVariable< double >( "node_scalars_double",
                                                     mesh_num_nodes, 1 );
  Field* node_field_i = new FieldVariable< int >( "node_scalars_int",
                                                  mesh_num_nodes, 1 );
  double* double_ptr = node_field_d->getDoublePtr();
  int* int_ptr = node_field_i->getIntPtr();
  for ( int idx = 0 ; idx < mesh_num_nodes ; ++idx )
  {
    double x = mesh->getMeshNodeCoordinate( idx, 0 );
    double y = 0;
    double z = 0;

    if ( mesh_dim > 1)
    {
      y = mesh->getMeshNodeCoordinate( idx, 1 );
    }
    if ( mesh_dim > 2)
    {
      z = mesh->getMeshNodeCoordinate( idx, 2 );
    }

    double r2 = (x * x) + (y * y) + (z * z);
    r2 = std::exp(-r2 / 100.0);
    double temp = (x * x) - 4 * y + (2 * y * z) - (5 * x * y * z + 1) * r2;
    double_ptr[ idx ] = temp;
    int_ptr[ idx ] = static_cast<int>(temp);
  }
  mesh->getNodeFieldData()->addField( node_field_d );
  mesh->getNodeFieldData()->addField( node_field_i );

  Field* cell_field_d = new FieldVariable< double >( "cell_scalars_double",
                                                     mesh_num_cells, 1 );
  Field* cell_field_i = new FieldVariable< int >( "cell_scalars_int",
                                                  mesh_num_cells, 1 );
  double_ptr = cell_field_d->getDoublePtr();
  int_ptr = cell_field_i->getIntPtr();
  for ( int idx = 0 ; idx < mesh_num_cells ; ++idx )
  {
    double_ptr[ idx ] = idx;
    int_ptr[ idx ] = idx;
  }
  mesh->getCellFieldData()->addField( cell_field_d );
  mesh->getCellFieldData()->addField( cell_field_i );
}

/*!
 * \brief Creates artificial 3 dimensional mesh data.
 * \param [in] mesh the mesh to populate.
 * \pre mesh != AXOM_NULLPTR
 */
void create_vector_data( Mesh* mesh )
{
  const int mesh_dim = mesh->getDimension();
  const int mesh_num_nodes = mesh->getMeshNumberOfNodes();
  const int mesh_num_cells = mesh->getMeshNumberOfCells();

  Field* node_field_3d = new FieldVariable< double >( "node_vectors_3double",
                                                      mesh_num_nodes, 3 );
  Field* node_field_3i = new FieldVariable< int >( "node_vectors_3int",
                                                   mesh_num_nodes, 3 );
  Field* node_field_2d = new FieldVariable< double >( "node_vectors_2double",
                                                      mesh_num_nodes, 2 );
  Field* node_field_2i = new FieldVariable< int >( "node_vectors_2int",
                                                   mesh_num_nodes, 2 );
  double* double_ptr3 = node_field_3d->getDoublePtr();
  int* int_ptr3 = node_field_3i->getIntPtr();
  double* double_ptr2 = node_field_2d->getDoublePtr();
  int* int_ptr2 = node_field_2i->getIntPtr();
  for ( int idx = 0 ; idx < mesh_num_nodes ; ++idx )
  {
    double x = mesh->getMeshNodeCoordinate( idx, 0 );
    double y = 0;
    double z = 0;

    if ( mesh_dim > 1)
    {
      y = mesh->getMeshNodeCoordinate( idx, 1 );
    }
    if ( mesh_dim > 2)
    {
      z = mesh->getMeshNodeCoordinate( idx, 2 );
    }

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
  mesh->getNodeFieldData()->addField( node_field_3d );
  mesh->getNodeFieldData()->addField( node_field_3i );
  mesh->getNodeFieldData()->addField( node_field_2d );
  mesh->getNodeFieldData()->addField( node_field_2i );

  Field* cell_field_3d = new FieldVariable< double >( "cell_vectors_3double",
                                                      mesh_num_cells, 3 );
  Field* cell_field_3i = new FieldVariable< int >( "cell_vectors_3int",
                                                   mesh_num_cells, 3 );
  Field* cell_field_2d = new FieldVariable< double >( "cell_vectors_2double",
                                                      mesh_num_cells, 2 );
  Field* cell_field_2i = new FieldVariable< int >( "cell_vectors_2int",
                                                   mesh_num_cells, 2 );
  double_ptr3 = cell_field_3d->getDoublePtr();
  int_ptr3 = cell_field_3i->getIntPtr();
  double_ptr2 = cell_field_2d->getDoublePtr();
  int_ptr2 = cell_field_2i->getIntPtr();
  for ( int idx = 0 ; idx < mesh_num_cells ; ++idx )
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
  mesh->getCellFieldData()->addField( cell_field_3d );
  mesh->getCellFieldData()->addField( cell_field_3i );
  mesh->getCellFieldData()->addField( cell_field_2d );
  mesh->getCellFieldData()->addField( cell_field_2i );
}

/*!
 * \brief Creates artificial 4 dimensional mesh data.
 * \param [in] mesh the mesh to populate.
 * \pre mesh != AXOM_NULLPTR
 */
void create_multidim_data( Mesh* mesh )
{
  const int mesh_dim = mesh->getDimension();
  const int mesh_num_nodes = mesh->getMeshNumberOfNodes();
  const int mesh_num_cells = mesh->getMeshNumberOfCells();

  Field* node_field_d = new FieldVariable< double >( "node_multidim_double",
                                                     mesh_num_nodes, 4);
  Field* node_field_i = new FieldVariable< int >( "node_multidim_int",
                                                  mesh_num_nodes, 4);
  double* double_ptr = node_field_d->getDoublePtr();
  int* int_ptr = node_field_i->getIntPtr();
  for ( int idx = 0 ; idx < mesh_num_nodes ; ++idx )
  {
    double x = mesh->getMeshNodeCoordinate( idx, 0 );
    double y = 0;
    double z = 0;

    if ( mesh_dim > 1)
    {
      y = mesh->getMeshNodeCoordinate( idx, 1 );
    }
    if ( mesh_dim > 2)
    {
      z = mesh->getMeshNodeCoordinate( idx, 2 );
    }

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
  mesh->getNodeFieldData()->addField( node_field_d );
  mesh->getNodeFieldData()->addField( node_field_i );

  Field* cell_field_d = new FieldVariable< double >( "cell_multidim_double",
                                                     mesh_num_cells, 4 );
  Field* cell_field_i = new FieldVariable< int >( "cell_multidim_int",
                                                  mesh_num_cells, 4 );
  double_ptr = cell_field_d->getDoublePtr();
  int_ptr = cell_field_i->getIntPtr();
  for ( int idx = 0 ; idx < mesh_num_cells ; ++idx )
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
  mesh->getCellFieldData()->addField( cell_field_d );
  mesh->getCellFieldData()->addField( cell_field_i );
}

/*!
 * \brief Creates artificial mesh data and then writes the mesh out to disk.
 * \param [in] mesh the mesh to write out.
 * \param [in] path the path of the file to be written.
 * \pre mesh != AXOM_NULLPTR
 */
void populate_and_write( Mesh* mesh, std::string const path )
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
void check_scalar( Field* const field, std::ifstream& file,
                   axom::common::uint32 offset = 0 )
{
  const int num_components = field->getNumComponents();
  const int num_values = field->getNumTuples();

  std::string buffer;
  file >> buffer;
  EXPECT_EQ(  buffer, "LOOKUP_TABLE" );
  file >> buffer;
  EXPECT_EQ(  buffer, "default" );

  if ( field->getType() == DOUBLE_FIELD_TYPE )
  {
    double* field_data = field->getDoublePtr();
    for ( int idx = 0 ; idx < num_values ; ++idx )
    {
      double temp;
      file >> temp;
      EXPECT_DOUBLE_EQ( temp, field_data[ num_components * idx + offset ] );
    }
  }
  else if ( field->getType() == INTEGER_FIELD_TYPE )
  {
    int* field_data = field->getIntPtr();
    for ( int idx = 0 ; idx < num_values ; ++idx )
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
void check_vector_data( Field* const field, std::ifstream& file )
{
  const int num_components = field->getNumComponents();
  const int num_values = field->getNumTuples();

  if ( field->getType() == DOUBLE_FIELD_TYPE )
  {
    double* field_data = field->getDoublePtr();
    double temp;
    for ( int idx = 0 ; idx < num_values ; ++idx )
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
  else if ( field->getType() == INTEGER_FIELD_TYPE )
  {
    int* field_data = field->getIntPtr();
    int temp;
    for ( int idx = 0 ; idx < num_values ; ++idx )
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
void check_multidim_data( Field* const field, std::ifstream& file )
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
      EXPECT_EQ( field_type, INTEGER_FIELD_TYPE );
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
void check_fieldData( FieldData* const field_data, std::ifstream& file )
{
  std::set< std::string > fields_read;
  std::string type, name, d_type;
  int cur_pos;
  while ( file.good() )
  {
    file >> type >> name >> d_type;

    if ( field_data->hasField( name ) )
    {
      fields_read.insert( name );
      Field* field = field_data->getField( name );

      if ( d_type == "double" )
      {
        EXPECT_EQ( field->getType(), DOUBLE_FIELD_TYPE );
      }
      else if ( d_type == "int" )
      {
        EXPECT_EQ( field->getType(), INTEGER_FIELD_TYPE );
      }
      else
      {
        EXPECT_TRUE(false) << "Unknown field data type: " << d_type;
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
      ASSERT_NE( underscore_pos, std::string::npos) << name;
      std::string true_name = name.substr( 0, underscore_pos );
      EXPECT_EQ( name.substr(underscore_pos), "_000" );
      ASSERT_TRUE( field_data->hasField( true_name ) ) << true_name;
      fields_read.insert( true_name );

      Field* field = field_data->getField( true_name );

      if ( d_type == "double" )
      {
        EXPECT_EQ( field->getType(), DOUBLE_FIELD_TYPE );
      }
      else if ( d_type == "int" )
      {
        EXPECT_EQ( field->getType(), INTEGER_FIELD_TYPE );
      }
      else
      {
        EXPECT_TRUE(false) << "Unknown field data type: " << d_type;
      }

      EXPECT_TRUE( field->getNumComponents() > 3 );
      check_multidim_data( field, file );
    }

    cur_pos = static_cast<int>(file.tellg());
    file >> type;
    file.seekg( cur_pos );
    if ( type == "CELL_DATA" || type == "POINT_DATA" )
    {
      break;
    }
  }

  EXPECT_EQ( static_cast< int >( fields_read.size() ),
             field_data->getNumberOfFields() );
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
void check_data(Mesh* const mesh, std::ifstream& file )
{
  std::string data_type;
  int data_size;
  while ( file.good() )
  {
    file >> data_type >> data_size;
    if ( data_type == "POINT_DATA" )
    {
      EXPECT_EQ( data_size, mesh->getMeshNumberOfNodes() );
      FieldData* node_data = mesh->getNodeFieldData();
      check_fieldData( node_data, file );
    }
    else if ( data_type == "CELL_DATA" )
    {
      EXPECT_EQ( data_size, mesh->getMeshNumberOfCells() );
      FieldData* cell_data = mesh->getCellFieldData();
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
void check_uniform_mesh( UniformMesh* const u_mesh, std::ifstream& file )
{
  std::string buffer;
  file >> buffer;
  EXPECT_EQ(  buffer, "DATASET" );
  file >> buffer;
  EXPECT_EQ(  buffer, "STRUCTURED_POINTS" );

  int ext_size[3];
  u_mesh->getExtentSize( ext_size );
  file >> buffer;
  EXPECT_EQ( buffer, "DIMENSIONS" );
  for ( int i = 0 ; i < 3 ; ++i )
  {
    int temp;
    file >> temp;
    EXPECT_EQ( temp, ext_size[ i ] );
  }

  double origin[3];
  u_mesh->getOrigin( origin );
  file >> buffer;
  EXPECT_EQ( buffer, "ORIGIN" );
  for ( int i = 0 ; i < 3 ; ++i )
  {
    double temp;
    file >> temp;
    EXPECT_DOUBLE_EQ( temp, origin[ i ] );
  }

  double spacing[3];
  u_mesh->getSpacing( spacing );
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
void check_rectilinear_mesh( RectilinearMesh* const r_mesh,
                             std::ifstream& file )
{
  std::string buffer;
  file >> buffer;
  EXPECT_EQ(  buffer, "DATASET" );
  file >> buffer;
  EXPECT_EQ(  buffer, "RECTILINEAR_GRID" );

  int ext_size[3];
  r_mesh->getExtentSize( ext_size );
  file >> buffer;
  EXPECT_EQ( buffer, "DIMENSIONS" );
  for ( int i = 0 ; i < 3 ; ++i )
  {
    int temp;
    file >> temp;
    EXPECT_EQ( temp, ext_size[ i ] );
  }

  std::string coord_names[3] = { "X_COORDINATES", "Y_COORDINATES",
                                 "Z_COORDINATES" };
  std::string extracted_name, extracted_type;
  int extracted_size;
  double extracted_coord;
  for ( int dim = 0 ; dim < r_mesh->getDimension() ; ++dim )
  {
    file >> extracted_name >> extracted_size >> extracted_type;
    EXPECT_EQ(  extracted_name, coord_names[ dim ] );
    EXPECT_EQ(  extracted_size, ext_size[ dim ] );
    EXPECT_EQ(  extracted_type, "double" );

    const double* coord_array = r_mesh->getCoordinateArray( dim );
    for (int i = 0 ; i < ext_size[ dim ] ; ++i )
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
void check_points( Mesh* const mesh, std::ifstream& file )
{
  const int num_nodes = mesh->getMeshNumberOfNodes();
  const int mesh_dim = mesh->getDimension();

  std::string extracted_name, extracted_type;
  int extracted_size;
  file >> extracted_name >> extracted_size >> extracted_type;
  EXPECT_EQ(  extracted_name, "POINTS" );
  EXPECT_EQ(  extracted_size, num_nodes );
  EXPECT_EQ(  extracted_type, "double" );

  double extracted_coord;
  for ( int nodeIdx = 0 ; nodeIdx < num_nodes ; ++nodeIdx )
  {
    file >> extracted_coord;
    EXPECT_EQ( extracted_coord, mesh->getMeshNodeCoordinate( nodeIdx, 0 ) );
    for ( int dim = 1 ; dim < mesh_dim ; ++dim )
    {
      file >> extracted_coord;
      EXPECT_EQ( extracted_coord, mesh->getMeshNodeCoordinate( nodeIdx, dim ) );
    }
    for ( int dim = 0 ; dim < 3 - mesh_dim ; ++dim )
    {
      file >> extracted_coord;
      EXPECT_EQ( extracted_coord, 0.0 );
    }
  }
}

/*!
 * \brief Checks that the mesh cells were written correctly.
 * \param [in] mesh the mesh to check against.
 * \param [in] file the file to parse.
 * \pre mesh != AXOM_NULLPTR
 */
void check_cells( Mesh* const mesh, std::ifstream& file )
{
  const int num_cells = mesh->getMeshNumberOfCells();

  /* First need to get total size of the connectivity array. */
  /* If the mesh only has one cell type we can calculate this directly. */
  int max_cell_nodes = mesh->getMeshNumberOfCellNodes( 0 );
  int total_size = ( max_cell_nodes + 1 ) * num_cells;

  /* If the mesh has mixed cells then we need to loop over the elements. */
  if ( mesh->getMeshType() == MINT_UNSTRUCTURED_MIXED_ELEMENT_MESH )
  {
    total_size = num_cells;
    for ( int cellIdx = 0 ; cellIdx < num_cells ; ++cellIdx )
    {
      const int num_cell_nodes = mesh->getMeshNumberOfCellNodes( cellIdx );
      max_cell_nodes = utilities::max(num_cell_nodes, max_cell_nodes);
      total_size += num_cell_nodes;
    }
  }

  std::string type;
  int extracted_cells, extracted_size;
  file >> type >> extracted_cells >> extracted_size;
  EXPECT_EQ(  type,             "CELLS" );
  EXPECT_EQ(  extracted_cells,  num_cells );
  EXPECT_EQ(  extracted_size,   total_size );

  /* Write out the mesh cell connectivity. */
  int temp;
  int* cell_nodes = new int[ max_cell_nodes ];
  for ( int cellIdx = 0 ; cellIdx < num_cells ; ++cellIdx )
  {
    const int num_cell_nodes = mesh->getMeshNumberOfCellNodes( cellIdx );
    mesh->getMeshCell( cellIdx, cell_nodes );

    file >> temp;
    EXPECT_EQ( temp, num_cell_nodes );
    for ( int i = 0 ; i < num_cell_nodes ; ++i )
    {
      file >> temp;
      EXPECT_EQ( temp, cell_nodes[ i ] );
    }
  }
  delete[] cell_nodes;

  /* Write out the mesh cell types. */
  file >> type >> extracted_cells;
  EXPECT_EQ(  type,             "CELL_TYPES");
  EXPECT_EQ(  extracted_cells,  num_cells );
  for ( int cellIdx = 0 ; cellIdx < num_cells ; ++cellIdx )
  {
    int cell_type = mesh->getMeshCellType( cellIdx );
    file >> temp;
    EXPECT_EQ( temp, cell::vtk_types[ cell_type ] );
  }
}

/*!
 * \brief Checks that the curvilinear mesh header was written correctly.
 * \param [in] mesh the mesh to check against.
 * \param [in] file the file to parse.
 * \pre mesh != AXOM_NULLPTR
 */
void check_curvilinear_mesh( CurvilinearMesh* const c_mesh,
                             std::ifstream& file )
{
  std::string buffer;
  file >> buffer;
  EXPECT_EQ(  buffer, "DATASET" );
  file >> buffer;
  EXPECT_EQ(  buffer, "STRUCTURED_GRID" );

  int ext_size[3];
  c_mesh->getExtentSize( ext_size );
  file >> buffer;
  EXPECT_EQ( buffer, "DIMENSIONS" );
  for ( int i = 0 ; i < 3 ; ++i )
  {
    int temp;
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
void check_unstructured_mesh( Mesh* const mesh, std::ifstream& file )
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
TEST( mint_write_vtk, UniformMesh3D )
{
  const std::string path = "uniformMesh3D.vtk";
  const int ext[6] = { 0, 10, 0, 10, 0, 10 };
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
TEST( mint_write_vtk, UniformMesh2D )
{
  const std::string path = "uniformMesh2D.vtk";
  const int ext[4] = { 0, 10, 0, 10 };
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
TEST( mint_write_vtk, UniformMesh1D )
{
  const std::string path = "uniformMesh1D.vtk";
  const int ext[2] = { 0, 10 };
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
TEST( mint_write_vtk, RectilinearMesh3D )
{
  const std::string path = "rectilinearMesh3D.vtk";
  int ext[6] = { 0, 10, 0, 11, 0, 12 };
  RectilinearMesh* r_mesh = new RectilinearMesh( 3, ext );

  int ext_size[3];
  r_mesh->getExtentSize( ext_size );
  for ( int dim = 0 ; dim < 3 ; ++dim )
  {
    for ( int i = 0 ; i < ext_size[ dim ] ; ++i )
    {
      r_mesh->setCoordinate( dim, i, i * i / 10.0);
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
TEST( mint_write_vtk, RectilinearMesh2D )
{
  const std::string path = "rectilinearMesh2D.vtk";
  int ext[4] = { 0, 10, 0, 11 };
  RectilinearMesh* r_mesh = new RectilinearMesh( 2, ext );

  int ext_size[3];
  r_mesh->getExtentSize( ext_size );
  for ( int dim = 0 ; dim < 2 ; ++dim )
  {
    for ( int i = 0 ; i < ext_size[ dim ] ; ++i )
    {
      r_mesh->setCoordinate( dim, i, i * i / 10.0);
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
TEST( mint_write_vtk, RectilinearMesh1D )
{
  const std::string path = "rectilinearMesh1D.vtk";
  int ext[2] = { 0, 10 };
  RectilinearMesh* r_mesh = new RectilinearMesh( 1, ext );

  int ext_size[3];
  r_mesh->getExtentSize( ext_size );
  for ( int dim = 0 ; dim < 1 ; ++dim )
  {
    for ( int i = 0 ; i < ext_size[ dim ] ; ++i )
    {
      r_mesh->setCoordinate( dim, i, i * i / 10.0);
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
TEST( mint_write_vtk, CurvilinearMesh3D )
{
  const std::string path = "curvilinearMesh3D.vtk";
  int ext[6] = { 0, 10, 0, 11, 0, 12 };
  CurvilinearMesh* c_mesh = new CurvilinearMesh( 3, ext );

  int ext_size[3];
  c_mesh->getExtentSize( ext_size );
  for ( int idx = 0 ; idx < ext_size[0] ; ++idx )
  {
    for ( int idy = 0 ; idy < ext_size[1] ; ++idy )
    {
      for ( int idz = 0 ; idz < ext_size[2] ; ++idz )
      {
        double x = idx + internal::randomD( -0.45, 0.45 );
        double y = idy + internal::randomD( -0.45, 0.45 );
        double z = idz + internal::randomD( -0.45, 0.45 );
        c_mesh->setNode( idx, idy, idz, x, y, z );
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
TEST( mint_write_vtk, CurvilinearMesh2D )
{
  const std::string path = "curvilinearMesh2D.vtk";
  int ext[6] = { 0, 2, 0, 2 };
  CurvilinearMesh* c_mesh = new CurvilinearMesh( 2, ext );

  int ext_size[3];
  c_mesh->getExtentSize( ext_size );
  for ( int idx = 0 ; idx < ext_size[0] ; ++idx )
  {
    for ( int idy = 0 ; idy < ext_size[1] ; ++idy )
    {
      double x = idx + internal::randomD( -0.45, 0.45 );
      double y = idy + internal::randomD( -0.45, 0.45 );
      c_mesh->setNode( idx, idy, x, y );
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
TEST( mint_write_vtk, CurvilinearMesh1D )
{
  const std::string path = "curvilinearMesh1D.vtk";
  int ext[6] = { 0, 1 };
  CurvilinearMesh* c_mesh = new CurvilinearMesh( 1, ext );

  int ext_size[3];
  c_mesh->getExtentSize( ext_size );
  for ( int idx = 0 ; idx < ext_size[0] ; ++idx )
  {
    double x = idx + internal::randomD( -0.45, 0.45 );
    c_mesh->setNode( idx, x );
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
TEST( mint_write_vtk, UnstructuredMesh3D )
{
  const std::string path = "unstructuredMesh3D.vtk";
  UnstructuredMesh< MINT_HEX >* u_mesh = new UnstructuredMesh< MINT_HEX >( 3 );

  int nx = 11;
  int ny = 12;
  int nz = 13;
  for ( int idx = 0 ; idx < nx ; ++idx )
  {
    for ( int idy = 0 ; idy < ny ; ++idy )
    {
      for ( int idz = 0 ; idz < nz ; ++idz )
      {
        double x = idx + internal::randomD( -0.45, 0.45 );
        double y = idy + internal::randomD( -0.45, 0.45 );
        double z = idz + internal::randomD( -0.45, 0.45 );
        u_mesh->insertNode( x, y, z );
      }
    }
  }

  int cell[8];
  for ( int idx = 0 ; idx < nx - 1 ; ++idx )
  {
    for ( int idy = 0 ; idy < ny - 1 ; ++idy )
    {
      for ( int idz = 0 ; idz < nz - 1 ; ++idz )
      {
        int index = idx * ny * nz + idy * nz + idz;
        cell[0] = index;
        cell[1] = index + ny * nz;
        cell[2] = index + nz + ny * nz;
        cell[3] = index + nz;
        cell[4] = cell[0] + 1;
        cell[5] = cell[1] + 1;
        cell[6] = cell[2] + 1;
        cell[7] = cell[3] + 1;

        u_mesh->insertCell( cell, MINT_HEX, 8 );
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
TEST( mint_write_vtk, UnstructuredMesh2D )
{
  const std::string path = "unstructuredMesh2D.vtk";
  UnstructuredMesh< MINT_QUAD >* u_mesh =
    new UnstructuredMesh< MINT_QUAD >( 2 );

  int nx = 11;
  int ny = 12;
  for ( int idx = 0 ; idx < nx ; ++idx )
  {
    for ( int idy = 0 ; idy < ny ; ++idy )
    {
      double x = idx + internal::randomD( -0.45, 0.45 );
      double y = idy + internal::randomD( -0.45, 0.45 );
      u_mesh->insertNode( x, y );
    }
  }

  int cell[4];
  for ( int idx = 0 ; idx < nx - 1 ; ++idx )
  {
    for ( int idy = 0 ; idy < ny - 1 ; ++idy )
    {
      int index = idx * ny + idy;
      cell[0] = index;
      cell[1] = index + ny;
      cell[2] = index + ny + 1;
      cell[3] = index + 1;

      u_mesh->insertCell( cell, MINT_QUAD, 8 );
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
TEST( mint_write_vtk, UnstructuredMesh1D )
{
  const std::string path = "unstructuredMesh1D.vtk";
  UnstructuredMesh< MINT_SEGMENT >* u_mesh =
    new UnstructuredMesh< MINT_SEGMENT >( 1 );

  int nx = 11;
  for ( int idx = 0 ; idx < nx ; ++idx )
  {
    double x = idx + internal::randomD( 0.45, 0.45 );
    u_mesh->insertNode( x );
  }

  int cell[2];
  for ( int idx = 0 ; idx < nx - 1 ; ++idx )
  {
    cell[0] = idx;
    cell[1] = idx + 1;
    u_mesh->insertCell( cell, MINT_SEGMENT, 2 );
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
TEST( mint_write_vtk, UnstructuredMixedMesh3D )
{
  const std::string path = "unstructuredMixedMesh3D.vtk";
  UnstructuredMesh< MINT_MIXED_CELL >* u_mesh = new UnstructuredMesh<
    MINT_MIXED_CELL >( 3 );

  /* Create the nodes for the hexahedron. */
  int nx = 2;
  int ny = 2;
  int nz = 2;
  for ( int idx = 0 ; idx < nx ; ++idx )
  {
    for ( int idy = 0 ; idy < ny ; ++idy )
    {
      for ( int idz = 0 ; idz < nz ; ++idz )
      {
        double x = idx + internal::randomD( -0.45, 0.45 );
        double y = idy + internal::randomD( -0.45, 0.45 );
        double z = idz + internal::randomD( -0.45, 0.45 );
        u_mesh->insertNode( x, y, z );
      }
    }
  }

  /* Create the nodes for the pyramids. */
  u_mesh->insertNode( -1.0, 0.5,  0.5 );
  u_mesh->insertNode( 2.0,  0.5,  0.5 );
  u_mesh->insertNode( 0.5,  -1.0, 0.5 );
  u_mesh->insertNode( 0.5,  2.0,  0.5 );
  u_mesh->insertNode( 0.5,  0.5,  -1.0 );
  u_mesh->insertNode( 0.5,  0.5,  2.0 );

  /* Create the hexahedron. */
  int hex[8];
  hex[0] = 0;
  hex[1] = 4;
  hex[2] = 6;
  hex[3] = 2;
  hex[4] = 1;
  hex[5] = 5;
  hex[6] = 7;
  hex[7] = 3;
  u_mesh->insertCell( hex, MINT_HEX, 8 );

  /* Create the pyramid for the -x face. */
  int pyramid[5];
  pyramid[0] = hex[0];
  pyramid[1] = hex[4];
  pyramid[2] = hex[7];
  pyramid[3] = hex[3];
  pyramid[4] = 8;
  u_mesh->insertCell( pyramid, MINT_PYRAMID, 5 );

  /* Create the pyramid for the +x face. */
  pyramid[0] = hex[1];
  pyramid[1] = hex[2];
  pyramid[2] = hex[6];
  pyramid[3] = hex[5];
  pyramid[4] = 9;
  u_mesh->insertCell( pyramid, MINT_PYRAMID, 5 );

  /* Create the pyramid for the -y face. */
  pyramid[0] = hex[0];
  pyramid[1] = hex[1];
  pyramid[2] = hex[5];
  pyramid[3] = hex[4];
  pyramid[4] = 10;
  u_mesh->insertCell( pyramid, MINT_PYRAMID, 5 );

  /* Create the pyramid for the +x face. */
  pyramid[0] = hex[2];
  pyramid[1] = hex[3];
  pyramid[2] = hex[7];
  pyramid[3] = hex[6];
  pyramid[4] = 11;
  u_mesh->insertCell( pyramid, MINT_PYRAMID, 5 );

  /* Create the pyramid for the -z face. */
  pyramid[0] = hex[3];
  pyramid[1] = hex[2];
  pyramid[2] = hex[1];
  pyramid[3] = hex[0];
  pyramid[4] = 12;
  u_mesh->insertCell( pyramid, MINT_PYRAMID, 5 );

  /* Create the pyramid for the +z face. */
  pyramid[0] = hex[4];
  pyramid[1] = hex[5];
  pyramid[2] = hex[6];
  pyramid[3] = hex[7];
  pyramid[4] = 13;
  u_mesh->insertCell( pyramid, MINT_PYRAMID, 5 );

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
TEST( mint_write_vtk, UnstructuredMixedMesh2D )
{
  const std::string path = "unstructuredMixedMesh2D.vtk";
  UnstructuredMesh< MINT_MIXED_CELL >* u_mesh = new UnstructuredMesh<
    MINT_MIXED_CELL >( 2 );

  /* Create the nodes for the hexahedron. */
  int nx = 2;
  int ny = 2;
  for ( int idx = 0 ; idx < nx ; ++idx )
  {
    for ( int idy = 0 ; idy < ny ; ++idy )
    {
      double x = idx + internal::randomD( -0.45, 0.45 );
      double y = idy + internal::randomD( -0.45, 0.45 );
      u_mesh->insertNode( x, y );
    }
  }

  /* Create the nodes for the pyramids. */
  u_mesh->insertNode( -1.0, 0.5 );
  u_mesh->insertNode( 2.0,  0.5 );
  u_mesh->insertNode( 0.5,  -1.0 );
  u_mesh->insertNode( 0.5,  2.0 );

  /* Create the hexahedron. */
  int quad[4];
  quad[0] = 0;
  quad[1] = 2;
  quad[2] = 3;
  quad[3] = 1;
  u_mesh->insertCell( quad, MINT_QUAD, 4 );

  /* Create the pyramid for the -x face. */
  int triangle[3];
  triangle[0] = quad[0];
  triangle[1] = quad[3];
  triangle[2] = 4;
  u_mesh->insertCell( triangle, MINT_TRIANGLE, 3 );

  /* Create the pyramid for the +x face. */
  triangle[0] = quad[2];
  triangle[1] = quad[1];
  triangle[2] = 5;
  u_mesh->insertCell( triangle, MINT_TRIANGLE, 3 );

  /* Create the pyramid for the -y face. */
  triangle[0] = quad[1];
  triangle[1] = quad[0];
  triangle[2] = 6;
  u_mesh->insertCell( triangle, MINT_TRIANGLE, 3 );

  /* Create the pyramid for the +y face. */
  triangle[0] = quad[3];
  triangle[1] = quad[2];
  triangle[2] = 7;
  u_mesh->insertCell( triangle, MINT_TRIANGLE, 3 );

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
TEST( mint_write_vtk, ParticleMesh3D )
{
  const std::string path = "particleMesh3D.vtk";
  ParticleMesh* p_mesh = new ParticleMesh( 3 );

  for ( int i = 0 ; i < 1000 ; ++i )
  {
    double x = 10.0 * internal::randomD( 0.0, 1.0 );
    double y = 10.0 * internal::randomD( 0.0, 1.0 );
    double z = 10.0 * internal::randomD( 0.0, 1.0 );
    p_mesh->insertParticle( x, y, z );
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
TEST( mint_write_vtk, ParticleMesh2D )
{
  const std::string path = "particleMesh2D.vtk";
  ParticleMesh* p_mesh = new ParticleMesh( 2 );

  for ( int i = 0 ; i < 1000 ; ++i )
  {
    double x = 10.0 * internal::randomD( 0.0, 1.0 );
    double y = 10.0 * internal::randomD( 0.0, 1.0 );
    p_mesh->insertParticle( x, y );
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
TEST( mint_write_vtk, ParticleMesh1D )
{
  const std::string path = "particleMesh1D.vtk";
  ParticleMesh* p_mesh = new ParticleMesh( 1 );

  for ( int i = 0 ; i < 1000 ; ++i )
  {
    double x = 10.0 * internal::randomD( 0.0, 1.0 );
    p_mesh->insertParticle( x );
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

} /* end namespace mint */
} /* end namespace axom */
