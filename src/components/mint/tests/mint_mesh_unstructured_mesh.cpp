/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017, Lawrence Livermore National Security, LLC.
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
#include "mint/config.hpp"              /* for mint defintions */
#include "mint/UnstructuredMesh.hpp"               /* for mint::Array */

#include "axom_utils/Utilities.hpp"     /* for utilities::max */

#include "slic/UnitTestLogger.hpp"      /* for UnitTestLogger */
#include "slic/slic.hpp"                /* for slic macros */

#include "gtest/gtest.h"                /* for TEST and EXPECT_* macros */

#ifdef MINT_USE_SIDRE
#include "sidre/sidre.hpp"
#endif

// C/C++ includes
#include <algorithm>                    /* for std::fill_n */
#include <string>

namespace axom
{
namespace mint
{

constexpr double PI = 3.14159265358979323846;
constexpr double E =  2.71828182845904523536;
const char IGNORE_OUTPUT[] = ".*";

namespace internal
{


UnstructuredMesh< Topology::SINGLE >* 
createExternalSingle( int ndims, CellType cell_type, IndexType cell_capacity, 
                      IndexType node_capacity )
{
  double* x = new double[ node_capacity ];
  double* y = AXOM_NULLPTR;
  double* z = AXOM_NULLPTR;
  if ( ndims > 1 )
  {
    y = new double[ node_capacity ];
  }
  if ( ndims > 2 )
  {
    z = new double[ node_capacity ];
  }

  IndexType connec_capacity = cell_info[ cell_type ].num_nodes * cell_capacity;
  IndexType* connectivity = new IndexType[ connec_capacity ];

  return new UnstructuredMesh< Topology::SINGLE >( ndims, cell_type, 0, 
                                                   cell_capacity, connectivity, 
                                                   0, node_capacity, x, y, z );
}


UnstructuredMesh< Topology::MIXED >* 
createExternalMixed( int ndims, IndexType cell_capacity, IndexType node_capacity,
                     IndexType connec_capacity )
{
  double* x = new double[ node_capacity ];
  double* y = AXOM_NULLPTR;
  double* z = AXOM_NULLPTR;
  if ( ndims > 1 )
  {
    y = new double[ node_capacity ];
  }
  if ( ndims > 2 )
  {
    z = new double[ node_capacity ];
  }

  IndexType* connectivity = new IndexType[ connec_capacity ];
  IndexType* offsets = new IndexType[ cell_capacity + 1 ];
  CellType* types = new CellType[ cell_capacity ];

  return new UnstructuredMesh< Topology::MIXED >( ndims, 0, cell_capacity, 
                                                  connec_capacity, connectivity,
                                                  offsets, types, 0, 
                                                  node_capacity, x, y, z );
}


void deleteExternalMesh( UnstructuredMesh< Topology::SINGLE >* mesh )
{
  ASSERT_TRUE( mesh->isExternal() );

  const int ndims = mesh->getDimension();
  double* x = mesh->getCoordinateArray( X_COORDINATE );
  double* y;
  double* z;
  if ( ndims > 1 )
  {
    y = mesh->getCoordinateArray( Y_COORDINATE );
  }
  if ( ndims > 2 )
  {
    z = mesh->getCoordinateArray( Z_COORDINATE );
  }

  IndexType* connectivity = mesh->getCellConnectivityArray();

  delete mesh;
  delete[] x;

  if ( ndims > 1 )
  {
    delete[] y;
  }
  if ( ndims > 2 )
  {
    delete[] z;
  }

  delete[] connectivity;
}


void deleteExternalMesh( UnstructuredMesh< Topology::MIXED >* mesh )
{
  ASSERT_TRUE( mesh->isExternal() );

  const int ndims = mesh->getDimension();
  double* x = mesh->getCoordinateArray( X_COORDINATE );
  double* y;
  double* z;
  if ( ndims > 1 )
  {
    y = mesh->getCoordinateArray( Y_COORDINATE );
  }
  if ( ndims > 2 )
  {
    z = mesh->getCoordinateArray( Z_COORDINATE );
  }

  IndexType* connectivity = mesh->getCellConnectivityArray();
  const IndexType* offsets = mesh->getOffsetsArray();
  const CellType* types = mesh->getTypesArray();

  delete mesh;
  delete[] x;

  if ( ndims > 1 )
  {
    delete[] y;
  }
  if ( ndims > 2 )
  {
    delete[] z;
  }

  delete[] connectivity;
  delete[] offsets;
  delete[] types;
}


void deleteAndDuplicateExternalMesh( UnstructuredMesh< Topology::SINGLE >*& mesh )
{
  ASSERT_TRUE( mesh->isExternal() );

  const int ndims = mesh->getDimension();
  const CellType cell_type = mesh->getCellType();
  const IndexType n_cells = mesh->getNumberOfCells();
  const IndexType cell_capacity = mesh->getCellCapacity();
  IndexType* connectivity = mesh->getCellConnectivityArray();
  const IndexType n_nodes = mesh->getNumberOfNodes();
  const IndexType node_capacity = mesh->getNodeCapacity();
  double* x = mesh->getCoordinateArray( X_COORDINATE );
  double* y = AXOM_NULLPTR;
  double* z = AXOM_NULLPTR;
  if ( ndims > 1 )
  {
    y = mesh->getCoordinateArray( Y_COORDINATE );
  }
  if ( ndims > 2 )
  {
    z = mesh->getCoordinateArray( Z_COORDINATE );
  }

  delete mesh;
  mesh = new UnstructuredMesh< Topology::SINGLE >( ndims, cell_type, n_cells,
                                                   cell_capacity, connectivity,
                                                   n_nodes, node_capacity,
                                                   x, y, z );

  EXPECT_EQ( ndims, mesh->getDimension() );
  EXPECT_EQ( cell_type, mesh->getCellType() );
  EXPECT_EQ( n_cells, mesh->getNumberOfCells() );
  EXPECT_EQ( cell_capacity, mesh->getCellCapacity() );
  EXPECT_EQ( connectivity, mesh->getCellConnectivityArray() );
  EXPECT_EQ( n_nodes, mesh->getNumberOfNodes() );
  EXPECT_EQ( node_capacity, mesh->getNodeCapacity() );
  EXPECT_EQ( x, mesh->getCoordinateArray( X_COORDINATE ) );
  if ( ndims > 1 )
  {
    EXPECT_EQ( y, mesh->getCoordinateArray( Y_COORDINATE ) );
  }
  if ( ndims > 2 )
  {
    EXPECT_EQ( z, mesh->getCoordinateArray( Z_COORDINATE ) );
  }
}


void deleteAndDuplicateExternalMesh( UnstructuredMesh< Topology::MIXED >*& mesh )
{
  ASSERT_TRUE( mesh->isExternal() );

  const int ndims = mesh->getDimension();
  const CellType cell_type = mesh->getCellType();
  const IndexType n_cells = mesh->getNumberOfCells();
  const IndexType cell_capacity = mesh->getCellCapacity();
  const IndexType connec_capacity = mesh->getCellConnectivityCapacity();
  IndexType* connectivity = mesh->getCellConnectivityArray();
  IndexType* offsets = const_cast< IndexType* >( mesh->getOffsetsArray() );
  CellType* types = const_cast< CellType* >( mesh->getTypesArray() );
  const IndexType n_nodes = mesh->getNumberOfNodes();
  const IndexType node_capacity = mesh->getNodeCapacity();
  double* x = mesh->getCoordinateArray( X_COORDINATE );
  double* y = AXOM_NULLPTR;
  double* z = AXOM_NULLPTR;
  if ( ndims > 1 )
  {
    y = mesh->getCoordinateArray( Y_COORDINATE );
  }
  if ( ndims > 2 )
  {
    z = mesh->getCoordinateArray( Z_COORDINATE );
  }

  delete mesh;
  mesh = new UnstructuredMesh< Topology::MIXED >( ndims, n_cells, cell_capacity, 
                                                  connec_capacity, connectivity,
                                                  offsets, types, n_nodes, 
                                                  node_capacity, x, y, z );

  EXPECT_EQ( ndims, mesh->getDimension() );
  EXPECT_EQ( cell_type, mesh->getCellType() );
  EXPECT_EQ( n_cells, mesh->getNumberOfCells() );
  EXPECT_EQ( cell_capacity, mesh->getCellCapacity() );
  EXPECT_EQ( connectivity, mesh->getCellConnectivityArray() );
  EXPECT_EQ( offsets, mesh->getOffsetsArray() );
  EXPECT_EQ( types, mesh->getTypesArray() );
  EXPECT_EQ( n_nodes, mesh->getNumberOfNodes() );
  EXPECT_EQ( node_capacity, mesh->getNodeCapacity() );
  EXPECT_EQ( x, mesh->getCoordinateArray( X_COORDINATE ) );
  if ( ndims > 1 )
  {
    EXPECT_EQ( y, mesh->getCoordinateArray( Y_COORDINATE ) );
  }
  if ( ndims > 2 )
  {
    EXPECT_EQ( z, mesh->getCoordinateArray( Z_COORDINATE ) );
  }
}


template < Topology TOPO >
void deleteAndDuplicateSidreMesh( UnstructuredMesh< TOPO >*& mesh )
{
  ASSERT_TRUE( mesh->isInSidre() );

  const int ndims = mesh->getDimension();
  const CellType cell_type = mesh->getCellType();
  const IndexType n_cells = mesh->getNumberOfCells();
  const IndexType cell_capacity = mesh->getCellCapacity();
  const IndexType n_nodes = mesh->getNumberOfNodes();
  const IndexType node_capacity = mesh->getNodeCapacity();

  sidre::Group* group = mesh->getSidreGroup();
  const std::string topo = mesh->getTopologyName();

  delete mesh;
  mesh = new UnstructuredMesh< TOPO >( group, topo );

  EXPECT_EQ( ndims, mesh->getDimension() );
  EXPECT_EQ( cell_type, mesh->getCellType() );
  EXPECT_EQ( n_cells, mesh->getNumberOfCells() );
  EXPECT_EQ( cell_capacity, mesh->getCellCapacity() );
  EXPECT_EQ( n_nodes, mesh->getNumberOfNodes() );
  EXPECT_EQ( node_capacity, mesh->getNodeCapacity() );
}


inline double getCoordValue( IndexType ndims, IndexType i, int dim )
{ return (ndims * i + dim) * PI; }


template < Topology TOPO >
void check_append_nodes( const UnstructuredMesh< TOPO >* mesh )
{
  IndexType n_nodes = mesh->getNumberOfNodes();
  
  const int ndims = mesh->getDimension();
  const double* x = mesh->getCoordinateArray( X_COORDINATE );

  /* Check using pointers */
  for ( IndexType i = 0; i < n_nodes; ++i )
  {
    EXPECT_EQ( x[ i ], getCoordValue( ndims, i, 0 ) );
  }
  if ( ndims > 1 )
  {
    const double* y = mesh->getCoordinateArray( Y_COORDINATE );
    for ( IndexType i = 0; i < n_nodes; ++i )
    {
      EXPECT_EQ( y[ i ], getCoordValue( ndims, i, 1 ) );
    }
  }
  if ( ndims > 2 )
  {
    const double* z = mesh->getCoordinateArray( Z_COORDINATE );
    for ( IndexType i = 0; i < n_nodes; ++i )
    {
      EXPECT_EQ( z[ i ], getCoordValue( ndims, i, 2 ) );
    }
  }

  /* Check using getNode */
  double coords[3];
  for ( IndexType i = 0; i < n_nodes; ++i )
  {
    mesh->getNode( i, coords );
    for ( int dim = 0; dim < ndims; ++dim )
    {
      EXPECT_EQ( coords[ dim ], getCoordValue( ndims, i, dim ) );
    }
  }

  /* Check using getNodeCoordinate */
  for ( IndexType i = 0; i < n_nodes; ++i )
  {
    for ( int dim = 0; dim < ndims; ++dim )
    {
      EXPECT_EQ( mesh->getNodeCoordinate( i, dim ), 
                 getCoordValue( ndims, i, dim ) );
    }
  }
}


template < Topology TOPO >
void append_node_single( UnstructuredMesh< TOPO >* mesh, IndexType n_nodes )
{
  const int ndims = mesh->getDimension();
  IndexType cur_n_nodes = mesh->getNumberOfNodes();

  if ( ndims == 1 )
  {
    for ( IndexType i = 0; i < n_nodes; ++i )
    {
      mesh->appendNode( getCoordValue( ndims, cur_n_nodes, 0 ) );
      EXPECT_EQ( ++cur_n_nodes, mesh->getNumberOfNodes() );
    }
  }
  else if ( ndims == 2 )
  {
    for ( IndexType i = 0; i < n_nodes; ++i )
    {
      mesh->appendNode( getCoordValue( ndims, cur_n_nodes, 0 ), 
                        getCoordValue( ndims, cur_n_nodes, 1 ) );
      EXPECT_EQ( ++cur_n_nodes, mesh->getNumberOfNodes() );
    }
  }
  else
  {
    for ( IndexType i = 0; i < n_nodes; ++i )
    {
      mesh->appendNode( getCoordValue( ndims, cur_n_nodes, 0 ), 
                        getCoordValue( ndims, cur_n_nodes, 1 ),
                        getCoordValue( ndims, cur_n_nodes, 2 ) );
      EXPECT_EQ( ++cur_n_nodes, mesh->getNumberOfNodes() );
    }
  }
}


template < Topology TOPO >
void append_node_structs( UnstructuredMesh< TOPO >* mesh, IndexType n_nodes,
                     double* coords=AXOM_NULLPTR )
{
  const int ndims = mesh->getDimension();
  IndexType cur_n_nodes = mesh->getNumberOfNodes();
  bool free_coords = false;
  if ( coords == AXOM_NULLPTR )
  {
    coords = new double[ ndims * n_nodes ];
    free_coords = true;
  }

  /* Initialize the coords in array of structs layout */
  for ( IndexType i = 0; i < n_nodes; ++i )
  {
    for ( IndexType dim = 0; dim < ndims; ++dim )
    {
      coords[ ndims * i + dim ] = getCoordValue( ndims, i + cur_n_nodes, dim );
    }
  }

  /* Append multiple nodes at once using the array of structs layout. */
  mesh->appendNodes( coords, n_nodes );
  cur_n_nodes += n_nodes;
  EXPECT_EQ( cur_n_nodes, mesh->getNumberOfNodes() );

  if ( free_coords )
  {
    delete[] coords;
    coords = AXOM_NULLPTR;
  }
}


template < Topology TOPO >
void append_node_arrays( UnstructuredMesh< TOPO >* mesh, IndexType n_nodes,
                    double* coords=AXOM_NULLPTR )
{
  const int ndims = mesh->getDimension();
  IndexType cur_n_nodes = mesh->getNumberOfNodes();
  bool free_coords = false;
  if ( coords == AXOM_NULLPTR )
  {
    coords = new double[ ndims * n_nodes ];
    free_coords = true;
  }

  /* Initialize the coords in struct of arrays layout */
  for ( IndexType dim = 0; dim < ndims; ++dim )
  {
    for ( IndexType i = 0; i < n_nodes; ++i )
    {
      coords[ dim * n_nodes + i ] = getCoordValue( ndims, i + cur_n_nodes, dim );
    }
  }

  /* Append multiple nodes at once using the struct of arrays layout. */
  if ( ndims == 1 )
  {
    mesh->appendNodes( coords, n_nodes );
    cur_n_nodes += n_nodes;
    EXPECT_EQ( cur_n_nodes, mesh->getNumberOfNodes() );
  }
  else if ( ndims == 2 )
  {
    const double* x = coords;
    const double* y = coords + n_nodes;
    mesh->appendNodes( x, y, n_nodes );
    cur_n_nodes += n_nodes;
    EXPECT_EQ( cur_n_nodes, mesh->getNumberOfNodes() );
  }
  else
  {
    const double* x = coords;
    const double* y = coords + n_nodes;
    const double* z = coords + 2 * n_nodes;
    mesh->appendNodes( x, y, z, n_nodes );
    cur_n_nodes += n_nodes;
    EXPECT_EQ( cur_n_nodes, mesh->getNumberOfNodes() );
  }

  if ( free_coords )
  {
    delete[] coords;
    coords = AXOM_NULLPTR;
  }
}


template < Topology TOPO >
void append_nodes( UnstructuredMesh< TOPO >* mesh, IndexType n_nodes )
{
  ASSERT_EQ( n_nodes % 3, 0 );

  const int ndims = mesh->getDimension();
  IndexType cur_n_nodes = mesh->getNumberOfNodes();
  ASSERT_EQ( cur_n_nodes, 0 );
  EXPECT_TRUE( mesh->empty() );

  /* Append one node at a time */
  append_node_single( mesh, n_nodes / 3 );
  cur_n_nodes += n_nodes / 3;
  EXPECT_EQ( cur_n_nodes, mesh->getNumberOfNodes() );

  /* Allocate the coords array of nodes to append */
  double* coords = new double[ ndims * n_nodes / 3 ];
  
  /* Append multiple nodes at once using the array of structs layout. */
  append_node_structs( mesh, n_nodes / 3, coords );
  cur_n_nodes += n_nodes / 3;
  EXPECT_EQ( cur_n_nodes, 2 * n_nodes / 3 );
  EXPECT_EQ( cur_n_nodes, mesh->getNumberOfNodes() );

  /* Append multiple nodes at once using the struct of arrays layout. */
  append_node_arrays( mesh, n_nodes / 3, coords );
  cur_n_nodes += n_nodes / 3;
  EXPECT_EQ( cur_n_nodes, n_nodes );
  EXPECT_EQ( cur_n_nodes, mesh->getNumberOfNodes() );

  check_append_nodes( mesh );
}


inline IndexType getCellConnecValue( IndexType cur_n_cells, IndexType cur_node )
{ return cur_n_cells * MAX_NUM_NODES + cur_node; }


inline CellType getCellType( IndexType cur_n_cells )
{
  CellType type = QUAD;
  if ( cur_n_cells % 2 != 0 )
  {
    type = TRIANGLE;
  }
  return type;
}

constexpr IndexType getConnectivitySize( IndexType n_cells )
{ return ( n_cells * ( 4 + 3 ) ) / 2; }


void check_append_cells( const UnstructuredMesh< Topology::SINGLE >* mesh )
{
  IndexType n_cells = mesh->getNumberOfCells();

  const CellType type = mesh->getCellType();
  const int nodes_per_cell = cell_info[ type ].num_nodes;

  /* Check using pointers */
  const IndexType* connectivity = mesh->getCellConnectivityArray();
  for ( IndexType i = 0; i < n_cells; ++i )
  {
    EXPECT_EQ( type, mesh->getCellType( i ) );
    EXPECT_EQ( nodes_per_cell, mesh->getNumberOfCellNodes( i ) );

    for ( IndexType j = 0; j < nodes_per_cell; ++j )
    {
      EXPECT_EQ( connectivity[ i * nodes_per_cell + j ], 
               getCellConnecValue( i, j ) );
    }
  }

  /* Check using getCell */
  IndexType cell[ nodes_per_cell ];
  for ( IndexType i = 0; i < n_cells; ++i )
  {
    mesh->getCell( i, cell );
    for ( IndexType j = 0; j < nodes_per_cell; ++j )
    {
      EXPECT_EQ( cell[ j ], getCellConnecValue( i, j ) );
    }
  }

  /* Check using the other getCell */
  for ( IndexType i = 0; i < n_cells; ++i )
  {
    const IndexType* cell = mesh->getCell( i );
    for ( IndexType j = 0; j < nodes_per_cell; ++j )
    {
      EXPECT_EQ( cell[ j ], getCellConnecValue( i, j ) );
    }
  }
}


void check_append_cells( const UnstructuredMesh< Topology::MIXED >* mesh )
{
  IndexType n_cells = mesh->getNumberOfCells();

  /* Check using pointers */
  const IndexType* connectivity = mesh->getCellConnectivityArray();
  const IndexType* offsets = mesh->getOffsetsArray();
  const CellType* types = mesh->getTypesArray();
  for ( IndexType i = 0; i < n_cells; ++i )
  {
    CellType expected_type = getCellType( i );
    IndexType expected_n_nodes = cell_info[ expected_type ].num_nodes;
    ASSERT_EQ( expected_type, types[ i ] );
    ASSERT_EQ( expected_type, mesh->getCellType( i ) );
    ASSERT_EQ( expected_n_nodes, offsets[ i + 1 ] - offsets [ i ] );
    ASSERT_EQ( expected_n_nodes, mesh->getNumberOfCellNodes( i ) );

    IndexType offset = offsets[ i ];
    for ( IndexType j = 0; j < expected_n_nodes; ++j )
    {
      EXPECT_EQ( connectivity[ offset + j ], 
                 getCellConnecValue( i, j ) );
    }
  }

  /* Check using getCell */
  IndexType cell[ MAX_NUM_NODES ];
  for ( IndexType i = 0; i < n_cells; ++i )
  {
    mesh->getCell( i, cell );
    IndexType num_nodes = mesh->getNumberOfCellNodes( i );
    for ( IndexType j = 0; j < num_nodes; ++j )
    {
      EXPECT_EQ( cell[ j ], getCellConnecValue( i, j ) );
    }
  }

  /* Check using the other getCell */
  for ( IndexType i = 0; i < n_cells; ++i )
  {
    const IndexType* cell = mesh->getCell( i );
    IndexType num_nodes = mesh->getNumberOfCellNodes( i );
    for ( IndexType j = 0; j < num_nodes; ++j )
    {
      EXPECT_EQ( cell[ j ], getCellConnecValue( i, j ) );
    }
  }
}

inline void getCellConnec( IndexType cur_n_cells, IndexType nodes_per_cell, 
                           IndexType* connec )
{
  for ( IndexType i = 0; i < nodes_per_cell; ++i )
  {
    connec[ i ] = getCellConnecValue( cur_n_cells, i);
  }
}


inline CellType getCellConnec( IndexType cur_n_cells, IndexType* connec )
{
  CellType type = getCellType( cur_n_cells );
  const IndexType nodes_per_cell = cell_info[ type ].num_nodes;
  for ( IndexType i = 0; i < nodes_per_cell; ++i )
  {
    connec[ i ] = getCellConnecValue( cur_n_cells, i);
  }

  return type;
}


void append_cell_single( UnstructuredMesh< Topology::SINGLE >* mesh, 
                         IndexType n_cells )
{
  IndexType cur_n_cells = mesh->getNumberOfCells();
  const IndexType nodes_per_cell = cell_info[ mesh->getCellType() ].num_nodes;

  IndexType connec[ nodes_per_cell ];
  for ( IndexType i = 0; i < n_cells; ++i )
  {
    getCellConnec( cur_n_cells, nodes_per_cell, connec );
    mesh->appendCell( connec );
    EXPECT_EQ( ++cur_n_cells, mesh->getNumberOfCells() );
    EXPECT_EQ( cur_n_cells * nodes_per_cell, 
               mesh->getCellConnectivitySize() );
  }
}

void append_cell_single( UnstructuredMesh< Topology::MIXED >* mesh, 
                         IndexType n_cells )
{
  IndexType cur_n_cells = mesh->getNumberOfCells();
  IndexType cur_connec_size = mesh->getCellConnectivitySize();

  IndexType connec[ MAX_NUM_NODES ];
  CellType type;
  for ( IndexType i = 0; i < n_cells; ++i )
  {
    type = getCellConnec( cur_n_cells, connec );
    cur_connec_size += cell_info[ type ].num_nodes;
    mesh->appendCell( connec, type );
    EXPECT_EQ( ++cur_n_cells, mesh->getNumberOfCells() );
    EXPECT_EQ( cur_connec_size, mesh->getCellConnectivitySize() );
  }
}


void append_cell_multiple( UnstructuredMesh< Topology::SINGLE >* mesh, 
                           IndexType n_cells )
{
  IndexType cur_n_cells = mesh->getNumberOfCells();
  IndexType cur_connec_size = mesh->getCellConnectivitySize();
  const IndexType nodes_per_cell = mesh->getNumberOfCellNodes();
  IndexType* connectivity = new IndexType[ nodes_per_cell * n_cells ];

  for ( IndexType i = 0; i < n_cells; ++i )
  {
    IndexType* cur_connec = connectivity + i * nodes_per_cell;
    getCellConnec( cur_n_cells + i, nodes_per_cell, cur_connec );
  }

  mesh->appendCells( connectivity, n_cells );
  cur_n_cells += n_cells;
  cur_connec_size += n_cells * nodes_per_cell;
  EXPECT_EQ( cur_n_cells, mesh->getNumberOfCells() );
  EXPECT_EQ( cur_connec_size, mesh->getCellConnectivitySize() );
}


void append_cell_multiple( UnstructuredMesh< Topology::MIXED >* mesh, 
                           IndexType n_cells )
{
  IndexType cur_n_cells = mesh->getNumberOfCells();
  IndexType cur_connec_size = mesh->getCellConnectivitySize();
  IndexType* connectivity = new IndexType[ getConnectivitySize( n_cells ) ];
  IndexType* offsets = new IndexType[ n_cells + 1 ];
  CellType* types = new CellType[ n_cells ];

  offsets[0] = 0;

  for ( IndexType i = 0; i < n_cells; ++i )
  {
    IndexType* cur_connec = connectivity + offsets[ i ];
    types[ i ] = getCellConnec( cur_n_cells + i, cur_connec );
    offsets[ i + 1 ] = offsets[ i ] + cell_info[ types[ i ] ].num_nodes;
  }

  mesh->appendCells( connectivity, n_cells, offsets, types );
  cur_n_cells += n_cells;
  cur_connec_size += getConnectivitySize( n_cells );
  EXPECT_EQ( cur_n_cells, mesh->getNumberOfCells() );
  EXPECT_EQ( cur_connec_size, mesh->getCellConnectivitySize() );
}


template < Topology TOPO >
void append_cells( UnstructuredMesh< TOPO >* mesh, IndexType n_cells )
{
  ASSERT_EQ( n_cells % 2, 0 );

  IndexType cur_n_cells = mesh->getNumberOfCells();
  ASSERT_EQ( cur_n_cells, 0 );

  /* Append cells one at a time */
  append_cell_single( mesh, n_cells / 2 );
  cur_n_cells += n_cells / 2;
  EXPECT_EQ( cur_n_cells, mesh->getNumberOfCells() );

  /* Append cells all at once */
  append_cell_multiple( mesh, n_cells / 2 );
  cur_n_cells += n_cells / 2;
  EXPECT_EQ( cur_n_cells, n_cells );
  EXPECT_EQ( cur_n_cells, mesh->getNumberOfCells() );
}


inline double getNewCoordValue( IndexType ndims, IndexType i, int dim )
{ return (ndims * i + dim) * E; }

template < Topology TOPO >
void check_set_nodes( const UnstructuredMesh< TOPO >* mesh )
{
  IndexType n_nodes = mesh->getNumberOfNodes();
  
  const int ndims = mesh->getDimension();
  const double* x = mesh->getCoordinateArray( X_COORDINATE );

  /* Check using pointers */
  for ( IndexType i = 0; i < n_nodes; ++i )
  {
    EXPECT_EQ( x[ i ], getNewCoordValue( ndims, i, 0 ) );
  }
  if ( ndims > 1 )
  {
    const double* y = mesh->getCoordinateArray( Y_COORDINATE );
    for ( IndexType i = 0; i < n_nodes; ++i )
    {
      EXPECT_EQ( y[ i ], getNewCoordValue( ndims, i, 1 ) );
    }
  }
  if ( ndims > 2 )
  {
    const double* z = mesh->getCoordinateArray( Z_COORDINATE );
    for ( IndexType i = 0; i < n_nodes; ++i )
    {
      EXPECT_EQ( z[ i ], getNewCoordValue( ndims, i, 2 ) );
    }
  }

  /* Check using getNode */
  double coords[3];
  for ( IndexType i = 0; i < n_nodes; ++i )
  {
    mesh->getNode( i, coords );
    for ( int dim = 0; dim < ndims; ++dim )
    {
      EXPECT_EQ( coords[ dim ], getNewCoordValue( ndims, i, dim ) );
    }
  }

  /* Check using getNodeCoordinate */
  for ( IndexType i = 0; i < n_nodes; ++i )
  {
    for ( int dim = 0; dim < ndims; ++dim )
    {
      EXPECT_EQ( mesh->getNodeCoordinate( i, dim ), 
                 getNewCoordValue( ndims, i, dim ) );
    }
  }
}


template < Topology TOPO >
void set_nodes( UnstructuredMesh< TOPO >* mesh )
{
  const int ndims = mesh->getDimension();
  const IndexType n_nodes = mesh->getNumberOfNodes();

  double* x = mesh->getCoordinateArray( X_COORDINATE );
  for ( IndexType i = 0; i < n_nodes; ++i )
  {
    x[ i ] = getNewCoordValue( ndims, i, 0 );
  }

  if ( ndims > 1 )
  {
    double* y = mesh->getCoordinateArray( Y_COORDINATE );
    for ( IndexType i = 0; i < n_nodes; ++i )
    {
      y[ i ] = getNewCoordValue( ndims, i, 1 );
    }
  }
  if ( ndims > 2 )
  {
    double* z = mesh->getCoordinateArray( Z_COORDINATE );
    for ( IndexType i = 0; i < n_nodes; ++i )
    {
      z[ i ] = getNewCoordValue( ndims, i, 2 );
    }
  }

  EXPECT_EQ( n_nodes, mesh->getNumberOfNodes() );
  check_set_nodes( mesh );
}


inline IndexType getNewCellConnecValue( IndexType cur_n_cells, 
                                        IndexType cur_node )
{ return (cur_n_cells * MAX_NUM_NODES + cur_node) * cur_node; }


void check_set_cells( const UnstructuredMesh< Topology::SINGLE >* mesh )
{
  IndexType n_cells = mesh->getNumberOfCells();

  const CellType type = mesh->getCellType();
  const int nodes_per_cell = cell_info[ type ].num_nodes;

  /* Check using pointers */
  const IndexType* connectivity = mesh->getCellConnectivityArray();
  for ( IndexType i = 0; i < n_cells; ++i )
  {
    EXPECT_EQ( type, mesh->getCellType( i ) );
    EXPECT_EQ( nodes_per_cell, mesh->getNumberOfCellNodes( i ) );

    for ( IndexType j = 0; j < nodes_per_cell; ++j )
    {
      EXPECT_EQ( connectivity[ i * nodes_per_cell + j ], 
               getNewCellConnecValue( i, j ) );
    }
  }

  /* Check using getCell */
  IndexType cell[ nodes_per_cell ];
  for ( IndexType i = 0; i < n_cells; ++i )
  {
    mesh->getCell( i, cell );
    for ( IndexType j = 0; j < nodes_per_cell; ++j )
    {
      EXPECT_EQ( cell[ j ], getNewCellConnecValue( i, j ) );
    }
  }

  /* Check using the other getCell */
  for ( IndexType i = 0; i < n_cells; ++i )
  {
    const IndexType* cell = mesh->getCell( i );
    for ( IndexType j = 0; j < nodes_per_cell; ++j )
    {
      EXPECT_EQ( cell[ j ], getNewCellConnecValue( i, j ) );
    }
  }
}


void check_set_cells( const UnstructuredMesh< Topology::MIXED >* mesh )
{
  IndexType n_cells = mesh->getNumberOfCells();

  /* Check using pointers */
  const IndexType* connectivity = mesh->getCellConnectivityArray();
  const IndexType* offsets = mesh->getOffsetsArray();
  const CellType* types = mesh->getTypesArray();
  for ( IndexType i = 0; i < n_cells; ++i )
  {
    CellType expected_type = getCellType( i );
    IndexType expected_n_nodes = cell_info[ expected_type ].num_nodes;
    ASSERT_EQ( expected_type, types[ i ] );
    ASSERT_EQ( expected_type, mesh->getCellType( i ) );
    ASSERT_EQ( expected_n_nodes, offsets[ i + 1 ] - offsets [ i ] );
    ASSERT_EQ( expected_n_nodes, mesh->getNumberOfCellNodes( i ) );

    IndexType offset = offsets[ i ];
    for ( IndexType j = 0; j < expected_n_nodes; ++j )
    {
      EXPECT_EQ( connectivity[ offset + j ], 
                 getNewCellConnecValue( i, j ) );
    }
  }

  /* Check using getCell */
  IndexType cell[ MAX_NUM_NODES ];
  for ( IndexType i = 0; i < n_cells; ++i )
  {
    mesh->getCell( i, cell );
    IndexType num_nodes = mesh->getNumberOfCellNodes( i );
    for ( IndexType j = 0; j < num_nodes; ++j )
    {
      EXPECT_EQ( cell[ j ], getNewCellConnecValue( i, j ) );
    }
  }

  /* Check using the other getCell */
  for ( IndexType i = 0; i < n_cells; ++i )
  {
    const IndexType* cell = mesh->getCell( i );
    IndexType num_nodes = mesh->getNumberOfCellNodes( i );
    for ( IndexType j = 0; j < num_nodes; ++j )
    {
      EXPECT_EQ( cell[ j ], getNewCellConnecValue( i, j ) );
    }
  }
}


template < Topology TOPO >
void set_cells( UnstructuredMesh< TOPO >* mesh )
{
  const IndexType n_cells = mesh->getNumberOfCells();

  for ( IndexType i = 0; i < n_cells; ++i )
  {
    IndexType* cur_cell = mesh->getCell( i );
    const IndexType n_nodes = mesh->getNumberOfCellNodes( i );
    for ( IndexType j = 0; j < n_nodes; ++j )
    {
      cur_cell[ j ] = getNewCellConnecValue( i, j );
    }
  }

  EXPECT_EQ( n_cells, mesh->getNumberOfCells() );
  check_set_cells( mesh );
}



}   /* end namespace internal */

//------------------------------------------------------------------------------
TEST( mint_mesh_particle_mesh, appendNodesSingle )
{
#ifdef MINT_USE_SIDRE
  sidre::DataStore ds;
  sidre::Group* root = ds.getRoot();
  constexpr int STRIDE = 3;
#else
  constexpr int STRIDE = 2;
#endif

  constexpr IndexType N_NODES = 1200;
  constexpr IndexType N_CELLS = 0;

  constexpr int N_MESHES = STRIDE * 3;
  UnstructuredMesh< Topology::SINGLE >* meshes[N_MESHES];

  /* Create the meshes. */
  int cur_mesh = 0;
  for ( int dim = 1; dim <= 3; ++dim )
  {
    meshes[ cur_mesh++ ] = new UnstructuredMesh< Topology::SINGLE >( dim, QUAD );
    meshes[ cur_mesh++ ] = internal::createExternalSingle( dim, QUAD, N_CELLS, 
                                                         N_NODES );

#ifdef MINT_USE_SIDRE
    const std::string topo = "t" + std::to_string( dim );
    const std::string coordset = "c" + std::to_string( dim );
    meshes[ cur_mesh++ ] = new UnstructuredMesh< Topology::SINGLE >( dim, QUAD, 
                                                                     root, topo,
                                                                     coordset );
#endif
  }

  /* Check that append functions properly. */
  for ( int i = 0; i < N_MESHES; ++i )
  {
    internal::append_nodes( meshes[ i ], N_NODES );
  }

  /* Check that the external meshes function properly. */
  for ( int i = 1; i < N_MESHES; i += STRIDE )
  {
    internal::deleteAndDuplicateExternalMesh( meshes[ i ] );
    internal::check_append_nodes( meshes[ i ] );
  }

  /* Check that the sidre meshes function properly. */
#ifdef MINT_USE_SIDRE
  for ( int i = 2; i < N_MESHES; i += STRIDE )
  {
    internal::deleteAndDuplicateSidreMesh( meshes[ i ] );
    internal::check_append_nodes( meshes[ i ] );
  }
#endif

  /* Delete the meshes. */
  cur_mesh = 0;
  for ( int dim = 0; dim < 3; ++dim )
  {
    delete meshes[ cur_mesh++ ];
    internal::deleteExternalMesh( meshes[ cur_mesh++ ] );

#ifdef MINT_USE_SIDRE
    delete meshes[ cur_mesh++ ];
#endif
  }
}

//------------------------------------------------------------------------------
TEST( mint_mesh_particle_mesh, appendNodesMixed )
{
#ifdef MINT_USE_SIDRE
  sidre::DataStore ds;
  sidre::Group* root = ds.getRoot();
  constexpr int STRIDE = 3;
#else
  constexpr int STRIDE = 2;
#endif

  constexpr IndexType N_NODES = 1200;
  constexpr IndexType N_CELLS = 0;
  constexpr IndexType CONNEC_SIZE = 0;

  constexpr int N_MESHES = STRIDE * 3;
  UnstructuredMesh< Topology::MIXED >* meshes[N_MESHES];

  /* Create the meshes. */
  int cur_mesh = 0;
  for ( int dim = 1; dim <= 3; ++dim )
  {
    meshes[ cur_mesh++ ] = new UnstructuredMesh< Topology::MIXED >( dim );
    meshes[ cur_mesh++ ] = internal::createExternalMixed( dim, N_CELLS, N_NODES, 
                                                         CONNEC_SIZE );

#ifdef MINT_USE_SIDRE
    const std::string topo = "t" + std::to_string( dim );
    const std::string coordset = "c" + std::to_string( dim );
    meshes[ cur_mesh++ ] = new UnstructuredMesh< Topology::MIXED >( dim, root, 
                                                                    topo,
                                                                    coordset );
#endif
  }

  /* Check that append functions properly. */
  for ( int i = 0; i < N_MESHES; ++i )
  {
    internal::append_nodes( meshes[ i ], N_NODES );
  }

  /* Check that the external meshes function properly. */
  for ( int i = 1; i < N_MESHES; i += STRIDE )
  {
    internal::deleteAndDuplicateExternalMesh( meshes[ i ] );
    internal::check_append_nodes( meshes[ i ] );
  }

  /* Check that the sidre meshes function properly. */
#ifdef MINT_USE_SIDRE
  for ( int i = 2; i < N_MESHES; i += STRIDE )
  {
    internal::deleteAndDuplicateSidreMesh( meshes[ i ] );
    internal::check_append_nodes( meshes[ i ] );
  }
#endif

  /* Delete the meshes. */
  cur_mesh = 0;
  for ( int dim = 0; dim < 3; ++dim )
  {
    delete meshes[ cur_mesh++ ];
    internal::deleteExternalMesh( meshes[ cur_mesh++ ] );

#ifdef MINT_USE_SIDRE
    delete meshes[ cur_mesh++ ];
#endif
  }
}

//------------------------------------------------------------------------------
TEST( mint_mesh_particle_mesh, appendCellsSingle )
{
#ifdef MINT_USE_SIDRE
  sidre::DataStore ds;
  sidre::Group* root = ds.getRoot();
  constexpr int STRIDE = 3;
#else
  constexpr int STRIDE = 2;
#endif

  constexpr IndexType N_NODES = 12 * 10 * 10;
  constexpr IndexType N_CELLS = 900;

  constexpr int N_MESHES = STRIDE * 3;
  UnstructuredMesh< Topology::SINGLE >* meshes[N_MESHES];

  /* Create the meshes. */
  int cur_mesh = 0;
  for ( int dim = 1; dim <= 3; ++dim )
  {
    meshes[ cur_mesh++ ] = new UnstructuredMesh< Topology::SINGLE >( dim, QUAD );
    meshes[ cur_mesh++ ] = internal::createExternalSingle( dim, QUAD, N_CELLS, 
                                                           N_NODES );

#ifdef MINT_USE_SIDRE
    const std::string topo = "t" + std::to_string( dim );
    const std::string coordset = "c" + std::to_string( dim );
    meshes[ cur_mesh++ ] = new UnstructuredMesh< Topology::SINGLE >( dim, QUAD, 
                                                                     root, topo,
                                                                     coordset );
#endif
  }

  /* Check that append functions properly. */
  for ( int i = 0; i < N_MESHES; ++i )
  {
    internal::append_nodes( meshes[ i ], N_NODES );
    internal::append_cells( meshes[ i ], N_CELLS );
  }

  /* Check that the external meshes function properly. */
  for ( int i = 1; i < N_MESHES; i += STRIDE )
  {
    internal::deleteAndDuplicateExternalMesh( meshes[ i ] );
    internal::check_append_nodes( meshes[ i ] );
    internal::check_append_cells( meshes[ i ] );
  }

  /* Check that the sidre meshes function properly. */
#ifdef MINT_USE_SIDRE
  for ( int i = 2; i < N_MESHES; i += STRIDE )
  {
    internal::deleteAndDuplicateSidreMesh( meshes[ i ] );
    internal::check_append_nodes( meshes[ i ] );
    internal::check_append_cells( meshes[ i ] );
  }
#endif

  /* Delete the meshes. */
  cur_mesh = 0;
  for ( int dim = 0; dim < 3; ++dim )
  {
    delete meshes[ cur_mesh++ ];
    internal::deleteExternalMesh( meshes[ cur_mesh++ ] );

#ifdef MINT_USE_SIDRE
    delete meshes[ cur_mesh++ ];
#endif
  }
}

//------------------------------------------------------------------------------
TEST( mint_mesh_particle_mesh, appendCellsMixed )
{
#ifdef MINT_USE_SIDRE
  sidre::DataStore ds;
  sidre::Group* root = ds.getRoot();
  constexpr int STRIDE = 3;
#else
  constexpr int STRIDE = 2;
#endif

  constexpr IndexType N_NODES = 1200;
  constexpr IndexType N_CELLS = 900;
  constexpr IndexType CONNEC_SIZE = internal::getConnectivitySize( N_CELLS );

  constexpr int N_MESHES = STRIDE * 3;
  UnstructuredMesh< Topology::MIXED >* meshes[N_MESHES];

  /* Create the meshes. */
  int cur_mesh = 0;
  for ( int dim = 1; dim <= 3; ++dim )
  {
    meshes[ cur_mesh++ ] = new UnstructuredMesh< Topology::MIXED >( dim );
    meshes[ cur_mesh++ ] = internal::createExternalMixed( dim, N_CELLS, N_NODES, 
                                                         CONNEC_SIZE );

#ifdef MINT_USE_SIDRE
    const std::string topo = "t" + std::to_string( dim );
    const std::string coordset = "c" + std::to_string( dim );
    meshes[ cur_mesh++ ] = new UnstructuredMesh< Topology::MIXED >( dim, root, 
                                                                    topo,
                                                                    coordset );
#endif
  }

  /* Check that append functions properly. */
  for ( int i = 0; i < N_MESHES; ++i )
  {
    internal::append_nodes( meshes[ i ], N_NODES );
    internal::append_cells( meshes[ i ], N_CELLS );
  }

  /* Check that the external meshes function properly. */
  for ( int i = 1; i < N_MESHES; i += STRIDE )
  {
    internal::deleteAndDuplicateExternalMesh( meshes[ i ] );
    internal::check_append_nodes( meshes[ i ] );
    internal::check_append_cells( meshes[ i ] );
  }

  /* Check that the sidre meshes function properly. */
#ifdef MINT_USE_SIDRE
  for ( int i = 2; i < N_MESHES; i += STRIDE )
  {
    internal::deleteAndDuplicateSidreMesh( meshes[ i ] );
    internal::check_append_nodes( meshes[ i ] );
    internal::check_append_cells( meshes[ i ] );
  }
#endif

  /* Delete the meshes. */
  cur_mesh = 0;
  for ( int dim = 0; dim < 3; ++dim )
  {
    delete meshes[ cur_mesh++ ];
    internal::deleteExternalMesh( meshes[ cur_mesh++ ] );

#ifdef MINT_USE_SIDRE
    delete meshes[ cur_mesh++ ];
#endif
  }
}

//------------------------------------------------------------------------------
TEST( mint_mesh_particle_mesh, setNodesSingle )
{
#ifdef MINT_USE_SIDRE
  sidre::DataStore ds;
  sidre::Group* root = ds.getRoot();
  constexpr int STRIDE = 3;
#else
  constexpr int STRIDE = 2;
#endif

  constexpr IndexType N_NODES = 1200;
  constexpr IndexType N_CELLS = 0;

  constexpr int N_MESHES = STRIDE * 3;
  UnstructuredMesh< Topology::SINGLE >* meshes[N_MESHES];

  /* Create the meshes. */
  int cur_mesh = 0;
  for ( int dim = 1; dim <= 3; ++dim )
  {
    meshes[ cur_mesh++ ] = new UnstructuredMesh< Topology::SINGLE >( dim, QUAD );
    meshes[ cur_mesh++ ] = internal::createExternalSingle( dim, QUAD, N_CELLS, 
                                                         N_NODES );

#ifdef MINT_USE_SIDRE
    const std::string topo = "t" + std::to_string( dim );
    const std::string coordset = "c" + std::to_string( dim );
    meshes[ cur_mesh++ ] = new UnstructuredMesh< Topology::SINGLE >( dim, QUAD, 
                                                                     root, topo,
                                                                     coordset );
#endif
  }

  /* Check that append functions properly. */
  for ( int i = 0; i < N_MESHES; ++i )
  {
    internal::append_nodes( meshes[ i ], N_NODES );
    internal::set_nodes( meshes[ i ] );
  }

  /* Check that the external meshes function properly. */
  for ( int i = 1; i < N_MESHES; i += STRIDE )
  {
    internal::deleteAndDuplicateExternalMesh( meshes[ i ] );
    internal::check_set_nodes( meshes[ i ] );
  }

  /* Check that the sidre meshes function properly. */
#ifdef MINT_USE_SIDRE
  for ( int i = 2; i < N_MESHES; i += STRIDE )
  {
    internal::deleteAndDuplicateSidreMesh( meshes[ i ] );
    internal::check_set_nodes( meshes[ i ] );
  }
#endif

  /* Delete the meshes. */
  cur_mesh = 0;
  for ( int dim = 0; dim < 3; ++dim )
  {
    delete meshes[ cur_mesh++ ];
    internal::deleteExternalMesh( meshes[ cur_mesh++ ] );

#ifdef MINT_USE_SIDRE
    delete meshes[ cur_mesh++ ];
#endif
  }
}

//------------------------------------------------------------------------------
TEST( mint_mesh_particle_mesh, setNodesMixed )
{
#ifdef MINT_USE_SIDRE
  sidre::DataStore ds;
  sidre::Group* root = ds.getRoot();
  constexpr int STRIDE = 3;
#else
  constexpr int STRIDE = 2;
#endif

  constexpr IndexType N_NODES = 1200;
  constexpr IndexType N_CELLS = 0;
  constexpr IndexType CONNEC_SIZE = 0;

  constexpr int N_MESHES = STRIDE * 3;
  UnstructuredMesh< Topology::MIXED >* meshes[N_MESHES];

  /* Create the meshes. */
  int cur_mesh = 0;
  for ( int dim = 1; dim <= 3; ++dim )
  {
    meshes[ cur_mesh++ ] = new UnstructuredMesh< Topology::MIXED >( dim );
    meshes[ cur_mesh++ ] = internal::createExternalMixed( dim, N_CELLS, N_NODES, 
                                                         CONNEC_SIZE );

#ifdef MINT_USE_SIDRE
    const std::string topo = "t" + std::to_string( dim );
    const std::string coordset = "c" + std::to_string( dim );
    meshes[ cur_mesh++ ] = new UnstructuredMesh< Topology::MIXED >( dim, root, 
                                                                    topo,
                                                                    coordset );
#endif
  }

  /* Check that append functions properly. */
  for ( int i = 0; i < N_MESHES; ++i )
  {
    internal::append_nodes( meshes[ i ], N_NODES );
    internal::set_nodes( meshes[ i ] );
  }

  /* Check that the external meshes function properly. */
  for ( int i = 1; i < N_MESHES; i += STRIDE )
  {
    internal::deleteAndDuplicateExternalMesh( meshes[ i ] );
    internal::check_set_nodes( meshes[ i ] );
  }

  /* Check that the sidre meshes function properly. */
#ifdef MINT_USE_SIDRE
  for ( int i = 2; i < N_MESHES; i += STRIDE )
  {
    internal::deleteAndDuplicateSidreMesh( meshes[ i ] );
    internal::check_set_nodes( meshes[ i ] );
  }
#endif

  /* Delete the meshes. */
  cur_mesh = 0;
  for ( int dim = 0; dim < 3; ++dim )
  {
    delete meshes[ cur_mesh++ ];
    internal::deleteExternalMesh( meshes[ cur_mesh++ ] );

#ifdef MINT_USE_SIDRE
    delete meshes[ cur_mesh++ ];
#endif
  }
}

//------------------------------------------------------------------------------
TEST( mint_mesh_particle_mesh, setCellsSingle )
{
#ifdef MINT_USE_SIDRE
  sidre::DataStore ds;
  sidre::Group* root = ds.getRoot();
  constexpr int STRIDE = 3;
#else
  constexpr int STRIDE = 2;
#endif

  constexpr IndexType N_NODES = 12 * 10 * 10;
  constexpr IndexType N_CELLS = 900;

  constexpr int N_MESHES = STRIDE * 3;
  UnstructuredMesh< Topology::SINGLE >* meshes[N_MESHES];

  /* Create the meshes. */
  int cur_mesh = 0;
  for ( int dim = 1; dim <= 3; ++dim )
  {
    meshes[ cur_mesh++ ] = new UnstructuredMesh< Topology::SINGLE >( dim, QUAD );
    meshes[ cur_mesh++ ] = internal::createExternalSingle( dim, QUAD, N_CELLS, 
                                                           N_NODES );

#ifdef MINT_USE_SIDRE
    const std::string topo = "t" + std::to_string( dim );
    const std::string coordset = "c" + std::to_string( dim );
    meshes[ cur_mesh++ ] = new UnstructuredMesh< Topology::SINGLE >( dim, QUAD, 
                                                                     root, topo,
                                                                     coordset );
#endif
  }

  /* Check that append functions properly. */
  for ( int i = 0; i < N_MESHES; ++i )
  {
    internal::append_nodes( meshes[ i ], N_NODES );
    internal::append_cells( meshes[ i ], N_CELLS );
    internal::set_cells( meshes[ i ] );
  }

  /* Check that the external meshes function properly. */
  for ( int i = 1; i < N_MESHES; i += STRIDE )
  {
    internal::deleteAndDuplicateExternalMesh( meshes[ i ] );
    internal::check_append_nodes( meshes[ i ] );
    internal::check_set_cells( meshes[ i ] );
  }

  /* Check that the sidre meshes function properly. */
#ifdef MINT_USE_SIDRE
  for ( int i = 2; i < N_MESHES; i += STRIDE )
  {
    internal::deleteAndDuplicateSidreMesh( meshes[ i ] );
    internal::check_append_nodes( meshes[ i ] );
    internal::check_set_cells( meshes[ i ] );
  }
#endif

  /* Delete the meshes. */
  cur_mesh = 0;
  for ( int dim = 0; dim < 3; ++dim )
  {
    delete meshes[ cur_mesh++ ];
    internal::deleteExternalMesh( meshes[ cur_mesh++ ] );

#ifdef MINT_USE_SIDRE
    delete meshes[ cur_mesh++ ];
#endif
  }
}

//------------------------------------------------------------------------------
TEST( mint_mesh_particle_mesh, setCellsMixed )
{
#ifdef MINT_USE_SIDRE
  sidre::DataStore ds;
  sidre::Group* root = ds.getRoot();
  constexpr int STRIDE = 3;
#else
  constexpr int STRIDE = 2;
#endif

  constexpr IndexType N_NODES = 1200;
  constexpr IndexType N_CELLS = 900;
  constexpr IndexType CONNEC_SIZE = internal::getConnectivitySize( N_CELLS );

  constexpr int N_MESHES = STRIDE * 3;
  UnstructuredMesh< Topology::MIXED >* meshes[N_MESHES];

  /* Create the meshes. */
  int cur_mesh = 0;
  for ( int dim = 1; dim <= 3; ++dim )
  {
    meshes[ cur_mesh++ ] = new UnstructuredMesh< Topology::MIXED >( dim );
    meshes[ cur_mesh++ ] = internal::createExternalMixed( dim, N_CELLS, N_NODES, 
                                                         CONNEC_SIZE );

#ifdef MINT_USE_SIDRE
    const std::string topo = "t" + std::to_string( dim );
    const std::string coordset = "c" + std::to_string( dim );
    meshes[ cur_mesh++ ] = new UnstructuredMesh< Topology::MIXED >( dim, root, 
                                                                    topo,
                                                                    coordset );
#endif
  }

  /* Check that append functions properly. */
  for ( int i = 0; i < N_MESHES; ++i )
  {
    internal::append_nodes( meshes[ i ], N_NODES );
    internal::append_cells( meshes[ i ], N_CELLS );
    internal::set_cells( meshes[ i ] );
  }

  /* Check that the external meshes function properly. */
  for ( int i = 1; i < N_MESHES; i += STRIDE )
  {
    internal::deleteAndDuplicateExternalMesh( meshes[ i ] );
    internal::check_append_nodes( meshes[ i ] );
    internal::check_set_cells( meshes[ i ] );
  }

  /* Check that the sidre meshes function properly. */
#ifdef MINT_USE_SIDRE
  for ( int i = 2; i < N_MESHES; i += STRIDE )
  {
    internal::deleteAndDuplicateSidreMesh( meshes[ i ] );
    internal::check_append_nodes( meshes[ i ] );
    internal::check_set_cells( meshes[ i ] );
  }
#endif

  /* Delete the meshes. */
  cur_mesh = 0;
  for ( int dim = 0; dim < 3; ++dim )
  {
    delete meshes[ cur_mesh++ ];
    internal::deleteExternalMesh( meshes[ cur_mesh++ ] );

#ifdef MINT_USE_SIDRE
    delete meshes[ cur_mesh++ ];
#endif
  }
}

} /* end namespace mint */
} /* end namespace axom */

//------------------------------------------------------------------------------
#include "slic/UnitTestLogger.hpp"
using axom::slic::UnitTestLogger;

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  UnitTestLogger logger;  // create & initialize test logger,

  // finalized when exiting main scope

  result = RUN_ALL_TESTS();

  return result;
}

