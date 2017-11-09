/*
 * Copyright (c) 2017, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */

#include "gtest/gtest.h"

#include "mint/CellConnectivity.hpp"
#include "mint/CellType.hpp"
#include "mint/DataTypes.hpp"
#include "slic/slic.hpp"
#include "slic/UnitTestLogger.hpp"      /* for UnitTestLogger */


#include <cmath>                        /* for std::ceil */


namespace axom {
namespace mint {

namespace internal {


template < int cell_type >
localIndex calc_new_capacity( CellConnectivity< cell_type > & c, 
                              localIndex increase )
{ 
  int tuple_size = 1;
  if ( cell_type != MINT_MIXED_CELL ) {
    tuple_size = cell::num_nodes[ cell_type ];
  }

  localIndex newSize = c.getSize() + increase;
  if ( newSize > c.getCapacity() ) {
    double n_tuples = newSize * c.getResizeRatio() / tuple_size;
    return std::ceil( n_tuples ) * tuple_size;
  }

  return c.getCapacity();
}


template < int cell_type >
void test_single_element_type( localIndex capacity, double resize_ratio )
{
  CellConnectivity< cell_type > connec( capacity, resize_ratio );
  localIndex size = 0;
  localIndex num_cells = 0;
  localIndex nodes_per_cell = cell::num_nodes[ cell_type ];
  capacity *= nodes_per_cell;

  EXPECT_TRUE( connec.empty() );
  EXPECT_FALSE( connec.hasMixedCellTypes() );
  EXPECT_EQ( connec.getSize(), size );
  EXPECT_EQ( connec.getCapacity(), capacity );
  EXPECT_EQ( connec.getResizeRatio(), resize_ratio );
  EXPECT_EQ( connec.getNumberOfCells(), num_cells );
  EXPECT_EQ( connec.getNumberOfNodes( 0 ), nodes_per_cell );
  EXPECT_EQ( connec.getCellType( 0 ), cell_type );

  localIndex num_values = 1000 * nodes_per_cell;
  localIndex values[ num_values ];
  for ( localIndex i = 0; i < num_values; ++i ) {
    values[ i ] = i * ( i * (i - 3) + 5 ) - 77;
  }

  for ( localIndex i = 0; i < num_values / nodes_per_cell; ++i ) {
    size += nodes_per_cell;
    num_cells++;
    capacity = calc_new_capacity( connec, nodes_per_cell );
    connec.addCell( values + i * nodes_per_cell, i ); 
  }

  EXPECT_FALSE( connec.empty() );
  EXPECT_FALSE( connec.hasMixedCellTypes() );
  EXPECT_EQ( connec.getSize(), size );
  EXPECT_EQ( connec.getCapacity(), capacity );
  EXPECT_EQ( connec.getResizeRatio(), resize_ratio );
  EXPECT_EQ( connec.getNumberOfCells(), num_cells );

  for ( localIndex i = 0; i < num_cells; ++i ) {
    EXPECT_EQ( connec.getNumberOfNodes( i ), nodes_per_cell );
    EXPECT_EQ( connec.getCellType( i ), cell_type );

    const localIndex* cell = connec[ i ];
    for ( localIndex j = 0; j < nodes_per_cell; ++j ) {
      EXPECT_EQ( cell[ j ], values[ i * nodes_per_cell + j ] );
    }
  }

  for ( localIndex i = 0; i < num_values; ++i ) {
    values[ i ] = i * ( 10536 - 11 * i );
  }

  for ( localIndex i = 0; i < num_cells; ++i ) {
    connec.setCell( i, values + i * nodes_per_cell );
  }

  EXPECT_FALSE( connec.empty() );
  EXPECT_FALSE( connec.hasMixedCellTypes() );
  EXPECT_EQ( connec.getSize(), size );
  EXPECT_EQ( connec.getCapacity(), capacity );
  EXPECT_EQ( connec.getResizeRatio(), resize_ratio );
  EXPECT_EQ( connec.getNumberOfCells(), num_cells );

  for ( localIndex i = 0; i < num_cells; ++i ) {
    EXPECT_EQ( connec.getNumberOfNodes( i ), nodes_per_cell );
    EXPECT_EQ( connec.getCellType( i ), cell_type );

    const localIndex* cell = connec[ i ];
    for ( localIndex j = 0; j < nodes_per_cell; ++j ) {
      EXPECT_EQ( cell[ j ], values[ i * nodes_per_cell + j ] );
    }
  }
}


template < int cell_type >
localIndex * create_connectivity( localIndex num_cells )
{
  localIndex nodes_per_cell = cell::num_nodes[ cell_type ];
  localIndex num_nodes = nodes_per_cell * num_cells;
  localIndex * connectivity_array = new localIndex[ num_nodes ];
  
  for ( localIndex i = 0; i < num_nodes; ++i ) {
    connectivity_array[ i ] = nodes_per_cell + (nodes_per_cell - 7) * i;
    connectivity_array[ i ] += nodes_per_cell * nodes_per_cell * i * i;
    connectivity_array[ i ] -= ( nodes_per_cell - 4 ) * i * i * i / 100.0;
  }

  return connectivity_array;
}


template < int cell_type >
void check_connectivity( CellConnectivity< MINT_MIXED_CELL > & connec,
                         localIndex num_cells, localIndex * connectivity_array, 
                         localIndex offset )
{
  localIndex nodes_per_cell = cell::num_nodes[ cell_type ];

  for ( localIndex i = 0; i < num_cells; ++i ) {
    EXPECT_EQ( connec.getNumberOfNodes( i + offset), nodes_per_cell );
    EXPECT_EQ( connec.getCellType( i + offset), cell_type );

    const localIndex* cell = connec[ i + offset ];
    for ( localIndex j = 0; j < nodes_per_cell; ++j ) {
      EXPECT_EQ( cell[ j ], connectivity_array[ i * nodes_per_cell + j ] );
    }
  }
}


template < int cell_type >
void add_connectivity( CellConnectivity< MINT_MIXED_CELL > & connec,
                       localIndex num_cells, localIndex * connectivity_array )
{
  localIndex nodes_per_cell = cell::num_nodes[ cell_type ];
  localIndex prev_size = connec.getSize();
  localIndex prev_total_cells = connec.getNumberOfCells();

  localIndex capacity;
  for ( localIndex i = 0; i < num_cells; ++i ) {
    capacity = calc_new_capacity( connec, nodes_per_cell );
    connec.addCell( connectivity_array + i * nodes_per_cell, cell_type );
  }

  EXPECT_FALSE( connec.empty() );
  EXPECT_TRUE( connec.hasMixedCellTypes() );
  EXPECT_EQ( connec.getSize(), prev_size + num_cells * nodes_per_cell );
  EXPECT_EQ( connec.getCapacity(), capacity );
  EXPECT_EQ( connec.getNumberOfCells(), prev_total_cells + num_cells );

  check_connectivity< cell_type >( connec, num_cells, connectivity_array, prev_total_cells );
}

template < int cell_type >
void update_connectivity( CellConnectivity< MINT_MIXED_CELL > & connec,
                          localIndex num_cells, localIndex * connectivity_array,
                          localIndex offset )
{
  localIndex nodes_per_cell = cell::num_nodes[ cell_type ];
  localIndex prev_size = connec.getSize();
  localIndex prev_capacity = connec.getCapacity();
  localIndex prev_total_cells = connec.getNumberOfCells();

  for ( localIndex i = 0; i < num_cells; ++i ) {
    for ( localIndex j = 0; j < nodes_per_cell; ++j ) {
      connectivity_array[ i * nodes_per_cell + j ] *= 10.0 / ( i + 3 );
    }

    connec.setCell( offset + i, connectivity_array + i * nodes_per_cell );
  }

  EXPECT_FALSE( connec.empty() );
  EXPECT_TRUE( connec.hasMixedCellTypes() );
  EXPECT_EQ( connec.getSize(), prev_size );
  EXPECT_EQ( connec.getCapacity(), prev_capacity );
  EXPECT_EQ( connec.getNumberOfCells(), prev_total_cells );

  check_connectivity< cell_type >( connec, num_cells, connectivity_array, offset );
}




void test_mixed_element_type( localIndex capacity, double resize_ratio )
{
  CellConnectivity< MINT_MIXED_CELL > connec( capacity, resize_ratio );
  localIndex size = 0;
  localIndex num_cells = 0;
  localIndex nodes_per_cell = cell::num_nodes[ MINT_MIXED_CELL ];
  capacity *= nodes_per_cell;

  EXPECT_TRUE( connec.empty() );
  EXPECT_TRUE( connec.hasMixedCellTypes() );
  EXPECT_EQ( connec.getSize(), size );
  EXPECT_EQ( connec.getCapacity(), capacity );
  EXPECT_EQ( connec.getResizeRatio(), resize_ratio );
  EXPECT_EQ( connec.getNumberOfCells(), num_cells );

  /* Create and add the connectivity arrays for each element type. */
  localIndex num_vertex_cells = 400;
  localIndex * vertex_connectivity = 
                      create_connectivity< MINT_VERTEX >( num_vertex_cells );
  add_connectivity< MINT_VERTEX >( connec, num_vertex_cells, 
                                   vertex_connectivity );

  localIndex num_segment_cells = 410;
  localIndex * segment_connectivity = 
                      create_connectivity< MINT_SEGMENT >( num_segment_cells );
  add_connectivity< MINT_SEGMENT >( connec, num_segment_cells, 
                                    segment_connectivity );

  localIndex num_triangle_cells = 420;
  localIndex * triangle_connectivity = 
                    create_connectivity< MINT_TRIANGLE >( num_triangle_cells );
  add_connectivity< MINT_TRIANGLE >( connec, num_triangle_cells, 
                                     triangle_connectivity );

  localIndex num_quad_cells = 430;
  localIndex * quad_connectivity = 
                            create_connectivity< MINT_QUAD >( num_quad_cells );
  add_connectivity< MINT_QUAD >( connec, num_quad_cells, quad_connectivity );

  localIndex num_tet_cells = 440;
  localIndex * tet_connectivity = 
                              create_connectivity< MINT_TET >( num_tet_cells );
  add_connectivity< MINT_TET >( connec, num_tet_cells, tet_connectivity );

  localIndex num_hex_cells = 450;
  localIndex * hex_connectivity = 
                              create_connectivity< MINT_HEX >( num_hex_cells );
  add_connectivity< MINT_HEX >( connec, num_hex_cells, hex_connectivity );

  localIndex num_prism_cells = 460;
  localIndex * prism_connectivity = 
                          create_connectivity< MINT_PRISM >( num_prism_cells );
  add_connectivity< MINT_PRISM >( connec, num_prism_cells, prism_connectivity );

  localIndex num_pyramid_cells = 470;
  localIndex * pyramid_connectivity = 
                      create_connectivity< MINT_PYRAMID >( num_pyramid_cells );
  add_connectivity< MINT_PYRAMID >( connec, num_pyramid_cells, 
                                    pyramid_connectivity );

  localIndex num_quad9_cells = 480;
  localIndex * quad9_connectivity = 
                          create_connectivity< MINT_QUAD9 >( num_quad9_cells );
  add_connectivity< MINT_QUAD9 >( connec, num_quad9_cells, quad9_connectivity );

  localIndex num_hex27_cells = 490;
  localIndex * hex27_connectivity = 
                          create_connectivity< MINT_HEX27 >( num_hex27_cells );
  add_connectivity< MINT_HEX27 >( connec, num_hex27_cells, hex27_connectivity );

  /* Check that the connectivity of each element type was added correctly */
  localIndex cell_offset = 0;
  check_connectivity< MINT_VERTEX >( connec, num_vertex_cells, 
                                     vertex_connectivity, cell_offset );

  cell_offset += num_vertex_cells;
  check_connectivity< MINT_SEGMENT >( connec, num_segment_cells, 
                                      segment_connectivity, cell_offset );

  cell_offset += num_segment_cells;
  check_connectivity< MINT_TRIANGLE >( connec, num_triangle_cells, 
                                     triangle_connectivity, cell_offset );

  cell_offset += num_triangle_cells;
  check_connectivity< MINT_QUAD >( connec, num_quad_cells, 
                                   quad_connectivity, cell_offset );

  cell_offset += num_quad_cells;
  check_connectivity< MINT_TET >( connec, num_tet_cells, 
                                  tet_connectivity, cell_offset );

  cell_offset += num_tet_cells;
  check_connectivity< MINT_HEX >( connec, num_hex_cells, 
                                  hex_connectivity, cell_offset );

  cell_offset += num_hex_cells;
  check_connectivity< MINT_PRISM >( connec, num_prism_cells, 
                                    prism_connectivity, cell_offset );

  cell_offset += num_prism_cells;
  check_connectivity< MINT_PYRAMID >( connec, num_pyramid_cells, 
                                      pyramid_connectivity, cell_offset );

  cell_offset += num_pyramid_cells;
  check_connectivity< MINT_QUAD9 >( connec, num_quad9_cells, 
                                    quad9_connectivity, cell_offset );

  cell_offset += num_quad9_cells;
  check_connectivity< MINT_HEX27 >( connec, num_hex27_cells, 
                                    hex27_connectivity, cell_offset );

  /* Update the connectivity of each element type */
  cell_offset = 0;
  update_connectivity< MINT_VERTEX >( connec, num_vertex_cells, 
                                      vertex_connectivity, cell_offset );

  cell_offset += num_vertex_cells;
  update_connectivity< MINT_SEGMENT >( connec, num_segment_cells, 
                                       segment_connectivity, cell_offset );

  cell_offset += num_segment_cells;
  update_connectivity< MINT_TRIANGLE >( connec, num_triangle_cells, 
                                        triangle_connectivity, cell_offset );

  cell_offset += num_triangle_cells;
  update_connectivity< MINT_QUAD >( connec, num_quad_cells, 
                                    quad_connectivity, cell_offset );

  cell_offset += num_quad_cells;
  update_connectivity< MINT_TET >( connec, num_tet_cells, 
                                   tet_connectivity, cell_offset );

  cell_offset += num_tet_cells;
  update_connectivity< MINT_HEX >( connec, num_hex_cells, 
                                   hex_connectivity, cell_offset );

  cell_offset += num_hex_cells;
  update_connectivity< MINT_PRISM >( connec, num_prism_cells, 
                                     prism_connectivity, cell_offset );

  cell_offset += num_prism_cells;
  update_connectivity< MINT_PYRAMID >( connec, num_pyramid_cells, 
                                       pyramid_connectivity, cell_offset );

  cell_offset += num_pyramid_cells;
  update_connectivity< MINT_QUAD9 >( connec, num_quad9_cells, 
                                     quad9_connectivity, cell_offset );

  cell_offset += num_quad9_cells;
  update_connectivity< MINT_HEX27 >( connec, num_hex27_cells, 
                                     hex27_connectivity, cell_offset );

  /* Check that the connectivity of each element type was updated correctly */
  cell_offset = 0;
  check_connectivity< MINT_VERTEX >( connec, num_vertex_cells, 
                                     vertex_connectivity, cell_offset );

  cell_offset += num_vertex_cells;
  check_connectivity< MINT_SEGMENT >( connec, num_segment_cells, 
                                      segment_connectivity, cell_offset );

  cell_offset += num_segment_cells;
  check_connectivity< MINT_TRIANGLE >( connec, num_triangle_cells, 
                                     triangle_connectivity, cell_offset );

  cell_offset += num_triangle_cells;
  check_connectivity< MINT_QUAD >( connec, num_quad_cells, 
                                   quad_connectivity, cell_offset );

  cell_offset += num_quad_cells;
  check_connectivity< MINT_TET >( connec, num_tet_cells, 
                                  tet_connectivity, cell_offset );

  cell_offset += num_tet_cells;
  check_connectivity< MINT_HEX >( connec, num_hex_cells, 
                                  hex_connectivity, cell_offset );

  cell_offset += num_hex_cells;
  check_connectivity< MINT_PRISM >( connec, num_prism_cells, 
                                    prism_connectivity, cell_offset );

  cell_offset += num_prism_cells;
  check_connectivity< MINT_PYRAMID >( connec, num_pyramid_cells, 
                                      pyramid_connectivity, cell_offset );

  cell_offset += num_pyramid_cells;
  check_connectivity< MINT_QUAD9 >( connec, num_quad9_cells, 
                                    quad9_connectivity, cell_offset );

  cell_offset += num_quad9_cells;
  check_connectivity< MINT_HEX27 >( connec, num_hex27_cells, 
                                    hex27_connectivity, cell_offset );

  /* Check that the number of nodes and cells is correct. */
  num_cells = 0;
  size = 0;

  num_cells += num_vertex_cells;
  size += num_vertex_cells * cell::num_nodes[ MINT_VERTEX ];

  num_cells += num_segment_cells;
  size += num_segment_cells * cell::num_nodes[ MINT_SEGMENT ];

  num_cells += num_triangle_cells;
  size += num_triangle_cells * cell::num_nodes[ MINT_TRIANGLE ];

  num_cells += num_quad_cells;
  size += num_quad_cells * cell::num_nodes[ MINT_QUAD ];

  num_cells += num_tet_cells;
  size += num_tet_cells * cell::num_nodes[ MINT_TET ];

  num_cells += num_hex_cells;
  size += num_hex_cells * cell::num_nodes[ MINT_HEX ];

  num_cells += num_prism_cells;
  size += num_prism_cells * cell::num_nodes[ MINT_PRISM ];

  num_cells += num_pyramid_cells;
  size += num_pyramid_cells * cell::num_nodes[ MINT_PYRAMID ];

  num_cells += num_quad9_cells;
  size += num_quad9_cells * cell::num_nodes[ MINT_QUAD9 ];

  num_cells += num_hex27_cells;
  size += num_hex27_cells * cell::num_nodes[ MINT_HEX27 ];

  EXPECT_FALSE( connec.empty() );
  EXPECT_TRUE( connec.hasMixedCellTypes() );
  EXPECT_EQ( connec.getSize(), size );
  EXPECT_EQ( connec.getResizeRatio(), resize_ratio );
  EXPECT_EQ( connec.getNumberOfCells(), num_cells );
}

}   /* end namespace internal */


TEST( mint_connectivity, SingleElementType )
{
  for ( localIndex capacity = 0; capacity <= 1000; capacity += 2 * capacity + 1) {
    for ( double resize_ratio = 1.0; resize_ratio <= 3.0; resize_ratio += 0.25 ) {
      internal::test_single_element_type< MINT_VERTEX >( capacity, resize_ratio );
      internal::test_single_element_type< MINT_SEGMENT >( capacity, resize_ratio );
      internal::test_single_element_type< MINT_TRIANGLE >( capacity, resize_ratio );
      internal::test_single_element_type< MINT_QUAD >( capacity, resize_ratio );
      internal::test_single_element_type< MINT_TET >( capacity, resize_ratio );
      internal::test_single_element_type< MINT_HEX >( capacity, resize_ratio );
      internal::test_single_element_type< MINT_PRISM >( capacity, resize_ratio );
      internal::test_single_element_type< MINT_QUAD9 >( capacity, resize_ratio );
      internal::test_single_element_type< MINT_HEX27 >( capacity, resize_ratio );
    }
  }
}


TEST( mint_connectivity, MixedElementType )
{
  for ( localIndex capacity = 0; capacity <= 1000; capacity += 2 * capacity + 1) {
    for ( double resize_ratio = 1.0; resize_ratio <= 3.0; resize_ratio += 0.25 ) {
      internal::test_mixed_element_type( capacity, resize_ratio );
    }
  }
}



}   /* end namespace mint */
}   /* end namespace axom */

//------------------------------------------------------------------------------
using axom::slic::UnitTestLogger;

int main( int argc, char * argv[] )
{
  int result = 0;
  ::testing::InitGoogleTest( &argc, argv );
  UnitTestLogger logger;
  result = RUN_ALL_TESTS();
  return result;
}