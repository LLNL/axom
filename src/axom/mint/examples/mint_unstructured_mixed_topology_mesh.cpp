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

/*!
 * \file
 *
 * \brief Illustrates how to construct and use a mixed topology
 *  UnstructuredMesh by creating a 2D structured grid with both quads
 *  and triangles.
 */

// Mint includes
#include "mint/config.hpp"
#include "mint/UnstructuredMesh.hpp"
#include "mint/vtk_utils.hpp"
#include "slic/UnitTestLogger.hpp"      /* for UnitTestLogger */

// C/C++ includes
#include <random>                       /* for random number generator */

using namespace axom;
using axom::slic::UnitTestLogger;


inline bool appendQuad( mint::IndexType i, mint::IndexType j )
{
  return (i % 2 == 0 && j % 2 == 0) || ( i % 2 == 1 && j % 2 == 1);
}


//------------------------------------------------------------------------------
int main( int AXOM_NOT_USED(argc), char** AXOM_NOT_USED(argv) )
{
  UnitTestLogger logger;  // create & initialize test logger,

  constexpr int DIMENSION = 2;
  constexpr mint::IndexType X_EXTENT = 11;
  constexpr mint::IndexType Y_EXTENT = 11;
  constexpr double SPACING = 1.0;
  constexpr mint::IndexType NUM_NODES = X_EXTENT * Y_EXTENT;
  constexpr mint::IndexType NUM_CELLS = (X_EXTENT - 1) * (Y_EXTENT - 1) * 3 / 2;

  constexpr double HI  = -10.0;
  constexpr double LO  = 10.0;
  constexpr double PLO = 0.0;
  constexpr double PHI = 1.0;

  /* STEP 0: create the UnstructuredMesh */
  mint::UnstructuredMesh< mint::MIXED_SHAPE > mesh( DIMENSION, NUM_NODES,
                                                    NUM_CELLS );

  /* STEP 1: Add fields to the nodes and cells.
   * Note that we can only use the pointers below because we specified
   * the node and cell capacities in the constructor. In general calling
   * append or insert can invalidate any associated pointers. */
  double* vx = mesh.createField< double >( "vx", mint::NODE_CENTERED );
  double* vy = mesh.createField< double >( "vy", mint::NODE_CENTERED );
  double* vz = mesh.createField< double >( "vz", mint::NODE_CENTERED );
  double* p = mesh.createField< double > ( "pressure", mint::CELL_CENTERED );
  double* p_avg = mesh.createField< double >( "average_pressure",
                                              mint::NODE_CENTERED );
  int* cells_per_node = mesh.createField< int >( "cells_per_node",
                                                 mint::NODE_CENTERED );

  /* STEP 3: Add the nodes */
  mint::IndexType node_ID = 0;
  for ( mint::IndexType j = 0 ; j < Y_EXTENT ; ++j )
  {
    for ( mint::IndexType i = 0 ; i < X_EXTENT ; ++i )
    {
      mesh.appendNode( i * SPACING, j * SPACING );

      vx[ node_ID ] = utilities::random_real( LO, HI );
      vy[ node_ID ] = utilities::random_real( LO, HI );
      vz[ node_ID ] = utilities::random_real( LO, HI );
      p_avg[ node_ID ] = 0.0;
      cells_per_node[ node_ID ] = 0.0;
      node_ID++;
    }
  }

  /* STEP 4: Add the cells */
  mint::IndexType cell_ID = 0;
  for ( mint::IndexType j = 0 ; j < Y_EXTENT - 1 ; ++j )
  {
    for ( mint::IndexType i = 0 ; i < X_EXTENT - 1 ; ++i )
    {
      const mint::IndexType bottom_left = j * X_EXTENT + i;
      const mint::IndexType bottom_right = bottom_left + 1;
      const mint::IndexType top_right = bottom_right + X_EXTENT;
      const mint::IndexType top_left = bottom_left + X_EXTENT;

      if ( appendQuad( i, j ) )
      {
        /* Append a quad. */
        const mint::IndexType quad[4] =
        { bottom_left, bottom_right, top_right, top_left };
        mesh.appendCell( quad, mint::QUAD);

        p[ cell_ID ] = utilities::random_real( PLO, PHI );
        cell_ID++;
      }
      else
      {
        /* Append two triangles. */
        const mint::IndexType tri0[3] =
        { bottom_left, bottom_right, top_right };
        mesh.appendCell( tri0, mint::TRIANGLE );
        p[ cell_ID ] = utilities::random_real( PLO, PHI );
        cell_ID++;

        const mint::IndexType tri1[3] =
        { top_right, top_left, bottom_left };
        mesh.appendCell( tri1, mint::TRIANGLE );
        p[ cell_ID ] = utilities::random_real( PLO, PHI );
        cell_ID++;
      }
    }
  }

  /* STEP 5: Calculate the average pressure at each of the nodes. */
  const mint::IndexType n_cells = mesh.getNumberOfCells();
  for ( cell_ID = 0 ; cell_ID < n_cells ; ++cell_ID )
  {
    const double cell_pressure = p[ cell_ID ];
    const int n_cell_nodes = mesh.getNumberOfCellNodes( cell_ID );
    const mint::IndexType* connec = mesh.getCell( cell_ID );
    for ( int i = 0 ; i < n_cell_nodes ; ++i )
    {
      node_ID = connec[ i ];
      p_avg[ node_ID ] += cell_pressure;
      cells_per_node[ node_ID ]++;
    }
  }

  const mint::IndexType n_nodes = mesh.getNumberOfNodes();
  for ( node_ID = 0 ; node_ID < n_nodes ; ++node_ID )
  {
    p_avg[ node_ID ] /= cells_per_node[ node_ID ];
  }

  /* STEP 6: write the particle mesh in VTK format for visualization. */
  mint::write_vtk( &mesh, "UnstructuredMesh_mixed_topology.vtk" );

  return 0;
}
