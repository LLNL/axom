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
 * \brief Illustrates how to construct and use a mint::UniformMesh object.
 */

// Mint includes
#include "mint/config.hpp"
#include "mint/UniformMesh.hpp"
#include "mint/vtk_utils.hpp"

// C/C++ includes
#include <cmath>

// namespace aliases
namespace mint  = axom::mint;
using IndexType = mint::IndexType;

//------------------------------------------------------------------------------
double himmelblaus_function( double x, double y )
{
  const double x_2 = x*x;
  const double y_2 = y*y;
  const double A   = x_2 + y - 11.0;
  const double B   = x   + y_2 - 7.0;
  const double f   = A*A + B*B;
  return f;
}

//------------------------------------------------------------------------------
int main ( int AXOM_NOT_USED(argc), char** AXOM_NOT_USED(argv) )
{
  const int dimension = 2;
  const double lo[]   = { -5.0, -5.0 };
  const double hi[]   = {  5.0,  5.0 };

  // STEP 0: construct a 100 x 100 grid
  mint::UniformMesh mesh( dimension, lo, hi, 100, 100 );

  // STEP 1: add a cell-centered and a node-centered field
  double* phi = mesh.createField< double >( "phi", mint::NODE_CENTERED );

  // STEP 2: loop over the nodes
  const IndexType Ni = mesh.getNumberOfNodesAlongDim( mint::I_DIRECTION );
  const IndexType Nj = mesh.getNumberOfNodesAlongDim( mint::J_DIRECTION );
  const IndexType jp = mesh.jp();

  for ( IndexType j=0 ; j < Nj ; ++j )
  {
    const IndexType offset = j*jp;
    for ( IndexType i=0 ; i < Ni ; ++i )
    {
      const IndexType idx = i + offset;
      const double x = mesh.evaluateCoordinate( i, mint::I_DIRECTION );
      const double y = mesh.evaluateCoordinate( j, mint::J_DIRECTION );

      phi[ idx ] = himmelblaus_function( x, y );
    } // END for all i
  } // END for all j

  // STEP 3: write the mesh in a VTK file for visualization
  mint::write_vtk( &mesh, "uniform_mesh.vtk" );

  return 0;
}
