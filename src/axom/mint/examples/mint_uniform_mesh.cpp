// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file
 *
 * \brief Illustrates how to construct and use a mint::UniformMesh object.
 */

// sphinx_tutorial_basic_example_start

// sphinx_tutorial_walkthrough_includes_start

#include "axom/mint.hpp"                  // for mint classes and functions
#include "axom/core/numerics/Matrix.hpp"  // for numerics::Matrix

// sphinx_tutorial_walkthrough_includes_end

// namespace aliases
namespace mint     = axom::mint;
namespace numerics = axom::numerics;
namespace xargs    = mint::xargs;

using IndexType    = axom::IndexType;

// compile-time switch for execution policy
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_CUDA)
  constexpr int NUM_BLOCKS = 512;
  using ExecPolicy = mint::policy::parallel_gpu< NUM_BLOCKS >;
#elif defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP)
  using ExecPolicy = mint::policy::parallel_cpu;
#else
  using ExecPolicy = mint::policy::serial;
#endif

constexpr IndexType NUM_COMPONENTS     = 2;
constexpr IndexType NUM_NODES_PER_CELL = 4;
constexpr double ONE_OVER_4 = 1. / static_cast< double >( NUM_NODES_PER_CELL );

//------------------------------------------------------------------------------
int main ( int AXOM_NOT_USED(argc), char** AXOM_NOT_USED(argv) )
{

// sphinx_tutorial_walkthrough_construct_mesh_start

  // construct a 100 x 100 grid within a domain defined in [-5.0, 5.0]
  const double lo[]   = { -5.0, -5.0 };
  const double hi[]   = {  5.0,  5.0 };
  mint::UniformMesh mesh( lo, hi, 100, 100 );

// sphinx_tutorial_walkthrough_construct_mesh_end

// sphinx_tutorial_walkthrough_add_fields_start

  // add a cell-centered and a node-centered field
  double* phi = mesh.createField< double >( "phi", mint::NODE_CENTERED );
  double* xc  = mesh.createField< double >( "xc", mint::CELL_CENTERED, 2 );

// sphinx_tutorial_walkthrough_add_fields_end

// sphinx_tutorial_walkthrough_compute_hf_start

  // loop over the nodes and evaluate Himmelblaus Function
  mint::for_all_nodes< ExecPolicy, xargs::xy >(
      &mesh, AXOM_LAMBDA( IndexType nodeIdx, double x, double y )
  {
    const double x_2 = x * x;
    const double y_2 = y * y;
    const double A   = x_2 + y   - 11.0;
    const double B   = x   + y_2 - 7.0;

    phi[ nodeIdx ] = A * A + B * B;
  } );

// sphinx_tutorial_walkthrough_compute_hf_end

// sphinx_tutorial_walkthrough_cell_centers_start

  // loop over cells and compute cell centers
  mint::for_all_cells< ExecPolicy, xargs::coords >(
      &mesh, AXOM_LAMBDA( IndexType cellIdx,
                          const numerics::Matrix< double >& coords,
                          const IndexType* AXOM_NOT_USED(nodeIds) )
  {
    // NOTE: A column vector of the coords matrix corresponds to a nodes coords

    // Sum the cell's nodal coordinates
    double xsum = 0.0;
    double ysum = 0.0;
    for ( IndexType inode=0; inode < NUM_NODES_PER_CELL; ++inode )
    {
      const double* node = coords.getColumn( inode );
      xsum += node[ mint::X_COORDINATE ];
      ysum += node[ mint::Y_COORDINATE ];
    } // END for all cell nodes

    // compute cell centroid by averaging the nodal coordinate sums
    const IndexType offset = cellIdx * NUM_COMPONENTS;
    xc[ offset   ] = xsum * ONE_OVER_4;
    xc[ offset+1 ] = ysum * ONE_OVER_4;

  } );

// sphinx_tutorial_walkthrough_cell_centers_end

// sphinx_tutorial_walkthrough_vtk_start

  // write the mesh in a VTK file for visualization
  mint::write_vtk( &mesh, "uniform_mesh.vtk" );

// sphinx_tutorial_walkthrough_vtk_end

  return 0;
}

// sphinx_tutorial_basic_example_end


