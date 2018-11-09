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

#ifndef MINT_EXECUTION_HELPERS_HPP_
#define MINT_EXECUTION_HELPERS_HPP_

// mint includes
#include "axom/mint/config.hpp"                 // for compile-time definitions
#include "axom/mint/mesh/Mesh.hpp"              // for Mesh

#include "axom/core/numerics/Matrix.hpp"        // for Matrix

#ifdef AXOM_USE_RAJA
#include "RAJA/RAJA.hpp"
#endif

namespace axom
{
namespace mint
{
namespace internal
{

//------------------------------------------------------------------------------
template < typename ExecPolicy, int NDIM, int NNODES, typename FOR_ALL,
           typename KernelType >
inline void for_all_coords( const FOR_ALL & for_all_nodes, const mint::Mesh* m,
                            KernelType&& kernel )
{
  AXOM_STATIC_ASSERT_MSG( NDIM >= 1 && NDIM <= 3, "NDIM must be a valid dimension." );
  AXOM_STATIC_ASSERT_MSG( NNODES > 0, "NNODES must be greater than zero." );
  SLIC_ASSERT( m != nullptr );
  SLIC_ASSERT( m->getDimension() == NDIM );

  const double * const coords[3] = { 
                               m->getCoordinateArray(X_COORDINATE),
                  (NDIM > 1) ? m->getCoordinateArray(Y_COORDINATE) : nullptr,
                  (NDIM > 2) ? m->getCoordinateArray(Z_COORDINATE) : nullptr };

  for_all_nodes( ExecPolicy(), m, 
    AXOM_LAMBDA( IndexType cellID, const IndexType * nodeIDs, 
                 IndexType numNodes )
    {
      AXOM_DEBUG_VAR(numNodes);
      SLIC_ASSERT( numNodes == NNODES );

      double localCoords[ NDIM * NNODES ];
      for ( int i = 0; i < NNODES; ++i )
      {
        for ( int dim = 0; dim < NDIM; ++dim )
        {
          localCoords[ NDIM * i + dim ] = coords[ dim ][ nodeIDs[ i ] ];
        }
      }

      numerics::Matrix<double> coordsMatrix( NDIM, NNODES, localCoords, true );
      kernel( cellID, coordsMatrix, nodeIDs );
    }
  );
}

} /* namespace internal */
} /* namespace mint     */
} /* namespace axom     */

#endif /* MINT_EXECUTION_HELPERS_HPP_ */