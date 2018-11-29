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

namespace axom
{
namespace mint
{
namespace internal
{

/*!
 * \brief Iterate over the objects (cells or faces) in a mesh and for each
 *  object construct a NDIM x NNODES matrix of the nodal coordinates of the
 *  object and call a kernel with the object ID, the coordinate matrix, and the
 *  node IDs. 
 *
 * \tparam ExecPolicy the execution policy
 * \tparam NDIM the number of coordinate dimensions.
 * \tparam NNODES the number of nodes per object.
 * \tparam FOR_ALL_FUNCTOR the type of the for_all_nodes functor.
 * \tparam KernelType the type of the supplied lambda expression.
 *
 * \param [in] for_all_nodes functor that iterates over the objects in the mesh
 *  and for each object calls a function with the object ID, the node IDs and
 *  the number of nodes.
 * \param [in] m the Mesh to iterate over.
 * \param [in] kernel the kernel to call on each object.
 */
template < typename ExecPolicy, int NDIM, int NNODES, typename FOR_ALL_FUNCTOR,
           typename KernelType >
inline void for_all_coords( const FOR_ALL_FUNCTOR & for_all_nodes, 
                            const mint::Mesh* m, KernelType&& kernel )
{
  AXOM_STATIC_ASSERT_MSG( NDIM >= 1 && NDIM <= 3,
                          "NDIM must be a valid dimension." );
  AXOM_STATIC_ASSERT_MSG( NNODES > 0, "NNODES must be greater than zero." );
  SLIC_ASSERT( m != nullptr );
  SLIC_ASSERT( m->getDimension() == NDIM );

  constexpr bool NO_COPY = true;

  const double * const coords[3] = { 
                               m->getCoordinateArray(X_COORDINATE),
                  (NDIM > 1) ? m->getCoordinateArray(Y_COORDINATE) : nullptr,
                  (NDIM > 2) ? m->getCoordinateArray(Z_COORDINATE) : nullptr };

  for_all_nodes( ExecPolicy(), m, 
    AXOM_LAMBDA( IndexType objectID, const IndexType * nodeIDs, 
                 IndexType numNodes )
    {
      AXOM_DEBUG_VAR(numNodes);
      SLIC_ASSERT( numNodes == NNODES );

      double localCoords[ NDIM * NNODES ];
      for ( int i = 0; i < NNODES; ++i )
      {
        const int i_offset = NDIM * i;
        for ( int dim = 0; dim < NDIM; ++dim )
        {
          localCoords[ i_offset + dim ] = coords[ dim ][ nodeIDs[ i ] ];
        }
      }

      numerics::Matrix<double> coordsMatrix( NDIM, NNODES, localCoords,
                                             NO_COPY );
      kernel( objectID, coordsMatrix, nodeIDs );
    }
  );
}

} /* namespace internal */
} /* namespace mint     */
} /* namespace axom     */

#endif /* MINT_EXECUTION_HELPERS_HPP_ */