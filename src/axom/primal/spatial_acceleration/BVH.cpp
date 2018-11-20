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

#include "axom/primal/spatial_acceleration/BVH.hpp"

#include "axom/slic/interface/slic.hpp" // for SLIC macros

#ifdef AXOM_USE_RAJA
#include "RAJA/RAJA.hpp"
#endif


namespace axom
{
namespace primal
{

BVH::BVH( int dimension, const double* boxes, IndexType numItems ) :
    m_dimension( dimension ),
    m_numItems( numItems ),
    m_boxes( boxes )
{
  SLIC_ASSERT( m_dimension >= 1 && m_dimension <= 3 );
  SLIC_ASSERT( m_boxes != nullptr );
  SLIC_ASSERT( m_numItems >= 1 );

  // TODO: implement this
}

//------------------------------------------------------------------------------
BVH::~BVH()
{
  // TODO: implement this
}

//------------------------------------------------------------------------------
int BVH::build( )
{
  // TODO: build the BVH
  return BVH_BUILD_OK;
}

//------------------------------------------------------------------------------
void BVH::find( const double* pt,
                IndexType*& candidates,
                IndexType& numCandidates )
{
  SLIC_ASSERT( pt != nullptr );
  SLIC_ASSERT( candidates == nullptr );

  numCandidates = 0;
  // TODO: implement this
}

//------------------------------------------------------------------------------
void BVH::find( IndexType* offsets,
                IndexType*& candidates,
                IndexType numPts,
                const double* x,
                const double* y,
                const double* z )
{
  // TODO: implement this
}

} /* namespace primal */
} /* namespace axom */
