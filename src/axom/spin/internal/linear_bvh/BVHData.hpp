// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_SPIN_BVHDATA_HPP_
#define AXOM_SPIN_BVHDATA_HPP_

// axom core includes
#include "axom/core/Types.hpp"               // for fixed bitwidth types
#include "axom/core/memory_management.hpp"   // for alloc()/free()

// spin includes
#include "axom/spin/internal/linear_bvh/vec.hpp"
#include "axom/spin/internal/linear_bvh/aabb.hpp"

namespace axom
{
namespace spin
{
namespace internal
{
namespace linear_bvh
{

/*!
 * \brief BVHData provides a data-structure that represent the internal data
 *  layout of the BVH.
 *
 * \note <b> Internal Data Layout </b>
 * \verbatim
 *                      |               |               |               |
 *             S0       |      S1       |       S2      |       S3      |
 *      +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
 * .... | 0 | 1 | 2 | 3 | 0 | 1 | 2 | 3 | 0 | 1 | 2 | 3 | 0 | 1 | X | X | .....
 *      +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
 *                      |               |               |
 *                      |               |               |
 *        S0.0 = L.xmin | S1.0 = L.ymax | S2.0 = R.zmin | S3.0 = left_child
 *        S0.1 = L.ymin | S1.1 = L.zmax | S2.1 = R.xmax | S3.1 = right_child
 *        S0.2 = L.zmin | S1.2 = R.xmin | S2.2 = R.ymax | S3.2 = padded
 *        S0.3 = L.xmax | S1.3 = R.ymin | S2.3 = R.zmax | S3.3 = padded
 *
 * \endverbatim
 *
 * \note The internal data layout is organized in a flat buffer of 4 segments,
 *  where each segment is a Vec< FloatType, 4 > type, that stores the left
 *  and right boxes of a given node, as well as, the IDs of the right and
 *  left children, as illustrated above.
 *
 * \note A Vec< FloatType,4 > is chosen b/c it fits into GPU texture memory and
 *  allows for additional performance optimizations down the road
 *
 */
template < typename FloatType, int NDIMS >
struct BVHData
{
  Vec< FloatType, 4 > *m_inner_nodes;  // BVH bins including leafs
  int32  *m_leaf_nodes;   // leaf data
  AABB< FloatType, NDIMS > m_bounds;

  BVHData() :
    m_inner_nodes( nullptr ),
    m_leaf_nodes( nullptr )
  {

  }

  void allocate( int32 size )
  {
    m_inner_nodes = axom::allocate< Vec< FloatType,4 > >( (size-1)*4 );
    m_leaf_nodes  = axom::allocate< int32 >( size );
  }

  void deallocate()
  {
    axom::deallocate( m_inner_nodes );
    axom::deallocate( m_leaf_nodes );
  }

  ~BVHData()
  {

  }

};


} /* namespace linear_bvh */
} /* namespace internal */
} /* namespace spin */
} /* namespace axom */
#endif /* AXOM_SPIN_BVHDATA_HPP_ */
