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
