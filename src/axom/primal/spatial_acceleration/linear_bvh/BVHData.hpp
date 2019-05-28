// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_PRIMAL_BVHDATA_HPP_
#define AXOM_PRIMAL_BVHDATA_HPP_

// axom core includes
#include "axom/core/Types.hpp"               // for fixed bitwidth types
#include "axom/core/memory_management.hpp"   // for alloc()/free()

// primal includes
#include "axom/primal/spatial_acceleration/linear_bvh/aabb.hpp"

namespace axom
{
namespace primal
{
namespace bvh
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

  void free()
  {
    axom::deallocate( m_inner_nodes );
    axom::deallocate( m_leaf_nodes );
  }

  ~BVHData()
  {

  }

};

} /* namespace bvh */
} /* namespace primal */
} /* namespace axom */

#endif /* AXOM_PRIMAL_BVHDATA_HPP_ */
