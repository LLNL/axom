// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)


#ifndef AXOM_PRIMAL_BVH_LINEAR_BUILDER_H_
#define AXOM_PRIMAL_BVH_LINEAR_BUILDER_H_

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
struct BVH
{
  Vec< FloatType, 4 > *m_inner_nodes;  // BVH bins including leafs
  int32  *m_leaf_nodes;   // leaf data
  AABB< FloatType, NDIMS > m_bounds;

  BVH() :
    m_inner_nodes( nullptr ),
    m_leaf_nodes( nullptr )
  {

  }

  void free()
  {
    axom::deallocate( m_inner_nodes );
    axom::deallocate( m_leaf_nodes );
  }

  ~BVH()
  {

  }

};

class LinearBVHBuilder
{

public:

  template < typename FloatType, int NDIMS >
  BVH< FloatType, NDIMS > construct(const FloatType *boxes, int size);

};

} /* namespace bvh    */
} /* namespace primal */
} /* namespace axom   */

//------------------------------------------------------------------------------
//  BVH BUILDINER IMPLEMENTATION
//------------------------------------------------------------------------------
#include "axom/primal/spatial_acceleration/linear_bvh/bvh_builder_impl.hpp"

#endif
