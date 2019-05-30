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
#include "axom/primal/spatial_acceleration/linear_bvh/BVHData.hpp"

namespace axom
{
namespace primal
{
namespace bvh
{



class LinearBVHBuilder
{

public:

  template < typename FloatType, int NDIMS >
  void construct( const FloatType *boxes,
                  int size,
                  BVHData< FloatType, NDIMS >& bvh_data );

};

} /* namespace bvh    */
} /* namespace primal */
} /* namespace axom   */

//------------------------------------------------------------------------------
//  BVH BUILDINER IMPLEMENTATION
//------------------------------------------------------------------------------
#include "axom/primal/spatial_acceleration/linear_bvh/bvh_builder_impl.hpp"

#endif
