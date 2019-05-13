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

using namespace axom::common;

namespace axom
{
namespace primal
{
namespace bvh
{

struct BVH
{
  Vec<float32, 4>  *m_inner_nodes;  // BVH bins including leafs
  int32            *m_leaf_nodes;   // leaf data
  AABB              m_bounds;

  BVH()
  {
    m_inner_nodes = nullptr;
    m_leaf_nodes = nullptr;
  }

  void free()
  {
    axom::free(m_inner_nodes);
    axom::free(m_leaf_nodes);
  }

  ~BVH()
  {
  }
};

class LinearBVHBuilder
{

public:
  BVH construct(const double *boxes, int size);

};

} /* namespace axom */
} /* namespace primal */
} /* namespace bvh */
#endif
