/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC.
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
  Vec<float32, 4>  *m_inner_nodes;
  int32            *m_leaf_nodes;
  AABB              m_bounds;

  BVH()
  {
    m_inner_nodes = nullptr;
    m_leaf_nodes = nullptr;
  }

  ~BVH()
  {
    axom::free(m_inner_nodes);
    axom::free(m_leaf_nodes);
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
