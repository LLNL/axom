// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_SPIN_BVHDATA_HPP_
#define AXOM_SPIN_BVHDATA_HPP_

// axom core includes
#include "axom/core/Types.hpp"              // for fixed bitwidth types
#include "axom/core/memory_management.hpp"  // for alloc()/free()

#include "axom/core/utilities/AnnotationMacros.hpp"  // for annotations

#include "axom/primal/geometry/BoundingBox.hpp"
#include "axom/primal/geometry/Vector.hpp"

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
 * \note The internal node data is organized into two arrays, one containing
 *  bounding boxes of the internal nodes' two child nodes, and the other
 *  containing indices of the two child nodes (scaled by two if inner node,
 *  ones-complement if leaf node).
 *
 */
template <typename FloatType, int NDIMS>
struct BVHData
{
  using BoundingBoxType = primal::BoundingBox<FloatType, NDIMS>;

  BoundingBoxType* m_inner_nodes;  // BVH bins including leafs
  int32* m_inner_node_children;
  int32* m_leaf_nodes;  // leaf data
  primal::BoundingBox<FloatType, NDIMS> m_bounds;

  BVHData() : m_inner_nodes(nullptr), m_leaf_nodes(nullptr) { }

  void allocate(int32 size, int allocID)
  {
    AXOM_PERF_MARK_FUNCTION("BVHData::allocate");
    m_inner_nodes = axom::allocate<BoundingBoxType>((size - 1) * 2, allocID);
    m_inner_node_children = axom::allocate<int32>((size - 1) * 2, allocID);
    m_leaf_nodes = axom::allocate<int32>(size, allocID);
  }

  void deallocate()
  {
    AXOM_PERF_MARK_FUNCTION("BVHData::deallocate");
    axom::deallocate(m_inner_nodes);
    axom::deallocate(m_inner_node_children);
    axom::deallocate(m_leaf_nodes);
  }

  ~BVHData() { }
};

} /* namespace linear_bvh */
} /* namespace internal */
} /* namespace spin */
} /* namespace axom */
#endif /* AXOM_SPIN_BVHDATA_HPP_ */
