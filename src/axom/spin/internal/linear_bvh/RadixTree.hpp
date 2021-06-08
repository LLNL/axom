// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_SPIN_RADIXTREE_HPP_
#define AXOM_SPIN_RADIXTREE_HPP_

#include "axom/core/memory_management.hpp"

#include "axom/core/utilities/AnnotationMacros.hpp"  // for annotations

#include "axom/primal/geometry/BoundingBox.hpp"

namespace axom
{
namespace spin
{
namespace internal
{
namespace linear_bvh
{
/*!
 * \brief RadixTree provides a binary radix tree representation that stores a
 *  list of axis-aligned bounding boxes sorted according to their Morton code.
 *
 * \note This data-structure provides an intermediate representation that serves
 *  as the building-block to construct a BVH in parallel.
 */
template <typename FloatType, int NDIMS>
struct RadixTree
{
  using BoxType = primal::BoundingBox<FloatType, NDIMS>;

  int32 m_size;
  int32 m_inner_size;

  int32* m_left_children;
  int32* m_right_children;
  int32* m_parents;
  BoxType* m_inner_aabbs;

  int32* m_leafs;
  uint32* m_mcodes;
  BoxType* m_leaf_aabbs;

  void allocate(int32 size, int allocID)
  {
    AXOM_PERF_MARK_FUNCTION("RadixTree::allocate");

    m_size = size;
    m_inner_size = m_size - 1;

    m_left_children = axom::allocate<int32>(m_inner_size, allocID);
    m_right_children = axom::allocate<int32>(m_inner_size, allocID);
    m_parents = axom::allocate<int32>((m_size + m_inner_size), allocID);
    m_inner_aabbs = axom::allocate<BoxType>(m_inner_size, allocID);

    m_leafs = axom::allocate<int32>(m_size, allocID);
    m_mcodes = axom::allocate<uint32>(m_size, allocID);
    m_leaf_aabbs = axom::allocate<BoxType>(m_size, allocID);
  }

  void deallocate()
  {
    AXOM_PERF_MARK_FUNCTION("RadixTree::deallocate");

    m_inner_size = 0;
    m_size = 0;

    axom::deallocate(m_left_children);
    axom::deallocate(m_right_children);
    axom::deallocate(m_parents);
    axom::deallocate(m_inner_aabbs);

    axom::deallocate(m_leafs);
    axom::deallocate(m_mcodes);
    axom::deallocate(m_leaf_aabbs);
  }
};

} /* namespace linear_bvh */
} /* namespace internal */
} /* namespace spin */
} /* namespace axom */

#endif /* AXOM_RADIXTREE_HPP_ */
