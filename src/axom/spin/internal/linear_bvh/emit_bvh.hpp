// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_SPIN_EMIT_BVH_H_
#define AXOM_SPIN_EMIT_BVH_H_

// axom core includes
#include "axom/core/Types.hpp"                       // for fixed bitwidth types
#include "axom/core/memory_management.hpp"           // for alloc()/free()
#include "axom/core/utilities/AnnotationMacros.hpp"  // for annotations
#include "axom/slic/interface/slic_macros.hpp"       // for SLIC_ASSERT()

#include "axom/core/execution/execution_space.hpp"
#include "axom/core/execution/for_all.hpp"

#include "axom/primal/geometry/BoundingBox.hpp"

#include "axom/spin/internal/linear_bvh/BVHData.hpp"
#include "axom/spin/internal/linear_bvh/RadixTree.hpp"

#include "RAJA/RAJA.hpp"

namespace axom
{
namespace spin
{
namespace internal
{
namespace linear_bvh
{
/*!
 * \brief Given a RadixTree, this method emits a corresponding BVH.
 *
 * \param [in]  data reference to the radix tree data
 * \param [out] bvh_data referene to the internal BVH data structure.
 *
 * \tparam FloatType the floating point precision, e.g., `double` or `float`
 *
 * \see BVH::build()
 */
template <typename ExecSpace, typename FloatType, int NDIMS>
void emit_bvh(RadixTree<FloatType, NDIMS>& data,
              BVHData<FloatType, NDIMS>& bvh_data)
{
  AXOM_PERF_MARK_FUNCTION("emit_bvh");

  using BoundingBoxType = primal::BoundingBox<FloatType, NDIMS>;

  const int32 size = data.m_size;
  const int32 inner_size = data.m_inner_size;
  SLIC_ASSERT(inner_size == size - 1);

  const int32* lchildren_ptr = data.m_left_children;
  const int32* rchildren_ptr = data.m_right_children;

  const BoundingBoxType* leaf_aabb_ptr = data.m_leaf_aabbs;
  const BoundingBoxType* inner_aabb_ptr = data.m_inner_aabbs;

  primal::BoundingBox<FloatType, NDIMS>* bvh_inner_nodes = bvh_data.m_inner_nodes;
  int32* bvh_inner_node_children = bvh_data.m_inner_node_children;

  AXOM_PERF_MARK_SECTION("emit_bvh_parents",
                         for_all<ExecSpace>(
                           inner_size,
                           AXOM_LAMBDA(int32 node) {
                             BoundingBoxType l_aabb, r_aabb;

                             int32 lchild = lchildren_ptr[node];
                             if(lchild >= inner_size)
                             {
                               l_aabb = leaf_aabb_ptr[lchild - inner_size];
                               lchild = -(lchild - inner_size + 1);
                             }
                             else
                             {
                               l_aabb = inner_aabb_ptr[lchild];
                               // do the offset now
                               lchild *= 2;
                             }

                             int32 rchild = rchildren_ptr[node];
                             if(rchild >= inner_size)
                             {
                               r_aabb = leaf_aabb_ptr[rchild - inner_size];
                               rchild = -(rchild - inner_size + 1);
                             }
                             else
                             {
                               r_aabb = inner_aabb_ptr[rchild];
                               // do the offset now
                               rchild *= 2;
                             }

                             const int32 out_offset = node * 2;
                             bvh_inner_nodes[out_offset + 0] = l_aabb;
                             bvh_inner_nodes[out_offset + 1] = r_aabb;

                             bvh_inner_node_children[out_offset + 0] = lchild;
                             bvh_inner_node_children[out_offset + 1] = rchild;
                           }););

  int32* radix_tree_leafs = data.m_leafs;
  int32* bvh_leafs = bvh_data.m_leaf_nodes;

  AXOM_PERF_MARK_SECTION(
    "emit_bvh_leafs",
    for_all<ExecSpace>(
      size,
      AXOM_LAMBDA(int32 i) { bvh_leafs[i] = radix_tree_leafs[i]; }););
}

} /* namespace linear_bvh */
} /* namespace internal */
} /* namespace spin */
} /* namespace axom */
#endif
