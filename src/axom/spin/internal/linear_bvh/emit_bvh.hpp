// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)


#ifndef AXOM_SPIN_EMIT_BVH_H_
#define AXOM_SPIN_EMIT_BVH_H_

// axom core includes
#include "axom/core/Types.hpp"                        // for fixed bitwidth types
#include "axom/core/memory_management.hpp"            // for alloc()/free()
#include "axom/core/utilities/AnnotationMacros.hpp"   // for annotations
#include "axom/slic/interface/slic_macros.hpp"        // for SLIC_ASSERT()

#include "axom/core/execution/execution_space.hpp"
#include "axom/core/execution/for_all.hpp"

#include "axom/spin/internal/linear_bvh/vec.hpp"
#include "axom/spin/internal/linear_bvh/aabb.hpp"
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
/// @{

template < typename ExecSpace, typename FloatType  >
void emit_bvh(RadixTree< FloatType,3 >& data, BVHData< FloatType,3 >& bvh_data);

template < typename ExecSpace, typename FloatType  >
void emit_bvh(RadixTree< FloatType,2 >& data, BVHData< FloatType,2 >& bvh_data);

/// @}

//------------------------------------------------------------------------------
//                        IMPLEMENTATION
//------------------------------------------------------------------------------
template < typename ExecSpace, typename FloatType  >
void emit_bvh( RadixTree<FloatType, 3>& data,
               BVHData< FloatType, 3 >& bvh_data )
{
  AXOM_PERF_MARK_FUNCTION( "emit_bvh3D" );

  const int32 size       = data.m_size;
  const int32 inner_size = data.m_inner_size;
  SLIC_ASSERT( inner_size == size-1 );

  const int32* lchildren_ptr = data.m_left_children;
  const int32* rchildren_ptr = data.m_right_children;

  const AABB<FloatType,3>* leaf_aabb_ptr  = data.m_leaf_aabbs;
  const AABB<FloatType,3>* inner_aabb_ptr = data.m_inner_aabbs;

  Vec<FloatType,4>* flat_ptr = bvh_data.m_inner_nodes;

  AXOM_PERF_MARK_SECTION( "emit_bvh_parents",
    for_all< ExecSpace >( inner_size, AXOM_LAMBDA (int32 node)
    {
      Vec<FloatType,4> vec1;
      Vec<FloatType,4> vec2;
      Vec<FloatType,4> vec3;
      Vec<FloatType,4> vec4;

      AABB<FloatType,3> l_aabb, r_aabb;

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
        lchild *= 4;
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
        rchild *= 4;
      }
      vec1[0] = l_aabb.m_x.min();
      vec1[1] = l_aabb.m_y.min();
      vec1[2] = l_aabb.m_z.min();

      vec1[3] = l_aabb.m_x.max();
      vec2[0] = l_aabb.m_y.max();
      vec2[1] = l_aabb.m_z.max();

      vec2[2] = r_aabb.m_x.min();
      vec2[3] = r_aabb.m_y.min();
      vec3[0] = r_aabb.m_z.min();

      vec3[1] = r_aabb.m_x.max();
      vec3[2] = r_aabb.m_y.max();
      vec3[3] = r_aabb.m_z.max();

      const int32 out_offset = node * 4;
      flat_ptr[out_offset + 0] = vec1;
      flat_ptr[out_offset + 1] = vec2;
      flat_ptr[out_offset + 2] = vec3;

      constexpr int32 isize = sizeof(int32);
      // memcopy so we do not truncate the ints
      memcpy(&vec4[0], &lchild, isize);
      memcpy(&vec4[1], &rchild, isize);
      flat_ptr[out_offset + 3] = vec4;
    } );
  );

  int32* radix_tree_leafs = data.m_leafs;
  int32* bvh_leafs        = bvh_data.m_leaf_nodes;

  AXOM_PERF_MARK_SECTION( "emit_bvh_leafs",
    for_all< ExecSpace >( size, AXOM_LAMBDA(int32 i)
    {
      bvh_leafs[ i ] = radix_tree_leafs[ i ];
    } );
  );

}

//------------------------------------------------------------------------------
template < typename ExecSpace, typename FloatType  >
void emit_bvh( RadixTree<FloatType, 2>& data,
               BVHData< FloatType, 2 >& bvh_data )
{
  AXOM_PERF_MARK_FUNCTION( "emit_bvh2D" );

  const int32 size       = data.m_size;
  const int32 inner_size = data.m_inner_size;
  SLIC_ASSERT( inner_size == size-1 );

  const int32* lchildren_ptr = data.m_left_children;
  const int32* rchildren_ptr = data.m_right_children;

  const AABB<FloatType,2>* leaf_aabb_ptr  = data.m_leaf_aabbs;
  const AABB<FloatType,2>* inner_aabb_ptr = data.m_inner_aabbs;


  Vec<FloatType,4>* flat_ptr = bvh_data.m_inner_nodes;

  AXOM_PERF_MARK_SECTION( "emit_bvh_parents",
    for_all< ExecSpace >( inner_size, AXOM_LAMBDA (int32 node)
    {
      Vec<FloatType,4> vec1;
      Vec<FloatType,4> vec2;
      Vec<FloatType,4> vec3;
      Vec<FloatType,4> vec4;

      AABB<FloatType,2> l_aabb, r_aabb;

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
        lchild *= 4;
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
        rchild *= 4;
      }
      vec1[0] = l_aabb.m_x.min();
      vec1[1] = l_aabb.m_y.min();
      vec1[2] = 0.0;

      vec1[3] = l_aabb.m_x.max();
      vec2[0] = l_aabb.m_y.max();
      vec2[1] = 0.0;

      vec2[2] = r_aabb.m_x.min();
      vec2[3] = r_aabb.m_y.min();
      vec3[0] = 0.0;

      vec3[1] = r_aabb.m_x.max();
      vec3[2] = r_aabb.m_y.max();
      vec3[3] = 0.0;

      const int32 out_offset = node * 4;
      flat_ptr[out_offset + 0] = vec1;
      flat_ptr[out_offset + 1] = vec2;
      flat_ptr[out_offset + 2] = vec3;

      constexpr int32 isize = sizeof(int32);
      // memcopy so we do not truncate the ints
      memcpy(&vec4[0], &lchild, isize);
      memcpy(&vec4[1], &rchild, isize);
      flat_ptr[out_offset + 3] = vec4;
    } );
  );

  int32* radix_tree_leafs = data.m_leafs;
  int32* bvh_leafs        = bvh_data.m_leaf_nodes;

  AXOM_PERF_MARK_SECTION( "emit_bvh_leafs",
    for_all< ExecSpace >( size, AXOM_LAMBDA(int32 i)
    {
      bvh_leafs[ i ] = radix_tree_leafs[ i ];
    } );
  );

}


} /* namespace linear_bvh */
} /* namespace internal */
} /* namespace spin */
} /* namespace axom */
#endif
