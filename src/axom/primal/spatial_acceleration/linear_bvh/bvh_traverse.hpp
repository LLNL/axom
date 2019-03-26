// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_PRIMAL_BVH_TRAVERSE_H_
#define AXOM_PRIMAL_BVH_TRAVERSE_H_

#include "axom/core/Types.hpp"
#include "axom/core/Macros.hpp"

#include "axom/primal/spatial_acceleration/linear_bvh/bvh_builder.hpp"

namespace axom
{
namespace primal
{

/*!
 * \brief Generic BVH traversal routine.
 *
 * \param [in] bvh the internal BVH data structure.
 * \param [in] inLeft lamda that determines if the left branch is traversed.
 * \param [in] inRight lambda that determines if the right branch is traversed.
 * \param [in] leafKernel lambda the defines what to when a leaf is reached.
 *
 * \note The supplied lambdas are intended to be thread local instances. They
 *  are created within a thread and therefore must be device decorated when
 *  called on an NVIDIA GPU device.
 */
template < int NDIMS, typename FloatType,
           typename InLeftType,
           typename InRightType,
           typename LeafKernelType >
AXOM_HOST_DEVICE void bvh_traverse( bvh::BVH< FloatType, NDIMS >& bvh,
                                    InLeftType&& inLeft,
                                    InRightType&& inRight,
                                    LeafKernelType&& leafKernel )
{
  constexpr int32 BARRIER    = -2000000000;
  constexpr int32 STACK_SIZE = 64;
  int32 current_node = 0;
  int32 stackptr     = 0;

  int32 todo[ STACK_SIZE ];
  todo[ stackptr ] = BARRIER;

  const bvh::Vec< FloatType, 4 >* inner_nodes = bvh.m_inner_nodes;

  while ( current_node != BARRIER )
  {

    if ( current_node <= -1 )
    {
      // leaf node
      current_node = -current_node - 1; // swap the neg address

      leafKernel( current_node, bvh );

      current_node = todo[ stackptr ];
      stackptr--;
    }
    else
    {
      // inner node
      const bvh::Vec< FloatType,4 > s1 = inner_nodes[ current_node     ];
      const bvh::Vec< FloatType,4 > s2 = inner_nodes[ current_node + 1 ];
      const bvh::Vec< FloatType,4 > s3 = inner_nodes[ current_node + 2 ];

      const bool in_left  = inLeft( s1, s2 );
      const bool in_right = inRight( s2, s3 );

      if ( !in_left && !in_right )
      {
        // pop the stack and continue
        current_node = todo[ stackptr ];
        stackptr--;
      }
      else
      {
        // traverse down the tree
        const bvh::Vec< FloatType,4 > children = inner_nodes[current_node + 3];

        constexpr int32 isize = sizeof(int32);
        int32 l_child;
        memcpy( &l_child, &children[0], isize );

        int32 r_child;
        memcpy( &r_child, &children[1], isize );

        current_node = ( in_left ) ? l_child : r_child;

        if ( in_left && in_right )
        {
          stackptr++;
          todo[ stackptr ] = r_child;
          // TODO: if we are in both children we could
          // go down the "closer" first by perhaps the distance
          // from the point to the center of the aabb
        }

      } // END else

    } // END else isLeaf()

  } // END while

}

} /* namespace primal */
} /* namespace axom */

#endif


