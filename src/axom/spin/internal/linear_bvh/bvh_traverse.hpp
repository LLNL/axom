// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_SPIN_BVH_TRAVERSE_HPP_
#define AXOM_SPIN_BVH_TRAVERSE_HPP_

#include "axom/config.hpp"       // compile-time definitions
#include "axom/core/Macros.hpp"  // for AXOM_HOST_DEVICE
#include "axom/core/Types.hpp"   // for axom types
#include "axom/slic.hpp"         // for SLIC macros

#include "axom/spin/internal/linear_bvh/BVHData.hpp"

namespace axom
{
namespace spin
{
namespace internal
{
namespace linear_bvh
{
/*!
 * \brief Checks if the node corresponding to the given node ID is a leaf node.
 * \param [in] nodeIdx index of the BVH node in query.
 * \return status true if the node is a leaf node, otherwise, false.
 */
AXOM_HOST_DEVICE
inline bool leaf_node(const int32& nodeIdx) { return (nodeIdx < 0); }

/*!
 * \brief Generic BVH traversal routine.
 *
 * \param [in] inner_nodes pointer to the BVH bins.
 * \param [in] inner_node_children pointer to pairs of child indices.
 * \param [in] leaf_nodes pointer to the leaf node IDs.
 * \param [in] p the primitive in query, e.g., a point, ray, etc.
 * \param [in] B functor that defines the check for the bins
 * \param [in] A functor that defines the leaf action
 *
 * \note The supplied functor `B` is expected to take the following two
 *  arguments:
 *    (1) The supplied primitive, p
 *    (2) a primal::BoundingBox< FloatType, NDIMS > of the BVH bin
 *
 * \see BVHData for the details on the internal data layout of the BVH.
 *
 * \note Moreover, the functor `B` returns a boolean status that indicates
 *  if the specified traversal predicate is satisfied.
 *
 */
template <int NDIMS, typename FloatType, typename PrimitiveType, typename InBinCheck, typename LeafAction>
AXOM_HOST_DEVICE inline void bvh_traverse(
  const primal::BoundingBox<FloatType, NDIMS>* inner_nodes,
  const int32* inner_node_children,
  const int32* leaf_nodes,
  const PrimitiveType& p,
  InBinCheck&& B,
  LeafAction&& A)
{
  using BBoxType = primal::BoundingBox<FloatType, NDIMS>;

  // setup stack
  constexpr int32 STACK_SIZE = 64;
  constexpr int32 BARRIER = -2000000000;
  int32 todo[STACK_SIZE];
  int32 stackptr = 0;
  todo[stackptr] = BARRIER;

  int32 current_node = 0;
  while(current_node != BARRIER)
  {
    if(!leaf_node(current_node))
    {
      BBoxType left_bin = inner_nodes[current_node + 0];
      BBoxType right_bin = inner_nodes[current_node + 1];
      const bool in_left = B(p, left_bin);
      const bool in_right = B(p, right_bin);
      int32 l_child = inner_node_children[current_node + 0];
      int32 r_child = inner_node_children[current_node + 1];

      if(!in_left && !in_right)
      {
        // pop the stack and continue
        current_node = todo[stackptr];
        stackptr--;
      }
      else
      {
        current_node = (in_left) ? l_child : r_child;

        if(in_left && in_right)
        {
          stackptr++;
          todo[stackptr] = r_child;
          // TODO: if we are in both children we could
          // go down the "closer" first by perhaps the distance
          // from the point to the center of the aabb
        }

      }  // END else

    }  // END if
    else
    {
      // compute leaf index
      current_node = -current_node - 1;  // swap the neg address

      // execute leaf action
      A(current_node, leaf_nodes);

      current_node = todo[stackptr];
      stackptr--;
    }  // END else

  }  // END while
}

} /* namespace linear_bvh */
} /* namespace internal */
} /* namespace spin */
} /* namespace axom */

#endif /* AXOM_SPIN_BVH_TRAVERSE_HPP_ */
