// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_SPIN_BVH_TRAVERSE_HPP_
#define AXOM_SPIN_BVH_TRAVERSE_HPP_

#include "axom/config.hpp"      // compile-time definitions
#include "axom/core/Macros.hpp" // for AXOM_HOST_DEVICE
#include "axom/core/Types.hpp"  // for axom types
#include "axom/slic.hpp"        // for SLIC macros

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
inline bool leaf_node( const int32& nodeIdx )
{
  return ( nodeIdx < 0 );
}

/*!
 * \brief Generic BVH traversal routine.
 *
 * \param [in] inner_nodes pointer to the BVH bins.
 * \param [in] leaf_nodes pointer to the leaf node IDs.
 * \param [in] p the primitive in query, e.g., a point, ray, etc.
 * \param [in] L functor that defines the check for the left bin
 * \param [in] R functor that defines the check for right bin
 * \param [in] A functor that defines the leaf action
 *
 * \note The supplied functors, `L`, `R`, are expected to take the following
 *  three arguments:
 *    (1) The supplied primitive, p
 *    (2) a vec4_t< FloatType > of the first segment that defines the BVH bin
 *    (3) a vec4_t< FloatType > of the second segment that defines the BVH bin
 *
 * \see BVHData for the details on the internal data layout of the BVH.
 *
 * \note Moreover, the functors `L`, `R` return a boolean status that indicates
 *  if the specified traversal predicate is satisfied.
 *
 * \see TraversalPredicates for the collection of defined predicates.
 */
template< typename FloatType,
          typename PrimitiveType,
          typename InLeftCheck,
          typename InRightCheck,
          typename LeafAction >
AXOM_HOST_DEVICE
inline void bvh_traverse( const vec4_t< FloatType >* inner_nodes,
                          const int32* leaf_nodes,
                          const PrimitiveType& p,
                          InLeftCheck&& L,
                          InRightCheck&& R,
                          LeafAction&& A )
{
  using VecType = vec4_t< FloatType >;

  // setup stack
  constexpr int32 ISIZE = sizeof(int32);
  constexpr int32 STACK_SIZE = 64;
  constexpr int32 BARRIER = -2000000000;
  int32 todo[ STACK_SIZE ];
  int32 stackptr   = 0;
  todo[ stackptr ] = BARRIER;

  int32 current_node      = 0;
  while ( current_node != BARRIER )
  {

    if ( !leaf_node( current_node ) )
    {
      const VecType first4  = inner_nodes[current_node + 0];
      const VecType second4 = inner_nodes[current_node + 1];
      const VecType third4  = inner_nodes[current_node + 2];

      const bool in_left  = L( p, first4, second4 );
      const bool in_right = R( p, second4, third4 );

      if ( !in_left && !in_right )
      {
        // pop the stack and continue
        current_node = todo[ stackptr ];
        stackptr--;
      }
      else
      {
        VecType children = inner_nodes[current_node + 3];

        // memcpy the int bits hidden in the floats
        int32 l_child;
        memcpy(&l_child, &children[0], ISIZE);
        int32 r_child;
        memcpy(&r_child, &children[1], ISIZE);

        current_node = (in_left) ? l_child : r_child;

        if (in_left && in_right)
        {
          stackptr++;
          todo[stackptr] = r_child;
          // TODO: if we are in both children we could
          // go down the "closer" first by perhaps the distance
          // from the point to the center of the aabb
        }

      } // END else

    } // END if
    else
    {
      // compute leaf index
      current_node = -current_node - 1; // swap the neg address

      // execute leaf action
      A( current_node, leaf_nodes );

      current_node = todo[stackptr];
      stackptr--;
    } // END else

  } // END while

}

} /* namespace linear_bvh */
} /* namespace internal */
} /* namespace spin */
} /* namespace axom */



#endif /* AXOM_SPIN_BVH_TRAVERSE_HPP_ */
