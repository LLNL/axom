// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file IndirectionSet.hpp
 *
 * \brief Defines some alias templates for OrderedSets with indirection
 */

#ifndef SLAM_INDIRECTION_SET_H_
#define SLAM_INDIRECTION_SET_H_

#include <cstddef>
#include <vector>

#include "axom/slam/OrderedSet.hpp"

namespace axom
{
namespace slam
{
/**
 * \brief Alias template for an OrderedSet with indirection over an array
 *
 * \tparam PosType The position type for indexing into the set
 * \tparam ElemType The type for the set's elements
 * \sa OrderedSet
 */
template <typename PosType = slam::DefaultPositionType,
          typename ElemType = slam::DefaultElementType>
using ArrayIndirectionSet =
  OrderedSet<PosType,
             ElemType,
             policies::RuntimeSize<PosType>,
             policies::ZeroOffset<PosType>,
             policies::StrideOne<PosType>,
             policies::ArrayIndirection<PosType, ElemType>>;

/**
 * \brief Alias template for an OrderedSet with indirection over an stl vector
 *
 * \tparam PosType The position type for indexing into the set
 * \tparam ElemType The type for the set's elements
 * \sa OrderedSet
 */
template <typename PosType = slam::DefaultPositionType,
          typename ElemType = slam::DefaultElementType>
using VectorIndirectionSet =
  OrderedSet<PosType,
             ElemType,
             policies::RuntimeSize<PosType>,
             policies::ZeroOffset<PosType>,
             policies::StrideOne<PosType>,
             policies::STLVectorIndirection<PosType, ElemType>>;

}  // end namespace slam
}  // end namespace axom

#endif  //  SLAM_INDIRECTION_SET_H_
