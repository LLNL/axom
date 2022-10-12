// Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef SLAM_SetIfacePolicies_HPP
#define SLAM_SetIfacePolicies_HPP

#include "axom/slam/Set.hpp"

namespace axom
{
namespace slam
{
namespace policies
{
template <typename PosType = slam::DefaultPositionType,
          typename ElemType = slam::DefaultElementType>
using VirtualSet = Set<PosType, ElemType>;

template <typename PosType = slam::DefaultPositionType,
          typename ElemType = slam::DefaultElementType>
struct ConcreteSet
{ };

}  // end namespace policies
}  // end namespace slam
}  // end namespace axom

#endif  // SLAM_SetIfacePolicies_HPP
