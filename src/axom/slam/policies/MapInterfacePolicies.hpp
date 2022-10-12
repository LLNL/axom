// Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef SLAM_MapIfacePolicies_HPP
#define SLAM_MapIfacePolicies_HPP

#include "axom/slam/Set.hpp"
#include "axom/slam/MapBase.hpp"

namespace axom
{
namespace slam
{
namespace policies
{
template <typename SetPositionType = slam::DefaultPositionType>
using VirtualMap = MapBase<SetPositionType>;

template <typename PosType = slam::DefaultPositionType,
          typename ElemType = slam::DefaultElementType>
struct ConcreteMap
{ };

}  // end namespace policies
}  // end namespace slam
}  // end namespace axom

#endif  // SLAM_MapIfacePolicies_HPP
