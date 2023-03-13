// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef SLAM_SetIfacePolicies_HPP
#define SLAM_SetIfacePolicies_HPP

#include <type_traits>

#include "axom/slam/Set.hpp"
#include "axom/slam/policies/InterfacePolicies.hpp"

namespace axom
{
namespace slam
{
namespace policies
{
namespace detail
{
template <typename InterfaceTag, typename PosType, typename ElemType>
struct SetInterfaceSelector
{
  static_assert(std::is_same<InterfaceTag, ConcreteInterface>::value,
                "InterfaceTag must be one of policies::ConcreteInterface or "
                "policies::VirtualInterface.");
  using Type = ConcreteInterface;
};

template <typename PosType, typename ElemType>
struct SetInterfaceSelector<VirtualInterface, PosType, ElemType>
{
  using Type = Set<PosType, ElemType>;
};
}  // namespace detail

template <typename InterfaceTag, typename PosType, typename ElemType>
using SetInterface =
  typename detail::SetInterfaceSelector<InterfaceTag, PosType, ElemType>::Type;

}  // end namespace policies
}  // end namespace slam
}  // end namespace axom

#endif  // SLAM_SetIfacePolicies_HPP
