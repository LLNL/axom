// Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef SLAM_BivarSetIfacePolicies_HPP
#define SLAM_BivarSetIfacePolicies_HPP

#include "axom/slam/BivariateSet.hpp"

namespace axom
{
namespace slam
{
namespace policies
{
template <typename Set1 = slam::Set<>, typename Set2 = slam::Set<>>
using VirtualBivariateSet = BivariateSet<Set1, Set2>;

template <typename Derived, typename Set1 = slam::Set<>, typename Set2 = slam::Set<>>
struct ConcreteBivariateSet
{
  using FirstSetType = Set1;
  using SecondSetType = Set2;

  using PositionType = typename FirstSetType::PositionType;
  using ElementType = typename FirstSetType::ElementType;

  using SubsetType = void;

  static constexpr PositionType INVALID_POS = PositionType(-1);

  /** \brief Size of the first set.   */
  inline PositionType firstSetSize() const
  {
    return getSize<FirstSetType>(getFirstSetDerived());
  }
  /** \brief Size of the second set.   */
  inline PositionType secondSetSize() const
  {
    return getSize<SecondSetType>(getSecondSetDerived());
  }

private:
  template <typename SetType>
  typename std::enable_if<std::is_abstract<SetType>::value, PositionType>::type
  getSize(const SetType* s) const
  {
    SLIC_ASSERT_MSG(s != nullptr, "nullptr in BivariateSet::getSize()");
    return s->size();
  }

  template <typename SetType>
  typename std::enable_if<!std::is_abstract<SetType>::value, PositionType>::type
  getSize(const SetType* s) const
  {
    if(policies::EmptySetTraits<SetType>::emptySet() == s)
    {
      return 0;
    }
    return static_cast<SetType>(*s).size();
  }

protected:
  const FirstSetType* getFirstSetDerived() const
  {
    return static_cast<const Derived*>(this)->getFirstSet();
  }
  const SecondSetType* getSecondSetDerived() const
  {
    return static_cast<const Derived*>(this)->getSecondSet();
  }
};

}  // end namespace policies
}  // end namespace slam
}  // end namespace axom

#endif  // SLAM_BivarSetIfacePolicies_HPP
