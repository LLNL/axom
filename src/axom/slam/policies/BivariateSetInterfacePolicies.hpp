// Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef SLAM_BivarSetIfacePolicies_HPP
#define SLAM_BivarSetIfacePolicies_HPP

#include "axom/slam/BivariateSet.hpp"
#include "axom/slam/policies/InterfacePolicies.hpp"

namespace axom
{
namespace slam
{
namespace policies
{
namespace detail
{
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

  bool isValid(bool verboseOutput = false) const
  {
    if(getFirstSetDerived() == nullptr || getSecondSetDerived() == nullptr)
    {
      if(verboseOutput)
      {
        SLIC_INFO("BivariateSet is not valid: "
                  << " Set pointers should not be null.");
      }
      return false;
    }
    return getFirstSetDerived()->isValid(verboseOutput) &&
      getSecondSetDerived()->isValid(verboseOutput);
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

template <typename InterfaceTag, typename Set1, typename Set2, typename Derived>
struct BSetInterfaceSelector
{
  static_assert(std::is_same<InterfaceTag, ConcreteInterface>::value,
                "InterfaceTag must be one of policies::ConcreteInterface or "
                "policies::VirtualInterface.");
  using Type = ConcreteBivariateSet<Derived, Set1, Set2>;
};

template <typename Set1, typename Set2, typename Derived>
struct BSetInterfaceSelector<VirtualInterface, Set1, Set2, Derived>
{
  using Type = BivariateSet<Set1, Set2>;
};
}  // namespace detail

template <typename InterfaceTag, typename FromSet, typename ToSet, typename Derived>
using BivariateSetInterface =
  typename detail::BSetInterfaceSelector<InterfaceTag, FromSet, ToSet, Derived>::Type;

}  // end namespace policies
}  // end namespace slam
}  // end namespace axom

#endif  // SLAM_BivarSetIfacePolicies_HPP
