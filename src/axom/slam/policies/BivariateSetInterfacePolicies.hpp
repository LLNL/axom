// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file BivariateSetInterfacePolicies.hpp
 * \brief Select virtual or concrete interfaces for bivariate sets.
 */

#pragma once

#include "axom/slam/BivariateSet.hpp"
#include "axom/slam/policies/InterfacePolicies.hpp"

#include <utility>

namespace axom::slam::policies
{
namespace detail
{
/*!
 * \brief The base class of concrete-interface bivariate sets.
 *
 * Stores the first and second set pointers and defines their position types.
 * Derived sets supply the pair access and traversal operations.
 */
template <typename Set1 = slam::Set<>,
          typename Set2 = slam::Set<>,
          typename Position = ::axom::slam::detail::default_flat_position_t<typename Set1::PositionType,
                                                                            typename Set2::PositionType>>
struct ConcreteBivariateSet
{
  using FirstSetType = Set1;
  using SecondSetType = Set2;

  using FirstPositionType = typename FirstSetType::PositionType;
  using SecondPositionType = typename SecondSetType::PositionType;
  using PositionType = Position;
  using ElementType = std::pair<FirstPositionType, SecondPositionType>;

  using SubsetType = void;
  using RangeSetType = typename RangeSet<PositionType, PositionType>::ConcreteSet;

  static constexpr PositionType INVALID_POS = PositionType(-1);

  /*!
   * \brief Bind the first and second sets, which must outlive this object.
   *
   * \param set1  Pointer to the first Set.
   * \param set2  Pointer to the second Set.
   */
  ConcreteBivariateSet(const Set1* set1 = policies::EmptySetTraits<Set1>::emptySet(),
                       const Set2* set2 = policies::EmptySetTraits<Set2>::emptySet())
    : m_set1(set1)
    , m_set2(set2)
  { }

  /// \brief Size of the first set.
  AXOM_HOST_DEVICE inline FirstPositionType firstSetSize() const
  {
    return getSize<FirstSetType>(m_set1);
  }
  /// \brief Size of the second set.
  AXOM_HOST_DEVICE inline SecondPositionType secondSetSize() const
  {
    return getSize<SecondSetType>(m_set2);
  }

  /// \brief Returns pointer to the first set.
  const FirstSetType* getFirstSet() const { return m_set1; }
  /// \brief Returns pointer to the second set.
  const SecondSetType* getSecondSet() const { return m_set2; }

  bool isValid(bool verboseOutput = false) const
  {
    if(getFirstSet() == nullptr || getSecondSet() == nullptr)
    {
      if(verboseOutput)
      {
        SLIC_INFO("BivariateSet is not valid: " << " Set pointers should not be null.");
      }
      return false;
    }
    if constexpr(Validatable<FirstSetType>)
    {
      if(!getFirstSet()->isValid(verboseOutput))
      {
        return false;
      }
    }
    if constexpr(Validatable<SecondSetType>)
    {
      if(!getSecondSet()->isValid(verboseOutput))
      {
        return false;
      }
    }
    return true;
  }

private:
  template <typename SetType>
  typename SetType::PositionType getSize(const SetType* s) const
    requires(std::is_abstract_v<SetType>)
  {
    SLIC_ASSERT_MSG(s != nullptr, "nullptr in BivariateSet::getSize()");
    return s->size();
  }

  template <typename SetType>
  AXOM_HOST_DEVICE typename SetType::PositionType getSize(const SetType* s) const
    requires(!std::is_abstract_v<SetType>)
  {
    SLIC_ASSERT_MSG(s != nullptr, "nullptr in BivariateSet::getSize()");
    return static_cast<SetType>(*s).size();
  }

protected:
  const FirstSetType* m_set1;
  const SecondSetType* m_set2;
};

/*!
 * \brief Select the base class for RelationSet and ProductSet.
 *
 *  - In the ConcreteInterface case, the base class ConcreteBivariateSet
 *    is used to avoid virtual function calls in the derived instance.
 *  - In the VirtualInterface case, the base class is BivariateSet, which
 *    contains a common virtual interface between Relation/ProductSets.
 *
 * \tparam InterfacePolicy the interface type to select (virtual or concrete)
 *
 * \see InterfacePolicies.hpp
 */
template <typename InterfacePolicy, typename Set1, typename Set2, typename Position>
struct BSetInterfaceSelector;

template <typename InterfacePolicy, typename Set1, typename Set2, typename Position>
struct BSetInterfaceSelector
{
  static_assert(std::is_same<InterfacePolicy, ConcreteInterface>::value,
                "InterfacePolicy must be one of policies::ConcreteInterface or "
                "policies::VirtualInterface.");
  using Type = ConcreteBivariateSet<Set1, Set2, Position>;
};

template <typename Set1, typename Set2, typename Position>
struct BSetInterfaceSelector<VirtualInterface, Set1, Set2, Position>
{
  using Type = BivariateSet<Set1, Set2, Position>;
};
}  // namespace detail

template <typename InterfacePolicy,
          typename FromSet,
          typename ToSet,
          typename Position = ::axom::slam::detail::
            default_flat_position_t<typename FromSet::PositionType, typename ToSet::PositionType>>
using BivariateSetInterface =
  typename detail::BSetInterfaceSelector<InterfacePolicy, FromSet, ToSet, Position>::Type;

}  // end namespace axom::slam::policies
