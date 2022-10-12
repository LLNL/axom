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

template <typename Set1 = slam::Set<>, typename Set2 = slam::Set<>>
struct ConcreteBivariateSet
{
  using FirstSetType = Set1;
  using SecondSetType = Set2;

  using PositionType = typename FirstSetType::PositionType;
  using ElementType = typename FirstSetType::ElementType;

  using SubsetType = void;

  static constexpr PositionType INVALID_POS = PositionType(-1);

  /**
   * \brief Constructor taking pointers to the two sets that defines the range
   *        of the indices of the BivariateSet.
   *
   * \param set1  Pointer to the first Set.
   * \param set2  Pointer to the second Set.
   */
  ConcreteBivariateSet(
    const Set1* set1 = policies::EmptySetTraits<Set1>::emptySet(),
    const Set2* set2 = policies::EmptySetTraits<Set2>::emptySet())
    : m_set1(set1)
    , m_set2(set2)
  { }

  /** \brief Size of the first set.   */
  inline PositionType firstSetSize() const
  {
    return getSize<FirstSetType>(m_set1);
  }
  /** \brief Size of the second set.   */
  inline PositionType secondSetSize() const
  {
    return getSize<SecondSetType>(m_set2);
  }

  /** \brief Returns pointer to the first set.   */
  const FirstSetType* getFirstSet() const { return m_set1; }
  /** \brief Returns pointer to the second set.   */
  const SecondSetType* getSecondSet() const { return m_set2; }

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
    SLIC_ASSERT_MSG(s != nullptr, "nullptr in BivariateSet::getSize()");
    return static_cast<SetType>(*s).size();
  }

protected:
  const FirstSetType* m_set1;
  const SecondSetType* m_set2;
};

}  // end namespace policies
}  // end namespace slam
}  // end namespace axom

#endif  // SLAM_BivarSetIfacePolicies_HPP
