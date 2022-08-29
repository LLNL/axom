// Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file ProductSet.hpp
 *
 * \brief Basic API for a SLAM Cartesian product set
 *
 */

#ifndef SLAM_PRODUCT_SET_H_
#define SLAM_PRODUCT_SET_H_

#include "axom/slam/BivariateSet.hpp"
#include "axom/slam/RangeSet.hpp"

#include <cassert>

namespace axom
{
namespace slam
{
/**
 * \class ProductSet
 *
 * \brief Models a set whose element is the Cartesian product of two sets. The
 *        number of elements in this set is the product of the sizes of the two
 *        input sets.
 *
 *        Users should refer to the BivariateSet documentation for descriptions
 *        of the different indexing names (SparseIndex, DenseIndex, FlatIndex).
 *
 * \see   BivariateSet
 */
template <typename SetType1 = slam::Set<>, typename SetType2 = slam::Set<>>
class ProductSet final
  : public BivariateSetBase<SetType1, SetType2, ProductSet<SetType1, SetType2>>,
    RangeSet<typename SetType1::PositionType, typename SetType1::ElementType>
{
public:
  using RangeSetType =
    RangeSet<typename SetType1::PositionType, typename SetType1::ElementType>;
  using FirstSetType = SetType1;
  using SecondSetType = SetType2;
  using BivariateSetType = BivariateSet<FirstSetType, SecondSetType>;
  using PositionType = typename BivariateSetType::PositionType;
  using ElementType = typename BivariateSetType::ElementType;
  using ProductSetType = ProductSet<SetType1, SetType2>;
  using RelationSubset =
    RangeSet<typename SetType1::PositionType, typename SetType2::PositionType>;

  using BaseClass = BivariateSetBase<SetType1, SetType2, ProductSet>;

  /** \brief Default constructor */
  ProductSet()
    : m_firstSet(policies::EmptySetTraits<FirstSetType>::emptySet())
    , m_secondSet(policies::EmptySetTraits<SecondSetType>::emptySet())
  { }

  /**
   * \brief Constructor taking in pointers of two Sets.
   *
   * \param set1  Pointer to the first Set.
   * \param set2  Pointer to the second Set.
   */
  template <typename UFirstSet, typename USecondSet>
  ProductSet(const UFirstSet* set1, const USecondSet* set2)
    : RangeSetType(set1->size() * set2->size())
    , m_firstSet(set1)
    , m_secondSet(set2)
  {
    static_assert(!std::is_abstract<UFirstSet>::value,
                  "Must pass in a concrete set instance.");
    static_assert(!std::is_abstract<USecondSet>::value,
                  "Must pass in a concrete set instance.");
  }

  /**
   * \brief Return the element SparseIndex. Since ProductSet is the full
   *        Cartesian product of the two sets, the SparseIndex is the same as
   *        its DenseIndex, which is the same as the \a pos2 parameter.
   *
   * \param pos1  The first set position.
   * \param pos2  The second set position.
   *
   * \return  The element's SparseIndex, which is the same as pos2
   */
  PositionType findElementIndex(PositionType pos1, PositionType pos2) const
  {
    verifyPositionImpl(pos1, pos2);
    return pos2;
  }

  /**
   * \brief Returns an element's FlatIndex given its DenseIndex. Since
   *        ProductSet is the full Cartesian product of the two sets, an
   *        element's FlatIndex is equal to `pos1*secondSetSize()+pos2`.
   *
   * \param pos1  The first set position.
   * \param pos2  The second set position.
   *
   * \return  The element's FlatIndex.
   */
  PositionType findElementFlatIndex(PositionType pos1, PositionType pos2) const
  {
    verifyPositionImpl(pos1, pos2);
    PositionType size2 = this->secondSetSize();
    PositionType pos = size2 * pos1 + pos2;

    return pos;
  }

  /**
   * \brief Returns the FlatIndex of the first element in the specified row.
   *        This is equal to `pos1*secondSetSize()`.
   *
   * \param pos1  The first set position that specifies the row.
   */
  PositionType findElementFlatIndex(PositionType pos1) const
  {
    return findElementFlatIndex(pos1, 0);
  }

  /**
   * \brief Return all elements from the second set associated with position
   *        \a pos1 in the first set.
   *
   * \param pos1   The first set position that specifies the row.
   *
   * \return  An OrderedSet of the elements in the row.
   */
  RelationSubset getElements(PositionType AXOM_DEBUG_PARAM(pos1)) const
  {
    SLIC_ASSERT(pos1 >= 0 && pos1 < this->firstSetSize());

    return RelationSubset(this->secondSetSize());
  }

  ElementType at(PositionType pos) const { return pos % this->secondSetSize(); }

  PositionType size() const
  {
    return this->firstSetSize() * this->secondSetSize();
  }

  PositionType size(PositionType) const { return this->secondSetSize(); }

  /** \brief Returns pointer to the first set.   */
  const FirstSetType* getFirstSet() const { return m_firstSet; }
  /** \brief Returns pointer to the second set.   */
  const SecondSetType* getSecondSet() const { return m_secondSet; }

  RangeSetType elementRangeSet(PositionType pos1) const
  {
    const auto sz = this->secondSetSize();
    return typename RangeSetType::SetBuilder().size(sz).offset(sz * pos1);
  }

  bool isValidIndex(PositionType s1, PositionType s2) const
  {
    PositionType size1 = this->firstSetSize();
    PositionType size2 = this->secondSetSize();
    return s1 >= 0 && s1 < size1 && s2 >= 0 && s2 < size2;
  }

  bool isValid(bool verboseOutput = false) const
  {
    return BaseClass::isValid(verboseOutput) &&
      RangeSetType::isValid(verboseOutput);
  }

  /** \brief verify the FlatIndex \a pos is within the valid range. */
  void verifyPosition(PositionType pos) const
  {  //from RangeSet, overloading to avoid warning in compiler
    verifyPositionImpl(pos);
  }

  /** \brief verify the SparseIndex (which is the same as its DenseIndex) is
   *         within the valid range. */
  void verifyPosition(PositionType pos1, PositionType pos2) const
  {
    verifyPositionImpl(pos1, pos2);
  }

private:
  /** \brief implementation for verifyPosition */
  inline void verifyPositionImpl(PositionType AXOM_DEBUG_PARAM(pos)) const
  {  //from RangeSet, overloading to avoid warning in compiler
    SLIC_ASSERT_MSG(
      pos >= 0 && pos < size(),
      "SLAM::ProductSet -- requested out-of-range element at position "
        << pos << ", but set only has " << size() << " elements.");
  }

  /** \brief verify the SparseIndex (which is the same as its DenseIndex) is
   *         within the valid range. */
  inline void verifyPositionImpl(PositionType AXOM_DEBUG_PARAM(pos1),
                                 PositionType AXOM_DEBUG_PARAM(pos2)) const
  {
    SLIC_ASSERT_MSG(
      isValidIndex(pos1, pos2),
      "SLAM::ProductSet -- requested out-of-range element at position ("
        << pos1 << "," << pos2 << "), but set only has " << this->firstSetSize()
        << "x" << this->secondSetSize() << " elements.");
  }

private:
  SetContainer<FirstSetType> m_firstSet;
  SetContainer<SecondSetType> m_secondSet;
};

}  // end namespace slam
}  // end namespace axom

#endif  //  SLAM_PRODUCT_SET_H
