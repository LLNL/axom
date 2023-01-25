// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
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
  : public BivariateSet<SetType1, SetType2>,
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
  using OrderedSetType = typename BivariateSetType::OrderedSetType;
  using ProductSetType = ProductSet<SetType1, SetType2>;

  /** \brief Default constructor */
  ProductSet() { }

  /**
   * \brief Constructor taking in pointers of two Sets.
   *
   * \param set1  Pointer to the first Set.
   * \param set2  Pointer to the second Set.
   */

  ProductSet(const FirstSetType* set1, const SecondSetType* set2)
    : BivariateSetType(set1, set2)
    , RangeSetType(set1->size() * set2->size())
  {
    //fill in the row data now for getElements(i) function,
    //since every row is the same, a call to getElements() returns the same set.
    //
    // HACK -- this should actually be returning a PositionSet since it always
    //         goes from 0 to secondSetSize()
    // This requires a change to the return type of BivariateSet::getElements()
    auto size2 = this->secondSetSize();
    m_rowSet_data.resize(size2);
    for(int s2 = 0; s2 < size2; ++s2)
    {
      m_rowSet_data[s2] = s2;
    }

    m_rowSet = typename OrderedSetType::SetBuilder().size(size2).offset(0).data(
      &m_rowSet_data);
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
  PositionType findElementIndex(PositionType pos1, PositionType pos2) const override
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
  PositionType findElementFlatIndex(PositionType pos1,
                                    PositionType pos2) const override
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
  PositionType findElementFlatIndex(PositionType pos1) const override
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
  const OrderedSetType getElements(PositionType AXOM_DEBUG_PARAM(pos1)) const override
  {
    SLIC_ASSERT(pos1 >= 0 && pos1 < this->firstSetSize());

    return m_rowSet;
  }

  ElementType at(PositionType pos) const override
  {
    return pos % this->secondSetSize();
  }

  PositionType size() const override
  {
    return this->firstSetSize() * this->secondSetSize();
  }

  PositionType size(PositionType) const override
  {
    return this->secondSetSize();
  }

  RangeSetType elementRangeSet(PositionType pos1) const override
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

  bool isValid(bool verboseOutput = false) const override
  {
    return BivariateSetType::isValid(verboseOutput) &&
      RangeSetType::isValid(verboseOutput);
  }

private:
  /** \brief verify the FlatIndex \a pos is within the valid range. */
  void verifyPosition(PositionType pos) const override
  {  //from RangeSet, overloading to avoid warning in compiler
    verifyPositionImpl(pos);
  }

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
  void verifyPosition(PositionType pos1, PositionType pos2) const override
  {
    verifyPositionImpl(pos1, pos2);
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
  std::vector<PositionType> m_rowSet_data;
  OrderedSetType m_rowSet;
};

}  // end namespace slam
}  // end namespace axom

#endif  //  SLAM_PRODUCT_SET_H
