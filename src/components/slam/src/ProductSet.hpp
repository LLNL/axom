/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-741217
 *
 * All rights reserved.
 *
 * This file is part of Axom.
 *
 * For details about use and distribution, please read axom/LICENSE.
 *
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * \file ProductSet.hpp
 *
 * \brief Basic API for a SLAM Cartesian product set
 *
 */

#ifndef SLAM_PRODUCT_SET_H_
#define SLAM_PRODUCT_SET_H_

#include "slam/BivariateSet.hpp"
#include "slam/RangeSet.hpp"

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
class ProductSet : public BivariateSet, RangeSet
{
  using SetType = Set;
public:
  using RangeSetType = RangeSet;
  using PositionType = BivariateSet::PositionType;
  using ElementType = BivariateSet::ElementType;
  using OrderedSetType = BivariateSet::OrderedSetType;

  /** \brief Default constructor */
  ProductSet() {}

  /**
   * \brief Constructor taking in pointers of two Sets.
   *
   * \param set1  Pointer to the first Set.
   * \param set2  Pointer to the second Set.
   */

  ProductSet(Set* set1, Set* set2) :
    BivariateSet(set1,set2), RangeSet(set1->size()*set2->size())
  {
    //fill in the row data now for getElements(i) function,
    //since every row is the same, a call to getElements() returns the same set.
    PositionType size2 = secondSetSize();
    m_rowSet_data.resize(size2);
    for (int s2 = 0 ; s2 < size2 ; s2++)
    {
      m_rowSet_data[s2] = s2;
    }
    m_rowSet = OrderedSetType::SetBuilder()
               .size(size2)
               .offset(0)
               .data(&m_rowSet_data);
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
  PositionType findElementIndex(PositionType pos1,
                                PositionType pos2) const override
  {
    isValidIndex(pos1, pos2);
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
    isValidIndex(pos1, pos2);
    PositionType size2 = secondSetSize();
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
  const OrderedSetType getElements(PositionType pos1) const override
  {
    SLIC_ASSERT(pos1 >= 0 && pos1 < firstSetSize());

    return m_rowSet;
  }

  ElementType at(PositionType pos) const override {
    return pos % firstSetSize();
  }

  PositionType size() const override
  {
    return firstSetSize()*secondSetSize();
  }

  PositionType size(PositionType) const override
  {
    return secondSetSize();
  }

  bool isValidIndex(PositionType s1, PositionType s2) const
  {
    PositionType size1 = firstSetSize();
    PositionType size2 = secondSetSize();
    return s1 >= 0 && s1 < size1 && s2 >= 0 && s2 < size2;
  }

  bool isValid(bool verboseOutput = false) const override
  {
    return BivariateSet::isValid(verboseOutput) &&
           RangeSet::isValid(verboseOutput);
  }

private:
  /** \brief verify the FlatIndex \a pos is within the valid range. */
  void verifyPosition(PositionType pos) const override
  { //from RangeSet, overloading to avoid warning in compiler
    SLIC_ASSERT_MSG(
      pos >= 0 && pos < size(),
      "SLAM::ProductSet -- requested out-of-range element at position "
      << pos << ", but set only has " << size() << " elements.");
  }

  /** \brief verify the SparseIndex (which is the same as its DenseIndex) is
   *         within the valid range. */
  void verifyPosition(PositionType pos1, PositionType pos2) const override
  {
    SLIC_ASSERT_MSG(
      isValidIndex(pos1,pos2),
      "SLAM::ProductSet -- requested out-of-range element at position ("
      << pos1 << "," << pos2 << "), but set only has "
      << firstSetSize() << "x" << secondSetSize() << " elements.");
  }

private:

  std::vector<PositionType> m_rowSet_data;
  OrderedSetType m_rowSet;

};



} // end namespace slam
} // end namespace axom

#endif //  SLAM_PRODUCT_SET_H
