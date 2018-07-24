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
 * \file BivariateSet.hpp
 *
 * \brief Contians the class BivariateSet and NullBivariateSet
 *
 */

#ifndef SLAM_BIVARIATE_SET_H_
#define SLAM_BIVARIATE_SET_H_

#include "slam/OrderedSet.hpp"
#include "slam/NullSet.hpp"

#include <cassert>

namespace axom
{
namespace slam
{

/**
 * \class BivariateSet
 *
 * \brief Abstract class that models a set whose elements are indexed by 2 indices
 *        (e.g. a matrix).
 * 
 * \detail A BivariateSet is like a Set, but uses 2 indices instead of 1. Its 
 *         constructor takes 2 Set pointers, and the size of the 2 Sets defines 
 *         the range of the 2 indices of the BivariateSet. 
 *         
 *  The class assumes that, for a sparse matrix representation, the second
 *  index is sparcified in the internal representation. The sparcified 
 *  second index is referred to as "SparseIndex" in the documentation, and the
 *  original dense index as "DenseIndex". To treat the elements as a contiguous 
 *  array and access via what's called the "FlatIndex".\n
 *  
 *  For example, a 2 x 4 sparse matrix below:
 *     \code
 *        0  1  2  3
 *     0  a  b
 *     1     c     d
 *     \endcode
 *   
 *   Access the elements using DenseIndex `(i,j)` would be...\n
 *   `(i = 1, j = 1) = c`\n
 *   `(i = 1, j = 3) = d`\n
 *  
 *   Using SparseIndex `(i,k)`...\n
 *   `(i = 1, k = 0) = c`\n
 *   `(i = 1, k = 1) = d`\n
 *  
 *   Using FlatIndex `[idx]`...\n
 *   `[idx = 0] = a`\n
 *   `[idx = 1] = b`\n
 *   `[idx = 2] = c`\n
 *   `[idx = 3] = d`\n
 *
 * 
 */

class BivariateSet
{
private:
  static const NullSet s_nullSet;

public:
  using PositionType = slam::Set::PositionType;
  using ElementType = slam::Set::ElementType;
  using OrderedSetType = OrderedSet<
      policies::RuntimeSize<Set::PositionType>,
      policies::RuntimeOffset<Set::PositionType>,
      policies::StrideOne<Set::PositionType>,
      policies::STLVectorIndirection<PositionType, PositionType>>;
  
  static const PositionType INVALID_POS = PositionType(-1);

public:

  /**
   * \brief Constructor taking pointers to the 2 sets that defines the range of the
   *        indices of the BivariateSet.
   *
   * \param set1  Pointer to the first Set.
   * \param set2  Pointer to the second Set.
   */
  BivariateSet(const Set* set1 = &s_nullSet, const Set* set2 = &s_nullSet) 
    : m_set1(set1), m_set2(set2)
  { }

  /**
   * \brief Searches for the SparseIndex of the element given its DenseIndex. 
   * \detail If the element (i,j) is the k<sup>th</sup> non-zero in the row, then 
   *         `findElementIndex(i,j)` returns `k`. If `element (i,j)` does not 
   *         exist (such as the case of a zero in a sparse matrix), then 
   *         `INVALID_POS` is returned.
   *
   * \param pos1  The first set position.
   * \param pos2  The second set position.
   * \return  The DenseIndex of the given element, or INVALID_POS if such element
   *          is missing from the set.
   * \pre   0 <= pos1 <= set1.size() && 0 <= pos2 <= size2.size()
   */
  virtual PositionType findElementIndex(PositionType pos1, 
                                        PositionType pos2) const = 0;

  /**
   * \brief Search for the FlatIndex of the element given its DenseIndex.
   *
   * \param pos1  The first set position.
   * \param pos2  The second set position.
   *
   * \return  The element's FlatIndex
   * \pre   0 <= pos1 <= set1.size() && 0 <= pos2 <= size2.size()
   */
  virtual PositionType findElementFlatIndex(PositionType pos1, 
                                            PositionType pos2) const = 0;

  /**
   * \brief Searches for the first existing element given the row index (first
   *        set position).
   *
   * \param pos1  The first set position.
   *
   * \return  The found element's FlatIndex.
   * \pre   0 <= pos1 <= set1.size()
   */
  virtual PositionType findElementFlatIndex(PositionType pos1) const = 0;

  /**
   * \brief Size of the BivariateSet, which is the number of non-zero entries
   *        in the BivariateSet.
   */
  virtual PositionType size() const = 0;

  /**
   * \brief Size of a row given the first set position. 
   * 
   * \pre  0 <= pos1 <= set1.size()
   */
  virtual PositionType size(PositionType pos1) const = 0; //size of a row

  /** \brief Size of the first set.   */
  PositionType firstSetSize() const { return m_set1 ? m_set1->size() : 0; }
  /** \brief Size of the second set.   */
  PositionType secondSetSize() const { return m_set2 ? m_set2->size() : 0; }

  /** \brief Returns pointer to the first set.   */
  virtual const Set* getFirstSet() const { return m_set1; }
  /** \brief Returns pointer to the second set.   */
  virtual const Set* getSecondSet() const { return m_set2; }

  /** \brief Value of the Set element given the FlatIndex of the element. */
  virtual ElementType at(PositionType pos) const = 0;

  /**
   * \brief A set of elements with the given first set index.
   *
   * \param s1  The first set index.
   * \return  An OrderedSet containing the elements in the row.
   * \pre  0 <= pos1 <= set1.size()
   */
  virtual const OrderedSetType getRow(PositionType s1) const = 0;

  virtual bool isValid(bool verboseOutput = false) const
  {

    if (m_set1 == nullptr || m_set2 == nullptr)
    {
      if (verboseOutput)
      {
        std::cout << "\n*** BivariateSet is not valid:\n"
          << "\t* Set pointers should not be null.\n"
          << std::endl;
      }
      return false;
    }
    return true;
  }

  virtual void verifyPosition(PositionType s1, PositionType s2) const = 0;
  
protected:
  const Set* m_set1;
  const Set* m_set2;
};




/**
* \class NullBivariateSet
*
* \brief A Null Bivariate class. Same as tehe NullSet for Set class. 
*/
class NullBivariateSet : public BivariateSet
{
public:
  NullBivariateSet() {}

  PositionType findElementIndex(PositionType pos1, PositionType pos2 = 0) const override
  {
    verifyPosition(pos1, pos2);
    return PositionType();
  }

  PositionType findElementFlatIndex(PositionType s1, PositionType s2) const override
  {
    verifyPosition(s1, s2);
    return PositionType();
  }

  PositionType findElementFlatIndex(PositionType s1) const override
  {
    return findElementFlatIndex(s1, 0);
  }

  ElementType at(PositionType) const override
  {
    return PositionType();
  }

  PositionType size() const override
  {
    return PositionType();
  }

  PositionType size(PositionType) const override
  {
    return PositionType();
  }

  const OrderedSetType getRow(PositionType) const override
  {
    return OrderedSetType::SetBuilder();
  }

private:
  void verifyPosition(PositionType AXOM_DEBUG_PARAM(pos1),
    PositionType AXOM_DEBUG_PARAM(pos2)) const override
  {
    SLIC_ASSERT_MSG(
      false,
      "Subscripting on NullSet is never valid."
      << "\n\tAttempted to access item at index " 
      << pos1 << "," << pos2 << ".");
  }
};


} // end namespace slam
} // end namespace axom

#endif //  SLAM_BIVARIATE_SET_H_
