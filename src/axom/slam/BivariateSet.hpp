// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file BivariateSet.hpp
 *
 * \brief Contains the class BivariateSet and NullBivariateSet
 *
 */

#ifndef SLAM_BIVARIATE_SET_H_
#define SLAM_BIVARIATE_SET_H_

#include "axom/slam/OrderedSet.hpp"
#include "axom/slam/NullSet.hpp"

#include <cassert>

namespace axom
{
namespace slam
{
/**
 * \class BivariateSet
 *
 * \brief Abstract class that models a set whose elements are indexed by two
 *        indices. Each element in a BivariateSet is equivalent to an ordered
 *        pair containing a row and column index, similar to indexing in a
 *        matrix.
 *
 * \detail BivariateSet models a subset of the Cartesian product of its two
 *         sets. Elements of a BivariateSet can be represented as an ordered
 *         pair of indices into the two sets.
 *
 *  For BivariateSets that do not model the entire Cartesian product, indices
 *  can be relative to the element positions in the original sets (in which
 *  case, we refer to them as a "DenseIndex"), or relative to the number of
 *  encoded indices, in which case we refer to them as a "SparseIndex".
 *  If we consider all the elements of a BivariateSet, we refer to this index
 *  space as the "FlatIndex". \n
 *
 *  For example, a 2 x 4 sparse matrix below:
 *     \code
 *        0  1  2  3
 *     0  a     b
 *     1     c     d
 *     \endcode
 *
 *   Access the elements using DenseIndex `(i,j)` would be...\n
 *   `(i = 0, j = 0) = a`\n
 *   `(i = 0, j = 2) = b`\n
 *   `(i = 1, j = 1) = c`\n
 *   `(i = 1, j = 3) = d`\n
 *
 *   Using SparseIndex `(i,k)`...\n
 *   `(i = 0, k = 0) = a`\n
 *   `(i = 0, k = 1) = b`\n
 *   `(i = 1, k = 0) = c`\n
 *   `(i = 1, k = 1) = d`\n
 *
 *   Using FlatIndex `[idx]`...\n
 *   `[idx = 0] = a`\n
 *   `[idx = 1] = b`\n
 *   `[idx = 2] = c`\n
 *   `[idx = 3] = d`\n
 *
 */

template <typename PosType = slam::DefaultPositionType,
          typename ElemType = slam::DefaultElementType>
class BivariateSet
{
public:
  using PositionType = PosType;
  using ElementType = ElemType;
  using SetType = Set<PositionType, ElementType>;
  using OrderedSetType =
    OrderedSet<PositionType,
               ElementType,
               policies::RuntimeSize<PositionType>,
               policies::RuntimeOffset<PositionType>,
               policies::StrideOne<PositionType>,
               policies::STLVectorIndirection<PositionType, ElementType>>;

  static const PositionType INVALID_POS = PositionType(-1);
  static const NullSet<PosType, ElemType> s_nullSet;

public:
  /**
   * \brief Constructor taking pointers to the two sets that defines the range
   *        of the indices of the BivariateSet.
   *
   * \param set1  Pointer to the first Set.
   * \param set2  Pointer to the second Set.
   */
  BivariateSet(const SetType* set1 = &s_nullSet, const SetType* set2 = &s_nullSet)
    : m_set1(set1)
    , m_set2(set2)
  { }

  /**
   * \brief Searches for the SparseIndex of the element given its DenseIndex.
   * \detail If the element (i,j) is the k<sup>th</sup> non-zero in the row,
   *         then `findElementIndex(i,j)` returns `k`. If `element (i,j)` does
   *         not exist (such as the case of a zero in a sparse matrix), then
   *         `INVALID_POS` is returned.
   *
   * \param pos1  The first set position.
   * \param pos2  The second set position.
   * \return  The DenseIndex of the given element, or INVALID_POS if such
   *          element is missing from the set.
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
   * \brief Number of elements of the BivariateSet whose first index is \a pos
   *
   * \pre  0 <= pos1 <= set1.size()
   */
  virtual PositionType size(PositionType pos1) const = 0;  //size of a row

  /** \brief Size of the first set.   */
  PositionType firstSetSize() const { return m_set1 ? m_set1->size() : 0; }
  /** \brief Size of the second set.   */
  PositionType secondSetSize() const { return m_set2 ? m_set2->size() : 0; }

  /** \brief Returns pointer to the first set.   */
  virtual const SetType* getFirstSet() const { return m_set1; }
  /** \brief Returns pointer to the second set.   */
  virtual const SetType* getSecondSet() const { return m_set2; }

  /** \brief Returns the element at the given FlatIndex \a pos */
  virtual ElementType at(PositionType pos) const = 0;

  /**
   * \brief A set of elements with the given first set index.
   *
   * \param s1  The first set index.
   * \return  An OrderedSet containing the elements
   * \pre  0 <= pos1 <= set1.size()
   */
  virtual const OrderedSetType getElements(PositionType s1) const = 0;

  virtual bool isValid(bool verboseOutput = false) const
  {
    if(m_set1 == nullptr || m_set2 == nullptr)
    {
      if(verboseOutput)
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
  const SetType* m_set1;
  const SetType* m_set2;
};

/**
 * \class NullBivariateSet
 *
 * \brief A Null BivariateSet class. Same as the NullSet for Set class.
 */
template <typename PosType = slam::DefaultPositionType,
          typename ElemType = slam::DefaultElementType>
class NullBivariateSet : public BivariateSet<PosType, ElemType>
{
public:
  using PositionType = PosType;
  using ElementType = ElemType;
  using OrderedSetType = typename BivariateSet<PosType, ElemType>::OrderedSetType;

public:
  NullBivariateSet() = default;

  PositionType findElementIndex(PositionType pos1,
                                PositionType pos2 = 0) const override
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

  ElementType at(PositionType) const override { return PositionType(); }

  PositionType size() const override { return PositionType(); }

  PositionType size(PositionType) const override { return PositionType(); }

  const OrderedSetType getElements(PositionType) const override
  {
    using OrderedSetBuilder = typename OrderedSetType::SetBuilder;
    return OrderedSetBuilder();
  }

private:
  void verifyPosition(PositionType AXOM_DEBUG_PARAM(pos1),
                      PositionType AXOM_DEBUG_PARAM(pos2)) const override
  {
    SLIC_ASSERT_MSG(false,
                    "Subscripting on NullSet is never valid."
                      << "\n\tAttempted to access item at index " << pos1 << ","
                      << pos2 << ".");
  }
};

template <typename PosType, typename ElemType>
const NullSet<PosType, ElemType> BivariateSet<PosType, ElemType>::s_nullSet;

}  // end namespace slam
}  // end namespace axom

#endif  //  SLAM_BIVARIATE_SET_H_
