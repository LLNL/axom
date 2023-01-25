// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
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

#include "axom/slic.hpp"

#include "axom/slam/Set.hpp"
#include "axom/slam/OrderedSet.hpp"
#include "axom/slam/NullSet.hpp"
#include "axom/slam/RangeSet.hpp"

#include <cassert>
#include <type_traits>

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
 *         0  1  2  3
 *         _  _  _  _
 *     0 | a     b
 *     1 |    c     d
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

template <typename Set1 = slam::Set<>, typename Set2 = slam::Set<>>
class BivariateSet
{
public:
  using FirstSetType = Set1;
  using SecondSetType = Set2;

  using PositionType = typename FirstSetType::PositionType;
  using ElementType = typename FirstSetType::ElementType;
  using NullSetType = NullSet<PositionType, ElementType>;

  using OrderedSetType =
    OrderedSet<PositionType,
               ElementType,
               policies::RuntimeSize<PositionType>,
               policies::RuntimeOffset<PositionType>,
               policies::StrideOne<PositionType>,
               policies::STLVectorIndirection<PositionType, ElementType>>;

  using RangeSetType = RangeSet<PositionType, ElementType>;

public:
  static const PositionType INVALID_POS = PositionType(-1);
  static const NullSetType s_nullSet;

public:
  /**
   * \brief Constructor taking pointers to the two sets that defines the range
   *        of the indices of the BivariateSet.
   *
   * \param set1  Pointer to the first Set.
   * \param set2  Pointer to the second Set.
   */
  BivariateSet(const Set1* set1 = policies::EmptySetTraits<Set1>::emptySet(),
               const Set2* set2 = policies::EmptySetTraits<Set2>::emptySet())
    : m_set1(set1)
    , m_set2(set2)
  { }

  /**
   * \brief Default virtual destructor
   *
   * \note BivariateSet does not own the two underlying sets
   */
  virtual ~BivariateSet() = default;

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
   * \brief Finds the range of indices of valid elements in the second set,
   *        given the index of an element in the first set.
   * \param Position of the element in the first set
   *
   * \return A range set of the positions in the second set
   */
  virtual RangeSetType elementRangeSet(PositionType pos1) const = 0;
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

  virtual bool isValid(bool verboseOutput = false) const;

  virtual void verifyPosition(PositionType s1, PositionType s2) const = 0;

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

template <typename Set1, typename Set2>
const typename BivariateSet<Set1, Set2>::NullSetType BivariateSet<Set1, Set2>::s_nullSet;

template <typename Set1, typename Set2>
bool BivariateSet<Set1, Set2>::isValid(bool verboseOutput) const
{
  if(m_set1 == nullptr || m_set2 == nullptr)
  {
    if(verboseOutput)
    {
      SLIC_INFO("BivariateSet is not valid: "
                << " Set pointers should not be null.");
    }
    return false;
  }
  return m_set1->isValid(verboseOutput) && m_set2->isValid(verboseOutput);
}

/**
 * \class NullBivariateSet
 *
 * \brief A Null BivariateSet class. Same as the NullSet for Set class.
 */
template <typename SetType1 = slam::Set<>, typename SetType2 = slam::Set<>>
class NullBivariateSet : public BivariateSet<SetType1, SetType2>
{
public:
  using FirstSetType = SetType1;
  using SecondSetType = SetType2;
  using BSet = BivariateSet<FirstSetType, SecondSetType>;
  using PositionType = typename BSet::PositionType;
  using ElementType = typename BSet::ElementType;
  using OrderedSetType = typename BSet::OrderedSetType;
  using RangeSetType = typename BSet::RangeSetType;

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

  RangeSetType elementRangeSet(PositionType) const override
  {
    return RangeSetType();
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

}  // end namespace slam
}  // end namespace axom

#endif  //  SLAM_BIVARIATE_SET_H_
