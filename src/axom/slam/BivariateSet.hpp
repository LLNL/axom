// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file BivariateSet.hpp
 *
 * \brief Contains the class BivariateSet and NullBivariateSet
 */

#pragma once

#include "axom/slic.hpp"

#include "axom/slam/Set.hpp"
#include "axom/slam/OrderedSet.hpp"
#include "axom/slam/NullSet.hpp"
#include "axom/slam/RangeSet.hpp"
#include "axom/slam/policies/PolicyTraits.hpp"

#include <cassert>
#include <optional>
#include <type_traits>
#include <utility>

namespace axom::slam
{
template <typename BivariateSetType>
struct BivariateSetIterator;

/**
 * \class BivariateSet
 *
 * \brief Abstract interface for a set of pairs of first- and second-set positions.
 *
 * A BivariateSet represents some or all pairs in the Cartesian product of two
 * sets. For example, a zone-node set contains a pair for each node of each zone.
 * getElements(i) returns the second-set positions associated with first-set position i.
 *
 * The API distinguishes three ways to identify an entry:
 * - DenseIndex uses positions in the two original sets, written as (i, j).
 * - SparseIndex uses i and the local position k within getElements(i).
 * - FlatIndex counts entries across all first-set positions in order.
 *
 * Suppose getElements(0) contains {0, 2} and getElements(1) contains {1, 3}.
 * The final entry has DenseIndex (1, 3), SparseIndex (1, 1), and FlatIndex 3.
 * at(3) returns the coordinate pair (1, 3), not an element stored in either
 * original set. Use getFirstSet()->at(1) and getSecondSet()->at(3) to access
 * those elements.
 *
 * \note The three position types are independent.
 *       \a FirstPositionType and \a SecondPositionType come from each of the sets,
 *       while the flat \a PositionType defaults to a type able to represent both.
 *       A set pair with 32-bit zones and 64-bit nodes therefore keeps each position at its
 *       own width while flattening into the wider one.
 */

template <typename Set1 = slam::Set<>,
          typename Set2 = slam::Set<>,
          typename Position =
            detail::default_flat_position_t<typename Set1::PositionType, typename Set2::PositionType>>
class BivariateSet
{
public:
  using FirstSetType = Set1;
  using SecondSetType = Set2;

  using FirstPositionType = typename FirstSetType::PositionType;
  using SecondPositionType = typename SecondSetType::PositionType;

  // PositionType indexes the flattened bivariate set.
  // Its elements are pairs of positions in the first and second sets.
  using PositionType = Position;
  using ElementType = std::pair<FirstPositionType, SecondPositionType>;
  using NullSetType = NullSet<PositionType, ElementType>;

  // A subset uses local positions to access positions in the second set.
  using SubsetType = OrderedSet<PositionType,
                                SecondPositionType,
                                policies::RuntimeSize<PositionType>,
                                policies::RuntimeOffset<PositionType>,
                                policies::StrideOne<PositionType>,
                                policies::ArrayViewIndirection<PositionType, SecondPositionType>>;

  // elementRangeSet() describes flat positions, not values in either constituent set.
  using RangeSetType = RangeSet<PositionType, PositionType>;
  using IteratorType = BivariateSetIterator<BivariateSet>;

public:
  static constexpr PositionType INVALID_POS = PositionType(-1);
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
   * If getElements(i) contains j at local position k, findElementIndex(i, j)
   * returns k. It returns INVALID_POS when that pair is absent.
   *
   * \param pos1  The first set position.
   * \param pos2  The second set position.
   * \return The local position in getElements(pos1), or INVALID_POS if absent.
   * \pre 0 <= pos1 < firstSetSize() && 0 <= pos2 < secondSetSize()
   */
  virtual PositionType findElementIndex(FirstPositionType pos1, SecondPositionType pos2) const = 0;

  /**
   * \brief Finds the SparseIndex of the element given its DenseIndex.
   *
   * \return An engaged `std::optional` containing the SparseIndex if the element exists,
   *         or an empty `std::optional` if the element does not exist.
   *
   * \note This is a convenience wrapper around `findElementIndex(...)` that avoids
   *       sentinel checks against `INVALID_POS`.
   */
  [[nodiscard]] std::optional<PositionType> findElementIndexOptional(FirstPositionType pos1,
                                                                     SecondPositionType pos2) const
  {
    const auto idx = findElementIndex(pos1, pos2);
    return idx != INVALID_POS ? std::optional<PositionType>(idx) : std::optional<PositionType> {};
  }

  /**
   * \brief Search for the FlatIndex of the element given its DenseIndex.
   *
   * \param pos1  The first set position.
   * \param pos2  The second set position.
   *
   * \return The element's FlatIndex, or INVALID_POS if the pair is absent.
   * \pre 0 <= pos1 < firstSetSize() && 0 <= pos2 < secondSetSize()
   */
  AXOM_HOST_DEVICE virtual PositionType findElementFlatIndex(FirstPositionType pos1,
                                                             SecondPositionType pos2) const = 0;

  /**
   * \brief Finds the FlatIndex of the element given its DenseIndex.
   *
   * \return An engaged `std::optional` containing the FlatIndex if the element exists,
   *         or an empty `std::optional` if the element does not exist.
   *
   * \note This is a convenience wrapper around `findElementFlatIndex(...)` that avoids
   *       sentinel checks against `INVALID_POS`.
   */
  [[nodiscard]] AXOM_HOST_DEVICE std::optional<PositionType> findElementFlatIndexOptional(
    FirstPositionType pos1,
    SecondPositionType pos2) const
  {
    const auto idx = findElementFlatIndex(pos1, pos2);
    return idx != INVALID_POS ? std::optional<PositionType>(idx) : std::optional<PositionType> {};
  }

  /**
   * \brief Find the first flat position associated with a first-set position.
   *
   * \param pos1  The first set position.
   *
   * \return The first FlatIndex, or INVALID_POS if the associated subset is empty.
   * \pre 0 <= pos1 < firstSetSize()
   */
  virtual PositionType findElementFlatIndex(FirstPositionType pos1) const = 0;

  /*!
   * \brief Find the first flat position associated with pos1, if one exists.
   *
   * \return The first FlatIndex, or an empty optional if the subset is empty.
   *
   * \note This is a convenience wrapper around `findElementFlatIndex(pos1)` that avoids
   *       sentinel checks against `INVALID_POS`.
   */
  [[nodiscard]] std::optional<PositionType> findElementFlatIndexOptional(FirstPositionType pos1) const
  {
    const auto idx = findElementFlatIndex(pos1);
    return idx != INVALID_POS ? std::optional<PositionType>(idx) : std::optional<PositionType> {};
  }

  /**
   * \brief Return the first-set position at the given flat position.
   *
   * \param flatIndex Position in the bivariate set.
   *
   * \pre 0 <= flatIndex < size()
   */
  AXOM_HOST_DEVICE virtual FirstPositionType flatToFirstIndex(PositionType flatIndex) const = 0;

  /**
   * \brief Return the second-set position at the given flat position.
   *
   * \param flatIndex Position in the bivariate set.
   *
   * \pre 0 <= flatIndex < size()
   */
  AXOM_HOST_DEVICE virtual SecondPositionType flatToSecondIndex(PositionType flatIndex) const = 0;

  /**
   * \brief Return the flat positions associated with a first-set position.
   * \param pos1 Position in the first set.
   * \return A RangeSet of flat positions, not second-set positions.
   */
  AXOM_HOST_DEVICE virtual RangeSetType elementRangeSet(FirstPositionType pos1) const = 0;

  /// \brief The number of coordinate pairs in the bivariate set.
  [[nodiscard]] AXOM_HOST_DEVICE virtual PositionType size() const = 0;

  /// \brief Checks if there are any elements in the set
  AXOM_SUPPRESS_HD_WARN
  [[nodiscard]] AXOM_HOST_DEVICE bool empty() const { return size() == PositionType {}; }

  /**
   * \brief Number of coordinate pairs whose first-set position is \a pos1.
   *
   * \pre 0 <= pos1 < firstSetSize()
   */
  virtual PositionType size(FirstPositionType pos1) const = 0;

  /// \brief Size of the first set.
  [[nodiscard]] AXOM_HOST_DEVICE inline FirstPositionType firstSetSize() const
  {
    return getSize<FirstSetType>(m_set1);
  }

  /// \brief Size of the second set.
  AXOM_SUPPRESS_HD_WARN
  [[nodiscard]] AXOM_HOST_DEVICE inline SecondPositionType secondSetSize() const
  {
    return getSize<SecondSetType>(m_set2);
  }

  /// \brief Returns pointer to the first set.
  const FirstSetType* getFirstSet() const { return m_set1; }

  /// \brief Returns pointer to the second set.
  const SecondSetType* getSecondSet() const { return m_set2; }

  /// \brief Returns the element at the given FlatIndex \a pos
  [[nodiscard]] AXOM_HOST_DEVICE virtual ElementType at(PositionType pos) const = 0;

  /**
   * \brief Return the second-set positions associated with s1.
   *
   * \param s1  The first set index.
   * \return An OrderedSet of second-set positions.
   * \pre 0 <= s1 < firstSetSize()
   */
  virtual SubsetType getElements(FirstPositionType s1) const = 0;

  /// \brief Return an iterator to the first pair of set positions.
  IteratorType begin() const { return IteratorType(this, 0); }

  /// \brief Return an iterator past the last pair of set positions.
  IteratorType end() const { return IteratorType(this, size()); }

  [[nodiscard]] virtual bool isValid(bool verboseOutput = false) const;

private:
  virtual void verifyPosition(FirstPositionType s1, SecondPositionType s2) const = 0;

  AXOM_SUPPRESS_HD_WARN
  template <typename SetType>
  AXOM_HOST_DEVICE typename SetType::PositionType getSize(const SetType* s) const
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

template <typename Set1, typename Set2, typename Position>
const typename BivariateSet<Set1, Set2, Position>::NullSetType BivariateSet<Set1, Set2, Position>::s_nullSet;

template <typename Set1, typename Set2, typename Position>
bool BivariateSet<Set1, Set2, Position>::isValid(bool verboseOutput) const
{
  if(m_set1 == nullptr || m_set2 == nullptr)
  {
    if(verboseOutput)
    {
      SLIC_INFO("BivariateSet is not valid: " << " Set pointers should not be null.");
    }
    return false;
  }
  if constexpr(Validatable<Set1>)
  {
    if(!m_set1->isValid(verboseOutput))
    {
      return false;
    }
  }
  if constexpr(Validatable<Set2>)
  {
    if(!m_set2->isValid(verboseOutput))
    {
      return false;
    }
  }
  return true;
}

/**
 * \class BivariateSetIterator
 *
 * \brief Implements a forward iterator concept on a BivariateSet type.
 */
template <typename BivariateSetType>
struct BivariateSetIterator
  : public IteratorBase<BivariateSetIterator<BivariateSetType>, typename BivariateSetType::PositionType>
{
public:
  using IndexType = typename BivariateSetType::PositionType;
  using FirstPositionType = typename BivariateSetType::FirstPositionType;
  using SecondPositionType = typename BivariateSetType::SecondPositionType;
  using BaseType = IteratorBase<BivariateSetIterator<BivariateSetType>, IndexType>;
  using difference_type = IndexType;
  using value_type = typename BivariateSetType::ElementType;
  using reference = value_type;
  using pointer = void;
  using iterator_concept = std::forward_iterator_tag;
  using iterator_category = std::forward_iterator_tag;

  BivariateSetIterator() = default;

  AXOM_HOST_DEVICE BivariateSetIterator(const BivariateSetType* bset, IndexType flatPos = 0)
    : BaseType(flatPos)
    , m_bset(bset)
  { }

  /*!
   * \brief Returns the (first, second) coordinate at this iterator's flat index.
   *
   * \note For ProductSet and RelationSet, secondIndex() takes O(1) time.
   *  firstIndex() also takes O(1) for ProductSet and MappedVariableCardinality.
   *  VariableCardinality uses a binary search that takes O(log(fromSetSize)).
   *  Dereferencing performs both lookups at each position. When the first-set
   *  position is already known, use secondIndex() to avoid looking it up again.
   */
  value_type operator*() const { return {firstIndex(), secondIndex()}; }

  /*!
   * \brief Return the first set index pointed to by this iterator.
   * \note See operator*() for the flat-index lookup cost.
   */
  FirstPositionType firstIndex() const { return m_bset->flatToFirstIndex(flatIndex()); }

  /// \brief Return the second-set position. O(1) for ProductSet and RelationSet.
  SecondPositionType secondIndex() const { return m_bset->flatToSecondIndex(flatIndex()); }

  /// \brief Return the flat iteration index of this iterator.
  AXOM_HOST_DEVICE IndexType flatIndex() const { return this->m_pos; }

protected:
  AXOM_HOST_DEVICE void advance(IndexType n) { this->m_pos += n; }

private:
  const BivariateSetType* m_bset {nullptr};
};

/**
 * \class NullBivariateSet
 *
 * \brief A Null BivariateSet class. Same as the NullSet for Set class.
 */
template <typename SetType1 = slam::Set<>,
          typename SetType2 = slam::Set<>,
          typename Position = detail::default_flat_position_t<typename SetType1::PositionType,
                                                              typename SetType2::PositionType>>
class NullBivariateSet : public BivariateSet<SetType1, SetType2, Position>
{
public:
  using FirstSetType = SetType1;
  using SecondSetType = SetType2;
  using BSet = BivariateSet<FirstSetType, SecondSetType, Position>;
  using FirstPositionType = typename BSet::FirstPositionType;
  using SecondPositionType = typename BSet::SecondPositionType;
  using PositionType = typename BSet::PositionType;
  using ElementType = typename BSet::ElementType;
  using SubsetType = typename BSet::SubsetType;
  using RangeSetType = typename BSet::RangeSetType;

public:
  NullBivariateSet() = default;

  PositionType findElementIndex(FirstPositionType pos1, SecondPositionType pos2 = 0) const override
  {
    verifyPosition(pos1, pos2);
    return PositionType();
  }

  AXOM_SUPPRESS_HD_WARN
  AXOM_HOST_DEVICE PositionType findElementFlatIndex(FirstPositionType s1,
                                                     SecondPositionType s2) const override
  {
    verifyPosition(s1, s2);
    return PositionType();
  }

  PositionType findElementFlatIndex(FirstPositionType s1) const override
  {
    return findElementFlatIndex(s1, 0);
  }

  AXOM_HOST_DEVICE FirstPositionType flatToFirstIndex(PositionType) const override
  {
    return FirstPositionType();
  }

  AXOM_HOST_DEVICE SecondPositionType flatToSecondIndex(PositionType) const override
  {
    return SecondPositionType();
  }

  AXOM_SUPPRESS_HD_WARN
  AXOM_HOST_DEVICE RangeSetType elementRangeSet(FirstPositionType) const override
  {
    return RangeSetType();
  }

  AXOM_HOST_DEVICE ElementType at(PositionType) const override { return ElementType {}; }

  AXOM_HOST_DEVICE PositionType size() const override { return PositionType(); }

  PositionType size(FirstPositionType) const override { return PositionType(); }

  SubsetType getElements(FirstPositionType) const override
  {
    using OrderedSetBuilder = typename SubsetType::SetBuilder;
    return OrderedSetBuilder();
  }

private:
  void verifyPosition(FirstPositionType AXOM_DEBUG_PARAM(pos1),
                      SecondPositionType AXOM_DEBUG_PARAM(pos2)) const override
  {
    SLIC_ASSERT_MSG(false,
                    "Subscripting on NullSet is never valid."
                      << "\n\tAttempted to access item at index " << pos1 << "," << pos2 << ".");
  }
};

}  // end namespace axom::slam
