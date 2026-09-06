// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#pragma once

#include "axom/slam/Concepts.hpp"
#include "axom/slam/RangeSet.hpp"
#include "axom/slam/BivariateSet.hpp"
#include "axom/slam/policies/BivariateSetInterfacePolicies.hpp"

#include <optional>

namespace axom::slam
{
namespace detail
{
/// Concrete adapters return the source row. Virtual adapters must convert it to
/// the row type required by the BivariateSet interface.
template <typename R, typename Interface>
concept RelationSetSourceFor = RelationSetSource<R> &&
  (std::same_as<Interface, policies::ConcreteInterface> ||
   (std::same_as<Interface, policies::VirtualInterface> &&
    std::convertible_to<relation_row_t<R>,
                        typename BivariateSet<typename R::FromSetType,
                                              typename R::ToSetType,
                                              typename R::FlatPositionType>::SubsetType>));
}  // namespace detail

/**
 * \class RelationSet
 *
 * \brief Models a Set whose elements are derived from a relation, one element
 *        per fromSet and toSet pair in the relation.
 *
 *  RelationSet models a subset of the Cartesian product of two sets. Users
 *  should refer to the BivariateSet documentation for descriptions of the
 *  different indexing names (SparseIndex, DenseIndex, FlatIndex).
 *
 * \tparam  Relation  The Relation type that this set uses.
 *
 * \see   BivariateSet
 */

template <typename Relation,
          typename SetType1 = typename Relation::FromSetType,
          typename SetType2 = typename Relation::ToSetType,
          typename InterfaceType = policies::VirtualInterface>
  requires(detail::RelationSetSourceFor<Relation, InterfaceType> &&
           std::is_same_v<typename Relation::FromSetType, SetType1> &&
           std::is_same_v<typename Relation::ToSetType, SetType2>)
class RelationSet final
  : public policies::BivariateSetInterface<InterfaceType, SetType1, SetType2, typename Relation::FlatPositionType>
{
public:
  using FirstSetType = typename Relation::FromSetType;
  using SecondSetType = typename Relation::ToSetType;

  using RelationType = Relation;

private:
  using BaseType =
    policies::BivariateSetInterface<InterfaceType, SetType1, SetType2, typename Relation::FlatPositionType>;
  using BaseSubsetType = typename BaseType::SubsetType;

public:
  /*!
   * \brief The type returned by elementRangeSet(): an ordered set of flat positions.
   *
   * \note Public, matching ProductSet and the BivariateSet interface policies.
   *  A private alias here would shadow the inherited one and leave callers
   *  unable to name the return type of a public method.
   */
  using RangeSetType = typename BaseType::RangeSetType;

  using FirstPositionType = typename BaseType::FirstPositionType;
  using SecondPositionType = typename BaseType::SecondPositionType;
  using PositionType = typename BaseType::PositionType;
  using ElementType = typename BaseType::ElementType;

  using RelationSubset = detail::relation_row_t<RelationType>;
  using SubsetType =
    std::conditional_t<std::is_same<void, BaseSubsetType>::value, RelationSubset, BaseSubsetType>;

  using BaseType::INVALID_POS;

  using IteratorType = BivariateSetIterator<RelationSet>;

private:
  // A concrete external row need not support the fixed virtual row interface.
  static auto otherInterfaceType()
  {
    using OtherInterface = std::conditional_t<std::same_as<InterfaceType, policies::VirtualInterface>,
                                              policies::ConcreteInterface,
                                              policies::VirtualInterface>;
    if constexpr(detail::RelationSetSourceFor<Relation, OtherInterface>)
    {
      return std::type_identity<RelationSet<Relation, SetType1, SetType2, OtherInterface>> {};
    }
    else
    {
      return std::type_identity<void> {};
    }
  }

public:
  /// The opposite interface type, or void when the source cannot provide its row type.
  using OtherSet = typename decltype(otherInterfaceType())::type;
  using ConcreteSet =
    std::conditional_t<std::same_as<InterfaceType, policies::ConcreteInterface>, RelationSet, OtherSet>;
  using VirtualSet =
    std::conditional_t<std::same_as<InterfaceType, policies::VirtualInterface>, RelationSet, OtherSet>;

  template <typename Other>
    requires std::same_as<Other, OtherSet>
  RelationSet(const Other& other)
    : BaseType(other.getFirstSet(), other.getSecondSet())
    , m_relation(other.getRelation())
  { }

public:
  RelationSet() = default;

  /**
   * \brief Constructor taking in the relation this BivariateSet is based on.
   * \pre relation pointer must not be a null pointer
   */
  RelationSet(RelationType* relation)
    : BaseType(relation ? relation->fromSet() : policies::EmptySetTraits<FirstSetType>::emptySet(),
               relation ? relation->toSet() : policies::EmptySetTraits<SecondSetType>::emptySet())
    , m_relation(relation)
  {
    SLIC_ASSERT(relation != nullptr);
  }

  /**
   * \brief Searches for the SparseIndex of the element given its DenseIndex.
   * \detail If the element (i,j) is the k<sup>th</sup> non-zero in the row,
   *         then `findElementIndex(i,j)` returns `k`. If `element(i,j)` does
   *         not exist (such as the case of a zero in a sparse matrix), then
   *         `INVALID_POS` is returned.
   *
   * \warning This function can be slow, since a linear search is performed on
   *          the row each time.
   *
   * \param pos1  The first set position.
   * \param pos2  The second set position.
   *
   * \return  The DenseIndex of the given element, or INVALID_POS if such
   *          element is missing from the set.
   * \pre   0 <= pos1 <= set1.size() && 0 <= pos2 <= size2.size()
   */

  AXOM_HOST_DEVICE PositionType findElementIndex(FirstPositionType pos1, SecondPositionType pos2) const
  {
    const PositionType begin = readRelation().offset(pos1);
    const PositionType count = size(pos1);
    for(PositionType i = 0; i < count; ++i)
    {
      if(static_cast<SecondPositionType>(readRelation().relationData()[begin + i]) == pos2)
      {
        return i;
      }
    }
    return BaseType::INVALID_POS;
  }

  /**
   * \brief Optional-returning wrapper for `findElementIndex`.
   *
   * \return An engaged `std::optional` with the SparseIndex if the element exists, else empty.
   */
  [[nodiscard]] std::optional<PositionType> findElementIndexOptional(FirstPositionType pos1,
                                                                     SecondPositionType pos2) const
  {
    const auto idx = findElementIndex(pos1, pos2);
    return idx != BaseType::INVALID_POS ? std::optional<PositionType>(idx)
                                        : std::optional<PositionType> {};
  }

  /**
   * \brief Search for the FlatIndex of the element given its DenseIndex.
   * \warning This function can be slow, since a linear search is performed on the row each time.
   *
   * \param pos1  The first set position.
   * \param pos2  The second set position.
   *
   * \return  The element's FlatIndex
   * \pre   0 <= pos1 <= set1.size() && 0 <= pos2 <= size2.size()
   */
  AXOM_HOST_DEVICE PositionType findElementFlatIndex(FirstPositionType s1, SecondPositionType s2) const
  {
    const PositionType index = findElementIndex(s1, s2);
    return index == BaseType::INVALID_POS
      ? index
      : static_cast<PositionType>(readRelation().offset(s1)) + index;
  }

  /**
   * \brief Optional-returning wrapper for `findElementFlatIndex(s1, s2)`.
   *
   * \return An engaged `std::optional` with the FlatIndex if the element exists, else empty.
   */
  [[nodiscard]] AXOM_HOST_DEVICE std::optional<PositionType> findElementFlatIndexOptional(
    FirstPositionType s1,
    SecondPositionType s2) const
  {
    const auto idx = findElementFlatIndex(s1, s2);
    return idx != BaseType::INVALID_POS ? std::optional<PositionType>(idx)
                                        : std::optional<PositionType> {};
  }

  /**
   * \brief Given the from-set index pos1, return the FlatIndex of the first
   *        existing to-set element in the relation pair, or `INVALID_POS` if
   *        this row contains no elements.
   *
   * \param pos1  Index into the from-set.
   * \param pos2  Index into the to-set.
   *
   * \return  The FlatIndex of the first existing to-set element.
   */
  PositionType findElementFlatIndex(FirstPositionType pos1) const
  {
    return size(pos1) == 0 ? BaseType::INVALID_POS
                           : static_cast<PositionType>(readRelation().offset(pos1));
  }

  /**
   * \brief Optional-returning wrapper for `findElementFlatIndex(pos1)`.
   *
   * \return An engaged `std::optional` with the FlatIndex if the row contains any elements, else empty.
   */
  [[nodiscard]] std::optional<PositionType> findElementFlatIndexOptional(FirstPositionType pos1) const
  {
    const auto idx = findElementFlatIndex(pos1);
    return idx != BaseType::INVALID_POS ? std::optional<PositionType>(idx)
                                        : std::optional<PositionType> {};
  }

  /**
   * \brief Given the flat index, return the associated to-set index in the relation pair.
   *
   * \param flatIndex The FlatIndex of the from-set/to-set pair.
   *
   * \return pos2  The to-set index.
   */
  AXOM_SUPPRESS_HD_WARN
  AXOM_HOST_DEVICE SecondPositionType flatToSecondIndex(PositionType flatIndex) const
  {
#ifndef AXOM_DEVICE_CODE
    SLIC_ASSERT_MSG(flatIndex >= 0 && flatIndex < size(),
                    "SLAM::RelationSet -- requested out-of-range flat index "
                      << flatIndex << "; set has " << size() << " elements.");
#endif
    return readRelation().relationData()[flatIndex];
  }

  /**
   * \brief Given the flat index, return the associated from-set index in the relation pair.
   *
   * \param flatIndex The FlatIndex of the from-set/to-set pair.
   *
   * \return pos1  The from-set index.
   */
  AXOM_HOST_DEVICE FirstPositionType flatToFirstIndex(PositionType flatIndex) const
  {
#ifndef AXOM_DEVICE_CODE
    SLIC_ASSERT_MSG(flatIndex >= 0 && flatIndex < size(),
                    "SLAM::RelationSet -- requested out-of-range flat index "
                      << flatIndex << "; set has " << size() << " elements.");
#endif
    return static_cast<FirstPositionType>(readRelation().firstIndex(flatIndex));
  }

  AXOM_HOST_DEVICE RangeSetType elementRangeSet(FirstPositionType pos1) const
  {
    return typename RangeSetType::SetBuilder().size(size(pos1)).offset(readRelation().offset(pos1));
  }

  /**
   * \brief A set of elements with the given first set index.
   *
   * \param s1  The first set index.
   * \return  The row of to-set positions associated with s1.
   * \pre  0 <= pos1 <= set1.size()
   */
  SubsetType getElements(FirstPositionType s1) const { return readRelation()[s1]; }

  AXOM_SUPPRESS_HD_WARN
  /*!
   * \brief Returns the (first, second) coordinate at flat index \a pos.
   *
   * \note Costs whatever the relation's cardinality policy charges to invert a
   *  flat index back to a row: O(1) for MappedVariableCardinality, which stores
   *  the inverse, and O(log(fromSetSize)) for plain VariableCardinality,
   *  which binary-searches the begins array. flatToSecondIndex() alone is O(1),
   *  so prefer it when the row is already known.
   */
  [[nodiscard]] AXOM_HOST_DEVICE ElementType at(PositionType pos) const
  {
#ifndef AXOM_DEVICE_CODE
    RelationSet::verifyPosition(pos);
#endif
    return {flatToFirstIndex(pos), flatToSecondIndex(pos)};
  }

  /// \brief Returns the relation pointer
  RelationType* getRelation() const { return m_relation; }

  RelationType* getRelation() { return m_relation; }

  /// \brief Return the size of the relation
  PositionType totalSize() const
  {
    return static_cast<PositionType>(readRelation().relationData().size());
  }

  /**
   * \brief Return the size of a row, which is the number of to-set
   *        elements associated with the given from-set index.
   *
   * \param pos The from-set position.
   */
  AXOM_HOST_DEVICE PositionType size(FirstPositionType pos) const
  {
    return static_cast<PositionType>(readRelation()[pos].size());
  }

  /// \brief Return an iterator to the first pair of set elements in the relation.
  IteratorType begin() const { return IteratorType(this, 0); }

  /// \brief Return an iterator to one past the last pair of set elements in the relation.
  IteratorType end() const { return IteratorType(this, totalSize()); }

  [[nodiscard]] bool isValid(bool verboseOutput = false) const
  {
    if(m_relation == nullptr)
    {
      if(verboseOutput)
      {
        std::cout << "\n*** RelationSet is not valid:\n"
                  << "\t* Relation pointer should not be null.\n"
                  << std::endl;
      }
      return false;
    }
    return readRelation().isValid(verboseOutput);
  }

public:
  //hiding size() from the Set base class, replaced with totalSize().
  //but still implemented due to the function being virtual
  //(and can be called from base ptr)
  // KW -- made this public to use from BivariateMap
  AXOM_SUPPRESS_HD_WARN
  [[nodiscard]] AXOM_HOST_DEVICE PositionType size() const
  {
    return static_cast<PositionType>(readRelation().relationData().size());
  }

  /// \brief Checks if there are any elements in the set
  AXOM_SUPPRESS_HD_WARN
  [[nodiscard]] AXOM_HOST_DEVICE bool empty() const { return size() == PositionType {}; }

private:
  AXOM_HOST_DEVICE const RelationType& readRelation() const { return *m_relation; }

  //range check only
  [[nodiscard]] bool isValidIndex(FirstPositionType s1, SecondPositionType s2) const
  {
    return s1 >= 0 && s1 < m_relation->fromSet()->size() && s2 >= 0 &&
      s2 < m_relation->toSet()->size();
  }

  void verifyPosition(PositionType AXOM_DEBUG_PARAM(sPos)) const
  {
    SLIC_ASSERT_MSG(sPos >= 0 && sPos < size(),
                    "SLAM::RelationSet -- requested out-of-range element at position "
                      << sPos << ", but set only has " << size() << " elements.");
  }

  void verifyPosition(FirstPositionType AXOM_DEBUG_PARAM(s1),
                      SecondPositionType AXOM_DEBUG_PARAM(s2)) const
  {
    SLIC_ASSERT_MSG(isValidIndex(s1, s2),
                    "SLAM::RelationSet -- requested out-of-range element at position ("
                      << s1 << "," << s2 << "), but set only has " << this->firstSetSize() << "x"
                      << this->secondSetSize() << " elements.");
  }

private:
  RelationType* m_relation {nullptr};  //the relation that this set is based off of
};

}  // end namespace axom::slam
