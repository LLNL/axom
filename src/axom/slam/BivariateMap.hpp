// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#pragma once

/**
 * \file BivariateMap.hpp
 *
 * \brief Contains the BivariateMap class, a map for a BivariateSet
 */

#include "axom/slam/Map.hpp"
#include "axom/slam/Relation.hpp"
#include "axom/slam/BivariateSet.hpp"
#include "axom/slam/SubMap.hpp"
#include "axom/slam/Set.hpp"
#include "axom/slam/policies/StridePolicies.hpp"
#include "axom/slam/policies/PolicyTraits.hpp"

#include <cassert>
#include <concepts>
#include <typeinfo>

namespace axom::slam
{
namespace detail
{
/// The extra operations BivariateMap consumes beyond the bivariate set contract.
template <typename T>
concept BivariateMapSet = BivariateSetLike<T> &&
  requires(const T& set,
           typename T::FirstSetType::PositionType first,
           typename T::SecondSetType::PositionType second) {
    { T::INVALID_POS } -> std::convertible_to<typename T::PositionType>;
    { set.findElementIndex(first, second) } -> PositionValueLike;
    { set.findElementFlatIndex(first, second) } -> std::convertible_to<typename T::PositionType>;
    {
      set.elementRangeSet(first)
    }
    -> std::convertible_to<typename RangeSet<typename T::PositionType, typename T::PositionType>::ConcreteSet>;
  };
template <typename Data, typename Set, typename Indirection, typename Stride>
concept BivariateMapParameters =
  BivariateMapSet<Set> && MapStridePolicyFor<Stride, typename Set::PositionType> &&
  MapIndirectionPolicyFor<Indirection, typename Set::PositionType, Data>;
}  // namespace detail

/**
 * \class BivariateMap
 * \brief Associate component values with each pair in a bivariate set.
 *
 * A cell-material map can store a volume fraction for each material present in
 * each cell. `map(i, j, c)` accesses component c for first-set position i and
 * second-set position j. `map(i)` returns a SubMap of all entries associated with i.
 * See BivariateSet for the DenseIndex, SparseIndex and FlatIndex conventions.
 *
 * ElementFlatIndex selects a coordinate pair in the bivariate set.
 * ComponentFlatIndex selects one component in the value buffer. With three
 * components per entry, `flatValue(1, 2)` and `map[5]` access the same value.
 * `index(1)` returns the pair of set positions associated with that entry.
 *
 * \note When \a IndPol is not specified, \c BivariateMap stores its values in an
 *       \c axom::Array via \c policies::ArrayIndirection, and manages that buffer itself.
 *       To refer to a buffer managed elsewhere, use \c policies::ArrayViewIndirection.
 *       For \c std::vector backing, specify \c policies::STLVectorIndirection explicitly.
 *
 * \note Pointer-bound sets are borrowed. Value-bound sets are copied, but their
 *       referenced sets and buffers remain borrowed. Keep these bindings valid
 *       while the map is used. Component shape and value constness follow the
 *       inner Map. A const map over ArrayView<T> can still return T&.
 *
 * \tparam T The type of each component value.
 * \tparam BSet The bivariate set type.
 * \tparam IndPol The value-buffer policy.
 * \tparam StrPol The component-count or shape policy.
 * \tparam IfacePol The concrete or virtual map interface.
 * \see BivariateSet, SubMap
 */

template <typename T,
          typename BSet = BivariateSet<>,
          typename IndPol = policies::ArrayIndirection<typename BSet::PositionType, T>,
          typename StrPol = policies::StrideOne<typename BSet::PositionType>,
          typename IfacePol = policies::ConcreteInterface>
  requires detail::BivariateMapParameters<T, BSet, IndPol, StrPol>
class BivariateMap : public policies::MapInterface<IfacePol, typename BSet::PositionType>,
                     public StrPol
{
public:
  using DataType = T;
  using BivariateSetType = BSet;
  /// The complete set bound by set(), not the internal flat-position set.
  using MappedSetType = BivariateSetType;
  using IndirectionPolicy = IndPol;
  using StridePolicyType = StrPol;

  using FirstPositionType = typename BSet::FirstSetType::PositionType;
  using SecondPositionType = typename BSet::SecondSetType::PositionType;
  using PositionType = typename BSet::PositionType;
  using SetElement = typename BSet::ElementType;

  using ElementShape = typename StridePolicyType::ShapeType;

  // The internal map is indexed by flat bivariate positions.
  // Its backing set is independent of the coordinate ElementType.
  using SetType = typename slam::RangeSet<PositionType, PositionType>::ConcreteSet;
  using MapType = Map<DataType, SetType, IndPol, StrPol, IfacePol>;
  using OrderedSetType = std::remove_cvref_t<decltype(std::declval<const BSet&>().getElements(
    std::declval<FirstPositionType>()))>;

  using ValueType = typename IndirectionPolicy::IndirectionResult;
  using ConstValueType = typename IndirectionPolicy::ConstIndirectionResult;

  using PointerType = std::remove_reference_t<ValueType>*;
  using ConstPointerType = std::remove_reference_t<ConstValueType>*;

  using BivariateMapType = BivariateMap<DataType, BSet, IndPol, StrPol, IfacePol>;

  template <bool Const>
  class FlatIterator;
  using iterator = FlatIterator<false>;
  using const_iterator = FlatIterator<true>;

  template <bool Const>
  class RangeIterator;
  using range_iterator = RangeIterator<false>;
  using const_range_iterator = RangeIterator<true>;

  using SubMapType = SubMap<BivariateMapType, SetType, IfacePol>;
  using ConstSubMapType = const SubMap<const BivariateMapType, SetType, IfacePol>;
  using SubMapIterator = detail::SubMapIterator<SubMapType>;
  using ConstSubMapIterator = detail::SubMapIterator<std::remove_const_t<ConstSubMapType>>;
  using SubMapRangeIterator = detail::SubMapRangeIterator<SubMapType>;
  using ConstSubMapRangeIterator = detail::SubMapRangeIterator<std::remove_const_t<ConstSubMapType>>;

  using NullBivariateSetType =
    NullBivariateSet<typename BSet::FirstSetType, typename BSet::SecondSetType, typename BSet::PositionType>;

private:
  static const NullBivariateSetType s_nullBiSet;

  template <typename USet = BivariateSetType, bool HasValue = !std::is_abstract<USet>::value>
  struct BSetContainer;

  template <typename USet>
  struct BSetContainer<USet, false>
  {
    BSetContainer(const USet* set) : m_pSet(set) { }

    AXOM_HOST_DEVICE const USet* get() const { return m_pSet; }

    const USet* m_pSet;
  };

  template <typename USet>
  struct BSetContainer<USet, true>
  {
    BSetContainer(const USet* set) : m_pSet(set) { }
    BSetContainer(const USet& set) : m_set(set) { }

    AXOM_HOST_DEVICE const USet* get() const
    {
      if(m_pSet)
      {
        return m_pSet;
      }
      else
      {
        return &m_set;
      }
    }

    const USet* m_pSet {nullptr};
    USet m_set;
  };

public:
  using ConcreteMap = BivariateMap<T, BSet, IndPol, StrPol, policies::ConcreteInterface>;
  using VirtualMap = BivariateMap<T, BSet, IndPol, StrPol, policies::VirtualInterface>;

public:
  /**
   * \brief Allocate values for a pointer-bound bivariate set.
   *
   * \param bSet The bivariate set, which must outlive the map.
   * \param defaultValue Initial value of every component.
   * \param shape Component count or multidimensional shape.
   * \param allocatorID Allocator used by the value-buffer policy.
   * \pre bSet is non-null and the shape agrees with any compile-time stride.
   */
  BivariateMap(const BivariateSetType* bSet = &s_nullBiSet,
               DataType defaultValue = DataType(),
               ElementShape shape = StridePolicyType::DefaultSize(),
               int allocatorID = axom::getDefaultAllocatorID())
    requires AllocatingMapIndirectionPolicyFor<IndirectionPolicy, PositionType, DataType>
    : StridePolicyType(shape)
    , m_bset(bSet)
    , m_map(SetType(bSet->size()), defaultValue, shape, allocatorID)
  { }

  /// \overload
  /// \note This value-storing overload accepts only the exact, non-abstract BivariateSetType.
  ///       Use the pointer overload for polymorphic sets.
  template <typename UBSet>
    requires(!std::is_abstract_v<BivariateSetType> && std::same_as<BivariateSetType, UBSet> &&
             AllocatingMapIndirectionPolicyFor<IndirectionPolicy, PositionType, DataType>)
  BivariateMap(const UBSet& bSet,
               DataType defaultValue = DataType(),
               ElementShape shape = StridePolicyType::DefaultSize(),
               int allocatorID = axom::getDefaultAllocatorID())
    : StridePolicyType(shape)
    , m_bset(bSet)
    , m_map(SetType(bSet.size()), defaultValue, shape, allocatorID)
  { }

  /**
   * \brief Bind a bivariate set by pointer and store the supplied value buffer.
   *
   * \param bSet The bivariate set, which must outlive the map.
   * \param data Value buffer passed by value. A view still borrows its allocation.
   * \param shape Component count or multidimensional shape.
   * \pre bSet is non-null. A non-resizable buffer has exactly size() * numComp() entries.
   */
  BivariateMap(const BivariateSetType* bSet,
               typename MapType::OrderedMap data,
               ElementShape shape = StridePolicyType::DefaultSize())
    : StridePolicyType(shape)
    , m_bset(bSet)
    , m_map(SetType(bSet->size()), data, shape)
  { }

  /**
   * \brief Copy a bivariate set and store the supplied value buffer.
   *
   * \param bSet The set to copy. Its referenced sets and buffers remain borrowed.
   * \param data Value buffer passed by value. A view still borrows its allocation.
   * \param shape Component count or multidimensional shape.
   * \pre A non-resizable buffer has exactly size() * numComp() entries.
   * \note This value-storing overload accepts only the exact, non-abstract BivariateSetType.
   *       Use the pointer overload for polymorphic sets.
   */
  template <typename UBSet>
    requires(!std::is_abstract_v<BivariateSetType> && std::same_as<BivariateSetType, UBSet>)
  BivariateMap(const UBSet& bSet,
               typename MapType::OrderedMap data,
               ElementShape shape = StridePolicyType::DefaultSize())
    : StridePolicyType(shape)
    , m_bset(bSet)
    , m_map(SetType(bSet.size()), data, shape)
  { }

  // (KW) Problem -- does not work with RelationSet
  template <typename BivariateSetRetType, typename RelType = void>
    requires(!traits::has_relation_ptr<BivariateSetRetType>::value)
  BivariateSetRetType getBivariateSet() const
  {
    using OuterSet = const typename BivariateSetRetType::FirstSetType;
    using InnerSet = const typename BivariateSetRetType::SecondSetType;
    OuterSet* outer = dynamic_cast<OuterSet*>(set()->getFirstSet());
    InnerSet* inner = dynamic_cast<InnerSet*>(set()->getSecondSet());

    return BivariateSetRetType(outer, inner);
  }

  template <typename BivariateSetRetType, typename RelType>
    requires traits::has_relation_ptr<BivariateSetRetType>::value
  BivariateSetRetType getBivariateSet() const
  {
    auto* rel = dynamic_cast<const RelType*>(m_bset)->getRelation();
    SLIC_ASSERT(rel != nullptr);

    return BivariateSetRetType(rel);
  }

  /// \name BivariateMap value access functions
  /// @{
  ///

  /**
   * \brief Access a value by ComponentFlatIndex.
   *
   * \return The value for the j<sup>th</sup> component of the i<sup>th</sup>
   *         element, where `setIndex = i * numComp() + j`.
   * \pre    0 <= setIndex < size() * numComp()
   */
  AXOM_HOST_DEVICE ConstValueType operator[](PositionType setIndex) const
  {
    return m_map[setIndex];
  }
  AXOM_HOST_DEVICE ValueType operator[](PositionType setIndex) { return m_map[setIndex]; }

public:
  /**
   * \brief Return the mapped subset associated with a first-set position.
   * \pre 0 <= firstIdx < firstSetSize()
   */
  AXOM_HOST_DEVICE ConstSubMapType operator()(FirstPositionType firstIdx) const
  {
#ifndef AXOM_DEVICE_CODE
    verifyFirstSetIndex(firstIdx);
#endif
    auto s = set()->elementRangeSet(firstIdx);
    return ConstSubMapType(this, s);
  }

  AXOM_SUPPRESS_HD_WARN
  AXOM_HOST_DEVICE SubMapType operator()(FirstPositionType firstIdx)
  {
#ifndef AXOM_DEVICE_CODE
    verifyFirstSetIndex(firstIdx);
#endif
    auto s = set()->elementRangeSet(firstIdx);
    return SubMapType(this, s);
  }

  /**
   * \brief Access the value associated with the given DenseIndex into the
   *        BivariateSet and the component index.
   *
   * \pre `0 <= s1 < firstSetSize()`
   * \pre `0 <= s2 < secondSetSize()`
   * \pre The pair (s1, s2) exists in the bivariate set.
   * \pre A single component index is in [0, numComp()). Shaped access takes one
   *      index per dimension, each in [0, shape()[dimension]).
   */
  template <typename... ComponentIndex>
  AXOM_HOST_DEVICE ConstValueType operator()(FirstPositionType s1,
                                             SecondPositionType s2,
                                             ComponentIndex... comp) const
  {
    auto idx = flatIndex(s1, s2);
    return flatValue(idx, comp...);
  }

  template <typename... ComponentIndex>
  AXOM_HOST_DEVICE ValueType operator()(FirstPositionType s1,
                                        SecondPositionType s2,
                                        ComponentIndex... comp)
  {
    auto idx = flatIndex(s1, s2);
    return flatValue(idx, comp...);
  }

  /**
   * \brief Access the value associated with the given FlatIndex into the
   *        BivariateSet and the component index.
   *
   * \pre `0 <= flatIndex < size()`
   * \pre A single component index is in [0, numComp()). Shaped access takes one
   *      index per dimension, each in [0, shape()[dimension]).
   */
  template <typename... ComponentIndex>
  AXOM_HOST_DEVICE ConstValueType flatValue(PositionType flatIndex, ComponentIndex... comp) const
  {
    return m_map(flatIndex, comp...);
  }

  template <typename... ComponentIndex>
  AXOM_HOST_DEVICE ValueType flatValue(PositionType flatIndex, ComponentIndex... comp)
  {
    return m_map(flatIndex, comp...);
  }

  /**
   * \brief Access the value associated with the given DenseIndex into the
   *        BivariateSet and the component index.
   *
   * \pre `0 <= s1 < firstSetSize()`
   * \pre `0 <= s2 < secondSetSize()`
   * \pre A single component index is in [0, numComp()). Shaped access takes one
   *      index per dimension, each in [0, shape()[dimension]).
   *
   * \return a DataType pointer to the value associated with the given index,
   *         or nullptr if there is no value for the given index.
   * \note A RelationSet searches the entries associated with s1.
   */
  template <typename... ComponentIndex>
  AXOM_HOST_DEVICE ConstPointerType findValue(FirstPositionType s1,
                                              SecondPositionType s2,
                                              ComponentIndex... comp) const
  {
    PositionType i = set()->findElementFlatIndex(s1, s2);
    if(i == BivariateSetType::INVALID_POS)
    {
      //the BivariateSet does not contain this index pair
      return nullptr;
    }
    return &(m_map(i, comp...));
  }

  template <typename... ComponentIndex>
  AXOM_HOST_DEVICE PointerType findValue(FirstPositionType s1,
                                         SecondPositionType s2,
                                         ComponentIndex... comp)
  {
    PositionType i = set()->findElementFlatIndex(s1, s2);
    if(i == BivariateSetType::INVALID_POS)
    {
      //the BivariateSet does not contain this index pair
      return nullptr;
    }
    return &(m_map(i, comp...));
  }

  /// @}

  /// \name BivariateMap index access functions
  /// @{
  ///

  /**
   * \brief Return the coordinate pair at the given ElementFlatIndex.
   *
   * Like Map::index(), this returns the set element at a flat position.
   * The pair contains positions in the first and second sets.
   * \pre 0 <= pos < size()
   */
  AXOM_SUPPRESS_HD_WARN
  AXOM_HOST_DEVICE SetElement index(PositionType pos) const { return set()->at(pos); }

  /**
   * \brief Returns the SparseIndex of the element given the DenseIndex
   *
   * \return The SparseIndex of the element, or BivariateSet::INVALID_POS if
   *         the set doesn't contain the given DenseIndex
   */
  auto index(FirstPositionType s1, SecondPositionType s2) const
  {
    return set()->findElementIndex(s1, s2);
  }

  /**
   * \brief Return the second-set positions associated with s1.
   *
   * \param s1 the first set index
   * \return The subset supplied by the bivariate set.
   */
  OrderedSetType indexSet(FirstPositionType s1) const { return set()->getElements(s1); }

  /// \brief Search for the FlatIndex of an element given its DenseIndex in the BivariateSet.
  AXOM_HOST_DEVICE inline PositionType flatIndex(FirstPositionType s1, SecondPositionType s2) const
  {
    return set()->findElementFlatIndex(s1, s2);
  }

  /// @}

public:
  /// BivariateMap iterator functions
  AXOM_HOST_DEVICE iterator begin() { return iterator(this, 0); }

  AXOM_HOST_DEVICE iterator end() { return iterator(this, totalSize() * numComp()); }

  AXOM_HOST_DEVICE const_iterator begin() const { return const_iterator(this, 0); }

  AXOM_HOST_DEVICE const_iterator end() const
  {
    return const_iterator(this, totalSize() * numComp());
  }

  AXOM_HOST_DEVICE range_iterator set_begin() { return range_iterator(this, 0); }

  AXOM_HOST_DEVICE range_iterator set_end() { return range_iterator(this, totalSize()); }

  AXOM_HOST_DEVICE const_range_iterator set_begin() const { return const_range_iterator(this, 0); }

  AXOM_HOST_DEVICE const_range_iterator set_end() const
  {
    return const_range_iterator(this, totalSize());
  }

  /// Iterator via Submap
  AXOM_HOST_DEVICE SubMapIterator begin(FirstPositionType i) { return (*this)(i).begin(); }

  AXOM_HOST_DEVICE SubMapIterator end(FirstPositionType i) { return (*this)(i).end(); }

  AXOM_HOST_DEVICE ConstSubMapIterator begin(FirstPositionType i) const
  {
    return (*this)(i).begin();
  }

  AXOM_HOST_DEVICE ConstSubMapIterator end(FirstPositionType i) const { return (*this)(i).end(); }

  AXOM_HOST_DEVICE SubMapRangeIterator set_begin(FirstPositionType i)
  {
    return (*this)(i).set_begin();
  }

  AXOM_HOST_DEVICE SubMapRangeIterator set_end(FirstPositionType i) { return (*this)(i).set_end(); }

  AXOM_HOST_DEVICE ConstSubMapRangeIterator set_begin(FirstPositionType i) const
  {
    return (*this)(i).set_begin();
  }

  AXOM_HOST_DEVICE ConstSubMapRangeIterator set_end(FirstPositionType i) const
  {
    return (*this)(i).set_end();
  }

public:
  AXOM_HOST_DEVICE const BivariateSetType* set() const { return m_bset.get(); }

  AXOM_HOST_DEVICE const MapType* getMap() const { return &m_map; }

  AXOM_HOST_DEVICE MapType* getMap() { return &m_map; }

  [[nodiscard]] bool isValid(bool verboseOutput = false) const
  {
    if(set() == nullptr || m_map.size() != set()->size() || !m_map.isValid(verboseOutput))
    {
      return false;
    }
    if constexpr(Validatable<BivariateSetType>)
    {
      return set()->isValid(verboseOutput);
    }
    return true;
  }

  /// \name BivariateMap cardinality functions
  /// @{
  ///

  /// \brief Returns the BivariateSet size.
  AXOM_HOST_DEVICE PositionType size() const { return set()->size(); }

  /// \brief Returns the BivariateSet size.
  AXOM_HOST_DEVICE PositionType totalSize() const { return set()->size(); }

  FirstPositionType firstSetSize() const { return set()->getFirstSet()->size(); }

  AXOM_HOST_DEVICE SecondPositionType secondSetSize() const
  {
    return set()->getSecondSet()->size();
  }

  /// \brief Returns the number of the BivariateSet ordered pairs with the given first set index.
  auto size(FirstPositionType s) const { return set()->getElements(s).size(); }

  /// \brief Return the number of components of the map
  AXOM_HOST_DEVICE PositionType numComp() const { return m_map.numComp(); }

  /// Component queries use the inner map's state. The inherited policy base is
  /// retained for compatibility and must not be used to configure this map.
  AXOM_HOST_DEVICE auto stride() const { return m_map.stride(); }
  AXOM_HOST_DEVICE ElementShape shape() const { return m_map.shape(); }
  AXOM_HOST_DEVICE ElementShape strides() const
    requires requires(const MapType& map) {
      { map.strides() } -> std::same_as<ElementShape>;
    }
  {
    return m_map.strides();
  }

  /// @}

  /**
   * \brief Given a DataType array of size `totalSize()*numComp()`, copy
   *        the data into the BivariateMap storage.
   *
   * \param data_arr The array of DataType that contains the data to be copied.
   */
  void copy(const DataType* data_arr)
    requires std::assignable_from<ValueType, const DataType&>
  {
    for(PositionType i = 0; i < m_map.size() * numComp(); i++)
    {
      m_map[i] = data_arr[i];
    }
  }

  /// \brief replace all elements in the Map with the default DataType
  void clear()
    requires requires(MapType& map) { map.clear(); }
  {
    m_map.clear();
  }

private:
  /// \brief Check the given ElementFlatIndex is valid.
  void verifyPosition(PositionType AXOM_DEBUG_PARAM(pos)) const
  {
    SLIC_ASSERT_MSG(pos >= 0 && pos < PositionType(m_map.size()),
                    "Attempted to access element " << pos << " but BivariateMap's data has size "
                                                   << m_map.size());
  }

  void verifyFirstSetIndex(FirstPositionType AXOM_DEBUG_PARAM(firstIdx)) const
  {
    SLIC_ASSERT_MSG(firstIdx >= 0 && firstIdx < firstSetSize(),
                    "Attempted to access elements with first set index "
                      << firstIdx << ", but BivariateMap's first set has size " << firstSetSize());
  }

private:
  BSetContainer<> m_bset;
  MapType m_map;
};  //end BivariateMap

template <typename T, typename BSet, typename IndPol, typename StrPol, typename IfacePol>
  requires detail::BivariateMapParameters<T, BSet, IndPol, StrPol>
typename BivariateMap<T, BSet, IndPol, StrPol, IfacePol>::NullBivariateSetType const
  BivariateMap<T, BSet, IndPol, StrPol, IfacePol>::s_nullBiSet;

/**
 * \class BivariateMap::FlatIterator
 * \brief Traverse individual component values by ComponentFlatIndex.
 */
template <typename T, typename BSet, typename IndPol, typename StrPol, typename IfacePol>
  requires detail::BivariateMapParameters<T, BSet, IndPol, StrPol>
template <bool Const>
class BivariateMap<T, BSet, IndPol, StrPol, IfacePol>::FlatIterator
  : public IteratorBase<FlatIterator<Const>, PositionType>
{
private:
  using IterBase = IteratorBase<FlatIterator<Const>, PositionType>;
  using iter = FlatIterator;

public:
  using DataRefType = std::conditional_t<Const, ConstValueType, ValueType>;
  using BivariateMapPtr = std::conditional_t<Const, const BivariateMap*, BivariateMap*>;

  using iterator_concept = std::random_access_iterator_tag;
  using iterator_category = std::random_access_iterator_tag;
  using value_type = std::remove_cvref_t<DataRefType>;
  using reference = DataRefType;
  using pointer = std::add_pointer_t<std::remove_reference_t<reference>>;
  using difference_type = PositionType;

  static constexpr PositionType INVALID_POS = -2;

public:
  FlatIterator() = default;

  /// \brief Construct an iterator at a ComponentFlatIndex.
  AXOM_HOST_DEVICE FlatIterator(BivariateMapPtr sMap, PositionType pos) : IterBase(pos), m_map(sMap)
  { }

  /// \brief Returns the current map element pointed to by the iterator.
  AXOM_SUPPRESS_HD_WARN
  AXOM_HOST_DEVICE DataRefType operator*() const
  {
    return m_map->flatValue(this->m_pos / numComp(), compIndex());
  }

  AXOM_HOST_DEVICE pointer operator->() const { return &this->operator*(); }

  /// \brief Returns the map value after advancing by \a n flat positions.
  AXOM_HOST_DEVICE reference operator[](PositionType n) const { return *(*this + n); }

  /// \brief return the current iterator's first index into the BivariateSet
  FirstPositionType firstIndex() const { return m_map->index(this->m_pos / numComp()).first; }

  /// \brief return the current iterator's second index (DenseIndex) into the BivariateSet
  SecondPositionType secondIndex() const { return m_map->index(this->m_pos / numComp()).second; }

  /// \brief return the current iterator's component index
  PositionType compIndex() const { return this->m_pos % numComp(); }

  /// \brief Returns the number of components per element in the map.
  AXOM_SUPPRESS_HD_WARN
  AXOM_HOST_DEVICE PositionType numComp() const { return m_map->numComp(); }

protected:
  AXOM_SUPPRESS_HD_WARN
  AXOM_HOST_DEVICE void advance(IndexType n) { this->m_pos += n; }

private:
  BivariateMapPtr m_map {nullptr};
};

/**
 * \class BivariateMap::RangeIterator
 *
 * \brief Traverse the component views of successive bivariate-set entries.
 *
 * Dereferencing returns a reference to the iterator's cached component view.
 * Copy the view to keep it after the iterator advances or is destroyed.
 * The map's value storage must remain valid while the view is used.
 */
template <typename T, typename BSet, typename IndPol, typename StrPol, typename IfacePol>
  requires detail::BivariateMapParameters<T, BSet, IndPol, StrPol>
template <bool Const>
class BivariateMap<T, BSet, IndPol, StrPol, IfacePol>::RangeIterator
  : public IteratorBase<RangeIterator<Const>, PositionType>
{
public:
  using IterBase = IteratorBase<RangeIterator<Const>, PositionType>;

  using MapIterator = typename MapType::template MapRangeIterator<Const>;

  // The underlying MapRangeIterator returns its cached view
  // by-reference from dereference and by-value from subscript.
  // Copies of the view retain access to values, but do not own their storage.
  using iterator_concept = std::bidirectional_iterator_tag;
  using iterator_category = std::bidirectional_iterator_tag;
  using value_type = typename MapIterator::value_type;
  using reference = typename MapIterator::reference;
  using pointer = typename MapIterator::pointer;
  using difference_type = PositionType;

public:
  using DataRefType = typename MapIterator::DataRefType;
  using BivariateMapPtr = std::conditional_t<Const, const BivariateMap*, BivariateMap*>;

  static constexpr PositionType INVALID_POS = -2;

public:
  RangeIterator() = default;

  /// \brief Construct a new BivariateMap Iterator given an ElementFlatIndex
  AXOM_HOST_DEVICE RangeIterator(BivariateMapPtr sMap, PositionType pos)
    : IterBase(pos)
    , m_map(sMap)
    , m_mapIterator(m_map->getMap()->set_begin() + pos)
  { }

  /// \brief Returns the range of elements pointed to by this iterator.
  AXOM_HOST_DEVICE reference operator*() const { return *m_mapIterator; }

  AXOM_HOST_DEVICE pointer operator->() const { return m_mapIterator.operator->(); }

  /**
   * \brief Returns the iterator's value at the given component index.
   *
   * \pre Supply one index in comp_idx per component-shape dimension.
   * \pre Each index is in [0, shape()[dimension]).
   */
  AXOM_SUPPRESS_HD_WARN
  template <typename... ComponentIndex>
  AXOM_HOST_DEVICE DataRefType operator()(ComponentIndex... comp_idx) const
  {
    return value(comp_idx...);
  }

  /// \brief Returns the map view after advancing by \a n set positions.
  AXOM_HOST_DEVICE value_type operator[](PositionType n) const { return *(*this + n); }

  /// \brief Return the value at the iterator's position for a given component index. Same as operator()
  template <typename... ComponentIndex>
  AXOM_HOST_DEVICE DataRefType value(ComponentIndex... comp) const
  {
    return m_mapIterator(comp...);
  }

  /// \brief return the current iterator's first index into the BivariateSet
  FirstPositionType firstIndex() const { return m_map->index(this->m_pos).first; }

  /// \brief return the current iterator's second index (DenseIndex) into the BivariateSet
  SecondPositionType secondIndex() const { return m_map->index(this->m_pos).second; }

  /// \brief Return the current iterator's flat bivariate index.
  AXOM_HOST_DEVICE PositionType flatIndex() const { return m_mapIterator.flatIndex(); }

  /// \brief Returns the number of components per element in the map.
  PositionType numComp() const { return m_map->numComp(); }

protected:
  AXOM_SUPPRESS_HD_WARN
  AXOM_HOST_DEVICE void advance(IndexType n)
  {
    this->m_pos += n;
    m_mapIterator += n;
  }

private:
  BivariateMapPtr m_map {nullptr};
  typename MapType::template MapRangeIterator<Const> m_mapIterator;
};

}  // end namespace axom::slam
