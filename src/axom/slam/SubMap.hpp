// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#pragma once

/**
 * \file SubMap.hpp
 *
 * \brief Contains SubMap, which is a subset of a Map
 */

#include "axom/slic.hpp"

#include "axom/slam/Map.hpp"
#include "axom/slam/OrderedSet.hpp"
#include "axom/slam/RangeSet.hpp"
#include "axom/core/IteratorBase.hpp"

#include <type_traits>
#include <utility>

namespace axom::slam
{
namespace detail
{
// Separate iterator templates let BivariateMap use them before its type is complete.
template <typename SubMapType>
class SubMapIterator;
template <typename SubMapType>
class SubMapRangeIterator;

/// \brief The parent operations used by a SubMap, independent of storage policies.
template <typename T>
concept SubMapSource = requires(T& map, const T& constMap, typename T::PositionType pos) {
  { constMap.size() } -> std::same_as<typename T::PositionType>;
  { constMap.numComp() } -> std::convertible_to<typename T::PositionType>;
  constMap.shape();
  constMap.index(pos);
  map[pos];
  requires std::is_lvalue_reference_v<decltype(map[pos])>;
  map.set_begin();
  { map.set_end() } -> std::same_as<decltype(map.set_begin())>;
};

}  // namespace detail

/**
 * \class SubMap
 * \brief The SubMap class provides an API to easily traverse a subset of a Map.
 *
 * A SubMap stores a pointer to its parent map and a set of positions selecting
 * entries in that map. It accesses the parent's values and component shape.
 * BivariateMap uses SubMap to return the values associated with one row.
 *
 * set()->at(i) is a position in the immediate parent. index(i) identifies the
 * selected set element, following index() through every parent SubMap.
 * For example, if the parent's set contains {10, 20, 30, 40} and the selected
 * positions are {3, 1}, set()->at(0) is 3 and index(0) is 40. A nested SubMap
 * selecting position 1 identifies element 20. A bivariate parent returns a pair
 * of positions in its first and second sets instead of a scalar set element.
 *
 * \tparam SuperMapType the type of SuperMap
 * \tparam SubsetType defines the indices in the super map. It cannot be abstract.
 *
 * \note Value access preserves the reference type returned by the super-map.
 *       A const SubMap wrapper does not add constness to the mapped values.
 *       To obtain deep-const access to an owning map, use a const SuperMapType.
 * \note The super-map and any storage referenced by the index set must outlive
 *       this SubMap and its iterators. The index set itself is stored by value.
 *       The selected positions must remain valid if the super-map is reassigned.
 *       Iterators must be recreated after changing the parent's value storage
 *       or component shape.
 *
 * \see Map, BivariateMap
 */

template <typename SuperMapType,
          typename SubsetType,  //= slam::RangeSet<PositionType, SetElement>
          typename InterfacePolicy = policies::ConcreteInterface>
class SubMap : public policies::MapInterface<InterfacePolicy, typename SubsetType::PositionType>,
               private detail::HostObjectView
{
public:
  static_assert(!std::is_abstract<SubsetType>::value, "SetType for slam::SubMap cannot be abstract");

  using ParentMapType = SuperMapType;
  /// The set of flat positions selecting entries in the parent map.
  using IndexSetType = SubsetType;
  /// The index set returned by set().
  using SetType = IndexSetType;
  using PositionType = typename SubsetType::PositionType;
  using SetElement = typename SubsetType::ElementType;
  using SuperPositionType = typename SuperMapType::PositionType;
  using ProjectedElement = std::remove_cvref_t<decltype(std::declval<const SuperMapType&>().index(
    std::declval<SuperPositionType>()))>;
  using ElementShape = std::remove_cvref_t<decltype(std::declval<const SuperMapType&>().shape())>;

  //iterator type aliases
  using Iterator = detail::SubMapIterator<SubMap>;
  using iterator = Iterator;
  using const_iterator = Iterator;
  using iterator_pair = std::pair<iterator, iterator>;

  using RangeIterator = detail::SubMapRangeIterator<SubMap>;
  using const_range_iterator = RangeIterator;
  using range_iterator = RangeIterator;

  /*!
   * \brief The reference type for a value access
   *
   * A SubMap is a view, so the constness comes from \a SuperMapType
   * rather than that of the SubMap object.
   */
  using reference = decltype(std::declval<SuperMapType&>()[std::declval<SuperPositionType>()]);
  using const_reference = reference;
  using DataType = std::remove_cvref_t<reference>;
  using DataRefType = reference;
  using ValueType = reference;
  using ConstValueType = const_reference;

public:
  /// Default Constructor
  SubMap() = default;

  /**
   * \brief Constructor for SubMap given the ElementFlatIndex into the SuperMap
   *
   * \param supermap The map that this SubMap is a subset of.
   * \param subset_idxset a Set of ElementFlatIndex into the SuperMap
   * \param indicesHaveIndirection Unused; retained for source compatibility.
   *
   * \note \a indicesHaveIndirection no longer selects between projecting a
   *       subset index through the SuperMap's set and returning it unchanged.
   *       index() now always projects. \see index()
   */
  AXOM_HOST_DEVICE SubMap(SuperMapType* supermap,
                          SubsetType subset_idxset,
                          bool AXOM_UNUSED_PARAM(indicesHaveIndirection) = true)
    : m_superMap(supermap)
    , m_subsetIdx(subset_idxset)
  {
    // Check the parent operations at construction, once the parent type is complete.
    static_assert(detail::SubMapSource<SuperMapType>,
                  "SubMap requires parent size, component, index, value, and range access");
    static_assert(FlatRangeOver<SubsetType, typename SuperMapType::PositionType>,
                  "SubMap requires an index set of flat positions into its super-map");
  }

  /// \name SubMap individual access functions
  /// @{
  ///

  /**
   * \brief Access the value in the SubMap given the ComponentFlatIndex
   *
   * \param idx the ComponentFlatIndex into the subset
   * \return The value for the j<sup>th</sup> component of the i<sup>th</sup>
   *         element, where `setIndex = i * numComp() + j`.
   * \pre    0 <= idx < size() * numComp()
   */
  AXOM_HOST_DEVICE DataRefType operator[](IndexType idx) const
  {
#ifndef AXOM_DEVICE_CODE
    verifyPositionImpl(idx);
#endif
    const SuperPositionType flat_idx = getMapCompFlatIndex(idx);
    return (*m_superMap)[flat_idx];
  }

  /**
   * \brief Access the value associated with the given position in the subset and the component index.
   *
   * \pre `0 <= idx < size()`
   * \pre `0 <= comp < numComp()`
   */
  template <typename... ComponentIndex>
  AXOM_HOST_DEVICE DataRefType operator()(IndexType idx, ComponentIndex... comp) const
  {
    return flatValue(idx, comp...);
  }

  /**
   * \brief Access the value associated with the given position in the subset and the component index.
   *
   * \pre `0 <= idx < size()`
   * \pre `0 <= comp < numComp()`
   */
  template <typename... ComponentIndex>
  AXOM_HOST_DEVICE DataRefType value(IndexType idx, ComponentIndex... comp) const
  {
    return flatValue(idx, comp...);
  }

  /// \brief Access components of the entry at the given subset position.
  template <typename... ComponentIndex>
  AXOM_HOST_DEVICE reference flatValue(PositionType idx, ComponentIndex... comp) const
  {
    static_assert((std::integral<ComponentIndex> && ...),
                  "SubMap component indices must be integral types");
    SLIC_ASSERT_MSG(m_superMap != nullptr, "SubMap's super-map was null.");
    SLIC_ASSERT_MSG(idx >= 0 && idx < size(), "SubMap position is outside the subset.");
    const SuperPositionType parentPos = getMapElemFlatIndex(idx);
    if constexpr(sizeof...(ComponentIndex) <= 1)
    {
      const SuperPositionType component = (SuperPositionType {} + ... + comp);
      SLIC_ASSERT_MSG(component >= 0 && component < numComp(),
                      "SubMap component is outside the element's component range.");
      return (*m_superMap)[parentPos * numComp() + component];
    }
    else
    {
      return m_superMap->flatValue(parentPos, comp...);
    }
  }

  /*!
   * \brief Return the set element selected by a position in this SubMap.
   *
   * Equivalent to parent.index(set()->at(idx)). For nested submaps this follows
   * the selections back to the original map. A bivariate map returns a pair of
   * positions in its first and second sets.
   * \pre 0 <= idx < size()
   */
  AXOM_SUPPRESS_HD_WARN
  AXOM_HOST_DEVICE ProjectedElement index(IndexType idx) const
  {
    const SuperMapType& parent = *m_superMap;
    return parent.index(getMapElemFlatIndex(idx));
  }

  /// @}

  /// \name SubMap cardinality functions
  /// @{
  ///

  /// \brief Return the index set selecting positions in the super-map.
  AXOM_HOST_DEVICE const SetType* set() const { return &m_subsetIdx; }

  /// \brief returns the size of the SubMap
  AXOM_HOST_DEVICE PositionType size() const { return m_subsetIdx.size(); }

  /// \brief returns the number of components (aka. stride) of the SubMap
  AXOM_HOST_DEVICE SuperPositionType numComp() const
  {
    return m_superMap == nullptr ? SuperPositionType {} : m_superMap->numComp();
  }

  /// \brief Return the parent's component count without exposing mutable policy state.
  AXOM_HOST_DEVICE SuperPositionType stride() const { return numComp(); }

  /// \brief Return the parent's component shape.
  AXOM_HOST_DEVICE ElementShape shape() const
  {
    return m_superMap == nullptr ? ElementShape {} : m_superMap->shape();
  }

  /// @}

  [[nodiscard]] bool isValid(bool verboseOutput = false) const
  {
    if constexpr(Validatable<IndexSetType>)
    {
      if(!m_subsetIdx.isValid(verboseOutput))
      {
        return false;
      }
    }
    if(size() < 0 || (m_superMap == nullptr && size() != 0))
    {
      return false;
    }
    for(PositionType pos = 0; pos < size(); ++pos)
    {
      const SuperPositionType parentPos = getMapElemFlatIndex(pos);
      if(parentPos < 0 || parentPos >= m_superMap->size())
      {
        if(verboseOutput)
        {
          SLIC_INFO("Subset index " << parentPos << " is outside the parent map of size "
                                    << m_superMap->size());
        }
        return false;
      }
    }
    return true;
  }

private:  //helper functions
  friend RangeIterator;
  /// \brief Get the ElementFlatIndex into the SuperMap given the subset's index.
  AXOM_HOST_DEVICE SuperPositionType getMapElemFlatIndex(PositionType idx) const
  {
    return set()->at(idx);
  }

  /**
   * \brief Get the ComponentFlatIndex into the SuperMap given the subset's
   * ComponentFlatIndex. This is used only with bracket [] access
   */
  AXOM_HOST_DEVICE SuperPositionType getMapCompFlatIndex(PositionType idx) const
  {
    const SuperPositionType comp = numComp();
    const SuperPositionType s = idx % comp;
    return getMapElemFlatIndex(idx / comp) * comp + s;
  }

  /// Checks the ComponentFlatIndex is valid
  void verifyPosition(PositionType idx) const { verifyPositionImpl(idx); }

  /// Checks the ComponentFlatIndex is valid
  void verifyPositionImpl(PositionType AXOM_DEBUG_PARAM(idx)) const
  {
    SLIC_ASSERT_MSG(idx >= 0 && idx < size() * numComp(),
                    "Attempted to access element " << idx << " but Submap's data has size "
                                                   << size() * numComp());
  }

public:  // Functions related to iteration
  AXOM_HOST_DEVICE iterator begin() const { return iterator(this, 0); }
  AXOM_HOST_DEVICE iterator end() const { return iterator(this, size() * numComp()); }
  AXOM_HOST_DEVICE range_iterator set_begin() const { return range_iterator(this, 0); }
  AXOM_HOST_DEVICE range_iterator set_end() const { return range_iterator(this, size()); }

protected:  //Member variables
  SuperMapType* m_superMap {nullptr};
  IndexSetType m_subsetIdx;

};  //end SubMap

namespace detail
{
/// \brief Scalar traversal of a SubMap, retaining its index metadata by value.
template <typename SubMapType>
class SubMapIterator
  : public IteratorBase<SubMapIterator<SubMapType>, typename SubMapType::PositionType>
{
public:
  using PositionType = typename SubMapType::PositionType;
  using DataType = typename SubMapType::DataType;
  using DataRefType = typename SubMapType::reference;
  using ProjectedElement = typename SubMapType::ProjectedElement;
  using iterator_concept = std::random_access_iterator_tag;
  using iterator_category = std::random_access_iterator_tag;
  using value_type = DataType;
  using reference = DataRefType;
  using pointer = std::add_pointer_t<std::remove_reference_t<reference>>;
  using difference_type = PositionType;

  using IterBase = IteratorBase<SubMapIterator, PositionType>;
  using IterBase::m_pos;
  using iter = SubMapIterator;

  SubMapIterator() = default;

  AXOM_HOST_DEVICE SubMapIterator(const SubMapType* sMap, PositionType pos)
    : IterBase(pos)
    , m_submap(*sMap)
  { }

  /// \brief Returns the current iterator value.
  AXOM_HOST_DEVICE DataRefType operator*() const { return m_submap[m_pos]; }

  AXOM_HOST_DEVICE pointer operator->() const { return &this->operator*(); }

  /// \brief Returns the first component value after n increments.
  DataRefType operator[](PositionType n) const { return *(*this + n); }

  /// \brief Returns the Set element at the iterator's position
  ProjectedElement index() const { return m_submap.index(m_pos / m_submap.numComp()); }

  /// \brief Returns the component index pointed to by this iterator.
  PositionType compIndex() const { return m_pos % m_submap.numComp(); }

  /// \brief Returns the flat index pointed to by this iterator.
  PositionType flatIndex() const { return this->m_pos; }

  /// \brief Returns the number of component per element in the SubMap.
  PositionType numComp() const { return m_submap.numComp(); }

protected:
  /// Implementation of advance() as required by IteratorBase
  AXOM_HOST_DEVICE void advance(PositionType pos) { m_pos += pos; }

private:
  SubMapType m_submap;
};

/// \brief Traversal of the component ranges selected by a SubMap.
template <typename SubMapType>
class SubMapRangeIterator
  : public IteratorBase<SubMapRangeIterator<SubMapType>, typename SubMapType::PositionType>
{
public:
  using PositionType = typename SubMapType::PositionType;
  using SuperPositionType = typename SubMapType::SuperPositionType;
  using SuperMapType = typename SubMapType::ParentMapType;
  using DataRefType = typename SubMapType::reference;
  using ProjectedElement = typename SubMapType::ProjectedElement;
  using IterBase = IteratorBase<SubMapRangeIterator, PositionType>;
  using IterBase::m_pos;
  using iter = SubMapRangeIterator;

private:
  using MapRangeIterator = decltype(std::declval<SuperMapType&>().set_begin());
  static_assert(
    std::default_initializable<MapRangeIterator> &&
      std::constructible_from<MapRangeIterator, SuperMapType*, SuperPositionType>,
    "SubMap range iteration requires a parent iterator constructible at a flat position");

public:
  // Dereference returns a reference to a cached ArrayView,
  // while subscript returns a value to avoid dangling from a temporary iterator.
  // \warning Inherits MapRangeIterator's multipass and dangling caveats;
  //  see the warning on Map::MapRangeIterator.
  using iterator_concept = std::bidirectional_iterator_tag;
  using iterator_category = std::bidirectional_iterator_tag;
  using reference = decltype(*std::declval<const MapRangeIterator&>());
  using value_type = std::remove_cvref_t<reference>;
  using pointer = std::add_pointer_t<std::remove_reference_t<reference>>;
  using difference_type = PositionType;

  SubMapRangeIterator() = default;

private:
  AXOM_HOST_DEVICE static MapRangeIterator makeParentIterator(const SubMapType& submap,
                                                              PositionType pos)
  {
    if(submap.m_superMap == nullptr)
    {
      return MapRangeIterator {};
    }

    // Use the parent's end position for every subset end, including empty subsets.
    const SuperPositionType parentPos =
      pos < submap.size() ? submap.getMapElemFlatIndex(pos) : submap.m_superMap->size();
    return MapRangeIterator(submap.m_superMap, parentPos);
  }

public:
  AXOM_HOST_DEVICE SubMapRangeIterator(const SubMapType* sMap, PositionType pos)
    : IterBase(pos)
    , m_submap(*sMap)
    , m_mapIter(makeParentIterator(*sMap, pos))
  { }

  /// \brief Returns the current iterator value.
  AXOM_HOST_DEVICE reference operator*() const { return (*m_mapIter); }

  AXOM_HOST_DEVICE pointer operator->() const { return m_mapIter.operator->(); }

  template <typename... ComponentIndex>
  AXOM_HOST_DEVICE DataRefType operator()(ComponentIndex... comp_idx) const
  {
    return m_mapIter(comp_idx...);
  }

  template <typename... ComponentIndex>
  AXOM_HOST_DEVICE DataRefType value(ComponentIndex... comp_idx) const
  {
    return m_mapIter.value(comp_idx...);
  }

  AXOM_HOST_DEVICE value_type operator[](PositionType n) const { return *(*this + n); }

  /// \brief Returns the set element mapped by this iterator.
  ProjectedElement index() const { return m_submap.index(this->m_pos); }

  /// \brief Returns the flat position in the original map.
  auto flatIndex() const { return m_mapIter.flatIndex(); }

  /// \brief Returns the index into the submap pointed to by this iterator.
  PositionType submapIndex() const { return this->m_pos; }

  /// \brief Returns the number of components per element in the Map.
  PositionType numComp() const { return m_mapIter.numComp(); }

protected:
  /// \brief Select the new parent position directly, also for nested or reordered subsets.
  AXOM_HOST_DEVICE void advance(PositionType n)
  {
    this->m_pos += n;
    m_mapIter = makeParentIterator(m_submap, this->m_pos);
  }

private:
  SubMapType m_submap;
  MapRangeIterator m_mapIter;
};
}  // namespace detail

}  // namespace axom::slam
