// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#pragma once

/**
 * \file Map.hpp
 *
 * \brief Values associated with set elements, with a fixed component count.
 */

#include <concepts>
#include <iostream>
#include <sstream>
#include <vector>

#include "axom/core.hpp"
#include "axom/slic.hpp"

#include "axom/slam/MapBase.hpp"
#include "axom/slam/Concepts.hpp"
#include "axom/slam/Set.hpp"
#include "axom/slam/NullSet.hpp"
#include "axom/slam/detail/SizeChecks.hpp"

#include "axom/core/IteratorBase.hpp"
#include "axom/core/RangeAdapter.hpp"

#include "axom/slam/policies/StridePolicies.hpp"
#include "axom/slam/policies/IndirectionPolicies.hpp"
#include "axom/slam/policies/PolicyTraits.hpp"
#include "axom/slam/policies/MapInterfacePolicies.hpp"

namespace axom::slam
{
// This class is missing some simplifying copy constructors
// -- or at least ways of interacting with the data store
// We should probably support copy on write shallow copies when possible...

/**
 * \class   Map
 *
 * \brief Associates a fixed number of component values with each set element.
 *
 * \tparam  T The data type of each value
 * \tparam  S The map's set type
 * \tparam  IndPol The map's indirection policy
 * \tparam  StrPol A policy class that determines how many values to
 *          associate with each element. There is a fixed \a stride between
 *          the data associated with each element of the set.
 * \tparam IfacePol Selects a virtual or concrete map interface.
 * \details A temperature map can store one value per cell, while a velocity map
 * stores several components. The stride policy sets the component count or shape
 * at compile time or runtime. Access uses set positions, not set element values.
 * `map(i, j)` and `map[i * numComp() + j]` access component j at set position i.
 * `index(i)` returns the associated set element.
 *
 * \note When \a IndPol is not specified, \c Map stores its values in an \c axom::Array
 *       via \c policies::ArrayIndirection, and manages that buffer itself.
 *       To refer to a buffer managed elsewhere, use \c policies::ArrayViewIndirection.
 *       For \c std::vector backing, specify \c policies::STLVectorIndirection explicitly.
 * \note Component counts must be positive and the storage size must fit both
 *       PositionType and axom::IndexType. Constructors check these conditions.
 *       Referenced sets must retain a size consistent with the value buffer.
 *       Construct or assign a map with the desired shape rather than changing
 *       its inherited stride policy.
 * \note A pointer-bound set must outlive the map. Value-bound sets are copied,
 *       but any storage they reference must remain valid. A const owning map
 *       returns const value references. ArrayView-backed maps preserve the
 *       view's constness, so a const map over ArrayView<T> can still return T&.
 */

template <typename T,
          typename S = Set<>,
          typename IndPol = policies::ArrayIndirection<typename S::PositionType, T>,
          typename StrPol = policies::StrideOne<typename S::PositionType>,
          typename IfacePol = policies::ConcreteInterface>
  requires detail::MapParameters<T, S, IndPol, StrPol>
class Map : public StrPol, public policies::MapInterface<IfacePol, typename S::PositionType>
{
public:
  using DataType = T;
  using SetType = S;
  using IndirectionPolicy = IndPol;
  using StridePolicyType = StrPol;

  using OrderedMap = typename IndirectionPolicy::IndirectionBufferType;

  using PositionType = typename SetType::PositionType;
  using SetElement = typename SetType::ElementType;
  /// The complete set bound by set(), as distinct from an indexing helper.
  using MappedSetType = SetType;
  static const NullSet<PositionType, SetElement> s_nullSet;

  using ElementShape = typename StridePolicyType::ShapeType;

  using ValueType = typename IndirectionPolicy::IndirectionResult;
  using ConstValueType = typename IndirectionPolicy::ConstIndirectionResult;

  class MapBuilder;

  // types for iterator
  template <bool Const>
  class MapIterator;
  using const_iterator = MapIterator<true>;
  using const_iterator_pair = std::pair<const_iterator, const_iterator>;
  using iterator = MapIterator<false>;
  using iterator_pair = std::pair<iterator, iterator>;

  template <bool Const>
  class MapRangeIterator;
  using const_range_iterator = MapRangeIterator<true>;
  using range_iterator = MapRangeIterator<false>;

public:
  using ConcreteMap = Map<T, S, IndPol, StrPol, policies::ConcreteInterface>;
  using VirtualMap = Map<T, S, IndPol, StrPol, policies::VirtualInterface>;

private:
  template <typename USet = SetType, bool HasValue = !std::is_abstract<USet>::value>
  struct SetContainer;

  template <typename USet>
  struct SetContainer<USet, false>
  {
    SetContainer(const USet* set) : m_pSet(set) { }

    AXOM_HOST_DEVICE const USet* get() const { return m_pSet; }

    const USet* m_pSet;
  };

  template <typename USet>
  struct SetContainer<USet, true>
  {
    SetContainer(const USet* set) : m_pSet(set) { }
    SetContainer(const USet& set) : m_set(set) { }

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
  /**
   * \brief Allocate values for a pointer-bound set.
   *
   * \param theSet The set, which must outlive the map.
   * \param defaultValue Initial value of every component.
   * \param shape Component count or multidimensional shape, as specified by StrPol.
   * \param allocatorID Allocator used by the buffer policy.
   * \pre The shape has positive dimensions and agrees with any compile-time stride.
   */

  Map(const SetType* theSet = policies::EmptySetTraits<SetType>::emptySet(),
      DataType defaultValue = DataType(),
      ElementShape shape = StridePolicyType::DefaultSize(),
      int allocatorID = axom::getDefaultAllocatorID())
    requires AllocatingMapIndirectionPolicyFor<IndirectionPolicy, PositionType, DataType>
    : StridePolicyType(shape)
    , m_set(theSet)
    , m_data(IndirectionPolicy::create(requiredStorageSize(), defaultValue, allocatorID))
  { }

  /**
   * \brief Constructor for Map from a Set pointer and an existing buffer.
   *
   * The buffer argument is moved into the map. An owning buffer retains ownership,
   * while a view continues to borrow its allocation.
   *
   * \param theSet The set, which must outlive the map.
   * \param data Value buffer. Any borrowed allocation must outlive the map's use of it.
   * \param shape Component count or multidimensional shape.
   * \pre A non-resizable buffer has exactly size() * numComp() entries.
   */
  Map(const SetType* theSet, OrderedMap data, ElementShape shape = StridePolicyType::DefaultSize())
    : StridePolicyType(shape)
    , m_set(theSet)
    , m_data(std::move(data))
  {
    checkBackingSize(std::integral_constant<bool, IndirectionPolicy::IsMutableBuffer> {});
  }

  /// \overload
  /// \note This value-storing overload accepts only the exact, non-abstract SetType.
  ///       Use the pointer overload for polymorphic sets.
  template <typename USet>
    requires(!std::is_abstract_v<SetType> && std::same_as<SetType, USet> &&
             AllocatingMapIndirectionPolicyFor<IndirectionPolicy, PositionType, DataType>)
  Map(const USet& theSet,
      DataType defaultValue = DataType(),
      ElementShape shape = StridePolicyType::DefaultSize(),
      int allocatorID = axom::getDefaultAllocatorID())
    : StridePolicyType(shape)
    , m_set(theSet)
    , m_data(IndirectionPolicy::create(requiredStorageSize(), defaultValue, allocatorID))
  { }

  /**
   * \brief Copy the set and store the supplied value buffer.
   *
   * \param theSet The set to copy. Any storage it references remains borrowed.
   * \param data Value buffer, moved into the map. A view still borrows its allocation.
   * \param shape Component count or multidimensional shape.
   * \pre A non-resizable buffer has exactly size() * numComp() entries.
   * \note This value-storing overload accepts only the exact, non-abstract SetType.
   *       Use the pointer overload for polymorphic sets.
   */
  template <typename USet>
    requires(!std::is_abstract_v<SetType> && std::same_as<SetType, USet>)
  Map(const USet& theSet, OrderedMap data, ElementShape shape = StridePolicyType::DefaultSize())
    : StridePolicyType(shape)
    , m_set(theSet)
    , m_data(std::move(data))
  {
    checkBackingSize(std::integral_constant<bool, IndirectionPolicy::IsMutableBuffer> {});
  }

  /// \brief Constructor for Map using a MapBuilder
  Map(const MapBuilder& builder)
    requires AllocatingMapIndirectionPolicyFor<IndirectionPolicy, PositionType, DataType> &&
    std::assignable_from<ValueType, DataType&>
    : Map(builder.m_set, builder.m_defaultValue, builder.m_stride.shape())
  {
    //copy the data if exists
    if(builder.m_data_ptr)
    {
      for(PositionType idx = PositionType(); idx < size() * numComp(); ++idx)
      {
        (*this)[idx] = builder.m_data_ptr[idx];
      }
    }
  }

  /// \brief Returns a pointer to the map's underlying set
  AXOM_HOST_DEVICE const SetType* set() const { return m_set.get(); }

  /// \name Map individual access functions
  /// @{
  ///

  /**
   * \brief  Access the value in the map using a flat index in the range of 0 to
   *         `size()*numComp()`
   *
   * \return The value for the j<sup>th</sup> component of the i<sup>th</sup>
   *         element, where `setIndex = i * numComp() + j`.
   * \pre    0 <= setIndex < size() * numComp()
   */
  AXOM_HOST_DEVICE ConstValueType operator[](PositionType setIndex) const
  {
#ifndef AXOM_DEVICE_CODE
    verifyPositionImpl(setIndex);
#endif
    return *IndirectionPolicy::getConstIndirection(m_data, setIndex);
  }

  AXOM_HOST_DEVICE ValueType operator[](PositionType setIndex)
  {
#ifndef AXOM_DEVICE_CODE
    verifyPositionImpl(setIndex);
#endif
    return *IndirectionPolicy::getIndirection(m_data, setIndex);
  }

  /// \brief Access the value associated with the given position in the set.
  AXOM_HOST_DEVICE ConstValueType operator()(PositionType setIdx) const
  {
    // TODO: validate that runtime stride is 1-D with value 1?
    return value(setIdx, 0);
  }

  /// \overload
  AXOM_HOST_DEVICE ValueType operator()(PositionType setIdx) { return value(setIdx, 0); }

  /**
   * \brief Access the value associated with the given position in the set and
   *        the component index.
   *
   * \pre `0 <= setIdx < size()`
   * A single component index is in [0, numComp()). Otherwise supply one index
   * per shape dimension, with each index in [0, shape()[idim]).
   */
  template <typename... ComponentPos>
  AXOM_HOST_DEVICE ConstValueType operator()(PositionType setIdx, ComponentPos... compIdx) const
  {
    return value(setIdx, compIdx...);
  }

  /// \overload
  template <typename... ComponentPos>
  AXOM_HOST_DEVICE ValueType operator()(PositionType setIdx, ComponentPos... compIdx)
  {
    return value(setIdx, compIdx...);
  }

  AXOM_HOST_DEVICE ConstValueType value(PositionType setIdx) const { return value(setIdx, 0); }

  AXOM_HOST_DEVICE ValueType value(PositionType setIdx) { return value(setIdx, 0); }

  /**
   * \brief Access the value associated with the given position in the set and
   *        the component index.
   *
   * \pre `0 <= setIdx < size()`
   * A single component index is in [0, numComp()). Otherwise supply one index
   * per shape dimension, with each index in [0, shape()[idim]).
   */
  template <typename... ComponentPos>
  AXOM_HOST_DEVICE ConstValueType value(PositionType setIdx, ComponentPos... compIdx) const
  {
    static_assert(sizeof...(ComponentPos) == 1 || sizeof...(ComponentPos) == StridePolicyType::NumDims,
                  "Invalid number of components provided for given Map's StridePolicy");
    static_assert(axom::detail::all_types_are_integral<ComponentPos...>::value,
                  "Map::value(...): index parameter pack must all be integral types.");
#ifndef AXOM_DEVICE_CODE
    verifyPositionImpl(setIdx, compIdx...);
#endif
    PositionType elemIndex = setIdx * StridePolicyType::stride();
    elemIndex += componentOffset(compIdx...);
    return *IndirectionPolicy::getConstIndirection(m_data, elemIndex);
  }

  /// \overload
  template <typename... ComponentPos>
  AXOM_HOST_DEVICE ValueType value(PositionType setIdx, ComponentPos... compIdx)
  {
    static_assert(sizeof...(ComponentPos) == 1 || sizeof...(ComponentPos) == StridePolicyType::NumDims,
                  "Invalid number of components provided for given Map's StridePolicy");
    static_assert(axom::detail::all_types_are_integral<ComponentPos...>::value,
                  "Map::value(...): index parameter pack must all be integral types.");
#ifndef AXOM_DEVICE_CODE
    verifyPositionImpl(setIdx, compIdx...);
#endif
    PositionType elemIndex = setIdx * StridePolicyType::stride();
    elemIndex += componentOffset(compIdx...);
    return *IndirectionPolicy::getIndirection(m_data, elemIndex);
  }

  /// \brief Return the set element at the given position, without a component offset.
  AXOM_SUPPRESS_HD_WARN
  AXOM_HOST_DEVICE SetElement index(IndexType idx) const { return set()->at(idx); }

  /// \brief Access components by set position. Equivalent to value(setIdx, compIdx...).
  template <typename... ComponentPos>
  AXOM_HOST_DEVICE ConstValueType flatValue(PositionType setIdx, ComponentPos... compIdx) const
  {
    return value(setIdx, compIdx...);
  }

  /// \overload
  template <typename... ComponentPos>
  AXOM_HOST_DEVICE ValueType flatValue(PositionType setIdx, ComponentPos... compIdx)
  {
    return value(setIdx, compIdx...);
  }

  /// @}

  /// \name Map cardinality functions
  /// @{

  /**
   * \brief Return the set size. Same as `set()->size()`.
   *
   * The total storage size for the map's values is `size() * numComp()`
   */
  AXOM_SUPPRESS_HD_WARN
  [[nodiscard]] AXOM_HOST_DEVICE PositionType size() const
  {
    return !policies::EmptySetTraits<SetType>::isEmpty(m_set.get())
      ? static_cast<PositionType>(m_set.get()->size())
      : PositionType(0);
  }

  /**
   * \brief  Gets the number of component values associated with each element.
   *         Equivalent to stride().
   */
  [[nodiscard]] AXOM_HOST_DEVICE PositionType numComp() const { return StridePolicyType::stride(); }

  /**
   * \brief Returns the shape of the component values associated with each element.
   *
   *  For one-dimensional strides, equivalent to stride(). Otherwise, returns
   *  an N-dimensional array with the number of values in each sub-component index.
   */
  AXOM_HOST_DEVICE ElementShape shape() const { return StridePolicyType::shape(); }

  /// @}

  /// \name Map modifying functions for all entries
  /// @{

  /** \brief replace all elements in the Map with the default DataType */
  void clear()
    requires std::assignable_from<ValueType, DataType>
  {
    fill();
  }

  /** Set each entry in the map to the given value  */
  void fill(DataType val = DataType())
    requires std::assignable_from<ValueType, DataType>
  {
    const PositionType sz = static_cast<PositionType>(m_data.size());

    for(PositionType idx = PositionType(); idx < sz; ++idx)
    {
      (*this)[idx] = val;
    }
  }

  /** \brief Element-wise copy of data from another map */
  void copy(const Map& other)
    requires std::assignable_from<ValueType, ConstValueType>
  {
    SLIC_ASSERT(other.size() == size());
    SLIC_ASSERT(other.stride() == StridePolicyType::stride());

    const PositionType sz = size() * StridePolicyType::stride();
    for(PositionType idx = PositionType(); idx < sz; ++idx)
    {
      (*this)[idx] = other[idx];
    }
  }

  ///@}

  /// \brief print information on the map, including every element inside Map
  void print() const;

  /// \brief returns true if the map is valid, false otherwise
  [[nodiscard]] bool isValid(bool verboseOutput = false) const;

public:
  /**
   * \class MapBuilder
   * \brief Helper class for constructing a Map
   **/
  class MapBuilder
  {
  public:
    friend class Map;

    MapBuilder() : m_set(policies::EmptySetTraits<SetType>::emptySet()) { }

    /// \brief Provide the Set to be used by the Map
    MapBuilder& set(const SetType* set)
    {
      m_set = set;
      return *this;
    }

    /// \brief Set the stride of the Map using StridePolicy
    MapBuilder& stride(PositionType str)
      requires std::constructible_from<StridePolicyType, PositionType> &&
      std::assignable_from<StridePolicyType&, StridePolicyType>
    {
      m_stride = StridePolicyType(str);
      return *this;
    }

    /// \brief Set the pointer to the array of data the Map will contain
    /// The array must contain size() * numComp() values for the configured map.
    MapBuilder& data(DataType* bufPtr)
    {
      m_data_ptr = bufPtr;
      return *this;
    }

  private:
    const SetType* m_set;
    StridePolicyType m_stride {StridePolicyType::DefaultSize()};
    DataType* m_data_ptr = nullptr;
    DataType m_defaultValue = DataType();
  };

public:
  /**
   * \class MapIterator
   * \brief An iterator type for a map. Each increment operation advances the
   *        iterator to the element at the next flat index.
   */
  template <bool Const>
  class MapIterator : public IteratorBase<MapIterator<Const>, PositionType>
  {
  public:
    using DataRefType = std::conditional_t<Const, ConstValueType, ValueType>;
    using DataType = std::remove_reference_t<DataRefType>;

    using iterator_concept = std::random_access_iterator_tag;
    using iterator_category = std::random_access_iterator_tag;
    using value_type = std::remove_cv_t<DataType>;
    using reference = DataRefType;
    using pointer = std::add_pointer_t<std::remove_reference_t<reference>>;
    using difference_type = PositionType;

    using IterBase = IteratorBase<MapIterator, PositionType>;
    using MapConstPtr = std::conditional_t<Const, const Map*, Map*>;
    using iter = MapIterator;
    using IterBase::m_pos;

  public:
    MapIterator() = default;

    MapIterator(MapConstPtr oMap, PositionType pos) : IterBase(pos), m_map(oMap) { }

    /// \brief Returns the current iterator value.
    AXOM_HOST_DEVICE reference operator*() const { return (*m_map)[m_pos]; }

    AXOM_HOST_DEVICE pointer operator->() const { return &this->operator*(); }

    /// \brief Returns the map value after advancing by \a n flat positions.
    AXOM_HOST_DEVICE reference operator[](PositionType n) const { return *(*this + n); }

    /// \brief Returns the set element mapped by this iterator.
    SetElement index() const { return m_map->index(this->m_pos / m_map->numComp()); }

    /// \brief Returns the component index pointed to by this iterator.
    PositionType compIndex() const { return m_pos % m_map->numComp(); }

    /// \brief Returns the flat index pointed to by this iterator.
    PositionType flatIndex() const { return this->m_pos; }

  protected:
    /// Implementation of advance() as required by IteratorBase
    AXOM_HOST_DEVICE void advance(PositionType n) { m_pos += n; }

  private:
    MapConstPtr m_map {nullptr};
  };

  /**
   * \class   MapRangeIterator
   * \brief   An iterator type for a map.
   *          Each increment operation advances the iterator to the next set element.
   *          To access the j<sup>th</sup> component values of the iterator's current element, use `iter(j)`.
   * `iter[off]` returns the component view at an offset of off set positions,
   * equivalent to `*(iter + off)`. `iter(j)` accesses component j of the current
   * set element for a one-dimensional shape. Multidimensional access takes one
   * index per dimension.
   *
   * Dereferencing returns a reference to the iterator's cached view. Copy the
   * view to keep it after the iterator advances or is destroyed. The map's
   * value storage must remain valid while either view is used.
   */
  template <bool Const>
  class MapRangeIterator : public IteratorBase<MapRangeIterator<Const>, PositionType>
  {
  public:
    using IterBase = IteratorBase<MapRangeIterator, PositionType>;
    using MapConstPtr = std::conditional_t<Const, const Map*, Map*>;

    using DataRefType = std::conditional_t<Const, ConstValueType, ValueType>;
    using DataType = std::remove_reference_t<DataRefType>;

    constexpr static int Dims = StridePolicyType::NumDims;

    // Dereference returns a reference to a cached ArrayView, while subscript
    // returns a value to avoid dangling from a temporary iterator.
    //
    // \warning Two consequences of that cached reference:
    //  (1) `*it` binds to a member of `it`, so equal iterators yield references to different objects.
    //      That satisfies std::bidirectional_iterator syntactically but violates the multipass guarantee in
    //      [forward.iterators]/6, which algorithms are entitled to rely on.
    //  (2) `const auto& v = *map.set_begin();` dangles, because the iterator temporary owns the referent;
    //      bind by value, or keep the iterator alive.
    //  Returning `value_type` by value from `operator*` would fix both and restore random access, at the cost of `operator->`.
    using iterator_concept = std::bidirectional_iterator_tag;
    using iterator_category = std::bidirectional_iterator_tag;
    using value_type = axom::ArrayView<DataType, Dims>;
    using reference = const value_type&;
    using pointer = const value_type*;
    using difference_type = PositionType;

  private:
    AXOM_HOST_DEVICE static value_type makeRange(MapConstPtr map, PositionType pos)
    {
      // An empty owning buffer has no element zero from which to obtain a pointer.
      auto* data = map->size() == 0 ? nullptr : map->data_ptr();
      if(data != nullptr)
      {
        data += pos * map->stride();
      }
      if constexpr(Dims == 1)
      {
        return value_type(data, static_cast<axom::IndexType>(map->shape()));
      }
      else
      {
        StackArray<axom::IndexType, Dims> shape;
        for(int dim = 0; dim < Dims; ++dim)
        {
          shape[dim] = map->shape()[dim];
        }
        return value_type(data, shape);
      }
    }

  public:
    MapRangeIterator() = default;

    AXOM_HOST_DEVICE MapRangeIterator(MapConstPtr oMap, PositionType pos)
      : IterBase(pos)
      , m_map(oMap)
      , m_currRange(makeRange(oMap, pos))
    { }

    /// \brief Returns the current iterator value.
    AXOM_HOST_DEVICE reference operator*() const { return m_currRange; }

    AXOM_HOST_DEVICE pointer operator->() const { return &m_currRange; }

    /**
     * \brief Returns the iterator's value at the specified component.
     * \param comp_idx Zero-based indices, one per component-shape dimension.
     */
    template <typename... ComponentIndex>
    AXOM_HOST_DEVICE DataRefType operator()(ComponentIndex... comp_idx) const
    {
      return value(comp_idx...);
    }
    template <typename ComponentIndex>
    AXOM_HOST_DEVICE DataRefType value(ComponentIndex comp_idx) const
    {
      static_assert(Dims == 1,
                    "Map::RangeIterator::value(): incorrect number of indexes "
                    "for the component dimensionality.");
      static_assert(std::is_integral<ComponentIndex>::value,
                    "Map::RangeIterator::value(): index must be an integral "
                    "type.");
      return m_currRange[comp_idx];
    }
    template <typename... ComponentIndex>
    AXOM_HOST_DEVICE DataRefType value(ComponentIndex... comp_idx) const
    {
      static_assert(sizeof...(ComponentIndex) == Dims,
                    "Map::RangeIterator::value(): incorrect number of indexes "
                    "for the component dimensionality.");
      static_assert(axom::detail::all_types_are_integral<ComponentIndex...>::value,
                    "Map::RangeIterator::value(...): index parameter pack must "
                    "all be integral types.");
      return m_currRange(comp_idx...);
    }
    AXOM_HOST_DEVICE value_type operator[](PositionType n) const { return *(*this + n); }

    /// \brief Returns the set element mapped by this iterator.
    SetElement index() const { return m_map->index(this->m_pos); }

    /// \brief Returns the flat index pointed to by this iterator.
    AXOM_HOST_DEVICE PositionType flatIndex() const { return this->m_pos; }

    /// \brief Returns the number of components per element in the Map.
    PositionType numComp() const { return m_map->stride(); }

  protected:
    /// Implementation of advance() as required by IteratorBase
    AXOM_HOST_DEVICE void advance(PositionType n)
    {
      this->m_pos += n;
      m_currRange = makeRange(m_map, this->m_pos);
    }

  private:
    MapConstPtr m_map {nullptr};
    value_type m_currRange;
  };

public:  // Functions related to iteration
  iterator begin() { return iterator(this, 0); }
  iterator end() { return iterator(this, size() * StridePolicyType::stride()); }
  const_iterator begin() const { return const_iterator(this, 0); }
  const_iterator end() const { return const_iterator(this, size() * StridePolicyType::stride()); }

  RangeAdapter<iterator> range() { return RangeAdapter<iterator> {begin(), end()}; }
  RangeAdapter<const_iterator> range() const
  {
    return RangeAdapter<const_iterator> {begin(), end()};
  }

  AXOM_HOST_DEVICE range_iterator set_begin() { return range_iterator(this, 0); }
  AXOM_HOST_DEVICE range_iterator set_end() { return range_iterator(this, size()); }
  AXOM_HOST_DEVICE const_range_iterator set_begin() const { return const_range_iterator(this, 0); }
  AXOM_HOST_DEVICE const_range_iterator set_end() const
  {
    return const_range_iterator(this, size());
  }
  RangeAdapter<range_iterator> set_elements()
  {
    return RangeAdapter<range_iterator> {set_begin(), set_end()};
  }
  RangeAdapter<const_range_iterator> set_elements() const
  {
    return RangeAdapter<const_range_iterator> {set_begin(), set_end()};
  }

public:
  /// \brief Returns a reference to the underlying map data
  OrderedMap& data() { return m_data; }
  const OrderedMap& data() const { return m_data; }

private:
  inline void verifyPosition(PositionType idx) const { verifyPositionImpl(idx); }

  inline void verifyPosition(PositionType setIdx, PositionType compIdx) const
  {
    verifyPositionImpl(setIdx, compIdx);
  }

  inline void verifyPositionImpl(PositionType AXOM_DEBUG_PARAM(idx)) const
  {
    SLIC_ASSERT_MSG(
      idx >= 0 && idx < PositionType(m_data.size()),
      "Attempted to access element " << idx << " but map's data has size " << m_data.size());
  }

  template <typename ComponentIndex>
  inline void verifyPositionImpl(PositionType AXOM_DEBUG_PARAM(setIdx),
                                 ComponentIndex AXOM_DEBUG_PARAM(compIdx)) const
  {
    SLIC_ASSERT_MSG(setIdx >= 0 && setIdx < size() && compIdx >= 0 && compIdx < numComp(),
                    "Attempted to access element at (" << setIdx << "," << compIdx
                                                       << ",) but map's set has size " << size()
                                                       << " with " << numComp() << " components.");
  }

  template <typename... ComponentIndex>
  inline void verifyPositionImpl(PositionType AXOM_DEBUG_PARAM(setIdx),
                                 ComponentIndex... AXOM_DEBUG_PARAM(compIdx)) const
  {
#ifdef AXOM_DEBUG
    const PositionType indexArray[] {static_cast<PositionType>(compIdx)...};
    PositionType shapeArray[StridePolicyType::NumDims];
    bool validIndexes = true;
    for(int dim = 0; dim < StridePolicyType::NumDims; dim++)
    {
      shapeArray[dim] = this->shape()[dim];
      validIndexes = validIndexes && (indexArray[dim] >= 0);
      validIndexes = validIndexes && (indexArray[dim] < this->shape()[dim]);
    }
    std::string invalid_message = fmt::format(
      "Attempted to access element at ({}, {}) but map's set has size {} with "
      "component shape ({})",
      setIdx,
      fmt::join(indexArray, ", "),
      size(),
      fmt::join(shapeArray, ", "));
    SLIC_ASSERT_MSG(setIdx >= 0 && setIdx < size() && validIndexes, invalid_message);
#endif
  }

  template <typename ComponentIndex>
  AXOM_HOST_DEVICE inline PositionType componentOffset(ComponentIndex componentIndex) const
  {
    return componentIndex;
  }

  template <typename... ComponentIndex>
  AXOM_HOST_DEVICE inline PositionType componentOffset(ComponentIndex... componentIndex) const
  {
    const PositionType indexArray[] {static_cast<PositionType>(componentIndex)...};
    ElementShape strides = StridePolicyType::strides();
    PositionType offset = 0;
    for(int dim = 0; dim < StridePolicyType::NumDims; dim++)
    {
      offset += indexArray[dim] * strides[dim];
    }
    return offset;
  }

  // setStride function should not be called after constructor is called.
  // This (should) override the StridePolicy setStride(s) function.
  void setStride(PositionType AXOM_UNUSED_PARAM(str))
  {
    SLIC_ASSERT_MSG(false, "Stride should not be changed after construction of map.");
  }

  PositionType requiredStorageSize() const
  {
    return detail::checkedMapStorageSize(size(), StridePolicyType::stride());
  }

  // If we can resize the underlying buffer, do so if the buffer is not large
  // enough to correspond to the size of the set.
  void checkBackingSize(std::true_type)
  {
    const PositionType neededSize = requiredStorageSize();
    SLIC_ERROR_IF(
      !std::in_range<PositionType>(m_data.size()) || !std::in_range<axom::IndexType>(m_data.size()),
      "SLAM map backing buffer size is not representable.");
    if(std::cmp_less(m_data.size(), neededSize))
    {
      m_data.resize(neededSize);
    }
  }

  void checkBackingSize(std::false_type)
  {
    const PositionType neededSize = requiredStorageSize();
    SLIC_ERROR_IF(!std::cmp_equal(m_data.size(), neededSize),
                  "SLAM map backing view must contain exactly " << neededSize << " elements.");
  }

  AXOM_HOST_DEVICE typename IndirectionPolicy::ConstResultPtr data_ptr() const
  {
    return IndirectionPolicy::getConstIndirection(m_data);
  }

  AXOM_HOST_DEVICE typename IndirectionPolicy::ResultPtr data_ptr()
  {
    return IndirectionPolicy::getIndirection(m_data);
  }

private:
  SetContainer<> m_set;
  OrderedMap m_data;
};

template <typename T, typename S, typename IndPol, typename StrPol, typename IfacePol>
  requires detail::MapParameters<T, S, IndPol, StrPol>
bool Map<T, S, IndPol, StrPol, IfacePol>::isValid(bool verboseOutput) const
{
  PositionType expected {};
  if(!detail::mapStorageSize(size(), StridePolicyType::stride(), expected))
  {
    SLIC_INFO_IF(verboseOutput,
                 "Map has an invalid component count or unrepresentable storage size.");
    return false;
  }
  bool bValid = true;

  std::stringstream errStr;

  if(policies::EmptySetTraits<S>::isEmpty(m_set.get()))
  {
    if(!m_data.empty())
    {
      if(verboseOutput)
      {
        errStr << "\n\t* the underlying set was never provided,"
               << " but its associated data is not empty"
               << " , data has size " << m_data.size();
      }

      bValid = false;
    }
  }
  else
  {
    if(!std::cmp_equal(m_data.size(), expected))
    {
      if(verboseOutput)
      {
        errStr << "\n\t* the underlying set and its associated mapped data"
               << " have different sizes"
               << " , underlying set has size " << m_set.get()->size() << " with stride "
               << StridePolicyType::stride() << " , data has size " << m_data.size();
      }

      bValid = false;
    }
  }

  if(verboseOutput && !bValid)
  {
    std::cout << "\n*** Detailed results of isValid on the map.\n"
              << "Map was NOT valid.\n"
              << errStr.str() << std::endl;
  }

  return bValid;
}

template <typename T, typename S, typename IndPol, typename StrPol, typename IfacePol>
  requires detail::MapParameters<T, S, IndPol, StrPol>
void Map<T, S, IndPol, StrPol, IfacePol>::print() const
{
  bool valid = isValid(true);
  std::stringstream sstr;

  if(valid)
  {
    if(!m_set.get())
    {
      sstr << "** map is empty.";
    }
    else
    {
      sstr << "** underlying set has size " << m_set.get()->size() << ": ";
      sstr << "\n** the stride of the map is " << StridePolicyType::stride() << ": ";

      sstr << "\n** Mapped data:";
      for(PositionType idx = 0; idx < this->size(); ++idx)
      {
        for(PositionType idx2 = 0; idx2 < StridePolicyType::stride(); ++idx2)
        {
          sstr << "\n\telt[" << idx << "," << idx2 << "]:\t"
               << (*this)[idx * StridePolicyType::stride() + idx2];
        }
      }
    }
  }

  std::cout << sstr.str() << std::endl;
}

}  // end namespace axom::slam
