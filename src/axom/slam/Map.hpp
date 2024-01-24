// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file Map.hpp
 *
 * \brief Basic API for a map from each element of a set to some domain
 *
 */

#ifndef SLAM_MAP_HPP_
#define SLAM_MAP_HPP_

#include <vector>
#include <sstream>
#include <iostream>

#include "axom/core.hpp"
#include "axom/slic.hpp"

#include "axom/slam/MapBase.hpp"
#include "axom/slam/Set.hpp"
#include "axom/slam/NullSet.hpp"

#include "axom/core/IteratorBase.hpp"
#include "axom/core/RangeAdapter.hpp"

#include "axom/slam/policies/StridePolicies.hpp"
#include "axom/slam/policies/IndirectionPolicies.hpp"
#include "axom/slam/policies/PolicyTraits.hpp"
#include "axom/slam/policies/MapInterfacePolicies.hpp"

namespace axom
{
namespace slam
{
// This class is missing some simplifying copy constructors
// -- or at least ways of interacting with the data store
// We should probably support copy on write shallow copies when possible...

/**
 * \class   Map
 *
 * \brief   A Map class that associates a constant number of values to every
 *          element in a set.
 *
 * \tparam  T The data type of each value
 * \tparam  S The map's set type
 * \tparam  IndPol The map's indirection policy
 * \tparam  StrPol A policy class that determines how many values to
 *          associate with each element. There is a fixed \a stride between
 *          the data associated with each element of the set.
 * \details The map class associates a fixed number of values, also referred
 *          to as \a components, with each element in its underlying \a Set
 *          instance. Depending on the \a StridePolicy, this can be fixed at
 *          compile time or at runtime.\n
 *          Access to the j<sup>th</sup> component of the i<sup>th</sup>
 *          element, can be obtained via the parenthesis operator
 *          ( `map(i,j)` ), or via the square bracket operator
 *          (i.e. `map[k]`, where `k = i * stride() + j` ).
 *
 */

template <typename T,
          typename S = Set<>,
          typename IndPol =
            policies::STLVectorIndirection<typename S::PositionType, T>,
          typename StrPol = policies::StrideOne<typename S::PositionType>,
          typename IfacePol = policies::VirtualInterface>
class Map : public StrPol,
            public policies::MapInterface<IfacePol, typename S::PositionType>
{
public:
  using DataType = T;
  using SetType = S;
  using IndirectionPolicy = IndPol;
  using StridePolicyType = StrPol;

  using OrderedMap = typename IndirectionPolicy::IndirectionBufferType;

  using SetPosition = typename SetType::PositionType;
  using SetElement = typename SetType::ElementType;
  static const NullSet<SetPosition, SetElement> s_nullSet;

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
    AXOM_HOST_DEVICE SetContainer(const USet* set) : m_pSet(set) { }

    AXOM_HOST_DEVICE const USet* get() const { return m_pSet; }

    const USet* m_pSet;
  };

  template <typename USet>
  struct SetContainer<USet, true>
  {
    AXOM_HOST_DEVICE SetContainer(const USet* set) : m_pSet(set) { }
    AXOM_HOST_DEVICE SetContainer(const USet& set) : m_set(set) { }

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
   * \brief Constructor for Map using a Set pointer
   *
   * \param theSet         (Optional) A pointer to the map's set
   * \param defaultValue   (Optional) If given, every entry in the map will be
   *                       initialized using defaultValue
   * \param shape   (Optional) The number of DataType that each element in the
   *                set will be mapped to.
   *                When using a \a RuntimeStridePolicy, the default is 1.
   * \note  When using a compile time StridePolicy, \a stride must be equal to
   *        \a stride(), when provided.
   */

  Map(const SetType* theSet = policies::EmptySetTraits<SetType>::emptySet(),
      DataType defaultValue = DataType(),
      ElementShape shape = StridePolicyType::DefaultSize(),
      int allocatorID = axom::getDefaultAllocatorID())
    : StridePolicyType(shape)
    , m_set(theSet)
  {
    m_data =
      IndirectionPolicy::create(size() * numComp(), defaultValue, allocatorID);
  }

  /// \overload
  template <typename USet,
            typename TSet = SetType,
            typename Enable = typename std::enable_if<
              !std::is_abstract<TSet>::value && std::is_base_of<TSet, USet>::value>::type>
  Map(const USet& theSet,
      DataType defaultValue = DataType(),
      ElementShape shape = StridePolicyType::DefaultSize(),
      int allocatorID = axom::getDefaultAllocatorID())
    : StridePolicyType(shape)
    , m_set(theSet)
  {
    static_assert(std::is_same<SetType, USet>::value,
                  "Argument set is of a more-derived type than the Map's set "
                  "type. This may lead to object slicing. Use Map's pointer "
                  "constructor instead to store polymorphic sets.");
    m_data =
      IndirectionPolicy::create(size() * numComp(), defaultValue, allocatorID);
  }

  /**
   * \brief Constructor for Map using a Set passed by-value and data passed in
   *        by-value.
   *
   * \param theSet  A reference to the map's set
   * \param data    Pointer to the externally-owned data
   * \param shape   (Optional) The number of DataType that each element in the
   *                set will be mapped to.
   *                When using a \a RuntimeStridePolicy, the default is 1.
   * \note  When using a compile time StridePolicy, \a stride must be equal to
   *        \a stride(), when provided.
   */
  template <typename USet,
            typename TSet = SetType,
            typename Enable = typename std::enable_if<
              !std::is_abstract<TSet>::value && std::is_base_of<TSet, USet>::value>::type>
  AXOM_HOST_DEVICE Map(const USet& theSet,
                       OrderedMap data,
                       ElementShape shape = StridePolicyType::DefaultSize())
    : StridePolicyType(shape)
    , m_set(theSet)
    , m_data(std::move(data))
  {
    static_assert(std::is_same<SetType, USet>::value,
                  "Argument set is of a more-derived type than the Map's set "
                  "type. This may lead to object slicing. Use Map's pointer "
                  "constructor instead to store polymorphic sets.");

#ifndef AXOM_DEVICE_CODE
    checkBackingSize(
      std::integral_constant<bool, IndirectionPolicy::IsMutableBuffer> {});
#endif
  }

  /**
   * \brief Constructor for Map using a MapBuilder
   */
  Map(const MapBuilder& builder)
    : Map(builder.m_set, builder.m_defaultValue, builder.m_stride.stride())
  {
    //copy the data if exists
    if(builder.m_data_ptr)
    {
      for(SetPosition idx = SetPosition(); idx < builder.m_set->size(); ++idx)
      {
        m_data[idx] = builder.m_data_ptr[idx];
      }
    }
  }

  /**
   * \brief Returns a pointer to the map's underlying set
   */
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
  AXOM_HOST_DEVICE ConstValueType operator[](SetPosition setIndex) const
  {
#ifndef AXOM_DEVICE_CODE
    verifyPositionImpl(setIndex);
#endif
    return IndirectionPolicy::getConstIndirection(m_data, setIndex);
  }

  AXOM_HOST_DEVICE ValueType operator[](SetPosition setIndex)
  {
#ifndef AXOM_DEVICE_CODE
    verifyPositionImpl(setIndex);
#endif
    return IndirectionPolicy::getIndirection(m_data, setIndex);
  }

  /**
   * \brief Access the value associated with the given position in the set.
   */
  AXOM_HOST_DEVICE ConstValueType operator()(SetPosition setIdx) const
  {
    // TODO: validate that runtime stride is 1-D with value 1?
    return value(setIdx, 0);
  }

  /// \overload
  AXOM_HOST_DEVICE ValueType operator()(SetPosition setIdx)
  {
    return value(setIdx, 0);
  }

  /**
   * \brief Access the value associated with the given position in the set and
   *        the component index.
   *
   * \pre `0 <= setIdx < size()`
   * \pre `0 <= compIdx[i] < shape()[i]`
   */
  template <typename... ComponentPos>
  AXOM_HOST_DEVICE ConstValueType operator()(SetPosition setIdx,
                                             ComponentPos... compIdx) const
  {
    return value(setIdx, compIdx...);
  }

  /// \overload
  template <typename... ComponentPos>
  AXOM_HOST_DEVICE ValueType operator()(SetPosition setIdx,
                                        ComponentPos... compIdx)
  {
    return value(setIdx, compIdx...);
  }

  AXOM_HOST_DEVICE ConstValueType value(SetPosition setIdx) const
  {
    return value(setIdx, 0);
  }

  AXOM_HOST_DEVICE ValueType value(SetPosition setIdx)
  {
    return value(setIdx, 0);
  }

  /**
   * \brief Access the value associated with the given position in the set and
   *        the component index.
   *
   * \pre `0 <= setIdx < size()`
   * \pre `0 <= compIdx[i] < shape()[i]`
   */
  template <typename... ComponentPos>
  AXOM_HOST_DEVICE ConstValueType value(SetPosition setIdx,
                                        ComponentPos... compIdx) const
  {
    static_assert(
      sizeof...(ComponentPos) == StridePolicyType::NumDims,
      "Invalid number of components provided for given Map's StridePolicy");
#ifndef AXOM_DEVICE_CODE
    verifyPositionImpl(setIdx, compIdx...);
#endif
    SetPosition elemIndex = setIdx * StridePolicyType::stride();
    elemIndex += componentOffset(compIdx...);
    return IndirectionPolicy::getConstIndirection(m_data, elemIndex);
  }

  /// \overload
  template <typename... ComponentPos>
  AXOM_HOST_DEVICE ValueType value(SetPosition setIdx, ComponentPos... compIdx)
  {
    static_assert(
      sizeof...(ComponentPos) == StridePolicyType::NumDims,
      "Invalid number of components provided for given Map's StridePolicy");
#ifndef AXOM_DEVICE_CODE
    verifyPositionImpl(setIdx, compIdx...);
#endif
    SetPosition elemIndex = setIdx * StridePolicyType::stride();
    elemIndex += componentOffset(compIdx...);
    return IndirectionPolicy::getIndirection(m_data, elemIndex);
  }

  AXOM_HOST_DEVICE SetElement index(IndexType idx) const
  {
    return set()->at(idx);
  }

  /// @}

  /// \name Map cardinality functions
  /// @{

  /**
   * \brief Return the set size. Same as `set()->size()`.
   *
   * The total storage size for the map's values is `size() * numComp()`
   */
  AXOM_HOST_DEVICE SetPosition size() const
  {
    return !policies::EmptySetTraits<SetType>::isEmpty(m_set.get())
      ? static_cast<SetPosition>(m_set.get()->size())
      : SetPosition(0);
  }

  /*
   * \brief  Gets the number of component values associated with each element.
   *         Equivalent to stride().
   */
  SetPosition numComp() const { return StridePolicyType::stride(); }

  /**
   * \brief Returns the shape of the component values associated with each
   *  element.
   *
   *  For one-dimensional strides, equivalent to stride(); otherwise, returns
   *  an N-dimensional array with the number of values in each sub-component
   *  index.
   */
  ElementShape shape() const { return StridePolicyType::shape(); }

  /// @}

  /// \name Map modifying functions for all entries
  /// @{

  /** \brief replace all elements in the Map with the default DataType */
  void clear() { fill(); }

  /** Set each entry in the map to the given value  */
  void fill(DataType val = DataType())
  {
    const SetPosition sz = static_cast<SetPosition>(m_data.size());

    for(SetPosition idx = SetPosition(); idx < sz; ++idx)
    {
      m_data[idx] = val;
    }
  }

  /** \brief Element-wise copy of data from another map */
  void copy(const Map& other)
  {
    SLIC_ASSERT(other.size() == size());
    SLIC_ASSERT(other.stride() == StridePolicyType::stride());

    const SetPosition sz = size() * StridePolicyType::stride();
    for(SetPosition idx = SetPosition(); idx < sz; ++idx)
    {
      m_data[idx] = other[idx];
    }
  }

  ///@}

  /** \brief print information on the map, including every element inside Map
   */
  void print() const;

  /** \brief returns true if the map is valid, false otherwise.  */
  bool isValid(bool verboseOutput = false) const;

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

    /** \brief Provide the Set to be used by the Map */
    MapBuilder& set(const SetType* set)
    {
      m_set = set;
      return *this;
    }

    /** \brief Set the stride of the Map using StridePolicy */
    MapBuilder& stride(SetPosition str)
    {
      m_stride = StridePolicyType(str);
      return *this;
    }

    /** \brief Set the pointer to the array of data the Map will contain
     *  (makes a copy of the array currently)
     */
    MapBuilder& data(DataType* bufPtr)
    {
      m_data_ptr = bufPtr;
      return *this;
    }

  private:
    const SetType* m_set;
    StridePolicyType m_stride;
    DataType* m_data_ptr = nullptr;
    DataType m_defaultValue = DataType();
  };

private:
  /*!
   * \class MapIteratorStorage
   *
   * \brief Helper class to handle storage for iterators.
   *
   *  When the Map is backed by a mutable buffer (std::vector or axom::Array),
   *  access to the underlying Map is provided via a pointer. When the map is
   *  backed by a view-like object (axom::ArrayView), we just keep a copy of the
   *  Map by-value.
   *
   *  This is necessary to allow BivariateMap to return subset iterators
   *  pointing to a temporary map. The original SubMap iterator would just keep
   *  the SubMap by-value.
   */
  template <bool Const, bool IsView>
  struct MapIteratorStorage;

  template <bool Const>
  struct MapIteratorStorage<Const, false>
  {
    using MapConstType = std::conditional_t<Const, const Map*, Map*>;
    using MapRefType = std::conditional_t<Const, const Map&, Map&>;
    AXOM_HOST_DEVICE MapIteratorStorage(MapConstType pMap) : m_map(pMap) { }

    AXOM_HOST_DEVICE MapRefType map() const { return *m_map; }

    MapConstType m_map;
  };

  template <bool Const>
  struct MapIteratorStorage<Const, true>
  {
    using MapConstType = std::conditional_t<Const, const Map*, Map*>;
    AXOM_HOST_DEVICE MapIteratorStorage(MapConstType pMap) : m_map(*pMap) { }

    AXOM_HOST_DEVICE const Map& map() const { return m_map; }

    Map m_map;
  };

public:
  /**
   * \class MapIterator
   * \brief An iterator type for a map. Each increment operation advances the
   *        iterator to the element at the next flat index.
   */
  template <bool Const>
  class MapIterator
    : public IteratorBase<MapIterator<Const>, SetPosition>,
      MapIteratorStorage<Const, !IndirectionPolicy::IsMutableBuffer>
  {
  public:
    using DataRefType = std::conditional_t<Const, ConstValueType, ValueType>;
    using DataType = std::remove_reference_t<DataRefType>;

    using iterator_category = std::random_access_iterator_tag;
    using value_type = DataType;
    using reference_type = DataRefType;
    using pointer_type = value_type*;
    using difference_type = SetPosition;

    using IterBase = IteratorBase<MapIterator, SetPosition>;
    using MapStorage =
      MapIteratorStorage<Const, !IndirectionPolicy::IsMutableBuffer>;
    using iter = MapIterator;
    using PositionType = SetPosition;
    using IterBase::m_pos;

  public:
    MapIterator(PositionType pos, Map* oMap) : IterBase(pos), MapStorage(oMap)
    { }

    /**
     * \brief Returns the current iterator value.
     */
    reference_type operator*() { return this->map()[m_pos]; }

    /// \brief Returns the set element mapped by this iterator.
    SetElement index() const
    {
      return this->map().index(this->m_pos / this->map().numComp());
    }

    /// \brief Returns the component index pointed to by this iterator.
    PositionType compIndex() const { return m_pos % this->map().numComp(); }

    /// \brief Returns the flat index pointed to by this iterator.
    SetPosition flatIndex() const { return this->m_pos; }

  protected:
    /** Implementation of advance() as required by IteratorBase */
    void advance(PositionType n) { m_pos += n; }
  };

  /**
   * \class   MapRangeIterator
   * \brief   An iterator type for a map.
   *          Each increment operation advances the iterator to the next set
   *          element.
   *          To access the j<sup>th</sup> component values of the iterator's
   *          current element, use `iter(j)`.
   * \warning Note the difference between the subscript operator ( `iter[off]` )
   *          and the parenthesis operator ( `iter(j)` ). \n
   *          `iter[off]` returns the value of the first component of the
   *          element at offset \a `off` from the currently pointed to
   *          element.\n
   *          And `iter(j)` returns the value of the j<sup>th</sup> component of
   *          the currently pointed to element (where 0 <= j < numComp()).\n
   *          For example: `iter[off]` is the same as `(iter+off)(0)`
   */
  template <bool Const>
  class MapRangeIterator
    : public IteratorBase<MapRangeIterator<Const>, SetPosition>,
      MapIteratorStorage<Const, !IndirectionPolicy::IsMutableBuffer>
  {
  public:
    using IterBase = IteratorBase<MapRangeIterator, SetPosition>;
    using MapStorage =
      MapIteratorStorage<Const, !IndirectionPolicy::IsMutableBuffer>;

    using DataRefType = std::conditional_t<Const, ConstValueType, ValueType>;
    using DataType = std::remove_reference_t<DataRefType>;

    using PositionType = SetPosition;
    using IterBase::m_pos;

    using iterator_category = std::random_access_iterator_tag;
    using value_type = axom::ArrayView<DataType>;
    using reference_type = value_type&;
    using pointer_type = value_type*;
    using difference_type = SetPosition;

  public:
    MapRangeIterator(PositionType pos, Map* oMap)
      : IterBase(pos)
      , MapStorage(oMap)
    { }

    /**
     * \brief Returns the current iterator value.
     */
    value_type operator*()
    {
      return value_type {this->map().data().data() + m_pos * this->map().stride(),
                         this->map().stride()};
    }

    /**
     * \brief Returns the iterator's value at the specified component.
     *        Returns the first component if comp_idx is not specified.
     * \param comp_idx  Zero-based index of the component.
     */
    DataRefType operator()(SetPosition comp_idx) const
    {
      return value(comp_idx);
    }
    DataRefType value(PositionType comp_idx) const
    {
      return this->map()(m_pos, comp_idx);
    }
    value_type operator[](PositionType n) const { return *(*this + n); }

    /// \brief Returns the set element mapped by this iterator.
    SetElement index() const { return this->map().index(this->m_pos); }

    /// \brief Returns the flat index pointed to by this iterator.
    SetPosition flatIndex() const { return this->m_pos; }

    /** \brief Returns the number of components per element in the Map. */
    PositionType numComp() const { return this->map().stride(); }

  protected:
    /** Implementation of advance() as required by IteratorBase */
    void advance(PositionType n) { m_pos += n; }
  };

public:  // Functions related to iteration
  iterator begin() { return iterator(0, this); }
  iterator end() { return iterator(size() * StridePolicyType::stride(), this); }
  const_iterator begin() const { return const_iterator(0, this); }
  const_iterator end() const
  {
    return const_iterator(size() * StridePolicyType::stride(), this);
  }

  RangeAdapter<iterator> range() const
  {
    return RangeAdapter<iterator> {begin(), end()};
  }

  range_iterator set_begin() { return range_iterator(0, this); }
  range_iterator set_end() { return range_iterator(size(), this); }
  const_range_iterator set_begin() const
  {
    return const_range_iterator(0, this);
  }
  const_range_iterator set_end() const
  {
    return const_range_iterator(size(), this);
  }
  RangeAdapter<range_iterator> set_elements()
  {
    return RangeAdapter<range_iterator> {set_begin(), set_end()};
  }

public:
  /**
   * \brief Returns a reference to the underlying map data
   */
  AXOM_HOST_DEVICE OrderedMap& data() { return m_data; }
  AXOM_HOST_DEVICE const OrderedMap& data() const { return m_data; }

private:
  inline void verifyPosition(SetPosition idx) const { verifyPositionImpl(idx); }

  inline void verifyPosition(SetPosition setIdx, SetPosition compIdx) const
  {
    verifyPositionImpl(setIdx, compIdx);
  }

  inline void verifyPositionImpl(SetPosition AXOM_DEBUG_PARAM(idx)) const
  {
    SLIC_ASSERT_MSG(idx >= 0 && idx < SetPosition(m_data.size()),
                    "Attempted to access element "
                      << idx << " but map's data has size " << m_data.size());
  }

  inline void verifyPositionImpl(SetPosition AXOM_DEBUG_PARAM(setIdx),
                                 SetPosition AXOM_DEBUG_PARAM(compIdx)) const
  {
    SLIC_ASSERT_MSG(
      setIdx >= 0 && setIdx < size() && compIdx >= 0 && compIdx < numComp(),
      "Attempted to access element at ("
        << setIdx << "," << compIdx << ",) but map's set has size " << size()
        << " with " << numComp() << " components.");
  }

  template <typename... ComponentIndex>
  inline void verifyPositionImpl(SetPosition AXOM_DEBUG_PARAM(setIdx),
                                 ComponentIndex... AXOM_DEBUG_PARAM(compIdx)) const
  {
#ifdef AXOM_DEBUG
    ElementShape indexArray {{compIdx...}};
    bool validIndexes = true;
    for(int dim = 0; dim < StridePolicyType::NumDims; dim++)
    {
      validIndexes = validIndexes && (indexArray[dim] >= 0);
      validIndexes = validIndexes && (indexArray[dim] < this->shape()[dim]);
    }
    std::string invalid_message = fmt::format(
      "Attempted to access element at ({}, {}) but map's set has size {} with "
      "component shape ({})",
      setIdx,
      fmt::join(indexArray, ", "),
      size(),
      fmt::join(this->shape(), ", "));
    SLIC_ASSERT_MSG(setIdx >= 0 && setIdx < size() && validIndexes,
                    invalid_message);
#endif
  }

  AXOM_HOST_DEVICE inline SetPosition componentOffset(SetPosition componentIndex) const
  {
    return componentIndex;
  }

  template <typename... ComponentIndex>
  AXOM_HOST_DEVICE inline SetPosition componentOffset(
    ComponentIndex... componentIndex) const
  {
    ElementShape indexArray {{componentIndex...}};
    ElementShape strides = StridePolicyType::strides();
    SetPosition offset = 0;
    for(int dim = 0; dim < StridePolicyType::NumDims; dim++)
    {
      offset += indexArray[dim] * strides[dim];
    }
    return offset;
  }

  // setStride function should not be called after constructor is called.
  // This (should) override the StridePolicy setStride(s) function.
  void setStride(SetPosition AXOM_UNUSED_PARAM(str))
  {
    SLIC_ASSERT_MSG(false,
                    "Stride should not be changed after construction of map.");
  }

  // If we can resize the underlying buffer, do so if the buffer is not large
  // enough to correspond to the size of the set.
  void checkBackingSize(std::true_type)
  {
    const IndexType neededSize = size() * numComp();
    if(m_data.size() < neededSize)
    {
      m_data.resize(neededSize);
    }
  }

  void checkBackingSize(std::false_type)
  {
    SLIC_ASSERT_MSG(m_data.size() == size() * numComp(),
                    "Not enough elements in buffer passed to Map constructor.");
  }

private:
  SetContainer<> m_set;
  OrderedMap m_data;
};

template <typename T, typename S, typename IndPol, typename StrPol, typename IfacePol>
bool Map<T, S, IndPol, StrPol, IfacePol>::isValid(bool verboseOutput) const
{
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
    if(static_cast<SetPosition>(m_data.size()) !=
       m_set.get()->size() * StridePolicyType::stride())
    {
      if(verboseOutput)
      {
        errStr << "\n\t* the underlying set and its associated mapped data"
               << " have different sizes"
               << " , underlying set has size " << m_set.get()->size()
               << " with stride " << StridePolicyType::stride()
               << " , data has size " << m_data.size();
      }

      bValid = false;
    }
  }

  if(verboseOutput)
  {
    std::stringstream sstr;

    sstr << "\n*** Detailed results of isValid on the map.\n";
    if(bValid)
    {
      sstr << "Map was valid." << std::endl;
    }
    else
    {
      sstr << "Map was NOT valid.\n" << sstr.str() << std::endl;
    }

    std::cout << sstr.str() << std::endl;
  }

  return bValid;
}

template <typename T, typename S, typename IndPol, typename StrPol, typename IfacePol>
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
      sstr << "\n** the stride of the map is " << StridePolicyType::stride()
           << ": ";

      sstr << "\n** Mapped data:";
      for(SetPosition idx = 0; idx < this->size(); ++idx)
      {
        for(SetPosition idx2 = 0; idx2 < StridePolicyType::stride(); ++idx2)
        {
          sstr << "\n\telt[" << idx << "," << idx2 << "]:\t"
               << (*this)[idx * StridePolicyType::stride() + idx2];
        }
      }
    }
  }

  std::cout << sstr.str() << std::endl;
}

}  // end namespace slam
}  // end namespace axom

#endif  // SLAM_MAP_HPP_
