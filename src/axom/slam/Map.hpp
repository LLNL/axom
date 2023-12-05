// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
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

  using ValueType = typename IndirectionPolicy::IndirectionResult;
  using ConstValueType = typename IndirectionPolicy::ConstIndirectionResult;

  class MapBuilder;

  // types for iterator
  class MapIterator;
  using const_iterator = MapIterator;
  using const_iterator_pair = std::pair<const_iterator, const_iterator>;
  using iterator = const_iterator;
  using iterator_pair = const_iterator_pair;

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
   * \brief Constructor for Map using a Set pointer
   *
   * \param theSet         (Optional) A pointer to the map's set
   * \param defaultValue   (Optional) If given, every entry in the map will be
   *                       initialized using defaultValue
   * \param stride  (Optional) The stride. The number of DataType that
   *                each element in the set will be mapped to.
   *                When using a \a RuntimeStridePolicy, the default is 1.
   * \note  When using a compile time StridePolicy, \a stride must be equal to
   *        \a stride(), when provided.
   */

  Map(const SetType* theSet = policies::EmptySetTraits<SetType>::emptySet(),
      DataType defaultValue = DataType(),
      SetPosition stride = StridePolicyType::DEFAULT_VALUE,
      int allocatorID = axom::getDefaultAllocatorID())
    : StridePolicyType(stride)
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
      SetPosition stride = StridePolicyType::DEFAULT_VALUE,
      int allocatorID = axom::getDefaultAllocatorID())
    : StridePolicyType(stride)
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
   * \param stride  (Optional) The stride. The number of DataType that
   *                each element in the set will be mapped to.
   *                When using a \a RuntimeStridePolicy, the default is 1.
   * \note  When using a compile time StridePolicy, \a stride must be equal to
   *        \a stride(), when provided.
   */
  template <typename USet,
            typename TSet = SetType,
            typename Enable = typename std::enable_if<
              !std::is_abstract<TSet>::value && std::is_base_of<TSet, USet>::value>::type>
  Map(const USet& theSet,
      OrderedMap data,
      SetPosition stride = StridePolicyType::DEFAULT_VALUE)
    : StridePolicyType(stride)
    , m_set(theSet)
    , m_data(std::move(data))
  {
    static_assert(std::is_same<SetType, USet>::value,
                  "Argument set is of a more-derived type than the Map's set "
                  "type. This may lead to object slicing. Use Map's pointer "
                  "constructor instead to store polymorphic sets.");

    checkBackingSize(
      std::integral_constant<bool, IndirectionPolicy::IsMutableBuffer> {});
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
  const SetType* set() const { return m_set.get(); }

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
   * \brief Access the value associated with the given position in the set and
   *        the component index.
   *
   * \pre `0 <= setIdx < size()`
   * \pre `0 <= comp < numComp()`
   */
  AXOM_HOST_DEVICE ConstValueType operator()(SetPosition setIdx,
                                             SetPosition comp = 0) const
  {
#ifndef AXOM_DEVICE_CODE
    verifyPositionImpl(setIdx, comp);
#endif
    SetPosition setIndex = setIdx * StridePolicyType::stride() + comp;
    return IndirectionPolicy::getConstIndirection(m_data, setIndex);
  }

  AXOM_HOST_DEVICE ValueType operator()(SetPosition setIdx, SetPosition comp = 0)
  {
#ifndef AXOM_DEVICE_CODE
    verifyPositionImpl(setIdx, comp);
#endif
    SetPosition setIndex = setIdx * StridePolicyType::stride() + comp;
    return IndirectionPolicy::getIndirection(m_data, setIndex);
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

  /**
   * \class   MapIterator
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
  class MapIterator : public IteratorBase<MapIterator, SetPosition>
  {
  public:
    using iterator_category = std::random_access_iterator_tag;
    using value_type = DataType;
    using difference_type = SetPosition;

    using IterBase = IteratorBase<MapIterator, SetPosition>;
    using iter = MapIterator;
    using PositionType = SetPosition;
    using IterBase::m_pos;

  public:
    MapIterator(PositionType pos, Map* oMap) : IterBase(pos), m_mapPtr(oMap) { }

    /**
     * \brief Returns the current iterator value. If the map has multiple
     *        components, this will return the first component. To access
     *        the other components, use iter(comp)
     */
    DataType& operator*() { return (*m_mapPtr)(m_pos, 0); }

    /**
     * \brief Returns the iterator's value at the specified component.
     *        Returns the first component if comp_idx is not specified.
     * \param comp_idx  (Optional) Zero-based index of the component.
     */
    DataType& operator()(SetPosition comp_idx = 0)
    {
      return (*m_mapPtr)(m_pos, comp_idx);
    }

    /** \brief Returns the first component value after n increments.  */
    const DataType& operator[](PositionType n) const { return *(*this + n); }

    DataType& operator[](PositionType n) { return *(*this + n); }

    /** \brief Returns the number of components per element in the Map. */
    PositionType numComp() const { return m_mapPtr->stride(); }

  protected:
    /** Implementation of advance() as required by IteratorBase */
    void advance(PositionType n) { m_pos += n; }

  protected:
    Map* m_mapPtr;
  };

public:  // Functions related to iteration
  MapIterator begin() { return MapIterator(0, this); }
  MapIterator end() { return MapIterator(size(), this); }
  const_iterator_pair range() const { return std::make_pair(begin(), end()); }

public:
  /**
   * \brief Returns a reference to the underlying map data
   */
  OrderedMap& data() { return m_data; }
  const OrderedMap& data() const { return m_data; }

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
