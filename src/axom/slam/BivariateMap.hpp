// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file BivariateMap.hpp
 *
 * \brief Contains the BivariateMap class, a map for a BivariateSet
 *
 */

#ifndef SLAM_BIVARIATE_MAP_HPP_
#define SLAM_BIVARIATE_MAP_HPP_

#include "axom/slam/Map.hpp"
#include "axom/slam/Relation.hpp"
#include "axom/slam/BivariateSet.hpp"
#include "axom/slam/SubMap.hpp"
#include "axom/slam/Set.hpp"
#include "axom/slam/policies/StridePolicies.hpp"
#include "axom/slam/policies/PolicyTraits.hpp"

#include <cassert>
#include <typeinfo>

namespace axom
{
namespace slam
{
/**
 * \class BivariateMap
 * \brief A Map for BivariateSet. It associates a constant number of values to
 *        every element in a BivariateSet (as determined by StridePolicy).
 *
 * \detail Like BivariateSet, every value in BivariateMap is indexed by two
 *         indices. BivariateMap's `operator(i)` returns a SubMap of all entries
 *         whose first index is `i` in the BivariateSet.
 *
 * The different indexing systems (DenseIndex, SparseIndex, FlatIndex) are
 * explained in BivariateSet. Because BivariateMap can have more than one
 * component, FlatIndex is further divided into ComponentFlatIndex, where each
 * component in each element is indexed separately, and ElementFlatIndex,
 * an index that disregards the individual components. Hence, to access
 * each component, one would need to provide a component index as well.
 *
 * Example:
 * For a 2 x 2 sparse matrix with 3 components below:
 *     \code
 *         0    1
 *     0  abc
 *     1       def
 *     \endcode
 *
 *   Access the elements using ElementFlatIndex `(e)` would be...\n
 *   `(e = 0) = abc`\n
 *   `(e = 1) = def`\n
 *   To access each component, provide a component index (c)
 *   `(e = 0, c = 0) = a`\n
 *   `(e = 0, c = 1) = b`\n
 *   `(e = 0, c = 2) = c`\n
 *   `(e = 1, c = 0) = d`\n
 *   `(e = 1, c = 1) = e`\n
 *   `(e = 1, c = 2) = f`\n
 *
 *   To access using ComponentFlatIndex `(idx)`...\n
 *   `(idx = 0) = a`\n
 *   `(idx = 1) = b`\n
 *   `(idx = 2) = c`\n
 *   `(idx = 3) = d`\n
 *   `(idx = 4) = e`\n
 *   `(idx = 5) = f`\n
 *
 * \tparam DataType the data type of each value
 * \tparam StridePolicy A policy class for configuring the number of components
 *         associate with each element. There is a fixed \a stride between
 *         the data associated with each element of the set.
 * \see BivariateSet, SubMap
 */

template <typename T,
          typename BSet = BivariateSet<>,
          typename IndPol =
            policies::STLVectorIndirection<typename BSet::PositionType, T>,
          typename StrPol = policies::StrideOne<typename BSet::PositionType>,
          typename IfacePol = policies::VirtualInterface>
class BivariateMap
  : public policies::MapInterface<IfacePol, typename BSet::PositionType>,
    public StrPol
{
public:
  using DataType = T;
  using BivariateSetType = BSet;
  using IndirectionPolicy = IndPol;
  using StridePolicyType = StrPol;

  using SetPosition = typename BSet::PositionType;
  using SetElement = typename BSet::ElementType;

  using SetType = typename slam::RangeSet<SetPosition, SetElement>::ConcreteSet;
  using MapType = Map<DataType, SetType, IndPol, StrPol, IfacePol>;
  using OrderedSetType = typename BSet::SubsetType;

  using ValueType = typename IndirectionPolicy::IndirectionResult;
  using ConstValueType = typename IndirectionPolicy::ConstIndirectionResult;
  using PointerType = std::remove_reference_t<ValueType>*;
  using ConstPointerType = std::remove_reference_t<ConstValueType>*;

  using BivariateMapType = BivariateMap<DataType, BSet, IndPol, StrPol, IfacePol>;

  template <bool Const>
  class BivariateMapIterator;
  using iterator = BivariateMapIterator<false>;
  using const_iterator = BivariateMapIterator<true>;

  using SubMapType = SubMap<BivariateMapType, SetType, IfacePol>;
  using ConstSubMapType = const SubMap<const BivariateMapType, SetType, IfacePol>;
  using SubMapIterator = typename SubMapType::iterator;
  using ConstSubMapIterator = typename ConstSubMapType::iterator;

  using NullBivariateSetType =
    NullBivariateSet<typename BSet::FirstSetType, typename BSet::SecondSetType>;

private:
  static const NullBivariateSetType s_nullBiSet;

  template <typename USet = BivariateSetType,
            bool HasValue = !std::is_abstract<USet>::value>
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
  using ConcreteMap =
    BivariateMap<T, BSet, IndPol, StrPol, policies::ConcreteInterface>;
  using VirtualMap =
    BivariateMap<T, BSet, IndPol, StrPol, policies::VirtualInterface>;

public:
  /**
   * \brief Constructor for a BivariateMap
   *
   * \param bSet          (Optional) Pointer to the BivariateSet.
   * \param defaultValue  (Optional) The default value used to initialize the
   *                      entries of the map.
   * \param stride        (Optional) The stride, or number of component, of the
   *                      map.
   *
   * \note  When using a compile time StridePolicy, \a stride must be equal to
   *        \a StridePolicy::stride(), when provided.
   */
  BivariateMap(const BivariateSetType* bSet = &s_nullBiSet,
               DataType defaultValue = DataType(),
               SetPosition stride = StridePolicyType::DEFAULT_VALUE,
               int allocatorID = axom::getDefaultAllocatorID())
    : StridePolicyType(stride)
    , m_bset(bSet)
    , m_map(SetType(bSet->size()), defaultValue, stride, allocatorID)
  { }

  /// \overload
  template <typename UBSet,
            typename Enable = typename std::enable_if<
              !std::is_abstract<BivariateSetType>::value &&
              std::is_base_of<BivariateSetType, UBSet>::value>::type>
  BivariateMap(const UBSet& bSet,
               DataType defaultValue = DataType(),
               SetPosition stride = StridePolicyType::DEFAULT_VALUE,
               int allocatorID = axom::getDefaultAllocatorID())
    : StridePolicyType(stride)
    , m_bset(bSet)
    , m_map(SetType(bSet->size()), defaultValue, stride, allocatorID)
  {
    static_assert(std::is_same<BivariateSetType, UBSet>::value,
                  "Argument set is of a more-derived type than the Map's set "
                  "type. This may lead to object slicing. Use Map's pointer "
                  "constructor instead to store polymorphic sets.");
  }

  /**
   * \brief Constructor for BivariateMap using a BivariateSet passed by-value
   *        and data passed in by-value.
   *
   * \param bSet    A reference to the map's associated bivariate set
   * \param data    The data buffer to set the map's data to.
   * \param stride  (Optional) The stride. The number of DataType that
   *                each element in the set will be mapped to.
   *                When using a \a RuntimeStridePolicy, the default is 1.
   * \note  When using a compile time StridePolicy, \a stride must be equal to
   *        \a stride(), when provided.
   */
  BivariateMap(const BivariateSetType* bSet,
               typename MapType::OrderedMap data,
               SetPosition stride = StridePolicyType::DEFAULT_VALUE)
    : StridePolicyType(stride)
    , m_bset(bSet)
    , m_map(SetType(bSet->size()), data, stride)
  { }

  /**
   * \brief Constructor for BivariateMap using a BivariateSet passed by-value
   *        and data passed in by-value.
   *
   * \param bSet    A reference to the map's associated bivariate set
   * \param data    The data buffer to set the map's data to.
   * \param stride  (Optional) The stride. The number of DataType that
   *                each element in the set will be mapped to.
   *                When using a \a RuntimeStridePolicy, the default is 1.
   * \note  When using a compile time StridePolicy, \a stride must be equal to
   *        \a stride(), when provided.
   */
  template <typename UBSet,
            typename TBSet = BivariateSetType,
            typename Enable =
              typename std::enable_if<!std::is_abstract<TBSet>::value &&
                                      std::is_base_of<TBSet, UBSet>::value>::type>
  BivariateMap(const UBSet& bSet,
               typename MapType::OrderedMap data,
               SetPosition stride = StridePolicyType::DEFAULT_VALUE)
    : StridePolicyType(stride)
    , m_bset(bSet)
    , m_map(SetType(bSet.size()), data, stride)
  {
    static_assert(std::is_same<BivariateSetType, UBSet>::value,
                  "Argument set is of a more-derived type than the Map's set "
                  "type. This may lead to object slicing. Use Map's pointer "
                  "constructor instead to store polymorphic sets.");
  }

  // (KW) Problem -- does not work with RelationSet
  template <typename BivariateSetRetType, typename RelType = void>
  typename std::enable_if<!traits::has_relation_ptr<BivariateSetRetType>::value,
                          BivariateSetRetType>::type
  getBivariateSet() const
  {
    using OuterSet = const typename BivariateSetRetType::FirstSetType;
    using InnerSet = const typename BivariateSetRetType::SecondSetType;
    OuterSet* outer = dynamic_cast<OuterSet*>(set()->getFirstSet());
    InnerSet* inner = dynamic_cast<InnerSet*>(set()->getSecondSet());

    return BivariateSetRetType(outer, inner);
  }

  template <typename BivariateSetRetType, typename RelType>
  typename std::enable_if<traits::has_relation_ptr<BivariateSetRetType>::value,
                          BivariateSetRetType>::type
  getBivariateSet() const
  {
    auto* rel = dynamic_cast<const RelType*>(m_bset)->getRelation();
    SLIC_ASSERT(rel != nullptr);

    return BivariateSetRetType(rel);
  }

  /// \name BivariateMap value access functions
  /// @{
  ///

  /**
   * \brief  Access the value in the map using a FlatIndex in the range of 0 to
   *         `size()*numComp()`
   *
   * \return The value for the j<sup>th</sup> component of the i<sup>th</sup>
   *         element, where `setIndex = i * numComp() + j`.
   * \pre    0 <= setIndex < size() * numComp()
   */
  AXOM_HOST_DEVICE ConstValueType operator[](SetPosition setIndex) const
  {
    return m_map[setIndex];
  }
  AXOM_HOST_DEVICE ValueType operator[](SetPosition setIndex)
  {
    return m_map[setIndex];
  }

public:
  /**
   * \brief Returns a SubMap containing the subset of the BivariateMap given the
   *        first set index
   * \pre 0 <= firstIdx < size(firstIdx)
   */
  AXOM_HOST_DEVICE ConstSubMapType operator()(SetPosition firstIdx) const
  {
#ifndef AXOM_DEVICE_CODE
    verifyFirstSetIndex(firstIdx);
#endif
    auto s = set()->elementRangeSet(firstIdx);
    const bool hasInd = submapIndicesHaveIndirection();
    return ConstSubMapType(this, s, hasInd);
  }

  AXOM_HOST_DEVICE SubMapType operator()(SetPosition firstIdx)
  {
#ifndef AXOM_DEVICE_CODE
    verifyFirstSetIndex(firstIdx);
#endif
    auto s = set()->elementRangeSet(firstIdx);
    const bool hasInd = submapIndicesHaveIndirection();
    return SubMapType(this, s, hasInd);
  }

  /**
   * \brief Access the value associated with the given SparseIndex into the
   *        BivariateSet and the component index.
   *
   * \pre `0 <= s1 < firstSetSize()`
   * \pre `0 <= s2 < size(s1)`
   * \pre `0 <= comp < numComp()`
   */
  AXOM_HOST_DEVICE ConstValueType operator()(SetPosition s1,
                                             SetPosition s2,
                                             SetPosition comp = 0) const
  {
    auto idx = flatIndex(s1, s2);
    return useCompIndexing() ? m_map(idx, comp) : m_map[idx];
  }

  AXOM_HOST_DEVICE ValueType operator()(SetPosition s1,
                                        SetPosition s2,
                                        SetPosition comp = 0)
  {
    auto idx = flatIndex(s1, s2);
    return useCompIndexing() ? m_map(idx, comp) : m_map[idx];
  }

  /**
   * \brief Access the value associated with the given DenseIndex into the
   *        BivariateSet and the component index.
   *
   * \pre `0 <= s1 < firstSetSize()`
   * \pre `0 <= s2 < secondSetSize()`
   * \pre `0 <= comp < numComp()`
   *
   * \return a DataType pointer to the value associated with the given index,
   *         or nullptr if there is no value for the given index.
   * \warning For sparse BivariateSet type, this function may have to do a
   *          linear search and can be slow.
   */
  AXOM_HOST_DEVICE ConstPointerType findValue(SetPosition s1,
                                              SetPosition s2,
                                              SetPosition comp = 0) const
  {
    SetPosition i = set()->findElementFlatIndex(s1, s2);
    if(i == BivariateSetType::INVALID_POS)
    {
      //the BivariateSet does not contain this index pair
      return nullptr;
    }
    return &(m_map(i, comp));
  }

  AXOM_HOST_DEVICE PointerType findValue(SetPosition s1,
                                         SetPosition s2,
                                         SetPosition comp = 0)
  {
    SetPosition i = set()->findElementFlatIndex(s1, s2);
    if(i == BivariateSetType::INVALID_POS)
    {
      //the BivariateSet does not contain this index pair
      return nullptr;
    }
    return &(m_map(i, comp));
  }

  /// @}

  /// \name BivariateMap index access functions
  /// @{
  ///

  /**
   * \brief Returns the SparseIndex of the element given the DenseIndex
   *
   * \return The SparseIndex of the element, or BivariateSet::INVALID_POS if
   *         the set doesn't contain the given DenseIndex
   */
  SetPosition index(SetPosition s1, SetPosition s2) const
  {
    return set()->findElementIndex(s1, s2);
  }

  /**
   * \brief Return a set of DenseIndex associated to the given first set index
   *
   * \param s1 the first set index
   * \return OrderedSet containing the elements
   */
  OrderedSetType indexSet(SetPosition s1) const
  {
    return set()->getElements(s1);
  }

  /**
   * \brief Search for the FlatIndex of an element given its DenseIndex in the
   *        BivariateSet.
   */
  AXOM_HOST_DEVICE inline SetPosition flatIndex(SetPosition s1, SetPosition s2) const
  {
    return set()->findElementFlatIndex(s1, s2);
  }

  /// @}

protected:
  /**
   * \brief Compile time predicate to check if we should use component indexing
   *
   * \note Intended to help optimize internal indexing
   */
  AXOM_HOST_DEVICE constexpr bool useCompIndexing() const
  {
    return !(StrPol::IS_COMPILE_TIME && StrPol::DEFAULT_VALUE == 1);
  }

  /**
   * \brief Utility function to determine if submaps should use indirection
   * when finding the set indices of their elements.
   *
   * This test distinguishes between ProductSet whose second set do not use
   * indirection and other BivariateSet types
   */
  AXOM_HOST_DEVICE constexpr bool submapIndicesHaveIndirection() const
  {
    return traits::indices_use_indirection<BivariateSetType>::value;
    //      || (set()->getSecondSet()->at(0) != 0);
  }

public:
  /**
   * \class BivariateMapIterator
   * \brief An iterator type for a BivariateMap, iterating via its
   *        ElementFlatIndex.
   *
   * This iterator class traverses the BivariateMap using its ElementFlatIndex.
   * In addition to m_pos from IteratorBase, this class also keeps track of the
   * iterator's first index (firstIdx), second dense index (secondIdx), and the
   * second sparse index (secondSparseIdx). The advance() function is
   * implemented to update those three additional indices.
   */
  template <bool Const>
  class BivariateMapIterator
    : public IteratorBase<BivariateMapIterator<Const>, SetPosition>
  {
  private:
    using iterator_category = std::random_access_iterator_tag;
    using value_type = DataType;
    using difference_type = SetPosition;

    using IterBase = IteratorBase<BivariateMapIterator, SetPosition>;
    using IterBase::m_pos;
    using iter = BivariateMapIterator;

  public:
    using DataRefType = std::conditional_t<Const, const DataType&, DataType&>;
    using BivariateMapPtr =
      std::conditional_t<Const, const BivariateMap*, BivariateMap*>;

    using PositionType = SetPosition;
    static constexpr PositionType INVALID_POS = -2;

  public:
    /**
     * \brief Construct a new BivariateMap Iterator given an ElementFlatIndex
     */
    BivariateMapIterator(BivariateMapPtr sMap, PositionType pos)
      : IterBase(pos)
      , m_map(sMap)
      , firstIdx(INVALID_POS)
      , secondIdx(INVALID_POS)
      , secondSparseIdx(INVALID_POS)
    {
      find_indices(pos);
    }

    bool operator==(const iter& other) const
    {
      return (m_map == other.m_map) && (m_pos == other.m_pos);
    }
    bool operator!=(const iter& other) const { return !operator==(other); }

    /**
     * \brief Returns the current iterator value. If the BivariateMap has
     *        multiple components, this will return the first component.
     *        To access the other components, use iter(comp)
     */
    DataRefType operator*() { return (*m_map)(firstIdx, secondIdx, 0); }

    /**
     * \brief Returns the iterator's value at the specified component.
     *        Returns the first component if comp_idx is not specified.
     * \param comp_idx  (Optional) Zero-based index of the component.
     */
    DataRefType operator()(PositionType comp_idx = 0)
    {
      return (*m_map)(firstIdx, secondIdx, comp_idx);
    }

    /** \brief Returns the first component value after n increments.  */
    DataRefType operator[](PositionType n) { return *(this->operator+(n)); }

    /**
     * \brief Return the value at the iterator's position. Same as operator()
     */
    DataRefType value(PositionType comp = 0)
    {
      return (*m_map)(firstIdx, secondIdx, comp);
    }

    /**
     * \brief return the current iterator's first index into the BivariateSet
     */
    PositionType firstIndex() const { return firstIdx; }

    /**
     * \brief return the current iterator's second index (DenseIndex)
     *        into the BivariateSet
     */
    PositionType secondIndex() const { return m_map->set()->at(m_pos); }

    /** \brief Returns the number of components per element in the map. */
    PositionType numComp() const { return m_map->numComp(); }

  private:
    /** Given the ElementFlatIndex, search for and update the other indices.
     *  This function does not depend on the three indices to be correct. */
    void find_indices(PositionType pos)
    {
      if(pos < 0 || pos > m_map->totalSize())
      {
        firstIdx = INVALID_POS;
        secondIdx = INVALID_POS;
        secondSparseIdx = INVALID_POS;
        return;
      }
      else if(pos == m_map->totalSize())
      {
        firstIdx = m_map->firstSetSize();
        secondIdx = 0;
        secondSparseIdx = 0;
        return;
      }

      firstIdx = 0;
      PositionType beginIdx = 0;
      while(beginIdx + m_map->set()->size(firstIdx) <= pos)
      {
        beginIdx += m_map->set()->size(firstIdx);
        firstIdx++;
      }

      SLIC_ASSERT(firstIdx < m_map->firstSetSize());
      secondIdx = m_map->set()->at(pos);
      secondSparseIdx = pos - beginIdx;
    }

    /* recursive helper function for advance(n) to update the three indices.
     * This is an updating function, which assumes the pre-advance state is
     * valid, i.e. the three indices were correct prior to advance(n).
     * It will recurse as many times as the change in firstIdx. */
    void advance_helper(PositionType n, PositionType idx1, PositionType idx2)
    {
      const BivariateSetType* set = m_map->set();
      if(idx2 + n < 0)
        advance_helper(n + (idx2 + 1), idx1 - 1, set->size(idx1 - 1) - 1);
      else if(idx2 + n >= set->size(idx1))
        advance_helper(n - (set->size(idx1) - idx2), idx1 + 1, 0);
      else
      {
        firstIdx = idx1;
        secondSparseIdx = idx2 + n;
        secondIdx = m_map->set()->at(m_pos);
      }
    }

  protected:
    /** Implementation of advance() as required by IteratorBase.
     *  It updates the three other indices as well. */
    void advance(PositionType n)
    {
      m_pos += n;
      PositionType size = m_map->totalSize();

      if(firstIdx == INVALID_POS)
      {  //iterator was in an invalid position. search for the indices.
        find_indices(m_pos);
      }
      else if(m_pos == size)
      {
        firstIdx = m_map->firstSetSize();
        secondIdx = 0;
        secondSparseIdx = 0;
      }
      else if(m_pos < 0 || m_pos > size)
      {
        firstIdx = INVALID_POS;
        secondIdx = INVALID_POS;
        secondSparseIdx = INVALID_POS;
      }
      else
      {
        advance_helper(n, firstIdx, secondSparseIdx);
      }
    }

  private:
    BivariateMapPtr m_map;
    PositionType firstIdx;
    PositionType secondIdx;
    PositionType secondSparseIdx;
  };

public:
  /** BivariateMap iterator functions */
  AXOM_HOST_DEVICE iterator begin() { return iterator(this, 0); }
  AXOM_HOST_DEVICE iterator end() { return iterator(this, totalSize()); }
  AXOM_HOST_DEVICE const_iterator begin() const
  {
    return const_iterator(this, 0);
  }
  AXOM_HOST_DEVICE const_iterator end() const
  {
    return const_iterator(this, totalSize());
  }

  /** Iterator via Submap */
  AXOM_HOST_DEVICE SubMapIterator begin(int i) { return (*this)(i).begin(); }
  AXOM_HOST_DEVICE SubMapIterator end(int i) { return (*this)(i).end(); }
  AXOM_HOST_DEVICE ConstSubMapIterator begin(int i) const
  {
    return (*this)(i).begin();
  }
  AXOM_HOST_DEVICE ConstSubMapIterator end(int i) const
  {
    return (*this)(i).end();
  }

public:
  AXOM_HOST_DEVICE const BivariateSetType* set() const { return m_bset.get(); }
  const MapType* getMap() const { return &m_map; }
  MapType* getMap() { return &m_map; }

  bool isValid(bool verboseOutput = false) const
  {
    return set()->isValid(verboseOutput) && m_map.isValid(verboseOutput);
  }

  /// \name BivariateMap cardinality functions
  /// @{
  ///

  /** \brief Returns the BivariateSet size. */
  AXOM_HOST_DEVICE SetPosition size() const { return set()->size(); }
  /** \brief Returns the BivariateSet size. */
  SetPosition totalSize() const { return set()->size(); }

  SetPosition firstSetSize() const { return set()->firstSetSize(); }
  AXOM_HOST_DEVICE SetPosition secondSetSize() const
  {
    return set()->secondSetSize();
  }
  /** \brief Returns the number of the BivariateSet ordered pairs with
   *         the given first set index. */
  SetPosition size(SetPosition s) const { return set()->size(s); }
  /** \brief Return the number of components of the map  */
  SetPosition numComp() const { return StrPol::stride(); }

  /// @}

  /**
   * \brief Given a DataType array of size `totalSize()*numComp()`, copy
   *        the data into the BivariateMap storage.
   *
   * \param data_arr The array of DataType that contains the data to be copied.
   */
  void copy(DataType* data_arr)
  {
    for(int i = 0; i < m_map.size() * StrPol::stride(); i++)
      m_map[i] = data_arr[i];
  }
  void copy(const DataType* data_arr)
  {
    for(int i = 0; i < m_map.size() * StrPol::stride(); i++)
      m_map[i] = data_arr[i];
  }

  /** \brief replace all elements in the Map with the default DataType */
  void clear() { m_map.clear(); }

private:
  /** \brief Check the indices (DenseIndex) are valid   */
  void verifyPosition(SetPosition s1, SetPosition s2) const
  {
    set()->verifyPosition(s1, s2);
  }

  /** \brief Check the given ElementFlatIndex is valid.  */
  void verifyPosition(SetPosition AXOM_DEBUG_PARAM(pos)) const
  {
    SLIC_ASSERT_MSG(pos >= 0 && pos < SetPosition(m_map.size()),
                    "Attempted to access element "
                      << pos << " but BivariateMap's data has size "
                      << m_map.size());
  }

  void verifyFirstSetIndex(SetPosition AXOM_DEBUG_PARAM(firstIdx)) const
  {
    SLIC_ASSERT_MSG(firstIdx >= 0 && firstIdx < firstSetSize(),
                    "Attempted to access elements with first set index "
                      << firstIdx << ", but BivariateMap's first set has size "
                      << firstSetSize());
  }

private:
  BSetContainer<> m_bset;
  MapType m_map;
};  //end BivariateMap

template <typename T, typename BSet, typename IndPol, typename StrPol, typename IfacePol>
typename BivariateMap<T, BSet, IndPol, StrPol, IfacePol>::NullBivariateSetType const
  BivariateMap<T, BSet, IndPol, StrPol, IfacePol>::s_nullBiSet;

}  // end namespace slam
}  // end namespace axom

#endif  // SLAM_BIVARIATE_MAP_HPP_
