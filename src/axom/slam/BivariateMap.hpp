// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
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

#include <cassert>

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

template <typename SetType,
          typename DataType,
          typename StridePolicy = policies::StrideOne<typename SetType::PositionType>>
class BivariateMap : public MapBase<typename SetType::PositionType>,
                     public StridePolicy
{
public:
  using SetPosition = typename SetType::PositionType;
  using SetElement = typename SetType::ElementType;

  using MapType = slam::Map<SetType, DataType, StridePolicy>;

  using BivariateSetType = BivariateSet<SetPosition, SetElement>;
  using OrderedSetType = typename BivariateSetType::OrderedSetType;

  using BivariateMapType = BivariateMap<SetType, DataType, StridePolicy>;
  using SubMapType = SubMap<SetType, DataType, BivariateMapType, StridePolicy>;
  using SubMapIterator = typename SubMapType::SubMapIterator;

private:
  static const NullBivariateSet<SetPosition, SetElement> s_nullBiSet;

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
               SetPosition stride = StridePolicy::DEFAULT_VALUE)
    : StridePolicy(stride)
    , m_set(bSet)
    , m_rangeSet(0, bSet->size())
    , m_map(&m_rangeSet, defaultValue, stride)
  { }

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
  const DataType& operator[](SetPosition setIndex) const
  {
    return m_map[setIndex];
  }
  DataType& operator[](SetPosition setIndex) { return m_map[setIndex]; }

private:
  //template for both const / non-const SubMap creator, given the firstIdx and
  //this pointer, which could be const or non-const.
  template <class constOrNonConstMap>
  SubMapType makeSubMap(SetPosition firstIdx, constOrNonConstMap* map_ptr) const
  {
    SLIC_ASSERT_MSG(firstIdx >= 0 && firstIdx < firstSetSize(),
                    "Attempted to access elements with first set index "
                      << firstIdx << ", but BivairateMap's first set has size "
                      << firstSetSize());

    SetPosition start_idx = m_set->findElementFlatIndex(firstIdx);
    SetPosition size = m_set->size(firstIdx);
    RangeSet<> rng_set(start_idx, start_idx + size);
    return SubMapType(map_ptr, rng_set);
  }

public:
  /**
   * \brief Returns a SubMap containing the subset of the BivariateMap given the
   *        first set index
   * \pre 0 <= firstIdx < size(firstIdx)
   */
  const SubMapType operator()(SetPosition firstIdx) const
  {
    return makeSubMap<const BivariateMapType>(firstIdx, this);
  }

  SubMapType operator()(SetPosition firstIdx)
  {
    return makeSubMap<BivariateMapType>(firstIdx, this);
  }

  /**
   * \brief Access the value associated with the given SparseIndex into the
   *        BivariateSet and the component index.
   *
   * \pre `0 <= s1 < firstSetSize()`
   * \pre `0 <= s2 < size(s1)`
   * \pre `0 <= comp < numComp()`
   */
  const DataType& operator()(SetPosition s1,
                             SetPosition s2,
                             SetPosition comp = 0) const
  {
    return m_map(m_set->findElementFlatIndex(s1, s2), comp);
  }

  DataType& operator()(SetPosition s1, SetPosition s2, SetPosition comp = 0)
  {
    const BivariateMap& constMe = *this;
    return const_cast<DataType&>(constMe(s1, s2, comp));
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
  const DataType* findValue(SetPosition s1, SetPosition s2, SetPosition comp = 0) const
  {
    SetPosition i = m_set->findElementFlatIndex(s1, s2);
    if(i == BivariateSetType::INVALID_POS)
    {
      //the BivariateSet does not contain this index pair
      return nullptr;
    }
    return &(m_map(i, comp));
  }

  DataType* findValue(SetPosition s1, SetPosition s2, SetPosition comp = 0)
  {
    const BivariateMap& constMe = *this;
    return const_cast<DataType*>(constMe.findValue(s1, s2, comp));
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
    return m_set->findElementIndex(s1, s2);
  }

  /**
   * \brief Return a set of DenseIndex associated to the given first set index
   *
   * \param s1 the first set index
   * \return OrderedSet containing the elements
   */
  OrderedSetType indexSet(SetPosition s1) const
  {
    return m_set->getElements(s1);
  }

  /**
   * \brief Search for the FlatIndex of an element given its DenseIndex in the
   *        BivariateSet.
   */
  SetPosition flatIndex(SetPosition s1, SetPosition s2) const
  {
    return m_set->findElementFlatIndex(s1, s2);
  }

  /// @}

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
  class BivariateMapIterator
    : public IteratorBase<BivariateMapIterator, SetPosition>
  {
  private:
    using iterator_category = std::random_access_iterator_tag;
    using value_type = DataType;
    using difference_type = SetPosition;

    using IterBase = IteratorBase<BivariateMapIterator, SetPosition>;
    using IterBase::m_pos;
    using iter = BivariateMapIterator;

  public:
    using PositionType = SetPosition;
    const PositionType INVALID_POS = -2;

  public:
    /**
     * \brief Construct a new BivariateMap Iterator given an ElementFlatIndex
     */
    BivariateMapIterator(BivariateMap* sMap, PositionType pos)
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
    DataType& operator*() { return (*m_map)(firstIdx, secondIdx, 0); }

    /**
     * \brief Returns the iterator's value at the specified component.
     *        Returns the first component if comp_idx is not specified.
     * \param comp_idx  (Optional) Zero-based index of the component.
     */
    DataType& operator()(PositionType comp_idx = 0)
    {
      return (*m_map)(firstIdx, secondIdx, comp_idx);
    }

    /** \brief Returns the first component value after n increments.  */
    DataType& operator[](PositionType n) { return *(this->operator+(n)); }

    /**
     * \brief Return the value at the iterator's position. Same as operator()
     */
    DataType& value(PositionType comp = 0)
    {
      return (*m_map)(firstIdx, secondIdx, comp);
    }

    /**
     * \brief return the current iterator's first index into the BivariateSet
     */
    PositionType firstIndex() { return firstIdx; }

    /**
     * \brief return the current iterator's second index (DenseIndex)
     *        into the BivariateSet
     */
    PositionType secondIndex() { return m_map->set()->at(m_pos); }

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
    BivariateMap* m_map;
    PositionType firstIdx;
    PositionType secondIdx;
    PositionType secondSparseIdx;
  };

public:
  /** BivariateMap iterator functions */
  BivariateMapIterator begin() { return BivariateMapIterator(this, 0); }
  BivariateMapIterator end() { return BivariateMapIterator(this, totalSize()); }

  /** Iterator via Submap */
  SubMapIterator begin(int i) { return (*this)(i).begin(); }
  SubMapIterator end(int i) { return (*this)(i).end(); }

public:
  const BivariateSetType* set() const { return m_set; }
  const MapType* getMap() const { return &m_map; }
  MapType* getMap() { return &m_map; }

  virtual bool isValid(bool verboseOutput = false) const override
  {
    return m_set->isValid(verboseOutput) && m_map.isValid(verboseOutput);
  }

  /// \name BivariateMap cardinality functions
  /// @{
  ///

  /** \brief Returns the BivariateSet size. */
  SetPosition size() const override { return m_set->size(); }
  /** \brief Returns the BivariateSet size. */
  SetPosition totalSize() const { return m_set->size(); }

  SetPosition firstSetSize() const { return m_set->firstSetSize(); }
  SetPosition secondSetSize() const { return m_set->secondSetSize(); }
  /** \brief Returns the number of the BivariateSet ordered pairs with
   *         the given first set index. */
  SetPosition size(SetPosition s) const { return m_set->size(s); }
  /** \brief Return the number of components of the map  */
  SetPosition numComp() const { return StridePolicy::stride(); }

  /// @}

  /**
   * \brief Given a DataType array of size `totalSize()*numComp()`, copy
   *        the data into the BivariateMap storage.
   *
   * \param data_arr The array of DataType that contains the data to be copied.
   */
  void copy(DataType* data_arr)
  {
    for(int i = 0; i < m_map.size() * StridePolicy::stride(); i++)
      m_map[i] = data_arr[i];
  }

private:
  /** \brief Check the indices (DenseIndex) are valid   */
  void verifyPosition(SetPosition s1, SetPosition s2) const
  {
    m_set->verifyPosition(s1, s2);
  }

  /** \brief Check the given ElementFlatIndex is valid.  */
  void verifyPosition(SetPosition AXOM_DEBUG_PARAM(pos)) const override
  {
    SLIC_ASSERT_MSG(pos >= 0 && pos < SetPosition(m_map.size()),
                    "Attempted to access element "
                      << pos << " but BivairateMap's data has size "
                      << m_map.size());
  }

private:
  const BivariateSetType* m_set;
  RangeSet<> m_rangeSet;  //for the map... since BivariateSet isn't a slam::Set
  MapType m_map;

};  //end BivariateMap

template <typename SetType, typename DataType, typename StridePolicy>
NullBivariateSet<typename SetType::PositionType, typename SetType::ElementType> const
  BivariateMap<SetType, DataType, StridePolicy>::s_nullBiSet;

}  // end namespace slam
}  // end namespace axom

#endif  // SLAM_BIVARIATE_MAP_HPP_
