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
 *
 */

#include "axom/slic.hpp"

#include "axom/slam/Map.hpp"
#include "axom/slam/OrderedSet.hpp"
#include "axom/slam/RangeSet.hpp"
#include "axom/core/IteratorBase.hpp"

#include <type_traits>
#include <sstream>

namespace axom::slam
{
/**
 * \class SubMap
 * \brief The SubMap class provides an API to easily traverse a subset of a Map.
 *
 * A SubMap is defined by a subset of the indices into a Map, which we refer to
 * as its SuperMap (of type SuperMapType). The indices are expressed as ElementFlatIndex.
 * See BivariateMap for an explanation of the various indexing schemes.
 *
 * SubMap is used by BivariateMap to return a set of values mapped to each item in its first set.
 *
 * \tparam SuperMapType the type of SuperMap
 * \tparam SetType defines the indices int the super map. It cannot be abstract.
 *
 * \warning SubMap constructor can take a const Map pointer or a non-const Map pointer.
 *        A non-const value access function in SubMap will fail
 *        if the SubMap is constructed using a const Map pointer.
 *
 * \see Map, BivariateMap
 */

template <typename SuperMapType,
          typename SubsetType,  //= slam::RangeSet<PositionType, SetElement>
          typename InterfacePolicy = policies::ConcreteInterface>
class SubMap : public policies::MapInterface<InterfacePolicy, typename SubsetType::PositionType>,
               public SuperMapType::StridePolicyType
{
public:
  static_assert(!std::is_abstract<SubsetType>::value, "SetType for slam::SubMap cannot be abstract");

  using DataType = typename SuperMapType::DataType;

  /// The subset carrying flat indices into the super-map's domain.
  using IndexSetType = SubsetType;
  /// The SubMap domain; retained as SetType for MapLike compatibility.
  using SetType = IndexSetType;
  using PositionType = typename SubsetType::PositionType;
  using SetElement = typename SubsetType::ElementType;
  /// The element obtained after projecting a subset index through the super-map's set.
  using ProjectedElement = typename SuperMapType::SetElement;

  using StridePolicyType = typename SuperMapType::StridePolicyType;
  using IndirectionPolicy = typename SuperMapType::IndirectionPolicy;

  using ElementShape = typename StridePolicyType::ShapeType;

  //iterator type aliases
  class Iterator;
  using iterator = Iterator;
  using const_iterator = Iterator;
  using iterator_pair = std::pair<iterator, iterator>;

  class RangeIterator;
  using const_range_iterator = RangeIterator;
  using range_iterator = RangeIterator;

  /*!
   * \brief The reference type for a value access
   *
   * A SubMap is a view, so the constness comes from \a SuperMapType
   * rather than that of the SubMap object.
   */
  using DataRefType = std::conditional_t<std::is_const<SuperMapType>::value,
                                         typename IndirectionPolicy::ConstIndirectionResult,
                                         typename IndirectionPolicy::IndirectionResult>;
  using ValueType = DataRefType;
  using ConstValueType = DataRefType;

public:
  /// Default Constructor
  SubMap() : m_superMap(nullptr) { }

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
    : StridePolicyType(*supermap)
    , m_superMap(supermap)
    , m_subsetIdx(subset_idxset)
  {
    // Checked here instead of in the class body since BivariateMap has a
    // `SubMapType::iterator` member aliases, and a class-body
    // constraint on SuperMapType would then depend on itself.
    static_assert(SubMappable<SuperMapType>,
                  "SubMap requires a super-map whose stride and indirection policies it can re-use");
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
    IndexType flat_idx = getMapCompFlatIndex(idx);
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
    static_assert(axom::detail::all_types_are_integral<ComponentIndex...>::value,
                  "SubMap::operator(): index parameter pack must all be integral types.");
#ifndef AXOM_DEVICE_CODE
    verifyPositionImpl(idx, comp...);
#endif
    SLIC_ASSERT_MSG(m_superMap != nullptr, "Submap's super map was null.");
    IndexType elemBegin = getMapElemFlatIndex(idx) * numComp();
    return (*m_superMap)[elemBegin + componentOffset(comp...)];
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
    return operator()(idx, comp...);
  }

  /*!
   * \brief Project a subset position to the corresponding super-map domain element.
   *
   * \note index() is not part of the MapLike contract. 
   * SetType describes the index-carrying subset returned by set(),
   * while this function returns a ProjectedElement from the super-map's semantic domain.
   * For a bivariate super-map, that projected element is a coordinate pair.
   */
  AXOM_SUPPRESS_HD_WARN
  AXOM_HOST_DEVICE ProjectedElement index(IndexType idx) const
  {
    return m_superMap->set()->at(m_subsetIdx[idx]);
  }

  /// @}

  /// \name SubMap cardinality functions
  /// @{
  ///

  /// \brief returns the subset that forms this SubMap's domain
  AXOM_HOST_DEVICE const SetType* set() const { return &m_subsetIdx; }

  /// \brief returns the size of the SubMap
  AXOM_HOST_DEVICE PositionType size() const { return m_subsetIdx.size(); }

  /// \brief returns the number of components (aka. stride) of the SubMap
  AXOM_HOST_DEVICE IndexType numComp() const { return StridePolicyType::stride(); }

  /// @}

  [[nodiscard]] bool isValid(bool VerboseOutput = false) const;

private:  //function inherit from StridePolicy that should not be accessible
  void setStride(IndexType)
  {
    SLIC_ASSERT_MSG(false, "Stride should not be changed after construction of SubMap.");
  }

private:  //helper functions
  friend class RangeIterator;
  /// \brief Get the ElementFlatIndex into the SuperMap given the subset's index.
  AXOM_HOST_DEVICE IndexType getMapElemFlatIndex(IndexType idx) const { return m_subsetIdx[idx]; }

  /**
   * \brief Get the ComponentFlatIndex into the SuperMap given the subset's
   * ComponentFlatIndex. This is used only with bracket [] access
   */
  AXOM_HOST_DEVICE IndexType getMapCompFlatIndex(IndexType idx) const
  {
    IndexType comp = numComp();
    IndexType s = idx % comp;
    return getMapElemFlatIndex(idx / comp) * comp + s;
  }

  /// Checks the ComponentFlatIndex is valid
  void verifyPosition(PositionType idx) const { verifyPositionImpl(idx); }

  /// Checks the ElementFlatIndex and the component index is valid
  void verifyPosition(PositionType idx, PositionType comp) const { verifyPositionImpl(idx, comp); }

  /// Checks the ComponentFlatIndex is valid
  void verifyPositionImpl(PositionType AXOM_DEBUG_PARAM(idx)) const
  {
    SLIC_ASSERT_MSG(idx >= 0 && idx < m_subsetIdx.size() * numComp(),
                    "Attempted to access element " << idx << " but Submap's data has size "
                                                   << m_subsetIdx.size() * numComp());
  }

  /// Checks the ElementFlatIndex and the component index is valid
  template <typename ComponentIndex>
  void verifyPositionImpl(PositionType AXOM_DEBUG_PARAM(idx), ComponentIndex AXOM_DEBUG_PARAM(comp)) const
  {
    SLIC_ASSERT_MSG(idx >= 0 && idx < m_subsetIdx.size() && comp >= 0 && comp < numComp(),
                    "Attempted to access element "
                      << idx << " component " << comp << ", but Submap's data has size "
                      << m_subsetIdx.size() << " with " << numComp() << " component");
  }

  /// Checks the ElementFlatIndex and the component index is valid
  template <typename... ComponentIndex>
  void verifyPositionImpl(PositionType AXOM_DEBUG_PARAM(idx),
                          ComponentIndex... AXOM_DEBUG_PARAM(comp)) const
  {
#ifdef AXOM_DEBUG
    ElementShape indexArray {{comp...}};
    bool validIndexes = true;
    for(int dim = 0; dim < StridePolicyType::NumDims; dim++)
    {
      validIndexes = validIndexes && (indexArray[dim] >= 0);
      validIndexes = validIndexes && (indexArray[dim] < m_superMap->shape()[dim]);
    }
    std::string invalid_message = fmt::format(
      "Attempted to access element {} component ({}), but SubMap's data has "
      "size {} with component shape ({})",
      idx,
      fmt::join(indexArray, ", "),
      m_subsetIdx.size(),
      fmt::join(m_superMap->shape(), ", "));
    SLIC_ASSERT_MSG(idx >= 0 && idx < m_subsetIdx.size() && validIndexes, invalid_message);
#endif
  }

  /// \brief Computes the flat indexing offset for a given component.
  AXOM_HOST_DEVICE inline PositionType componentOffset() const { return PositionType {}; }
  template <typename ComponentIndex>
  AXOM_HOST_DEVICE inline PositionType componentOffset(ComponentIndex componentIndex) const
  {
    return componentIndex;
  }

  template <typename... ComponentIndex>
  AXOM_HOST_DEVICE inline PositionType componentOffset(ComponentIndex... componentIndex) const
  {
    ElementShape indexArray {{componentIndex...}};
    ElementShape strides = StridePolicyType::strides();
    PositionType offset = 0;
    for(int dim = 0; dim < StridePolicyType::NumDims; dim++)
    {
      offset += indexArray[dim] * strides[dim];
    }
    return offset;
  }

public:  // Functions related to iteration
  AXOM_HOST_DEVICE iterator begin() const { return iterator(this, 0); }
  AXOM_HOST_DEVICE iterator end() const { return iterator(this, m_subsetIdx.size() * numComp()); }
  AXOM_HOST_DEVICE range_iterator set_begin() const { return range_iterator(this, 0); }
  AXOM_HOST_DEVICE range_iterator set_end() const
  {
    return range_iterator(this, m_subsetIdx.size());
  }

protected:  //Member variables
  SuperMapType* m_superMap;
  SubsetType m_subsetIdx;

};  //end SubMap

template <typename SuperMapType, typename SetType, typename InterfacePolicy>
bool SubMap<SuperMapType, SetType, InterfacePolicy>::isValid(bool verboseOutput) const
{
  bool isValid = true;
  std::stringstream errStr;

  if(m_superMap == nullptr)
  {
    if(m_subsetIdx.size() > 0)
    {
      isValid = false;
      if(verboseOutput)
      {
        errStr << "\n\t*SuperMap pointer is null, "
               << "but the subset index is non-empty";
      }
    }
  }
  else
  {
    int map_size = m_superMap->size();
    //Check all indices is inside the SuperMap range
    for(int i = 0; i < m_subsetIdx.size(); i++)
    {
      PositionType pos = m_subsetIdx[i];
      if(pos < 0 || pos >= map_size)
      {
        isValid = false;
        if(verboseOutput)
        {
          errStr << "\n\t* The given subset index " << pos
                 << "is outside of the SuperMap range of 0 to " << map_size;
        }
      }
    }
  }

  if(verboseOutput)
  {
    SLIC_INFO("Detailed results of isValid on the SubMap.\n"
              << "SubMap was " << (isValid ? "valid" : "NOT valid") << "\n"
              << errStr.str());
  }

  return isValid;
}

/**
 * \class SubMap::Iterator
 * \brief An iterator for SubMap, based on MapIterator
 *
 * \see MapIterator
 */
template <typename SuperMapType, typename SubsetType, typename InterfacePolicy>
class SubMap<SuperMapType, SubsetType, InterfacePolicy>::Iterator
  : public IteratorBase<Iterator, PositionType>
{
public:
  using iterator_concept = std::random_access_iterator_tag;
  using iterator_category = std::random_access_iterator_tag;
  using value_type = DataType;
  using reference = DataRefType;
  using pointer = std::add_pointer_t<std::remove_reference_t<reference>>;
  using difference_type = PositionType;

  using IterBase = IteratorBase<Iterator, PositionType>;
  using IterBase::m_pos;
  using iter = Iterator;

  Iterator() = default;

  AXOM_HOST_DEVICE Iterator(const SubMap* sMap, PositionType pos)
    : IterBase(pos)
    , m_submap(*sMap) { }

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
  SubMap m_submap;
};

/**
 * \class SubMap::RangeIterator
 * \brief An iterator for SubMap, based on MapIterator
 *
 * \see MapIterator
 */
template <typename SuperMapType, typename SubsetType, typename InterfacePolicy>
class SubMap<SuperMapType, SubsetType, InterfacePolicy>::RangeIterator
  : public IteratorBase<RangeIterator, PositionType>
{
private:
public:
  using IterBase = IteratorBase<RangeIterator, PositionType>;
  using IterBase::m_pos;
  using iter = Iterator;

private:
  using MapRangeIterator = std::conditional_t<std::is_const<SuperMapType>::value,
                                              typename SuperMapType::const_range_iterator,
                                              typename SuperMapType::range_iterator>;

public:
  // Dereference returns a reference to a cached ArrayView,
  // while subscript returns a value to avoid dangling from a temporary iterator.
  // \warning Inherits MapRangeIterator's multipass and dangling caveats;
  //  see the warning on Map::MapRangeIterator.
  using iterator_concept = std::bidirectional_iterator_tag;
  using iterator_category = std::bidirectional_iterator_tag;
  using value_type = typename MapRangeIterator::value_type;
  using reference = typename MapRangeIterator::reference;
  using pointer = typename MapRangeIterator::pointer;
  using difference_type = PositionType;

  RangeIterator() = default;

  AXOM_HOST_DEVICE PositionType getParentPosition(PositionType subset_pos)
  {
    PositionType subsetEnd = m_submap.m_subsetIdx.size() - 1;
    // End element is one past the last subset element.
    PositionType parentIndex = m_submap.m_subsetIdx[subsetEnd] + 1;
    if(subset_pos < m_submap.m_subsetIdx.size())
    {
      parentIndex = m_submap.m_subsetIdx[subset_pos];
    }
    return parentIndex;
  }

public:
  AXOM_HOST_DEVICE RangeIterator(const SubMap* sMap, PositionType pos)
    : IterBase(pos)
    , m_submap(*sMap)
    , m_mapIter(m_submap.m_superMap, getParentPosition(pos))
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

  /// \brief Returns the flat index in the original map pointed to by thisiterator.
  PositionType flatIndex() const { return m_mapIter.flatIndex(); }

  /// \brief Returns the index into the submap pointed to by this iterator.
  PositionType submapIndex() const { return this->m_pos; }

  /// \brief Returns the number of components per element in the Map.
  PositionType numComp() const { return m_mapIter.numComp(); }

protected:
  /*!
   * \brief Implementation of advance() as required by IteratorBase
   *
   * Class invariant: m_mapIter sits at getParentPosition(m_pos) in the super-map's index space,
   * so the step is the difference between the parent positions of the old and new subset positions.
   *
   * \note Does not ask m_mapIter where it is since flatIndex() reports a position
   *  in the original map, which coincides with the super-map's index space only
   *  when the super-map is a Map or BivariateMap. When the super-map is itself a SubMap,
   *  flatIndex() skips a level and the difference would be taken across two different index spaces.
   */
  AXOM_HOST_DEVICE void advance(PositionType n)
  {
    const PositionType currIndex = getParentPosition(this->m_pos);
    const PositionType nextIndex = getParentPosition(this->m_pos + n);
    m_mapIter += (nextIndex - currIndex);

    this->m_pos += n;
  }

private:
  SubMap m_submap;
  MapRangeIterator m_mapIter;
};

}  // end namespace axom::slam
