// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file SubMap.hpp
 *
 * \brief Contains SubMap, which is a subset of a Map
 *
 */

#ifndef SLAM_SUBMAP_HPP_
#define SLAM_SUBMAP_HPP_

#include "axom/slic.hpp"

#include "axom/slam/Map.hpp"
#include "axom/slam/OrderedSet.hpp"
#include "axom/slam/RangeSet.hpp"
#include "axom/slam/IteratorBase.hpp"

#include <type_traits>
#include <sstream>

namespace axom
{
namespace slam
{

/**
 * \class SubMap
 * \brief The SubMap class provides an API to easily traverse a subset of a Map.
 *
 * A SubMap is defined by a subset of the indices into a Map, which we refer to
 * as its SuperMap (of type SuperMapType). The indices are expressed as
 * ElementFlatIndex.\n
 * Please see BivariateMap for an explanation of the various indexing schemes.
 *
 * SubMap is used by BivariateMap to return a set of values mapped to each item
 * in its first set.
 *
 * \tparam SuperMapType the type of SuperMap
 * \tparam SetType defines the Position and Element types of the underlying set
 *
 * \warning SubMap constructor can take a const Map pointer or a non-const Map
 *        pointer. A non-const value access function in SubMap will fail if the
 *        Submap is constructed using a const Map pointer.
 *
 * \see Map, BivariateMap
 */

/// TODO: Modify code so that default SetType uses an indirection set, as it currently is
///       But when the user passes in a custom set type (e.g. RangeSet) -- use that.

template<
  typename SuperMapType,
  typename SetType
  >
class SubMap : public MapBase, public SuperMapType::StridePolicyType
{
public:
  using DataType = typename SuperMapType::DataType;

  using SetPosition = typename SetType::PositionType;
  using SetElement = typename SetType::ElementType;

  using StridePolicyType = typename SuperMapType::StridePolicyType;

private:
  using RangeSetType = RangeSet<SetPosition, SetElement>;

  template<typename T>  // TODO: Generalize this
  using IndPol = policies::STLVectorIndirection<SetPosition, T>;

  using MapType = Map<SetType, DataType, IndPol<DataType>, StridePolicyType>;

  using IndirSet = OrderedSet<SetPosition,SetElement,
                              policies::RuntimeSize<SetPosition>,
                              policies::ZeroOffset<SetPosition>,
                              policies::StrideOne<SetPosition>,
                              IndPol<SetElement> >;

  enum { USE_RANGE_SET = std::is_same<SetType, RangeSetType>::value };

public:
  using SubsetType = typename std::conditional<USE_RANGE_SET,
                                               RangeSetType,
                                               IndirSet
                                              >::type;

  using SubsetBuilder = typename SubsetType::SetBuilder;

  //iterator type aliases
  class SubMapIterator;
  using const_iterator = SubMapIterator;
  using const_iterator_pair = std::pair<const_iterator, const_iterator>;
  using iterator = const_iterator;
  using iterator_pair = const_iterator_pair;

private:

//  template<typename S>
//  typename std::enable_if< USE_RANGE_SET, void>::type
//  initSubset(const S& subset_idx)
//  {
//     m_subsetIdx = subset_idx;
//  }
//
//  template<typename S>
//  typename std::enable_if< ! USE_RANGE_SET, void>::type
//  initSubset(const S& subset_idx)
//  {
//     m_subsetIdx_data.resize(subset_idx.size());
//     m_subsetIdx = SubsetBuilder()
//                         .size(subset_idx.size())
//                         .offset()
//                         .data(&m_subsetIdx_data);
//
//     //copy the elements in the SuperMap's Set
//     for (int i = 0 ; i < subset_idx.size() ; i++)
//     {
//       m_subsetIdx_data[i] = subset_idx.at(i);
//     }
//  }
//
//  template<typename S>
//  typename std::enable_if< USE_RANGE_SET, void>::type
//  copySubset(const SubMap& other)
//  {
//     m_subsetIdx = other.m_subsetIdx;
//  }
//
//  template<typename S>
//  typename std::enable_if< ! USE_RANGE_SET, void>::type
//  copySubset(const SubMap& other)
//  {
//     m_subsetIdx_data = other.m_subsetIdx_data;
//     m_subsetIdx = SubsetBuilder()
//                   .size(m_subsetIdx_data.size())
//                   .data(&m_subsetIdx_data);
//  }
public:
  /** Default Constructor */
  SubMap() : m_superMap_constptr(nullptr), m_superMap_ptr(nullptr) {}

  /**
   * \brief Constructor for SubMap given the ElementFlatIndex into the SuperMap
   *
   * \param supermap The map that this SubMap is a subset of.
   * \param subset_idx a Set of ElementFlatIndex into the SuperMap
   */
  template<typename TheSetType,
           typename std::enable_if<std::is_same<TheSetType, RangeSetType>::value,int>::type = 0>
  SubMap(const SuperMapType* supermap,  TheSetType& subset_idx)
    : StridePolicyType(supermap->stride())
    , m_superMap_constptr(supermap), m_superMap_ptr(nullptr)
    , m_subsetIdx(subset_idx)
  {
  }

  template<typename TheSetType,
           typename std::enable_if<!std::is_same<TheSetType, RangeSetType>::value,int>::type = 0>
  SubMap(const SuperMapType* supermap,  TheSetType& subset_idx)
    : StridePolicyType(supermap->stride())
    , m_superMap_constptr(supermap), m_superMap_ptr(nullptr)
    , m_subsetIdx_data( subset_idx.size() )
    , m_subsetIdx ( SubsetBuilder()
                        .size(subset_idx.size())
                        .data(&m_subsetIdx_data) )
  {
     //copy the elements in the SuperMap's Set
     for (int i = 0 ; i < subset_idx.size() ; i++)
     {
       m_subsetIdx_data[i] = subset_idx.at(i);
     }
  }

  //non-const version of above constructor
  template<typename TheSetType,
           typename std::enable_if<std::is_same<TheSetType, RangeSetType>::value,int>::type = 0>
  SubMap(SuperMapType* supermap,  TheSetType& subset_idx)
    : StridePolicyType(supermap->stride())
    , m_superMap_constptr(supermap), m_superMap_ptr(supermap)
    , m_subsetIdx(subset_idx)
  {
  }

  template<typename TheSetType,
           typename std::enable_if<!std::is_same<TheSetType, RangeSetType>::value,int>::type = 0>
  SubMap(SuperMapType* supermap,  TheSetType& subset_idx)
    : StridePolicyType(supermap->stride())
    , m_superMap_constptr(supermap), m_superMap_ptr(supermap)
    , m_subsetIdx_data( subset_idx.size() )
    , m_subsetIdx ( SubsetBuilder()
                        .size(subset_idx.size())
                        .data(&m_subsetIdx_data) )
  {
     //copy the elements in the SuperMap's Set
     for (int i = 0 ; i < subset_idx.size() ; i++)
     {
       m_subsetIdx_data[i] = subset_idx.at(i);
     }
  }
  /** Copy Constructor */
  SubMap(const SubMap& otherMap)
    : StridePolicyType(otherMap)
    , m_superMap_constptr(otherMap.m_superMap_constptr)
    , m_superMap_ptr(otherMap.m_superMap_ptr)
  {
     //copySubset<SetType>(otherMap);
  }

  /** Assignment Operator */
  SubMap& operator=(const SubMap& otherMap)
  {
    StridePolicyType::operator=(otherMap);
    m_superMap_constptr = otherMap.m_superMap_constptr;
    m_superMap_ptr = otherMap.m_superMap_ptr;
    //copySubset<SetType>(otherMap);

    return *this;
  }

//  /** Move Constructor */
//  SubMap(SubMap&& otherMap)
//    : StridePolicyType(otherMap),
//    , m_superMap_constptr(otherMap.m_superMap_constptr),
//    , m_superMap_ptr(otherMap.m_superMap_ptr),
//    , m_subsetIdx_data(std::move(otherMap.m_subsetIdx_data)),
//    , m_subsetIdx( SubsetBuilder()
//                 .size(m_subsetIdx_data.size())
//                 //.data(&m_subsetIdx_data)
//                 )
//  { }
//
//  /** Move Assignment Operator */
//  SubMap& operator=(SubMap&& otherMap)
//  {
//    StridePolicyType::operator=(otherMap);
//    m_superMap_constptr = otherMap.m_superMap_constptr;
//    m_superMap_ptr = otherMap.m_superMap_ptr,
//    m_subsetIdx_data = std::move(otherMap.m_subsetIdx_data);
//    m_subsetIdx = SubsetBuilder()
//                  .size(m_subsetIdx_data.size())
//                  //.data(&m_subsetIdx_data)
//                  ;
//    return *this;
//  }


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
  const DataType & operator[](IndexType idx) const
  {
    verifyPositionImpl(idx);
    IndexType flat_idx = getMapCompFlatIndex(idx);
    return (*m_superMap_constptr)[flat_idx];
  }

  //non-const version
  DataType & operator[](IndexType idx)
  {
    verifyPositionImpl(idx);
    SLIC_ASSERT_MSG(m_superMap_ptr != nullptr,
                    "Submap was constructed with const Map pointer, "
                    << "non-const functions should not be called");

    IndexType flat_idx = getMapCompFlatIndex(idx);
    return (*m_superMap_ptr)[flat_idx];
  }

  /**
   * \brief Access the value associated with the given position in the subset
   *        and the component index.
   *
   * \pre `0 <= idx < size()`
   * \pre `0 <= comp < numComp()`
   */
  const DataType & operator()(IndexType idx, IndexType comp = 0) const
  {
    verifyPositionImpl(idx, comp);
    return (*m_superMap_constptr)[getMapElemFlatIndex(idx)*numComp() + comp];
  }

  //non-const version
  DataType & operator()(IndexType idx, IndexType comp = 0)
  {
    verifyPositionImpl(idx, comp);
    SLIC_ASSERT_MSG(m_superMap_ptr != nullptr,
                    "Submap was constructed with const Map pointer, "
                    << "non-const functions should not be called");

    return (*m_superMap_ptr)[getMapElemFlatIndex(idx)*numComp() + comp];
  }

  /**
   * \brief Access the value associated with the given position in the subset
   *        and the component index.
   *
   * \pre `0 <= idx < size()`
   * \pre `0 <= comp < numComp()`
   */
  const DataType & value(IndexType idx, IndexType comp = 0) const
  {
    return operator()(idx, comp);
  }

  DataType & value(IndexType idx, IndexType comp = 0)
  {
    return operator()(idx, comp);
  }

  /**
   * \brief Return the set element in the SuperMap at the given subset index
   */
  IndexType index(IndexType idx) const
  {
    return m_superMap_constptr->set()->at( m_subsetIdx[idx] );
  }

  /// @}


  /// \name SubMap cardinality functions
  /// @{
  ///

  /** \brief returns the size of the SubMap  */
  IndexType size() const override
  {
    return m_subsetIdx.size();
  }

  /** \brief returns the number of components (aka. stride) of the SubMap  */
  IndexType numComp() const { return StridePolicyType::stride(); }

  /// @}


  bool isValid(bool VerboseOutput = false) const override;

private: //function inherit from StridePolicy that should not be accessible
  void setStride(IndexType)
  {
    SLIC_ASSERT_MSG(false,
                    "Stride should not be changed after construction of SubMap.");
  }

private: //helper functions

  /**
   * \brief Get the ElementFlatIndex into the SuperMap given the subset's index.
   */
  IndexType getMapElemFlatIndex(IndexType idx) const
  {
    return m_subsetIdx[idx];
  }

  /*
   * \brief Get the ComponentFlatIndex into the SuperMap given the subset's
   * ComponentFlatIndex. This is used only with bracket [] access
   */
  IndexType getMapCompFlatIndex(IndexType idx) const
  {
    IndexType comp = numComp();
    IndexType s = idx % comp;
    return getMapElemFlatIndex(idx / comp) * comp + s;
  }

  /** Checks the ComponentFlatIndex is valid */
  void verifyPosition(SetPosition idx ) const override
  {
    verifyPositionImpl(idx);
  }

  /** Checks the ElementFlatIndex and the component index is valid */
  void verifyPosition(SetPosition idx,
                      SetPosition comp ) const
  {
    verifyPositionImpl(idx, comp);
  }

  /** Checks the ComponentFlatIndex is valid */
  void verifyPositionImpl(SetPosition AXOM_DEBUG_PARAM(idx) ) const
  {
    SLIC_ASSERT_MSG(idx >= 0 && idx < m_subsetIdx.size()*numComp(),
                    "Attempted to access element "
                    << idx << " but Submap's data has size "
                    << m_subsetIdx.size() * numComp() );
  }

  /** Checks the ElementFlatIndex and the component index is valid */
  void verifyPositionImpl(SetPosition AXOM_DEBUG_PARAM(idx),
                          SetPosition AXOM_DEBUG_PARAM(comp) ) const
  {
    SLIC_ASSERT_MSG(idx >= 0 && idx < m_subsetIdx.size() &&
                    comp >= 0 && comp < numComp(),
                    "Attempted to access element "
                    << idx << " component " << comp
                    << ", but Submap's data has size " << m_subsetIdx.size()
                    << " with " << numComp() << " component");
  }

public:

  /**
   * \class SubMapIterator
   * \brief An iterator for SubMap, based on MapIterator
   *
   * \see MapIterator
   */
  class SubMapIterator
    : public IteratorBase<SubMapIterator, SetPosition>
  {
  public:
    using iterator_category = std::random_access_iterator_tag;
    using value_type = DataType;
    using difference_type = SetPosition;

    using IterBase = IteratorBase<SubMapIterator, SetPosition>;
    using IterBase::m_pos;
    using iter = SubMapIterator;
    using PositionType = SetPosition;

    SubMapIterator(PositionType pos, SubMap* sMap)
      : IterBase(pos), m_submap(*sMap)
    {}

    /**
     * \brief Returns the current iterator value. If the SubMap has multiple
     *        components, this will return the first component. To access
     *        the other components, use iter(comp)
     */
    DataType & operator*()
    {
      return m_submap(m_pos, 0);
    }

    /**
     * \brief Returns the iterator's value at the specified component.
     *        Returns the first component if comp_idx is not specified.
     * \param comp_idx  (Optional) Zero-based index of the component.
     */
    DataType & operator()(IndexType comp= 0)
    {
      return m_submap(m_pos, comp);
    }

    /** \brief Returns the first component value after n increments.  */
    DataType & operator[](PositionType n)
    {
      return *(*this+n);
    }

    /** \brief Same as operator() */
    DataType& value(IndexType comp = 0)
    {
      return m_submap(m_pos, comp);
    }

    /** \brief Returns the Set element at the iterator's position */
    IndexType index()
    {
      return m_submap.index(m_pos);
    }

    /** \brief Returns the number of component per element in the SubMap. */
    PositionType numComp() const { return m_submap.numComp(); }

  protected:
    /* Implementation of advance() as required by IteratorBase */
    void advance(PositionType pos)
    {
      m_pos += pos;
    }

  private:
    SubMap m_submap;

  };

public:     // Functions related to iteration
  SubMapIterator begin() { return SubMapIterator(0, this); }
  SubMapIterator end() { return SubMapIterator(m_subsetIdx.size(), this); }


protected:     //Member variables
  const SuperMapType* m_superMap_constptr;
  SuperMapType* m_superMap_ptr;

  //Stores the ElementFlatIndex into the SuperMap
  std::vector<IndexType> m_subsetIdx_data;
  SubsetType m_subsetIdx;

}; //end SubMap




template<typename SuperMapType,typename SetType>
bool SubMap<SuperMapType,SetType>::isValid(bool verboseOutput) const
{
  bool isValid = true;
  std::stringstream errStr;

  if (m_superMap_constptr == nullptr)
  {
    if (m_subsetIdx.size() > 0)
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
    int map_size = ((const MapBase*)m_superMap_constptr)->size();
    //Check all indices is inside the SuperMap range
    for (int i = 0 ; i < m_subsetIdx.size() ; i++)
    {
      SetPosition pos = m_subsetIdx[i];
      if (pos < 0 || pos >= map_size)
      {
        isValid = false;
        if(verboseOutput)
        {
           errStr << "\n\t* The given subset index " << pos
                  << "is outside of the SuperMap range of 0 to "
                  << map_size;
        }
      }
    }
  }

  if(verboseOutput)
  {
    SLIC_INFO("Detailed results of isValid on the SubMap.\n"
       << "SubMap was " << (isValid ? "valid" : "NOT valid") << "\n"
       << errStr.str() );
  }

  return isValid;
}

} // end namespace slam
} // end namespace axom



#endif // SLAM_SUBMAP_HPP_
