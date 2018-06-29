/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-741217
 *
 * All rights reserved.
 *
 * This file is part of Axom.
 *
 * For details about use and distribution, please read axom/LICENSE.
 *
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

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

#include "axom/Macros.hpp"
#include "axom/Types.hpp"
#include "slic/slic.hpp"

#include "slam/MapBase.hpp"
#include "slam/Set.hpp"
#include "slam/NullSet.hpp"

#include "slam/StridePolicies.hpp"

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
 * \brief   A Map class that associates a value to every element in a set
 *          
 * \detail  The Map will take a Set pointer, and maintain a value of DataType 
 *          for each element in the Set. If the stride of the map is larger
 *          than 1, then each element of the set will be mapped to stride 
 *          number of DataType
 *          
 *          Accessing the element with map(i,j) will give you the ith element,
 *          jth component. Accessing the map with map[k] will give you the kth 
 *          element in the underhood array s.t. i*stride+j=k
 *
 */

template<
  typename DataType, 
  typename StridePolicy = policies::StrideOne<Set::IndexType>
>
class Map : public MapBase, public StridePolicy
{
public:

  using OrderedMap = std::vector<DataType>;
  using StridePolicyType = StridePolicy;

  static const NullSet s_nullSet;

  class MapBuilder;
  
#ifdef AXOM_USE_CXX11
  class MapIterator;
  using const_iterator = MapIterator;
  using const_iterator_pair = std::pair<const_iterator, const_iterator>;
  using iterator = const_iterator;
  using iterator_pair = const_iterator_pair;
#endif

public:

  /**
   * \brief Constructor for Map using a Set pointer
   *
   * \param theSet         (Optional) A pointer to the map's set
   * \param defaultValue   (Optional) If given, every entry in the map will be 
   *                       initialized using defaultValue
   * \param stride         (Optional) The stride. The number of DataType that 
   *                       each element in the set will be mapped to. 
   *                       Default is 1.
   */
  
  Map(const Set* theSet = &s_nullSet, DataType defaultValue = DataType(),
    Set::IndexType stride = StridePolicyType::DEFAULT_VALUE ) :
    StridePolicyType(stride), m_set(theSet)
  {
    m_data.resize( m_set->size() * StridePolicy::stride(), defaultValue);
  }

  /**
   * Copy constructor from another map
   */
  Map(const Map& otherMap) : 
    StridePolicyType(otherMap.StridePolicyType::stride()),
    m_set(otherMap.m_set)
  {
    m_data.resize( otherMap.m_data.size() );
    copy( otherMap );
  }

  /**
   * \brief Constructor for Map using a MapBuilder
   */
  Map(const MapBuilder& builder)
    : Map(builder.m_set, builder.m_defaultValue, builder.m_stride.stride())
  {
    //copy the data if exists
    if (builder.m_data_ptr) 
    {
      for (SetIndex idx = SetIndex(); idx < builder.m_set->size(); ++idx)
      {
        m_data[idx] = builder.m_data_ptr[idx];
      }
    }
  }

  /** 
   * \brief Assignment operator for Map     
   */
  Map& operator=(const Map& otherMap)
  {
    if(this != &otherMap)
    {
      m_set = otherMap.m_set;
      m_data = otherMap.m_data;
      StridePolicy::operator=(otherMap);
    }
    return *this;
  }

  //~Map(){}

  /**
   * \brief Returns a pointer to the map's underlying set
   */
  const Set* set() const { return m_set; }

  /// \name Map individual access functions
  /// @{
  ///

  /**
   * \brief Access the element given a flat index into the size()*numComp() range of the map
   */
  const DataType & operator[](SetPosition setIndex) const
  {
    verifyPosition(setIndex);
    return m_data[setIndex];
  }

  DataType & operator[](SetPosition setIndex)
  {
    verifyPosition(setIndex);
    return m_data[setIndex];
  }

  /**
   * \brief Access the element given a position in the set and a component index.   
   * 
   * \pre setIdx must be between 0 and size()
   * \pre comp must be between 0 and numComp()
   */
  const DataType & operator()(SetPosition setIdx, SetPosition comp = 0) const
  {
    SetPosition setIndex = setIdx * StridePolicyType::stride() + comp;
    verifyPosition(setIndex);
    return m_data[setIndex];
  }
  
  DataType & operator()(SetPosition setIdx, SetPosition comp = 0)
  {
    SetPosition setIndex = setIdx * StridePolicyType::stride() + comp;
    verifyPosition(setIndex);
    return m_data[setIndex];
  }

  /// @}


  /// \name Map cardinality functions
  /// @{

  /**
   * \brief Return the set size. Same as set()->size(). 
   * 
   * The total storage size would be size() * numComp()
   */
  SetPosition size() const { return SetPosition(m_set->size()); }

  /* 
   * Return the stride, equivalent to the number of components of the map
   */
  SetPosition numComp() const { return StridePolicyType::stride(); }

  /// @}


  /// \name Map modifying functions for all entries
  /// @{

  /** \brief replace all elements in the Map with the default DataType */
  void        clear() { fill(); }

  /** Set each entry in the map to the given value  */
  void        fill(DataType val = DataType())
  {
    const SetIndex sz = m_data.size();

    for(SetIndex idx = SetIndex() ; idx < sz ; ++idx)
    {
      m_data[idx] = val;
    }
  }

  /** \brief Element-wise copy of data from another map */
  void copy(const Map& other)
  {
    SLIC_ASSERT( other.size() == size() );
    SLIC_ASSERT( other.stride() == StridePolicyType::stride() );

    const SetIndex sz = size() * StridePolicyType::stride();
    for(SetIndex idx = SetIndex() ; idx < sz ; ++idx)
    {
      m_data[idx] = other[idx];
    }
  }

  ///@}
  


  /** \brief print information on the map, including every element inside Map  */
  void        print() const;

  /** \brief returns true if the map is valid, false otherwise.  */
  bool        isValid(bool verboseOutput = false) const;

public:
  /**
   * \class MapBuilder
   * \brief Helper class for constructing a Map
   **/
  class MapBuilder
  {
  public:
    friend class Map;

    MapBuilder(): m_set(&s_nullSet) {}

    /** \brief Provide the Set to be used by the Map */
    MapBuilder& set(const Set* set)
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
    const Set* m_set = &s_nullSet;
    StridePolicyType m_stride;
    DataType* m_data_ptr = nullptr;
    DataType m_defaultValue = DataType();
  };



#ifdef AXOM_USE_CXX11
  /**
  * \class MapIterator
  * \brief An iterator type for a map. It will increment/decrement by the 
  *        number of stride. To access the components, do iter(c) where c 
  *        is the component you want to access.
  */
  class MapIterator : public std::iterator<std::random_access_iterator_tag, DataType>
  {
  public:
    using iter = MapIterator;
    using PositionType = SetPosition;

  public:

    MapIterator(PositionType pos, Map* oMap)
      : m_pos(pos), m_mapPtr(oMap) {}

    /**
     * \brief Returns the current iterator value. If the map has multiple 
     *        components, this will return the first component. To access
     *        the other components, use iter(comp)
     */
    DataType & operator*()    
    {
      return (*m_mapPtr)(m_pos, 0);
    }

    /**
     * \brief Returns the iterator's value at the specified component. 
     *        Returns the first component if comp_idx is not specified.
     * \param comp_idx  (Optional) Zero-based index of the component.
     */
    DataType & operator()(SetPosition comp_idx = 0)
    {
      return (*m_mapPtr)(m_pos, comp_idx);
    }

    bool operator==(const iter& other) const
    {
      return (m_mapPtr == other.m_mapPtr) && (m_pos == other.m_pos);
    }

    bool operator!=(const iter& other) const { return !operator==(other); }
    bool operator<(const iter& other) const { return m_pos < other.m_pos; }

    iter& operator++() { advance(1); return *this; }
    iter operator++(int) { iter ret = *this; advance(1); return ret; }
    iter& operator--() { advance(-1); return *this; }
    iter operator--(int) { iter ret = *this; advance(-1); return ret; }

    iter& operator+=(PositionType n) { advance(n); return *this; }
    iter& operator-=(PositionType n) { advance(-n); return *this; }

    iter operator+(PositionType n) const
    {
      iter ret = *this; ret.advance(n);
      return ret;
    }

    iter operator-(PositionType n) const
    {
      iter ret = *this; ret.advance(-n);
      return ret;
    }

    /** \brief Returns the first component value after n increments.  */
    const DataType & operator[](PositionType n) const
    {
      return *(this->operator+(n));
    }

    DataType & operator[](PositionType n) 
    {
      return *(this->operator+(n));
    }

    friend PositionType operator-(const iter& a, const iter& b)
    {
      return (a.m_pos - b.m_pos);
    }

    /** \brief Returns the number of component the map has. */
    PositionType numComp() const { return m_mapPtr->stride(); }

  private:
    void advance(PositionType n) { m_pos += n; }

  protected:
    PositionType m_pos;
    Map* const m_mapPtr;
  };

public:     // Functions related to iteration

  MapIterator         begin()   {    return MapIterator(0, this);  }
  MapIterator         end()     {    return MapIterator(size(), this); }
  const_iterator_pair range() const { return std::make_pair(begin(), end()); }

#endif //AXOM_USE_CXX11




public:
  /**
   * \name DirectDataAccess
   * \brief Accessor functions to get the underlying map data

   * \note We will have to figure out a good way to limit
   * this access to situations where it makes sense.
   */

  /// \{

  //* Placeholder for function that returns the (pointer to) underlying data **/
  OrderedMap &        data()       { return m_data; }
  //* Placeholder for function that returns the (const pointer to) underlying
  // data **/
  const OrderedMap &  data() const { return m_data; }

  /// \}

private:
  inline void verifyPosition(SetPosition AXOM_DEBUG_PARAM(setIndex))      const
  {
    SLIC_ASSERT_MSG(
      setIndex >= 0 && setIndex < SetPosition( m_data.size()),
      "Attempted to access element "
      << setIndex << " but map's data has size "  << m_data.size() );
  }

  // setStride function should not be called after constructor is called.
  // This (should) override the StridePolicy setStride(s) function.
  void setStride(SetPosition str)
  {
    SLIC_ASSERT_MSG(false, "Stride should not be changed after construction of map.");
  }

private:
  const Set* m_set;
  OrderedMap m_data;
};



/**
 * \brief Definition of static instance of nullSet for all maps
 * \note Should this be a singleton or a global object?  Should the scope be
 * public?
 */
template<typename DataType, typename StridePolicy>
NullSet const Map<DataType, StridePolicy>::s_nullSet;


template<typename DataType, typename StridePolicy>
bool Map<DataType, StridePolicy>::isValid(bool verboseOutput) const
{
  bool bValid = true;

  std::stringstream errStr;

  if(*m_set == s_nullSet)
  {
    if(!m_data.empty() )
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
    if( static_cast<SetPosition>(m_data.size()) != m_set->size() * StridePolicyType::stride())
    {
      if(verboseOutput)
      {
        errStr << "\n\t* the underlying set and its associated mapped data"
               << " have different sizes"
               << " , underlying set has size " << m_set->size()
               << " with stride " << StridePolicy::stride()
               << " , data has size " << m_data.size();
      }

      bValid = false;
    }
  }


  if(verboseOutput)
  {
    std::stringstream sstr;

    if(bValid)
    {
      sstr << "Map was valid." << std::endl;
    }
    else
    {
      sstr << "Map was NOT valid.\n"
           << sstr.str()
           << std::endl;
    }

    std::cout << sstr.str() << std::endl;
  }

  return bValid;
}


template<typename DataType, typename StridePolicy>
void Map<DataType, StridePolicy>::print() const
{
  bool valid = isValid(true);
  std::stringstream sstr;

  if (valid)
  {
    if (!m_set)
    {
      sstr << "** map is empty.";
    }
    else
    {
      sstr << "** underlying set has size " << m_set->size() << ": ";
      sstr << "\n** the stride of the map is " << StridePolicy::stride() << ": ";

      sstr << "\n** Mapped data:";
      for (SetPosition idx = 0; idx < this->size(); ++idx)
      {
        for (SetPosition idx2 = 0; idx2 < StridePolicy::stride(); ++idx2)
        {
          sstr << "\n\telt[" << idx << "," << idx2 << "]:\t"
            << (*this)[idx*StridePolicyType::stride() + idx2];
        }
      }
    }
  }

  std::cout << sstr.str() << std::endl;
}


} // end namespace slam
} // end namespace axom



#endif // SLAM_MAP_HPP_
