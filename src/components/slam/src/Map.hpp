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

#include "slam/Set.hpp"
#include "slam/NullSet.hpp"

namespace axom
{
namespace slam
{

class NullSet;


// This class is missing some simplifying copy constructors
// -- or at least ways of interacting with the data store
// We should probably support copy on write shallow copies when possible...


template<typename DataType>
class Map
{
public:
  typedef Set::IndexType SetIndex;
  typedef Set::PositionType SetPosition;

  typedef std::vector<DataType> OrderedMap;

  static const NullSet s_nullSet;

public:
  Map(const Set* theSet = &s_nullSet) : m_set(theSet)
  {
    m_data.resize( m_set->size());
  }

  Map(const Set* theSet, DataType defaultValue) : m_set(theSet)
  {
    m_data.resize( m_set->size(), defaultValue );
  }

  /**
   * Copy constructor from another map
   */
  Map(const Map& otherMap) : m_set(otherMap.m_set)
  {
    m_data.resize( m_set->size());
    copy( otherMap);
  }

  /** Assignment operator for Map     */
  Map& operator=(const Map& otherMap)
  {
    if(this != &otherMap)
    {
      m_set = otherMap.m_set;
      m_data = otherMap.m_data;
    }

    return *this;
  }

  ~Map(){}

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


  Set const* set() const { return m_set; }


  SetPosition size() const { return m_set->size(); }

  bool        isValid(bool verboseOutput = false) const;


  void        clear() { fill(); }

  /**
   * Set each entry in the map to the given value
   */
  void        fill(DataType val = DataType())
  {
    const SetIndex sz = size();

    for(SetIndex idx = SetIndex() ; idx < sz ; ++idx)
    {
      m_data[idx] = val;
    }
  }

  /**
   * Element-wise copy of data from another map
   */
  void copy(const Map& other)
  {
    SLIC_ASSERT( other.size() == size() );

    const SetIndex sz = size();
    for(SetIndex idx = SetIndex() ; idx < sz ; ++idx)
    {
      m_data[idx] = other[idx];
    }
  }


public:
  /**
   * \name DirectDataAccess
   * \brief Accessor functions to get the underlying map data

   * \note We will have to figure out a good way to limit
   * this access to situations where it makes sense.
   */

  /// \{

  //* Placeholder for function that returns the (pointer to) underlying data **/
  OrderedMap &        data()        { return m_data; }
  //* Placeholder for function that returns the (const pointer to) underlying
  // data **/
  const OrderedMap &  data() const { return m_data; }

  /// \}

private:
  inline void verifyPosition(SetPosition AXOM_DEBUG_PARAM(setIndex))       const
  {
    SLIC_ASSERT_MSG(
      setIndex >= 0 && setIndex < m_set->size(),
      "Attempted to access element "
      << setIndex << " but map's set has size "  << m_set->size() );
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
template<typename DataType>
NullSet const Map<DataType>::s_nullSet;

template<typename DataType>
bool Map<DataType>::isValid(bool verboseOutput) const
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
    if( static_cast<SetPosition>(m_data.size()) != m_set->size())
    {
      if(verboseOutput)
      {
        errStr << "\n\t* the underlying set and its associated mapped data"
               << " have different sizes"
               << " , underlying set has size " << m_set->size()
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
      sstr << "Map was NOT valid.\n"
           << sstr.str()
           << std::endl;
    }

    if(!m_set)
    {
      sstr << "\n** map is empty.";
    }
    else
    {
      sstr << "\n** underlying set has size " << m_set->size() << ": ";

      sstr << "\n** Mapped data:";
      for(SetPosition idx = 0 ; idx < this->size() ; ++idx)
      {
        sstr << "\n\telt[" << idx << "]:\t" << (*this)[idx];
      }
    }
    std::cout << sstr.str() << std::endl;
  }

  return bValid;
}




} // end namespace slam
} // end namespace axom



#endif // SLAM_MAP_HPP_
