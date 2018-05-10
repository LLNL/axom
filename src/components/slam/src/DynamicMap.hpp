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

#ifndef SLAM_DYNAMIC_MAP_HPP_
#define SLAM_DYNAMIC_MAP_HPP_

#include <vector>
#include <sstream>
#include <iostream>

#include "axom/Macros.hpp"
#include "axom/Types.hpp"
#include "slic/slic.hpp"

#include "slam/DynamicSet.hpp"
#include "slam/Map.hpp"

namespace axom
{
namespace slam
{


/**
 * \brief A map that supports adding and removing entries.
 *
 * \detail An entry in the map is considered valid if
 * its corresponding set's entry is valid
 */
template<typename DataType>
class DynamicMap : public Map<DataType>
{

public:
  typedef Set::IndexType SetIndex;
  typedef Set::PositionType SetPosition;

  typedef std::vector<DataType> OrderedMap;

public:
  DynamicMap() : m_set(AXOM_NULLPTR){}

  DynamicMap(DynamicSet<>* theSet ) : m_set(theSet)
  {
    m_data.resize( m_set->size());
  }

  DynamicMap(DynamicSet<>* theSet, DataType defaultValue) : m_set(theSet)
  {
    m_data.resize( m_set->size(), defaultValue );
  }

  ~DynamicMap(){}

  const DataType & operator[](SetPosition setIndex) const
  {
    verifyPosition(setIndex);
    return m_data[setIndex];
  }

  SetPosition numberOfValidEntries() const
  {
    return m_set->numberOfValidEntries();
  }

  bool isValidEntry(SetPosition pos) const
  {
    if( m_set == AXOM_NULLPTR )
      return false;

    return m_set->isValidEntry( pos );
  }

  SetPosition size() const {
    return m_data.size();
  }

  OrderedMap & data(){
    return m_data;
  }

  const OrderedMap & data() const {
    return m_data;
  }

  bool isValid(bool verboseOutput = false) const;

private:

  inline void verifyPosition(SetPosition AXOM_DEBUG_PARAM(setIndex))       const
  {
    SLIC_ASSERT_MSG(
      setIndex >= 0 && setIndex < (int)m_data.size(),
      "Attempted to access entry "
      << setIndex << " but map's set has size " << m_data.size() );

    SLIC_ASSERT_MSG(
      m_set != AXOM_NULLPTR && m_set->isValidEntry( setIndex),
      "Attempted to access entry "
      << setIndex << " but the set entry is invalid" );
  }

public:
  /* Modifying functions */
  DataType & operator[](SetPosition position)
  {
    if((int)m_data.size() < position+1 )
    {
      resize( position + 1 );
    }

    return m_data[position];
  }

  void insert(SetPosition position, DataType value){
    operator[]( position ) = value;
  }

  void resize(SetPosition s){
    SLIC_ASSERT_MSG( s >= 0,
                     "Attempted to resize vector with a negative size " << s );
    m_data.resize(s);
  }

private:
  DynamicSet<>* m_set;
  OrderedMap m_data;

};


template< typename T>
bool DynamicMap<T>::isValid(bool verboseOutput) const
{
  bool bValid = true;

  std::stringstream errStr;

  if(m_set == AXOM_NULLPTR)
  {
    if(!m_data.empty() )
    {
      if(verboseOutput)
      {
        errStr
          << "\n\t* the underlying set was never provided,"
          << " but its associated data is not empty"
          << " , data has size " << m_data.size();
      }
      bValid = false;
    }
  }
  else
  {

    // Check the data array and set data have equal size
    if( static_cast<SetPosition>(m_data.size()) != m_set->size())
    {
      if(verboseOutput)
      {
        errStr
          <<"\n\t* the underlying set and its associated mapped data"
          <<" have different sizes, underlying set has size "
          << m_set->size() << " , data has size " << m_data.size();
        ;
      }

      bValid = false;
    }

  }

  if(verboseOutput)
  {
    if(bValid)
    {
      SLIC_DEBUG( "Map was valid." );
    }
    else
    {
      SLIC_DEBUG( "Map was not valid. " << errStr.str() );
    }
  }

  return bValid;
}


} // end namespace slam
} // end namespace axom

#endif // SLAM_DYNAMIC_MAP_HPP_
