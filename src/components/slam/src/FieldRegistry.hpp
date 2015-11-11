/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */


#ifndef SLAM_FIELD_REGISTRY_H_
#define SLAM_FIELD_REGISTRY_H_


#include <sstream>

#include "slic/slic.hpp"

#include "slam/Utilities.hpp"
#include "slam/Set.hpp"
#include "slam/Map.hpp"

namespace asctoolkit {
namespace slam {



/**
 * \brief Very simple container for fields of a given type DataType with minimal error checking.
 * \note
 *         We are using concrete instances for int and double in the code below.
 *         This should eventually be replaced with the sidre datastore.
 */
  template<typename TheDataType>
  class FieldRegistry
  {
  public:
    typedef TheDataType                     DataType;
    typedef std::string                     KeyType;
    typedef asctoolkit::slam::Map<DataType> MapType;
    typedef typename MapType::OrderedMap    BufferType;

    typedef std::map<KeyType, MapType>      DataVecMap;
    typedef std::map<KeyType, DataType>     DataAttrMap;

  public:
    MapType&  addField(KeyType key, Set const* theSet) { return m_dataVecs[key] = MapType(theSet); }
    DataType& addScalar(KeyType key, DataType val)     { return m_dataScalars[key] = val; }

    MapType&  addNamelessField(Set const* theSet)
    {
      static int cnt = 0;
      std::stringstream key;

      key << "__field_" << cnt++;
      return m_dataVecs[key.str()] = MapType(theSet);
    }


    MapType& getField(KeyType key)
    {
      verifyFieldsKey(key);
      return m_dataVecs[key];
    }
    MapType const& getField(KeyType key) const
    {
      verifyFieldsKey(key);
      return m_dataVecs[key];
    }

    DataType& getScalar(KeyType key)
    {
      verifyScalarsKey(key);
      return m_dataScalars[key];
    }
    DataType const& getScalar(KeyType key) const
    {
      verifyScalarsKey(key);
      return m_dataScalars[key];
    }

  private:
    inline void verifyFieldsKey(KeyType key){
      SLIC_ASSERT_MSG( m_dataVecs.find(key) != m_dataVecs.end()
          , "Didn't find " << asctoolkit::slam::util::TypeToString<DataType>::to_string() << " field named " << key );
    }
    inline void verifyScalarsKey(KeyType key){
      SLIC_ASSERT_MSG( m_dataScalars.find(key) != m_dataScalars.end()
          , "Didn't find " << asctoolkit::slam::util::TypeToString<DataType>::to_string() << " scalar named " << key );
    }
  private:
    DataVecMap m_dataVecs;
    DataAttrMap m_dataScalars;
  };
} // end namespace slam
} // end namespace asctoolkit

#endif // SLAM_FIELD_REGISTRY_H_
