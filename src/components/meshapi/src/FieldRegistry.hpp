#ifndef MESHAPI_FIELD_REGISTRY_H_
#define MESHAPI_FIELD_REGISTRY_H_


#include <sstream>

#include "slic/slic.hpp"

#include "meshapi/Utilities.hpp"
#include "meshapi/Set.hpp"
#include "meshapi/Map.hpp"

namespace asctoolkit {
namespace meshapi {



/**
 * \brief Very simple container for fields of a given type DataType with minimal error checking.
 * \note
 *         We are using concrete instances for int and double in the code below.
 *         This should eventually be replaced with the sidre datastore.
 */
  template<typename DataType>
  class FieldRegistry
  {
  public:
    typedef std::string                         KeyType;
    typedef asctoolkit::meshapi::Map<DataType>  MapType;

    typedef std::map<KeyType, MapType>          DataVecMap;
    typedef std::map<KeyType, DataType>         DataAttrMap;

  public:
    MapType&  addField(KeyType key, Set const* theSet) { return m_dataVecs[key] = MapType(theSet); }
    DataType& addScalar(KeyType key, DataType val)     { return m_dataScalars[key] = val; }

    MapType&  getField(KeyType key)
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
          , "Didn't find " << asctoolkit::meshapi::util::TypeToString<DataType>::to_string() << " field named " << key );
    }
    inline void verifyScalarsKey(KeyType key){
      SLIC_ASSERT_MSG( m_dataScalars.find(key) != m_dataScalars.end()
          , "Didn't find " << asctoolkit::meshapi::util::TypeToString<DataType>::to_string() << " scalar named " << key );
    }
  private:
    DataVecMap m_dataVecs;
    DataAttrMap m_dataScalars;
  };
} // end namespace meshapi
} // end namespace asctoolkit

#endif // MESHAPI_FIELD_REGISTRY_H_
