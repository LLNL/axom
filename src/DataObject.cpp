/*
 * DataSet.cpp
 *
 *  Created on: Dec 1, 2014
 *      Author: settgast
 */

#include "DataObject.hpp"
#include "DataGroup.hpp"

namespace DataStoreNS
{

DataObject::DataObject(const std::string& name,
                       const DataGroup* const parent) :
    m_parent(parent),
    m_name(name),
    m_data(nullptr),
    m_dataShape(),
    m_dataType(rtTypes::TypeID::undefined),
    m_groups()
//m_dataType(typeid(void*))
{
  if( m_parent!=nullptr )
  {
    m_dataShape = m_parent->GetDataShape();
  }

}

DataObject::DataObject(const DataObject& source) :
    m_parent(source.m_parent),
    m_name(source.m_name),
    m_data(source.m_data),
    m_dataShape(source.m_dataShape),
    m_dataType(source.m_dataType),
    m_groups(source.m_groups)
{

}

DataObject::DataObject( DataObject&& source) :
    m_parent(std::move(source.m_parent)),
    m_name(std::move(source.m_name)),
    m_data(source.m_data),
    m_dataShape(source.m_dataShape),
    m_dataType(std::move(source.m_dataType)),
    m_groups(std::move(source.m_groups))
{

}

DataObject::DataObject::~DataObject()
{
}
/*
DataObject* DataObject::SetAttribute(const Attribute& newAttribute)
{
  m_attributes.insert(std::make_pair(newAttribute.Name(), newAttribute));
  return this;
}

bool DataObject::HasAttribute(const std::string& attributeKey) const
{
  return m_attributes.count(attributeKey);
}

bool DataObject::DeleteAttribute(const std::string& attributeKey)
{
  bool rval;
  if (HasAttribute(attributeKey))
  {
    rval = true;
    m_attributes.erase(attributeKey);
  }
  else
  {
    rval = false;
  }
  return rval;
}
*/

DataObject* DataObject::Allocate()
{
  if ( m_dataShape.m_dimensions != nullptr && m_dataType!=rtTypes::TypeID::undefined )
  {
    std::size_t size = 1;
    for (int dim = 0; dim < m_dataShape.m_numDimensions; ++dim)
    {
      size *= m_dataShape.m_dimensions[dim];
    }
    m_memblob.resize( size * rtTypes::sizeofType(m_dataType) );
    m_data = m_memblob.data();
  }
  else
  {
    throw std::exception();
  }

  return this;
}

DataObject* DataObject::SetLength(const std::size_t newsize)
{
  if (m_dataShape.m_dimensions != nullptr)
  {
    m_dataShape.m_dimensions[0] = newsize;
    Allocate();
  }
  else
  {
    throw std::exception();
  }
  return this;
}

} /* namespace Datastore */
