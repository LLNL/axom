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


DataObject::DataObject( const IDType uid ) :
    m_uid(uid),
    m_stringDescriptor(),
    m_GroupSet(),
    m_data(nullptr),
    m_dataShape(),
    m_dataType(rtTypes::undefined),
    m_memblob()
{}

DataObject::DataObject( const IDType uid,
                        const std::string& stringDescriptor ) :
    m_uid(uid),
    m_stringDescriptor(stringDescriptor),
    m_GroupSet(),
    m_data(nullptr),
    m_dataShape(),
    m_dataType(rtTypes::undefined),
    m_memblob()
{}

DataObject::DataObject(const DataObject& source ) :
    m_uid(source.m_uid),
    m_stringDescriptor(source.m_stringDescriptor),
    m_GroupSet(source.m_GroupSet),
    m_data(source.m_data),
    m_dataShape(source.m_dataShape),
    m_dataType(source.m_dataType),
    m_memblob(source.m_memblob)
{
}


DataObject::~DataObject()
{
}



DataObject* DataObject::Allocate()
{
  if ( m_dataShape.m_dimensions != nullptr && m_dataType!=rtTypes::undefined )
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
