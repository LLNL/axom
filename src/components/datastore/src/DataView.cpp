/*
 * DataSet.cpp
 *
 *  Created on: Dec 1, 2014
 *      Author: settgast
 */

#include "DataView.hpp"
#include "DataGroup.hpp"

namespace DataStoreNS
{


DataView::DataView( const std::string& name,
                    DataGroup* const parentGroup,
                    DataBuffer* const dataBuffer ):
    m_name(name),
    m_parentGroup(parentGroup),
    m_dataBuffer(dataBuffer),
    m_viewStart(nullptr),
    m_dataShape(),
    m_dataType(rtTypes::undefined)
{}

DataView::DataView( const std::string& name,
                    DataGroup* const parentGroup,
                    DataStore* const dataStore ) :
    m_uid(uid),
    m_GroupSet(),
    m_viewStart(nullptr),
    m_dataShape(),
    m_dataType(rtTypes::undefined)
{}

DataView::DataView(const DataView& source ) :
    m_uid(source.m_uid),
    m_GroupSet(source.m_GroupSet),
    m_viewStart(source.m_viewStart),
    m_dataShape(source.m_dataShape),
    m_dataType(source.m_dataType)
{
}


DataView::~DataView()
{
}



DataView* DataView::Allocate()
{
  if ( m_dataShape.m_dimensions != nullptr && m_dataType!=rtTypes::undefined )
  {
    std::size_t size = 1;
    for (int dim = 0; dim < m_dataShape.m_numDimensions; ++dim)
    {
      size *= m_dataShape.m_dimensions[dim];
    }
//    m_memblob.resize( size * rtTypes::sizeofType(m_dataType) );
//    m_data = m_memblob.data();
  }
  else
  {
    throw std::exception();
  }

  return this;
}

DataView* DataView::SetLength(const std::size_t newsize)
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
