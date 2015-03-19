/*
 * DataSet.cpp
 *
 *  Created on: Dec 1, 2014
 *      Author: settgast
 */

#include "DataView.hpp"
#include "DataGroup.hpp"
#include "DataStore.hpp"
#include "DataBuffer.hpp"

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
{

}

DataView::DataView( const std::string& name,
                    DataGroup* const parentGroup,
                    DataStore* const dataStore ) :
  m_name(name),
  m_parentGroup(parentGroup),
  m_dataBuffer(nullptr),
  m_viewStart(nullptr),
  m_dataShape(),
  m_dataType(rtTypes::undefined)
{
  m_dataBuffer = dataStore->CreateDataBuffer();
}

DataView::DataView(const DataView& source ) :
    m_name(source.m_name),
    m_parentGroup(source.m_parentGroup),
    m_dataBuffer(source.m_dataBuffer),
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

  m_dataBuffer->Allocate();
  return this;
}


DataView* DataView::SetDataShape( const DataShape& dataShape )
{
  // HACK
  m_dataBuffer->SetDataShape( dataShape );

  // check to see what conditions m_dataDescriptor can be set.
  m_dataShape = dataShape;
  m_viewStart = dataShape.m_dataPtr;
  return this;
}

void DataView::ReconcileWithBuffer()
{
  m_dataShape = m_dataBuffer->GetDataShape();
}

} /* namespace Datastore */
