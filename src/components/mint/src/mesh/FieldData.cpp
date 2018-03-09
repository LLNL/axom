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

#include "FieldData.hpp"

// axom includes
#include "mint/Field.hpp"
#include "slic/slic.hpp"


namespace axom
{
namespace mint
{

FieldData::FieldData()
{}

//------------------------------------------------------------------------------
FieldData::~FieldData()
{
  this->clear();
}

//------------------------------------------------------------------------------
bool FieldData::hasField( const std::string& name ) const
{
  bool status = false;

  if ( m_container.find( name ) != m_container.end() )
  {
    status = true;
  }

  return status;
}

//------------------------------------------------------------------------------
void FieldData::addField( Field* f )
{
  SLIC_ASSERT(  f != AXOM_NULLPTR );
  SLIC_ASSERT(  this->hasField( f->getName() )==false );

  m_fields.push_back( f->getName() );
  m_container[ f->getName() ] = f;

  SLIC_ASSERT( m_fields.size() == m_container.size() );
}

//------------------------------------------------------------------------------
int FieldData::getNumberOfFields() const
{
  SLIC_ASSERT( m_fields.size() == m_container.size() );
  return static_cast< int >( m_fields.size() );
}

//------------------------------------------------------------------------------
Field* FieldData::getField( int i )
{
  SLIC_ASSERT( i >= 0 && i < this->getNumberOfFields() );

  if ( i < 0 || i >= this->getNumberOfFields() )
  {
    return AXOM_NULLPTR;
  }

  return this->getField( m_fields[ i ] );
}

//------------------------------------------------------------------------------
const Field* FieldData::getField( int i ) const
{
  return const_cast< const Field* >(
    const_cast< FieldData* >( this )->getField( i ) );
}

//------------------------------------------------------------------------------
Field* FieldData::getField( const std::string& name )
{
  SLIC_ASSERT( this->hasField( name ) );
  return m_container[ name ];
}

//------------------------------------------------------------------------------
const Field* FieldData::getField( const std::string& name ) const
{
  return const_cast< const Field* >(
    const_cast< FieldData* >( this )->getField( name ) );
}

//------------------------------------------------------------------------------
void FieldData::clear()
{
  std::map< std::string, Field* >::iterator iter = m_container.begin();
  for ( ; iter != m_container.end() ; ++iter )
  {
    delete iter->second;
    iter->second = AXOM_NULLPTR;
  }
  m_container.clear();
}

//------------------------------------------------------------------------------
bool FieldData::empty() const
{
  return( m_container.empty( ) );
}

} /* namespace mint */
} /* namespace axom */
