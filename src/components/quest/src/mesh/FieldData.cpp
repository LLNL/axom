/*
 * $Id$
 */

/*!
 *******************************************************************************
 * \file FieldData.cpp
 *
 * \date Sep 19, 2015
 * \author George Zagaris (zagaris2@llnl.gov)
 *******************************************************************************
 */

#include "FieldData.hpp"

// ATK includes
#include "quest/Field.hpp"
#include "slic/slic.hpp"

// C/C++ includes

namespace meshtk {


FieldData::FieldData()
{

}

//------------------------------------------------------------------------------
FieldData::~FieldData()
{
   this->clear();
}

//------------------------------------------------------------------------------
bool FieldData::hasField( const std::string& name )
{
   bool status = false;

   if ( m_container.find( name ) != m_container.end() ) {
       status = true;
   }

   return status;
}

//------------------------------------------------------------------------------
void FieldData::addField( Field* f )
{
   SLIC_ASSERT( f != ATK_NULLPTR );
   SLIC_ASSERT( this->hasField( f->getName() )==false );

   m_fields.push_back( f->getName() );
   m_container[ f->getName() ] = f;

   SLIC_ASSERT( m_fields.size()==m_container.size() );
}

//------------------------------------------------------------------------------
int FieldData::getNumberOfFields() const
{
   SLIC_ASSERT( m_fields.size()==m_container.size() );
   return static_cast< int >( m_fields.size() );
}

//------------------------------------------------------------------------------
Field* FieldData::getField( int i )
{
   SLIC_ASSERT( i >= 0 && i < this->getNumberOfFields() );

   if ( i < 0 || i >= this->getNumberOfFields() ) {
     return ATK_NULLPTR;
   }

   return this->getField( m_fields[ i ] );
}

//------------------------------------------------------------------------------
Field* FieldData::getField( const std::string& name )
{
   SLIC_ASSERT( this->hasField( name ) );
   return m_container[ name ];
}

//------------------------------------------------------------------------------
void FieldData::clear()
{
   std::map< std::string, Field* >::iterator iter = m_container.begin();
   for ( ; iter != m_container.end(); ++iter ) {
      delete iter->second;
      iter->second = ATK_NULLPTR;
   }
   m_container.clear();
}

//------------------------------------------------------------------------------
bool FieldData::empty() const
{
   return( m_container.empty( ) );
}

} /* namespace meshtk */
