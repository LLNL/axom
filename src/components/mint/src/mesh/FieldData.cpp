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

#include "axom/config.hpp"     // for axom compile-time definitions
#include "axom/Types.hpp"      // for AXOM_NULLPTR

#include "mint/config.hpp"     // for mint compile-time type definitions
#include "mint/FieldData.hpp"

#ifdef MINT_USE_SIDRE
#include "sidre/sidre.hpp"  // for sidre::Group, sidre::View
#endif

namespace axom
{
namespace mint
{

//------------------------------------------------------------------------------
// INTERNAL HELPER METHODS
//------------------------------------------------------------------------------
namespace internal
{

#ifdef MINT_USE_SIDRE
mint::Field* getFieldFromView( const std::string& name, sidre::View* view )
{
  SLIC_ASSERT( view != AXOM_NULLPTR );
  SLIC_ASSERT( !view->isEmpty() );

  using int32 = axom::common::int32;
  using int64 = axom::common::int64;

  mint::Field* f = AXOM_NULLPTR;

  switch( view->getTypeID() )
  {
  case sidre::INT32_ID:
    f = new mint::FieldVariable< int32 >( name, view );
    break;
  case sidre::INT64_ID:
    f = new mint::FieldVariable< int64 >( name, view );
    break;
  case sidre::FLOAT64_ID:
    f = new mint::FieldVariable< double >( name, view );
    break;
  case sidre::FLOAT32_ID:
    f= new mint::FieldVariable< float >( name, view );
    break;
  default:
    SLIC_ERROR( "Encountered unsupported type [" << view->getTypeID() << "]" );
  } // END switch

  SLIC_ERROR_IF( (f==AXOM_NULLPTR), "null field!" );
  return ( f );
}

//------------------------------------------------------------------------------
void removeFromSidre( sidre::Group* grp, const std::string& name )
{
  SLIC_ASSERT( grp != AXOM_NULLPTR );
  SLIC_ASSERT( grp->hasChildGroup( name ) );

  grp->destroyGroup( name );
}

#endif

}

//------------------------------------------------------------------------------
//  FIELDDATA IMPLEMENTATION
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
FieldData::FieldData( int association ) :
    m_association( association ),
    m_fields(),
    m_nextidx(),
    m_name2idx()
{
#ifdef MINT_USE_SIDRE
  m_fields_group = AXOM_NULLPTR;
#endif

  SLIC_ERROR_IF(
      (m_association < 0) || (m_association >= NUM_FIELD_ASSOCIATIONS),
      "Invalid field association!" );
}

//------------------------------------------------------------------------------
#ifdef MINT_USE_SIDRE
FieldData::FieldData( int association, sidre::Group* fields_group, const std::string& topo ) :
    m_association( association ),
    m_fields(),
    m_nextidx(),
    m_name2idx(),
    m_fields_group( fields_group ),
    m_topology( topo )
{
  SLIC_ERROR_IF(
      (m_association < 0) || (m_association >= NUM_FIELD_ASSOCIATIONS),
      "Invalid field association!" );

  SLIC_ERROR_IF( m_fields_group==AXOM_NULLPTR, "NULL sidre group!" );

  size_t numGroups = m_fields_group->getNumGroups( );
  for ( size_t i=0; i < numGroups; ++i )
  {
    sidre::Group* gp = m_fields_group->getGroup( i );
    SLIC_ERROR_IF( gp == AXOM_NULLPTR, "Encountered a NULL group" );

    SLIC_ERROR_IF( !gp->hasChildView( "topology" ),
        "field [" << gp->getName() << "] does not conform to blueprint!" <<
        " Missing 'topology' view" );
    SLIC_ERROR_IF( !gp->getView( "topology" )->isString(), 
        "topology view needs to hold a string." );
    if ( gp->getView( "topology" )->getString() != m_topology )
    {
      continue;
    }

    SLIC_ERROR_IF( !gp->hasChildView( "association" ),
        "field [" << gp->getName() << "] does not conform to blueprint!" <<
        " Missing 'association' view" );
    SLIC_ERROR_IF( !gp->getView( "association" )->isString(), 
        "association view needs to hold a string." );

    SLIC_ERROR_IF( !gp->hasChildView( "volume_dependent" ),
        "field [" << gp->getName() << "] does not conform to blueprint!" <<
        " Missing 'volume_dependent' view" );
    SLIC_ERROR_IF( !gp->getView( "volume_dependent" )->isString(), 
        "volume_dependent view needs to hold a string." );

    SLIC_ERROR_IF( !gp->hasChildView( "values" ),
        "field [" << gp->getName() << "] does not conform to blueprint!" <<
        " Missing 'values' view" );

    // NOTE: currently the blue-print supports
    const char* assoc    = gp->getView( "association" )->getString( );
    const bool isVertex  = (strcmp( assoc, "vertex" ) == 0);
    const bool isElement = (strcmp( assoc, "element" ) == 0);
    SLIC_ERROR_IF( (!isVertex && !isElement),
        "field [" << gp->getName() << "] has invalid association!" <<
        " => association= " << assoc );

    int centering = ( isVertex )? NODE_CENTERED : CELL_CENTERED;

    IndexType num_tuples = -1;
    if ( centering == m_association )
    {
      const std::string name = gp->getName();
      SLIC_ERROR_IF( hasField( name ), "Encountered a duplicate field!" );

      sidre::View* view = gp->getView( "values" );
      const int fldIdx = getNewFieldIndex();
      SLIC_ASSERT( ( fldIdx >= 0 ) &&
                  ( fldIdx < static_cast< int >( m_fields.size() ) ) );

      Field* field = internal::getFieldFromView( name, view );
      if ( num_tuples == -1 )
      {
        num_tuples = field->getNumTuples();
      }

      SLIC_ERROR_IF( field->getNumTuples() != num_tuples, "Inconsistent number of tuples" );

      m_fields[ fldIdx ] = field;
      m_name2idx[ name ] = fldIdx;
    } // END if centering
  } // END for all fields
}
#endif

//------------------------------------------------------------------------------
IndexType FieldData::getNumTuples() const
{
  const int numFields = getNumFields( );
  if ( numFields == 0 )
  {
    return 0;
  }

  const Field* first_field = m_fields[ m_name2idx.begin()->second ];
  SLIC_ASSERT( first_field != AXOM_NULLPTR );
  IndexType num_tuples = first_field->getNumTuples();
  
  for_all_fields( [num_tuples]( const Field* field )
  {
    SLIC_ERROR_IF( field->getNumTuples() != num_tuples, "Inconsistent number of tuples" );
  });

  return num_tuples;
}

//------------------------------------------------------------------------------
IndexType FieldData::getCapacity() const
{
  IndexType min_capacity = std::numeric_limits< IndexType >::max(); 
  for_all_fields( [&min_capacity]( const Field* field )
  {
    if ( field->getCapacity() < min_capacity )
    {
      min_capacity = field->getCapacity();
    }
  });

  return min_capacity;
}


//------------------------------------------------------------------------------
void FieldData::clear()
{
  for ( size_t i=0; i < m_fields.size(); ++i )
  {
     if ( m_fields[ i ] != AXOM_NULLPTR )
     {
       delete m_fields[ i ] ;
       m_fields[ i ] = AXOM_NULLPTR;
     }

  } // END for all fields

  m_fields.clear( );
  m_nextidx.clear( );
  m_name2idx.clear( );
}

//------------------------------------------------------------------------------
void FieldData::resize( IndexType newNumTuples )
{
  for_all_fields( [ newNumTuples ]( Field* field )
  {
    field->resize( newNumTuples );
  });
}

//------------------------------------------------------------------------------
void FieldData::reserveForInsert( IndexType pos, IndexType num_tuples )
{
  for_all_fields( [ pos, num_tuples ]( Field* field )
  {
    field->reserveForInsert( pos, num_tuples );
  });
}

//------------------------------------------------------------------------------
void FieldData::reserve( IndexType newCapacity )
{
  for_all_fields( [ newCapacity ]( Field* field )
  {
    field->reserve( newCapacity );
  });
}

//------------------------------------------------------------------------------
void FieldData::shrink()
{
  for_all_fields( []( Field* field )
  {
    field->shrink();
  });
}

//------------------------------------------------------------------------------
void FieldData::setResizeRatio( double ratio )
{
  for_all_fields( [ ratio ]( Field* field )
  {
    field->setResizeRatio( ratio );
  });
}


//------------------------------------------------------------------------------
void FieldData::removeField( const std::string& name )
{
  SLIC_ERROR_IF( name.empty(), "field name cannot be empty!" );
  SLIC_ERROR_IF( !hasField( name ),
                "cannot remove field [" << name << "]; field does not exist!");

  const int fldIdx = this->getFieldIndex( name );
  SLIC_ERROR_IF( fldIdx==INVALID_FIELD_INDEX,
                "cannot remove field [" << name << "]; invalid field index!");

  removeFieldAt( fldIdx );
  m_name2idx.erase( name );

#ifdef MINT_USE_SIDRE
  if ( hasSidreGroup() && m_fields_group->hasChildGroup(name) )
  {
    internal::removeFromSidre( m_fields_group, name );
    SLIC_ASSERT( ! m_fields_group->hasChildGroup( name ) );
  }
#endif
}

//------------------------------------------------------------------------------
void FieldData::removeField( int i )
{
  SLIC_ERROR_IF( (i < 0) && (i >= getNumFields() ),
                "specified field index is out-of-bounds!" );

  mint::Field* f = getField( i );
  SLIC_ASSERT( f != AXOM_NULLPTR );

  removeField( f->getName() );
  f = AXOM_NULLPTR;
}

//------------------------------------------------------------------------------
void FieldData::removeFieldAt( int i )
{
  SLIC_ASSERT( ( i >= 0 ) && ( i < static_cast< int >(m_fields.size() ) ) );
  SLIC_ASSERT( m_fields[ i ] != AXOM_NULLPTR );

  delete m_fields[ i ];
  m_fields[ i ] = AXOM_NULLPTR;
  m_nextidx.push_back( i );
}

} /* namespace mint */
} /* namespace axom */
