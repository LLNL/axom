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

  SLIC_ERROR_IF( f == AXOM_NULLPTR, "null field!" );
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
    m_resize_ratio( Array< double >::DEFAULT_RESIZE_RATIO ),
    m_fields()
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
    m_resize_ratio( Array< double >::DEFAULT_RESIZE_RATIO ),
    m_fields(),
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
      Field* field = internal::getFieldFromView( name, view );
      if ( num_tuples == -1 )
      {
        num_tuples = field->getNumTuples();
      }

      SLIC_ERROR_IF( field->getNumTuples() != num_tuples, "Inconsistent number of tuples" );

      m_fields[ name ] = field;
    } // END if centering
  } // END for all fields
}
#endif

//------------------------------------------------------------------------------
IndexType FieldData::getNumTuples() const
{
  const int numFields = getNumFields();
  if ( numFields == 0 )
  {
    return 0;
  }

  IndexType num_tuples = getField(0)->getNumTuples();
  bool status = true;

  for ( int i = 1; i < numFields; ++i )
  {
    const Field* f = getField( i );
    status &= f->getNumTuples() == num_tuples; 
  }

  SLIC_WARNING_IF( !status, "Inconsistent number of tuples." ); 
  return status ? num_tuples : -1;
}

//------------------------------------------------------------------------------
IndexType FieldData::getCapacity() const
{
  const int numFields = getNumFields();
  if ( numFields == 0 )
  {
    return 0;
  }

  IndexType min_capacity = std::numeric_limits< IndexType >::max(); 
  for ( int i = 0; i < numFields; ++i )
  {
    const Field* f = getField( i );
    const IndexType f_capacity = f->getCapacity();
    min_capacity = (f_capacity < min_capacity) ? f_capacity : min_capacity;
  }

  return min_capacity;
}

//------------------------------------------------------------------------------
void FieldData::clear()
{
  const int numFields = getNumFields();
  for ( int i = 0; i < numFields; ++i )
  {
    delete getField( i );
  }

  m_fields.clear( );
}

//------------------------------------------------------------------------------
void FieldData::resize( IndexType newNumTuples )
{
  const IndexType numFields = getNumFields();
  for ( int i = 0; i < numFields; ++i )
  {
    getField( i )->resize( newNumTuples );
  };
}

//------------------------------------------------------------------------------
void FieldData::emplace( IndexType pos, IndexType num_tuples )
{
  const IndexType numFields = getNumFields();
  for ( int i = 0; i < numFields; ++i )
  {
    getField( i )->emplace( pos, num_tuples );
  };
}

//------------------------------------------------------------------------------
void FieldData::reserve( IndexType newCapacity )
{
  const IndexType numFields = getNumFields();
  for ( int i = 0; i < numFields; ++i )
  {
    getField( i )->reserve( newCapacity );
  };
}

//------------------------------------------------------------------------------
void FieldData::shrink()
{
  const IndexType numFields = getNumFields();
  for ( int i = 0; i < numFields; ++i )
  {
    getField( i )->shrink();
  };
}

//------------------------------------------------------------------------------
void FieldData::setResizeRatio( double ratio )
{
  m_resize_ratio = ratio;
  const IndexType numFields = getNumFields();
  for ( int i = 0; i < numFields; ++i )
  {
    getField( i )->setResizeRatio( ratio );
  };
}

//------------------------------------------------------------------------------
void FieldData::removeField( const std::string& name )
{
  mint::Field* f = getField( name );
  SLIC_ASSERT( f != AXOM_NULLPTR );
  m_fields.erase( name );
  delete f;

#ifdef MINT_USE_SIDRE
  if ( hasSidreGroup() && m_fields_group->hasChildGroup(name) )
  {
    internal::removeFromSidre( m_fields_group, name );
    SLIC_ASSERT( !m_fields_group->hasChildGroup( name ) );
  }
#endif
}

//------------------------------------------------------------------------------
void FieldData::removeField( int i )
{
  mint::Field* f = getField( i );
  removeField( f->getName() );
  f = AXOM_NULLPTR;
}

} /* namespace mint */
} /* namespace axom */
