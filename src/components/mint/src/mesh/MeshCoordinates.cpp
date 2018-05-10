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

#include "mint/MeshCoordinates.hpp"  // for mint::MeshCoordinates

// Axom includes
#include "axom_utils/Utilities.hpp"  // for utilities::max()
#include "mint/config.hpp"           // for IndexType
#include "mint/Array.hpp"            // for mint::Array
#include "slic/slic.hpp"             // for slic macros

#ifdef MINT_USE_SIDRE
#include "sidre/sidre.hpp"          // for sidre::Group, sidre::View
#endif

// C/C++ includes
#include <cstring>                   // for std::strcmp()

namespace axom
{
namespace mint
{

constexpr IndexType DEFAULT_CAPACITY = 100;

//------------------------------------------------------------------------------
MeshCoordinates::MeshCoordinates( int dimension,
                                  IndexType numNodes,
                                  IndexType capacity ) :
#ifdef MINT_USE_SIDRE
  m_group( AXOM_NULLPTR ),
#endif
  m_ndims( dimension )
{
  SLIC_ERROR_IF( this->invalidDimension(), "invalid dimension" );

  IndexType max_capacity = -1;
  if ( capacity==USE_DEFAULT )
  {
    const double ratio = mint::Array< double >::DEFAULT_RESIZE_RATIO;
    max_capacity = utilities::max(
             DEFAULT_CAPACITY, static_cast< IndexType >( numNodes*ratio+0.5 ) );
  }
  else
  {
    max_capacity = capacity;
  }

  SLIC_ERROR_IF( numNodes > max_capacity, "numNodes > capacity!" );
  this->initialize( numNodes, max_capacity );
}

//------------------------------------------------------------------------------
MeshCoordinates::MeshCoordinates( IndexType numNodes, IndexType capacity, 
                                  double* x, double* y, double* z ) :
#ifdef MINT_USE_SIDRE
  m_group( AXOM_NULLPTR ),
#endif
  m_ndims( 0 )
{

  m_ndims = ( z != AXOM_NULLPTR ) ? 3 : ( (y != AXOM_NULLPTR ) ? 2 : 1 ) ;
  SLIC_ERROR_IF( invalidDimension(), "invalid dimension" );
  SLIC_ERROR_IF( capacity < 1, "capacity < 1" );

  double* ptrs[3];
  ptrs[ 0 ] = x;
  ptrs[ 1 ] = y;
  ptrs[ 2 ] = z;

  for ( int i=0; i < m_ndims; ++i )
  {
    SLIC_ERROR_IF( ptrs[ i ]==AXOM_NULLPTR,
                  "encountered null coordinate array for i=" << i );

    m_coordinates[ i ] = new Array< double >( ptrs[i], numNodes, 1, capacity );
  }

  SLIC_ASSERT( consistencyCheck() );
}

#ifdef MINT_USE_SIDRE
//------------------------------------------------------------------------------
MeshCoordinates::MeshCoordinates( sidre::Group* group ) :
  m_group( group ),
  m_ndims( 0 )
{
  SLIC_ERROR_IF( m_group == AXOM_NULLPTR, "null sidre::Group" );
  SLIC_ERROR_IF( !m_group->hasChildView( "type" ),
                 "sidre::Group does not conform to mesh blueprint" );

  sidre::View* type_view = m_group->getView( "type" );
  SLIC_ERROR_IF( !type_view->isString(),
                 "sidre::Group does not conform to mesh blueprint" );

  SLIC_ERROR_IF( std::strcmp( type_view->getString(), "explicit") != 0,
                "sidre::Group does not conform to mesh blueprint" );

  SLIC_ERROR_IF( !m_group->hasChildGroup("values"),
                 "sidre::Group does not conform to mesh blueprint" );

  // NOTE: here we should support cylindrical and spherical coordinates
  sidre::Group* values_group = m_group->getGroup( "values" );
  SLIC_ERROR_IF( !values_group->hasChildView( "x" ),
                 "sidre::Group does not conform to mesh blueprint" );

  const bool hasZ = values_group->hasChildView( "z" );
  const bool hasY = values_group->hasChildView( "y" );

  m_ndims = ( hasZ ) ? 3 : ( ( hasY ) ? 2 : 1 );
  SLIC_ERROR_IF( this->invalidDimension(), "invalid dimension" );

  const char* coord_names[3] = { "x", "y", "z" };
  for ( int i=0 ; i < m_ndims ; ++i )
  {
    const char* coord_name = coord_names[ i ];
    SLIC_ASSERT( values_group->hasView( std::string( coord_name ) ) );

    sidre::View* coord_view = values_group->getView( coord_name );
    SLIC_ASSERT( coord_view != AXOM_NULLPTR );
    SLIC_ERROR_IF( coord_view->getNumDimensions() != 2,
                   "view has invalid dimensions" );

    sidre::SidreLength dims[ 2 ];
    coord_view->getShape( 2, dims );
    SLIC_ERROR_IF( dims[1] != 1, "number of components is expected to be 1" );

    m_coordinates[ i ] = new Array< double >( coord_view );

  } // END for all dimensions

}

//------------------------------------------------------------------------------
MeshCoordinates::MeshCoordinates( sidre::Group* group, int dimension,
                                  IndexType numNodes, IndexType capacity ) :
  m_group( group ),
  m_ndims( dimension )
{
  SLIC_ERROR_IF( m_group==AXOM_NULLPTR, "null sidre::Group" );
  SLIC_ERROR_IF( m_group->getNumGroups() != 0, "sidre::Group is not empty!" );
  SLIC_ERROR_IF( m_group->getNumViews() !=0, "sidre::Group is not empty!" );
  SLIC_ERROR_IF( (capacity != USE_DEFAULT) && (numNodes > capacity),
                "numNodes < capacity pre-condition violated!" );

  m_group->createView( "type" )->setString( "explicit" );

  sidre::Group* values = m_group->createGroup( "values" );
  SLIC_ASSERT( values != AXOM_NULLPTR );

  const char* coord_names[3] = { "x", "y", "z" };

  for ( int dim = 0 ; dim < m_ndims ; ++dim )
  {
    const char* coord_name = coord_names[ dim ];
    sidre::View* coord_view = values->createView( coord_name );
    m_coordinates[ dim ] =
        new Array< double > ( coord_view, numNodes, 1, capacity );
  }

}

#endif

//------------------------------------------------------------------------------
MeshCoordinates::~MeshCoordinates()
{

  for ( int dim = 0 ; dim < m_ndims ; ++dim )
  {
    SLIC_ASSERT( m_coordinates[ dim ] != AXOM_NULLPTR );

    delete m_coordinates[ dim ];
    m_coordinates[ dim ] = AXOM_NULLPTR;
  }

}

//------------------------------------------------------------------------------
bool MeshCoordinates::consistencyCheck() const
{
  bool status = true;

  SLIC_ASSERT( !invalidDimension( ) );
  SLIC_ASSERT( m_coordinates[ 0 ] != AXOM_NULLPTR );

  const IndexType NUM_COMPONENTS     = 1;
  const IndexType expected_size      = m_coordinates[ 0 ]->size();
  const IndexType expected_capacity  = m_coordinates[ 0 ]->capacity();
  const double expected_resize_ratio = m_coordinates[ 0 ]->getResizeRatio();
  const bool expected_is_external    = m_coordinates[ 0 ]->isExternal();

  for ( int i = 1; i < m_ndims; ++i )
  {
    const IndexType actual_size       = m_coordinates[ i ]->size();
    const IndexType actual_components = m_coordinates[ i ]->numComponents();
    const IndexType actual_capacity   = m_coordinates[ i ]->capacity();
    const double actual_resize_ratio  = m_coordinates[ i ]->getResizeRatio();

    const bool size_mismatch      = ( actual_size != expected_size );
    const bool component_mismatch = ( actual_components != NUM_COMPONENTS );
    const bool capacity_mismatch  = ( actual_capacity != expected_capacity );
    const bool ratio_mismatch     =
        !utilities::isNearlyEqual( actual_resize_ratio, expected_resize_ratio );

    SLIC_WARNING_IF( size_mismatch, "coordinate array size mismatch!" );
    SLIC_WARNING_IF( component_mismatch, 
                                 "coordinate array number of components != 1" );
    SLIC_WARNING_IF( capacity_mismatch, "coordinate array capacity mismatch!" );
    SLIC_WARNING_IF( ratio_mismatch, "coordinate array ratio mismatch!" );

    if ( size_mismatch || capacity_mismatch || ratio_mismatch )
    {
      status = false;
      break;
    }

    if ( expected_is_external != m_coordinates[ i ]->isExternal( ) )
    {
      SLIC_WARNING( "external propery mismatch!" );
      status = false;
      break;
    }

  } // END for all dimensions

  return status;
}

}   /* end namespace mint */
}   /* end namespace axom */
