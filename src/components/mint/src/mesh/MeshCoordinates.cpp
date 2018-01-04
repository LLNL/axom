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

#include "mint/MeshCoordinates.hpp"     /* for MeshCoordinates */

#include "mint/config.hpp"           /* for IndexType */
#include "mint/Array.hpp"               /* for Array */
#include "slic/slic.hpp"                /* for slic macros */

#ifdef MINT_USE_SIDRE
#include "sidre/sidre.hpp"
#endif

#include <cstring>                      /* for strcmp */

namespace axom
{
namespace mint
{

//------------------------------------------------------------------------------
MeshCoordinates::MeshCoordinates( int dimension ) :
  m_ndims( dimension )
{
  // TODO: implement this
}


//------------------------------------------------------------------------------
MeshCoordinates::MeshCoordinates( int dimension, IndexType capacity,
                                  IndexType size, double resize_ratio ) :
  m_ndims( dimension )
{
  SLIC_ERROR_IF( m_ndims < 0 || m_ndims > 3, "invalid dimension" );

  for ( int dim = 0 ; dim < m_ndims ; ++dim )
  {
    Array< double >* coord_array;
    coord_array = new Array< double >( size, capacity );
    coord_array->setResizeRatio( resize_ratio );
    m_coordinates[ dim ] = coord_array;
  }

  for (int dim = m_ndims ; dim < 3 ; ++dim )
  {
    m_coordinates[ dim ] = AXOM_NULLPTR;
  }

}

#ifdef MINT_USE_SIDRE
//------------------------------------------------------------------------------
MeshCoordinates::MeshCoordinates( sidre::Group* group, int dimension,
                                  IndexType size, double resize_ratio ) :
  m_ndims(0)
{
  SLIC_ERROR_IF( group == AXOM_NULLPTR, "" );
  SLIC_ERROR_IF( !group->hasChildView( "type" ), "" );

  sidre::View* type_view = group->getView( "type" );
  SLIC_ERROR_IF( !type_view->isString(), "" );
  SLIC_ERROR_IF( std::strcmp( type_view->getString(), "explicit") != 0, "" );

  SLIC_ERROR_IF( !group->hasChildGroup("values"), "" );
  sidre::Group* values_group = group->getGroup( "values" );
  SLIC_ERROR_IF( !values_group->hasChildView( "x" ), "" );

  const char coord_names[3][2] = { "x", "y", "z" };
  for ( ; m_ndims < 3 ; ++m_ndims )
  {
    const char* coord_name = coord_names[ m_ndims ];

    if ( !values_group->hasChildView( coord_name ) )
    {
      break;
    }

    sidre::View* coord_view = values_group->getView( coord_name );
    SLIC_ERROR_IF( coord_view->getNumDimensions() != 2, "" );

    sidre::SidreLength dims[ 2 ];
    coord_view->getShape( 2, dims );
    SLIC_ERROR_IF( dims[1] != 1, "" );

    m_coordinates[ m_ndims ] =
      new Array< double >( coord_view, size );
    m_coordinates[ m_ndims ]->setResizeRatio( resize_ratio );
  }

  SLIC_ERROR_IF( m_ndims != dimension, "");

  for ( int i = m_ndims ; i < 3 ; ++i )
  {
    m_coordinates[ i ] = AXOM_NULLPTR;
  }
}

//------------------------------------------------------------------------------
MeshCoordinates::MeshCoordinates( sidre::Group* group, int dimension,
                                  IndexType capacity, IndexType size,
                                  double resize_ratio ) :
  m_ndims( dimension )
{
  SLIC_ERROR_IF( m_ndims < 0 || m_ndims > 3, "" );
  SLIC_ERROR_IF( group == AXOM_NULLPTR, "" );
  SLIC_ERROR_IF( group->getNumGroups() != 0, "" );
  SLIC_ERROR_IF( group->getNumViews() != 0, "" );

  group->createView( "type" )->setString( "explicit" );
  sidre::Group* values_group = group->createGroup( "values" );

  const char coord_names[3][2] = { "x", "y", "z" };

  for ( int dim = 0 ; dim < m_ndims ; ++dim )
  {
    const char* coord_name = coord_names[ dim ];
    sidre::View* coord_view = values_group->createView( coord_name );
    m_coordinates[ dim ] =
      new Array< double > ( coord_view, capacity, size, 1 );
    m_coordinates[ dim ]->setResizeRatio( resize_ratio );
  }

  for ( int i = m_ndims ; i < 3 ; ++i )
  {
    m_coordinates[ i ] = AXOM_NULLPTR;
  }
}

#endif

//------------------------------------------------------------------------------
MeshCoordinates::~MeshCoordinates()
{
  for ( int dim = 0 ; dim < 3 ; ++dim )
  {
    if ( m_coordinates[ dim ] != AXOM_NULLPTR )
    {
      delete m_coordinates[ dim ];
    }

    m_coordinates[ dim ] = AXOM_NULLPTR;
  }
}

}   /* end namespace mint */
}   /* end namespace axom */
