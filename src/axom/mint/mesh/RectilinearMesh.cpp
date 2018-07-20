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
#include "axom/mint/mesh/RectilinearMesh.hpp"

#include "axom/mint/core/Array.hpp"          // for mint::Array
#include "axom/mint/mesh/blueprint.hpp"      // for blueprint functions
#include "axom/mint/config.hpp"         // for compile-time definitions
#include "axom/mint/mesh/MeshTypes.hpp"      // for STRUCTURED_RECTILINEAR_MESH

#include "axom/mint/mesh/internal/MeshHelpers.hpp"    // for internal helper

namespace axom
{
namespace mint
{


//------------------------------------------------------------------------------
// RECTILINEAR MESH IMPLEMENTATION
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
RectilinearMesh::RectilinearMesh( int dimension, const IndexType* ext ) :
  StructuredMesh( STRUCTURED_RECTILINEAR_MESH, dimension, ext )
{
  initialize();
  allocateCoords();
}

//------------------------------------------------------------------------------
RectilinearMesh::RectilinearMesh( IndexType Ni, IndexType Nj, IndexType Nk ) :
  StructuredMesh( STRUCTURED_RECTILINEAR_MESH, Ni, Nj, Nk )
{
  initialize();
  allocateCoords();
}

//------------------------------------------------------------------------------
RectilinearMesh::RectilinearMesh( const IndexType* ext,
                                  double* x,
                                  double* y,
                                  double* z ) :
  StructuredMesh( STRUCTURED_RECTILINEAR_MESH, internal::dim(x, y, z), ext )
{
  initialize( );

  double* ptrs[3];
  ptrs[ 0 ] = x;
  ptrs[ 1 ] = y;
  ptrs[ 2 ] = z;

  for ( int dim = 0; dim < m_ndims; ++dim )
  {
    SLIC_ERROR_IF( ptrs[ dim ] == AXOM_NULLPTR,
                   "encountered null coordinate array for dim=" << dim );
    const IndexType N = getNodeExtent( dim );
    m_coordinates[ dim ] = new Array< double >( ptrs[ dim ], N, 1, N );
  }
}

#ifdef MINT_USE_SIDRE

//------------------------------------------------------------------------------
RectilinearMesh::RectilinearMesh( sidre::Group* group, 
                                  const std::string& topo ) :
  StructuredMesh( group, topo )
{
  SLIC_ERROR_IF( m_type != STRUCTURED_RECTILINEAR_MESH,
                 "supplied Sidre group does not correspond to a RectilinearMesh" );

  initialize();

  sidre::Group* c = getCoordsetGroup();
  SLIC_ERROR_IF( !blueprint::isValidCoordsetGroup( c ), "invalid coordset!" );

  const char* coords[] = { "values/x", "values/y", "values/z" };

  // initialize coordinates
  for ( int dim = 0; dim < m_ndims; ++dim )
  {
    m_coordinates[ dim ] = new Array< double >( c->getView( coords[ dim ] ) );
    SLIC_ERROR_IF( getNodeExtent( dim ) != m_coordinates[ dim ]->size(),
                   "coordinates size does not match rectilinear mesh extent" );
  }

}

//------------------------------------------------------------------------------
RectilinearMesh::RectilinearMesh( int dimension, const IndexType* ext, 
                                  sidre::Group* group, const std::string& topo,
                                  const std::string& coordset ) :
  StructuredMesh( STRUCTURED_RECTILINEAR_MESH, dimension, ext, group, topo, 
                  coordset )
{
  initialize();
  allocateCoordsOnSidre( );
}

//------------------------------------------------------------------------------
RectilinearMesh::RectilinearMesh( sidre::Group* group, const std::string& topo,
                                  const std::string& coordset, IndexType Ni,
                                  IndexType Nj, IndexType Nk ) :
  StructuredMesh( STRUCTURED_RECTILINEAR_MESH, Ni, Nj, Nk, group, topo, 
                  coordset )
{
  initialize();
  allocateCoordsOnSidre( );
}

//------------------------------------------------------------------------------
void RectilinearMesh::allocateCoordsOnSidre()
{
<<<<<<< 34872490e5ae665038fbc51bba86b61e5fd77246:src/axom/mint/mesh/RectilinearMesh.cpp
  SLIC_ASSERT( m_extent != nullptr );

  sidre::Group* coordsgrp = getCoordsetGroup();
  SLIC_ERROR_IF( coordsgrp==nullptr, "coordset group is null!" );
  SLIC_ERROR_IF( coordsgrp->getNumViews() != 0,  "coordset group is not empty");
  SLIC_ERROR_IF( coordsgrp->getNumGroups() != 0, "coordset group is not empty");
=======
  sidre::Group* coordsgrp = getCoordsetGroup();
  SLIC_ERROR_IF( coordsgrp == AXOM_NULLPTR, "coordset group is null!" );
>>>>>>> ENH: Removed Extent class.:src/components/mint/src/mesh/RectilinearMesh.cpp

  coordsgrp->createView( "type" )->setString( "rectilinear" );

  const char* coords[] = { "values/x", "values/y", "values/z" };

  for ( int dim=0 ; dim < m_ndims ; ++dim )
  {
    IndexType N          = getNodeExtent( dim );
    sidre::View* view    = coordsgrp->createView( coords[ dim ] );
    m_coordinates[ dim ] = new Array< double >( view, N, 1, N );
    m_coordinates[ dim ]->setResizeRatio( 0.0 );
  }

  SLIC_ERROR_IF( !blueprint::isValidCoordsetGroup( getCoordsetGroup() ),
                 "invalid coordset group!" );
}

#endif /* MINT_USE_SIDRE */

//------------------------------------------------------------------------------
RectilinearMesh::~RectilinearMesh()
{
  for ( int dim = 0; dim < 3; ++dim )
  {
    if ( m_coordinates[ dim ] != nullptr )
    {
      delete m_coordinates[ dim ];
      m_coordinates[ dim ] = nullptr;
    }
  }
}

//------------------------------------------------------------------------------
void RectilinearMesh::initialize()
{
<<<<<<< 34872490e5ae665038fbc51bba86b61e5fd77246:src/axom/mint/mesh/RectilinearMesh.cpp
  SLIC_ASSERT( m_extent != nullptr );

=======
>>>>>>> ENH: Removed Extent class.:src/components/mint/src/mesh/RectilinearMesh.cpp
  m_explicit_coords       = true;
  m_explicit_connectivity = false;
  m_has_mixed_topology    = false;
}

//------------------------------------------------------------------------------
void RectilinearMesh::allocateCoords()
{
<<<<<<< 34872490e5ae665038fbc51bba86b61e5fd77246:src/axom/mint/mesh/RectilinearMesh.cpp
  SLIC_ASSERT( m_extent != nullptr );
=======
>>>>>>> ENH: Removed Extent class.:src/components/mint/src/mesh/RectilinearMesh.cpp
  SLIC_ASSERT( (m_ndims >= 1) && (m_ndims <= 3) );

  for ( int dim = 0; dim < m_ndims; ++dim )
  {
    const IndexType N     = getNodeExtent( dim );
    m_coordinates[ dim ] = new Array< double >( N, 1, N );
    m_coordinates[ dim ]->setResizeRatio( 0.0 );
  } // END for all dimensions

}

//------------------------------------------------------------------------------
double* RectilinearMesh::getCoordinateArray( int dim )
{
  SLIC_ASSERT( 0 <= dim && dim < getDimension() );
  SLIC_ASSERT( m_coordinates[ dim ] != nullptr );
  return m_coordinates[ dim ]->getData();
}

//------------------------------------------------------------------------------
const double* RectilinearMesh::getCoordinateArray( int dim ) const
{
  SLIC_ASSERT( 0 <= dim && dim < getDimension() );
  SLIC_ASSERT( m_coordinates[ dim ] != nullptr );
  return m_coordinates[ dim ]->getData();
}

//------------------------------------------------------------------------------
void RectilinearMesh::getNode( IndexType nodeID, double* node ) const
{
  SLIC_ASSERT( 0 <= nodeID && nodeID < getNumberOfNodes() );
  SLIC_ASSERT( node != nullptr );

  IndexType i = -1;
  IndexType j = -1;
  IndexType k = -1;

  switch ( m_ndims )
  {
  case 1:
    node[ 0 ] = getCoordinateArray( X_COORDINATE )[ nodeID ];
    break;
  case 2:
    getNodeGridIndex( nodeID, i, j );
    node[ 0 ] = getCoordinateArray( X_COORDINATE )[ i ];
    node[ 1 ] = getCoordinateArray( Y_COORDINATE )[ j ];
    break;
  default:
    SLIC_ASSERT( m_ndims==3 );
    getNodeGridIndex( nodeID, i, j, k );
    node[ 0 ] = getCoordinateArray( X_COORDINATE )[ i ];
    node[ 1 ] = getCoordinateArray( Y_COORDINATE )[ j ];
    node[ 2 ] = getCoordinateArray( Z_COORDINATE )[ k ];
  } // END switch

}

} /* namespace mint */
} /* namespace axom */
