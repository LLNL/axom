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
#include "mint/RectilinearMesh.hpp"

#include "mint/Array.hpp"          // for mint::Array
#include "mint/blueprint.hpp"      // for blueprint functions
#include "mint/config.hpp"         // for compile-time definitions
#include "mint/MeshTypes.hpp"      // for STRUCTURED_RECTILINEAR_MESH

namespace axom
{
namespace mint
{

//------------------------------------------------------------------------------
// HELPER METHODS
//------------------------------------------------------------------------------
namespace
{

inline int dim( const double* AXOM_NOT_USED(x),
                const double* y,
                const double* z )
{
  return ( ( z != AXOM_NULLPTR ) ? 3 : ( (y != AXOM_NULLPTR ) ? 2 : 1 ) );
}

//------------------------------------------------------------------------------
inline int dim ( const IndexType& AXOM_NOT_USED( Ni ),
                 const IndexType& Nj,
                 const IndexType& Nk  )
{
  return ( (Nk >= 1)? 3 : ( (Nj >= 1) ? 2 : 1 ) );
}

} /* anonymous namespace */

//------------------------------------------------------------------------------
// RECTILINEAR MESH IMPLEMENTATION
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
RectilinearMesh::RectilinearMesh( int dimension, const int64* ext ) :
  StructuredMesh( STRUCTURED_RECTILINEAR_MESH, dimension, ext )
{
  initialize();
  allocateCoords();
}

//------------------------------------------------------------------------------
RectilinearMesh::RectilinearMesh( IndexType Ni, IndexType Nj, IndexType Nk ) :
  StructuredMesh( STRUCTURED_RECTILINEAR_MESH, dim(Ni,Nj,Nk) )
{
  int64 extent[ ] = { 0, Ni-1, 0, Nj-1, 0, Nk-1 };
  m_extent        = new mint::Extent( m_ndims, extent );

  initialize( );
  allocateCoords();
}

//------------------------------------------------------------------------------
RectilinearMesh::RectilinearMesh( const int64* ext,
                                  double* x,
                                  double* y,
                                  double* z ) :
  StructuredMesh( STRUCTURED_RECTILINEAR_MESH, dim(x,y,z), ext )
{
  initialize( );

  double* ptrs[3];
  ptrs[ 0 ] = x;
  ptrs[ 1 ] = y;
  ptrs[ 2 ] = z;

  for ( int i=0; i < m_ndims; ++i )
  {
    SLIC_ERROR_IF( ptrs[ i ]==AXOM_NULLPTR,
                  "encountered null coordinate array for i=" << i );
    const IndexType N  = m_extent->size( i );
    m_coordinates[ i ] = new Array< double >( ptrs[i], N, 1, N );
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

  int64 extent[ 6 ];
  blueprint::getRectilinearMeshExtent(
      m_ndims, getCoordsetGroup(), getTopologyGroup(), extent );
  m_extent = new mint::Extent( m_ndims, extent );

  initialize( );

  sidre::Group* c = getCoordsetGroup();
  SLIC_ERROR_IF( !blueprint::validCoordsetGroup( c ), "invalid coordset!" );

  const char* coords[] = { "values/x", "values/y", "values/z" };

  // initialize coordinates
  for ( int i=0; i < m_ndims; ++i )
  {
    m_coordinates[ i ] = new Array< double >( c->getView( coords[ i ] ) );
    SLIC_ERROR_IF( m_extent->size( i ) != m_coordinates[ i ]->size(),
        "coordinates size does not match rectilinear mesh extent" );
  }

}

//------------------------------------------------------------------------------
RectilinearMesh::RectilinearMesh( int dimension,
                                  const int64* ext,
                                  sidre::Group* group,
                                  const std::string& topo,
                                  const std::string& coordset ) :
  StructuredMesh( STRUCTURED_RECTILINEAR_MESH, dimension,group,topo,coordset )
{
  blueprint::initializeTopologyGroup( m_group, m_topology, m_coordset,
                                     "rectilinear" );
  SLIC_ERROR_IF( !blueprint::validTopologyGroup( getTopologyGroup() ),
                "invalid topology group!" );

  m_extent = new mint::Extent( m_ndims, ext );
  blueprint::setRectilinearMeshExtent( m_ndims, m_extent, getTopologyGroup() );

  initialize( );

  allocateCoordsOnSidre( );
}

//------------------------------------------------------------------------------
RectilinearMesh::RectilinearMesh( sidre::Group* group,
                                  const std::string& topo,
                                  const std::string& coordset,
                                  IndexType Ni,
                                  IndexType Nj,
                                  IndexType Nk ) :
  StructuredMesh( STRUCTURED_RECTILINEAR_MESH,dim(Ni,Nj,Nk),group,topo,coordset)
{
  blueprint::initializeTopologyGroup( m_group, m_topology, m_coordset,
                                      "rectilinear" );
  SLIC_ERROR_IF( !blueprint::validTopologyGroup( getTopologyGroup() ),
                 "invalid topology group!" );

  int64 extent[ ] = { 0, Ni-1, 0, Nj-1, 0, Nk-1 };
  m_extent        = new mint::Extent( m_ndims, extent );
  blueprint::setRectilinearMeshExtent( m_ndims, m_extent, getTopologyGroup() );

  initialize( );

  allocateCoordsOnSidre( );
}

//------------------------------------------------------------------------------
RectilinearMesh::RectilinearMesh( sidre::Group* group,
                                  IndexType Ni,
                                  IndexType Nj,
                                  IndexType Nk ) :
  RectilinearMesh( group, "", "", Ni, Nj, Nk )
{ }

//------------------------------------------------------------------------------
void RectilinearMesh::allocateCoordsOnSidre()
{
  SLIC_ASSERT( m_extent != AXOM_NULLPTR );

  sidre::Group* coordsgrp = getCoordsetGroup();
  SLIC_ERROR_IF( coordsgrp==AXOM_NULLPTR, "coordset group is null!" );
  SLIC_ERROR_IF( coordsgrp->getNumViews() != 0,  "coordset group is not empty");
  SLIC_ERROR_IF( coordsgrp->getNumGroups() != 0, "coordset group is not empty");

  coordsgrp->createView( "type" )->setString( "rectilinear" );

  const char* coords[] = { "values/x", "values/y", "values/z" };

  for ( int idim=0; idim < m_ndims; ++idim )
  {
    IndexType N           = m_extent->size( idim );
    sidre::View* view     = coordsgrp->createView( coords[ idim ] );
    m_coordinates[ idim ] = new Array< double >( view, N, 1, N );
    m_coordinates[ idim ]->setResizeRatio( 0.0 );
  }

  SLIC_ERROR_IF( !blueprint::validCoordsetGroup( getCoordsetGroup() ),
                 "invalid coordset group!" );
}

#endif /* MINT_USE_SIDRE */

//------------------------------------------------------------------------------
RectilinearMesh::~RectilinearMesh()
{
  for ( int dim = 0 ; dim < 3 ; ++dim )
  {
    if ( m_coordinates[ dim ] != AXOM_NULLPTR )
    {
      delete m_coordinates[ dim ];
      m_coordinates[ dim ] = AXOM_NULLPTR;
    }
  }
}

//------------------------------------------------------------------------------
void RectilinearMesh::initialize()
{
  SLIC_ASSERT( m_extent != AXOM_NULLPTR );

  m_explicit_coords       = true;
  m_explicit_connectivity = false;
  m_has_mixed_topology    = false;

  initializeFields();
}

//------------------------------------------------------------------------------
void RectilinearMesh::allocateCoords()
{
  SLIC_ASSERT( m_extent != AXOM_NULLPTR );
  SLIC_ASSERT( (m_ndims >= 1) && (m_ndims <= 3) );
  SLIC_ASSERT( m_extent->getDimension() == m_ndims );

  for ( int idim=0; idim < m_ndims; ++idim )
  {
    const IndexType N     = m_extent->size( idim );
    m_coordinates[ idim ] = new Array< double >( N, 1, N );
    m_coordinates[ idim ]->setResizeRatio( 0.0 );
  } // END for all dimensions

}

//------------------------------------------------------------------------------
double* RectilinearMesh::getCoordinateArray( int dim )
{
  SLIC_ASSERT( 0 <= dim && dim < getDimension() );
  SLIC_ASSERT( m_coordinates[ dim ] != AXOM_NULLPTR );
  return m_coordinates[ dim ]->getData();
}

//------------------------------------------------------------------------------
const double* RectilinearMesh::getCoordinateArray( int dim ) const
{
  SLIC_ASSERT( 0 <= dim && dim < getDimension() );
  SLIC_ASSERT( m_coordinates[ dim ] != AXOM_NULLPTR );
  return m_coordinates[ dim ]->getData();
}

//------------------------------------------------------------------------------
void RectilinearMesh::getNode( IndexType nodeID, double* node ) const
{
  SLIC_ASSERT( 0 <= nodeID && nodeID < getNumberOfNodes() );
  SLIC_ASSERT( node != AXOM_NULLPTR );

  IndexType i = -1;
  IndexType j = -1;
  IndexType k = -1;

  switch ( m_ndims )
  {
  case 1:
    node[ 0 ] = getCoordinateArray( X_COORDINATE )[ nodeID ];
    break;
  case 2:
    m_extent->getGridIndex( nodeID, i, j );
    node[ 0 ] = getCoordinateArray( X_COORDINATE )[ i ];
    node[ 1 ] = getCoordinateArray( Y_COORDINATE )[ j ];
    break;
  default:
    SLIC_ASSERT( m_ndims==3 );
    m_extent->getGridIndex( nodeID, i, j, k );
    node[ 0 ] = getCoordinateArray( X_COORDINATE )[ i ];
    node[ 1 ] = getCoordinateArray( Y_COORDINATE )[ j ];
    node[ 2 ] = getCoordinateArray( Z_COORDINATE )[ k ];
  } // END switch

}

} /* namespace mint */
} /* namespace axom */
