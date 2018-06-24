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
#include "mint/CurvilinearMesh.hpp"

// mint includes
#include "mint/blueprint.hpp"        // for blueprint functions
#include "mint/MeshTypes.hpp"        // for STRUCTURED_CURVILINEAR_MESH
#include "mint/MeshCoordinates.hpp"  // for MeshCoordinates class definition

#include "mint/MeshHelpers.hpp"      // for internal helper methods

// slic includes
#include "slic/slic.hpp"             // for SLIC macros

namespace axom
{
namespace mint
{

//------------------------------------------------------------------------------
// CURVILINEAR MESH IMPLEMENTATION
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
CurvilinearMesh::CurvilinearMesh( int ndims, const int64* ext ) :
  StructuredMesh( STRUCTURED_CURVILINEAR_MESH, ndims, ext ),
  m_coordinates( new mint::MeshCoordinates( ndims, m_extent->getNumNodes() ) )
{
  initialize( );

  // sanity checks
  SLIC_ASSERT( m_coordinates != AXOM_NULLPTR );
  SLIC_ASSERT( m_extent != AXOM_NULLPTR );
  SLIC_ASSERT( m_extent->getNumNodes()==m_coordinates->numNodes() );
}

//------------------------------------------------------------------------------
CurvilinearMesh::CurvilinearMesh( IndexType Ni,
                                  IndexType Nj,
                                  IndexType Nk ) :
  StructuredMesh( STRUCTURED_CURVILINEAR_MESH, internal::dim( Ni,Nj,Nk ) )
{
  int64 extent[]  = { 0,Ni-1, 0, Nj-1, 0,Nk-1 };
  m_extent        = new mint::Extent( m_ndims, extent );
  m_coordinates   = new mint::MeshCoordinates( m_ndims,getNumberOfNodes() );

  initialize( );

  // sanity checks
  SLIC_ASSERT( m_coordinates != AXOM_NULLPTR );
  SLIC_ASSERT( m_extent != AXOM_NULLPTR );
  SLIC_ASSERT( m_extent->getNumNodes()==m_coordinates->numNodes() );
}

//------------------------------------------------------------------------------
CurvilinearMesh::CurvilinearMesh( const int64* ext,
                                  double* x,
                                  double* y,
                                  double* z ) :
  StructuredMesh( STRUCTURED_CURVILINEAR_MESH, internal::dim(x,y,z), ext ),
  m_coordinates( new mint::MeshCoordinates( m_extent->getNumNodes(), x,y,z ) )
{
  initialize( );

  // sanity checks
  SLIC_ASSERT( m_coordinates != AXOM_NULLPTR );
  SLIC_ASSERT( m_extent != AXOM_NULLPTR );
  SLIC_ASSERT( m_extent->getNumNodes()==m_coordinates->numNodes() );
}

#ifdef MINT_USE_SIDRE

//------------------------------------------------------------------------------
CurvilinearMesh::CurvilinearMesh( sidre::Group* group,
                                  const std::string& topo ) :
  StructuredMesh( group, topo ),
  m_coordinates( new MeshCoordinates( getCoordsetGroup() ) )
{
  SLIC_ERROR_IF( m_type != STRUCTURED_CURVILINEAR_MESH,
                 "supplied Sidre group does not correspond to a CurvilinearMesh" );

  int64 extent[ 6 ];
  blueprint::getCurvilinearMeshExtent( m_ndims, getTopologyGroup(), extent );

  m_extent = new mint::Extent( m_ndims, extent );

  initialize( );

  // sanity checks
  SLIC_ASSERT( m_coordinates != AXOM_NULLPTR );
  SLIC_ASSERT( m_extent != AXOM_NULLPTR );
  SLIC_ASSERT( m_extent->getNumNodes()==m_coordinates->numNodes() );
  SLIC_ASSERT( m_coordinates->dimension()==m_ndims );
}

//------------------------------------------------------------------------------
CurvilinearMesh::CurvilinearMesh( int dimension,
                                  const int64* ext,
                                  sidre::Group* group,
                                  const std::string& topo,
                                  const std::string& coordset ) :
  StructuredMesh( STRUCTURED_CURVILINEAR_MESH, dimension,group,topo,coordset)
{

  blueprint::initializeTopologyGroup( m_group, m_topology, m_coordset,
                                      "structured" );
  SLIC_ERROR_IF( !blueprint::isValidTopologyGroup( getTopologyGroup() ),
                 "invalid topology group!" );

  m_extent        = new mint::Extent( m_ndims, ext );
  m_coordinates   = new mint::MeshCoordinates( getCoordsetGroup(),
                                               m_ndims,
                                               getNumberOfNodes(),
                                               getNumberOfNodes() );

  blueprint::setCurvilinearMeshExtent( m_ndims, m_extent, getTopologyGroup() );

  initialize( );

  // sanity checks
  SLIC_ASSERT( m_coordinates != AXOM_NULLPTR );
  SLIC_ASSERT( m_extent != AXOM_NULLPTR );
  SLIC_ASSERT( m_extent->getNumNodes()==m_coordinates->numNodes() );
  SLIC_ASSERT( m_coordinates->dimension()==m_ndims );
}

//------------------------------------------------------------------------------
CurvilinearMesh::CurvilinearMesh( sidre::Group* group,
                                  const std::string& topo,
                                  const std::string& coordset,
                                  IndexType Ni,
                                  IndexType Nj,
                                  IndexType Nk  ) :
  StructuredMesh( STRUCTURED_CURVILINEAR_MESH, internal::dim(Ni,Nj,Nk),
                  group,topo,coordset)
{
  blueprint::initializeTopologyGroup( m_group, m_topology, m_coordset,
                                      "structured" );
  SLIC_ERROR_IF( !blueprint::isValidTopologyGroup( getTopologyGroup() ),
                 "invalid topology group!" );

  int64 extent[]  = { 0,Ni-1, 0, Nj-1, 0,Nk-1 };
  m_extent        = new mint::Extent( m_ndims, extent );
  m_coordinates   = new mint::MeshCoordinates( getCoordsetGroup(),
                                               m_ndims,
                                               getNumberOfNodes(),
                                               getNumberOfNodes() );

  blueprint::setCurvilinearMeshExtent( m_ndims, m_extent, getTopologyGroup() );

  initialize( );

  // sanity checks
  SLIC_ASSERT( m_coordinates != AXOM_NULLPTR );
  SLIC_ASSERT( m_extent != AXOM_NULLPTR );
  SLIC_ASSERT( m_extent->getNumNodes()==m_coordinates->numNodes() );
  SLIC_ASSERT( m_coordinates->dimension()==m_ndims );
}

//------------------------------------------------------------------------------
CurvilinearMesh::CurvilinearMesh( sidre::Group* group,
                                  IndexType Ni,
                                  IndexType Nj,
                                  IndexType Nk  ) :
  CurvilinearMesh( group, "", "", Ni, Nj, Nk )
{ }

#endif

//------------------------------------------------------------------------------
CurvilinearMesh::~CurvilinearMesh()
{
  delete m_coordinates;
  m_coordinates = AXOM_NULLPTR;
}

//------------------------------------------------------------------------------
void CurvilinearMesh::initialize()
{
  m_explicit_coords       = true;
  m_explicit_connectivity = false;
  m_has_mixed_topology    = false;

  initializeFields( );
}

} /* namespace mint */
} /* namespace axom */
