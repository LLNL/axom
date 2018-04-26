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
#include "mint/UniformMesh.hpp"

// mint includes
#include "mint/blueprint.hpp"       // for blueprint functions
#include "mint/config.hpp"          // for compile-time definitions
#include "mint/MeshTypes.hpp"       // for STRUCTURED_UNIFORM_MESH

// axom_utils includes
#include "axom_utils/Utilities.hpp" // for utilities::isNearlyEqual

namespace axom
{
namespace mint
{

//------------------------------------------------------------------------------
// HELPER METHODS
//------------------------------------------------------------------------------
namespace
{

void set_spacing_and_origin( int ndims,
                             const mint::Extent* extent,
                             const double* lo,
                             const double* hi,
                             double* h,
                             double* origin )
{
  SLIC_ASSERT( ndims >= 1 && ndims <= 3 );
  SLIC_ASSERT( extent != AXOM_NULLPTR );
  SLIC_ASSERT( lo != AXOM_NULLPTR );
  SLIC_ASSERT( hi != AXOM_NULLPTR );
  SLIC_ASSERT( h != AXOM_NULLPTR );
  SLIC_ASSERT( origin != AXOM_NULLPTR );

  for ( int i=0; i < ndims; ++i )
  {
    origin[ i ] = lo[ i ];
    double dx   = hi[ i ] - lo[ i ];
    SLIC_ERROR_IF( utilities::isNearlyEqual( dx, 0.0 ) || dx < 0.0,
                   "supplied invalid bounds!" );
    h[ i ] = dx / static_cast< double >( ( extent->size( i ) - 1.0 ) );
  }

}

} // end anonymous namespace

//------------------------------------------------------------------------------
// UNIFORM MESH IMPLEMENTATION
//------------------------------------------------------------------------------
UniformMesh::UniformMesh( int dimension,
                          const double* origin,
                          const double* h,
                          const int64* ext ) :
  StructuredMesh( STRUCTURED_UNIFORM_MESH, dimension, ext )
{
  SLIC_ERROR_IF( origin==AXOM_NULLPTR, "supplied origin buffer is null" );
  SLIC_ERROR_IF( h==AXOM_NULLPTR, "supplied spacing buffer is null" );
  SLIC_ERROR_IF( ext==AXOM_NULLPTR, "supplied extent buffer is null" );

  const size_t bytesize = dimension * sizeof( double );
  memcpy( m_origin, origin, bytesize );
  memcpy( m_h,      h,      bytesize );
}

//------------------------------------------------------------------------------
UniformMesh::UniformMesh( int dimension,
                          const int64* ext,
                          const double* lower_bound,
                          const double* upper_bound ) :
  StructuredMesh( STRUCTURED_UNIFORM_MESH, dimension, ext )

{
  SLIC_ERROR_IF( ext==AXOM_NULLPTR, "supplied extent buffer is null" );
  SLIC_ERROR_IF( lower_bound==AXOM_NULLPTR, "supplied null for lower_bound" );
  SLIC_ERROR_IF( upper_bound==AXOM_NULLPTR, "supplied null for upper_bound" );

  set_spacing_and_origin( m_ndims, m_extent,
                          lower_bound, upper_bound, m_h, m_origin );
}

//------------------------------------------------------------------------------
UniformMesh::UniformMesh( int dimension,
                          const double* lower_bound,
                          const double* upper_bound,
                          IndexType Ni,
                          IndexType Nj,
                          IndexType Nk ) :
       StructuredMesh( STRUCTURED_UNIFORM_MESH, dimension )
{
  SLIC_ERROR_IF( lower_bound==AXOM_NULLPTR, "supplied null for lower_bound" );
  SLIC_ERROR_IF( upper_bound==AXOM_NULLPTR, "supplied null for upper_bound" );
  SLIC_ERROR_IF( Ni <= 0, "Ni must be greater or equal to 1" );
  SLIC_ERROR_IF( (dimension >= 2 ) && (Nj <= 0),
                 "Nj must be greater or equal to 1" );
  SLIC_ERROR_IF( (dimension==3) && (Nk <= 0), ""
                 "Nk must be greater or equal to 1" );

  int64 extent[ ] = { 0,Ni-1, 0, Nj-1, 0,Nk-1 };
  m_extent        = new mint::Extent( dimension, extent );
  set_spacing_and_origin( m_ndims, m_extent,
                          lower_bound, upper_bound, m_h, m_origin );
}

//------------------------------------------------------------------------------
void UniformMesh::getNode( IndexType nodeID, double* node ) const
{
  SLIC_ASSERT( 0 <= nodeID && nodeID < getNumberOfNodes() );
  SLIC_ASSERT( node != AXOM_NULLPTR );

  IndexType i=-1;
  IndexType j=-1;
  IndexType k=-1;
  switch ( m_ndims )
  {
  case 1:
    node[ 0 ] =evaluateCoordinate( nodeID, I_DIRECTION );
    break;
  case 2:
    m_extent->getGridIndex( nodeID, i, j );
    node[ 0 ] = evaluateCoordinate( i, I_DIRECTION );
    node[ 1 ] = evaluateCoordinate( j, J_DIRECTION );
    break;
  default:
    SLIC_ASSERT( m_ndims==3 );
    m_extent->getGridIndex( nodeID, i, j, k );
    node[ 0 ] = evaluateCoordinate( i, I_DIRECTION );
    node[ 1 ] = evaluateCoordinate( j, J_DIRECTION );
    node[ 2 ] = evaluateCoordinate( k, K_DIRECTION );
  } // END switch()

}

#ifdef MINT_USE_SIDRE

//------------------------------------------------------------------------------
UniformMesh::UniformMesh( sidre::Group* group,
                          const std::string& topo ) :
     StructuredMesh( group, topo )
{
  SLIC_ERROR_IF( m_type != STRUCTURED_UNIFORM_MESH,
                 "supplied Sidre group does not corresond to a UniformMesh!" );

  int64 extent[ 6 ];
  blueprint::getUniformMesh( m_ndims,
                             getCoordsetGroup(), getTopologyGroup(),
                             m_origin,
                             m_h,
                             extent );

  m_extent = new mint::Extent( m_ndims, extent );
}

//------------------------------------------------------------------------------
UniformMesh::UniformMesh( sidre::Group* group ) : UniformMesh( group, "" )
{ }

//------------------------------------------------------------------------------
UniformMesh::UniformMesh( int dimension,
                          const double* lower_bound,
                          const double* upper_bound,
                          const int64* extent,
                          sidre::Group* group,
                          const std::string& topo,
                          const std::string& coordset ) :
   StructuredMesh( STRUCTURED_UNIFORM_MESH, dimension, group, topo, coordset )
{
  SLIC_ERROR_IF( extent==AXOM_NULLPTR, "supplied extent buffer is null" );
  SLIC_ERROR_IF( lower_bound==AXOM_NULLPTR, "supplied null for lower_bound" );
  SLIC_ERROR_IF( upper_bound==AXOM_NULLPTR, "supplied null for upper_bound" );

  // STEP 0: initialize mesh
  m_extent = new mint::Extent( dimension, extent );
  set_spacing_and_origin( m_ndims, m_extent,
                            lower_bound, upper_bound, m_h, m_origin );

  // STEP 1: populate sidre
  blueprint::initializeTopologyGroup(  m_group, m_topology, m_coordset,
                                       "uniform" );

  blueprint::setUniformMesh( dimension, m_origin, m_h, m_extent,
                             getCoordsetGroup(),
                             getTopologyGroup() );
}

//------------------------------------------------------------------------------
UniformMesh::UniformMesh( int dimension,
                          const double* lower_bound,
                          const double* upper_bound,
                          const int64* extent,
                          sidre::Group* group ) :
  UniformMesh( dimension, lower_bound, upper_bound, extent, group, "", "" )
{ }

//------------------------------------------------------------------------------
UniformMesh::UniformMesh( int dimension,
                          const double* lower_bound,
                          const double* upper_bound,
                          sidre::Group* group,
                          const std::string& topo,
                          const std::string& coordset,
                          IndexType Ni,
                          IndexType Nj,
                          IndexType Nk ) :
  StructuredMesh( STRUCTURED_UNIFORM_MESH, dimension, group, topo, coordset )
{
  SLIC_ERROR_IF( lower_bound==AXOM_NULLPTR, "supplied null for lower_bound" );
  SLIC_ERROR_IF( upper_bound==AXOM_NULLPTR, "supplied null for upper_bound" );
  SLIC_ERROR_IF( Ni <= 0, "Ni must be greater or equal to 1" );
  SLIC_ERROR_IF( (dimension >= 2 ) && (Nj <= 0),
                 "Nj must be greater or equal to 1" );
  SLIC_ERROR_IF( (dimension==3) && (Nk <= 0), ""
                 "Nk must be greater or equal to 1" );

  // STEP 0: initialize mesh
  int64 extent[ ] = { 0,Ni-1, 0, Nj-1, 0,Nk-1 };
  m_extent        = new mint::Extent( dimension, extent );
  set_spacing_and_origin( m_ndims, m_extent,
                          lower_bound, upper_bound, m_h, m_origin );

  // STEP 1: populate sidre
  blueprint::initializeTopologyGroup(  m_group, m_topology, m_coordset,
                                       "uniform" );

  blueprint::setUniformMesh( dimension, m_origin, m_h, m_extent,
                             getCoordsetGroup(),
                             getTopologyGroup() );
}

//------------------------------------------------------------------------------
UniformMesh::UniformMesh( int dimension,
                          const double* lower_bound,
                          const double* upper_bound,
                          sidre::Group* group,
                          IndexType Ni,
                          IndexType Nj,
                          IndexType Nk ) :
  UniformMesh( dimension, lower_bound, upper_bound, group, "", "", Ni, Nj, Nk )
{ }

#endif

} /* namespace mint */
} /* namespace axom */
