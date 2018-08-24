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
#include "axom/mint/mesh/UniformMesh.hpp"

// mint includes
#include "axom/mint/mesh/blueprint.hpp"       // for blueprint functions
#include "axom/mint/config.hpp"          // for compile-time definitions
#include "axom/mint/mesh/MeshTypes.hpp"       // for STRUCTURED_UNIFORM_MESH

#include "axom/mint/mesh/internal/MeshHelpers.hpp"     // for internal helper

// core includes
#include "axom/core/utilities/Utilities.hpp" // for utilities::isNearlyEqual

namespace axom
{
namespace mint
{

//------------------------------------------------------------------------------
// UNIFORM MESH IMPLEMENTATION
//------------------------------------------------------------------------------
UniformMesh::UniformMesh( int dimension, const double* origin, const double* h,
                          const IndexType* ext ) :
  StructuredMesh( STRUCTURED_UNIFORM_MESH, dimension, ext )
{
  SLIC_ERROR_IF( origin==nullptr, "supplied origin buffer is null" );
  SLIC_ERROR_IF( h==nullptr, "supplied spacing buffer is null" );
  SLIC_ERROR_IF( ext==nullptr, "supplied extent buffer is null" );

  const size_t bytesize = dimension * sizeof( double );
  memcpy( m_origin, origin, bytesize );
  memcpy( m_h,      h,      bytesize );
}

//------------------------------------------------------------------------------
UniformMesh::UniformMesh( int dimension,
                          const IndexType* ext, 
                          const double* lower_bound, 
                          const double* upper_bound ) :
  StructuredMesh( STRUCTURED_UNIFORM_MESH, dimension, ext )
{
  SLIC_ERROR_IF( ext==nullptr, "supplied extent buffer is null" );
  SLIC_ERROR_IF( lower_bound==nullptr, "supplied null for lower_bound" );
  SLIC_ERROR_IF( upper_bound==nullptr, "supplied null for upper_bound" );

  setSpacingAndOrigin( lower_bound, upper_bound );
}

//------------------------------------------------------------------------------
UniformMesh::UniformMesh( const double* lower_bound, const double* upper_bound,
                          IndexType Ni, IndexType Nj, IndexType Nk ) :
  StructuredMesh( STRUCTURED_UNIFORM_MESH, Ni, Nj, Nk )
{
  SLIC_ERROR_IF( lower_bound==nullptr, "supplied null for lower_bound" );
  SLIC_ERROR_IF( upper_bound==nullptr, "supplied null for upper_bound" );

  setSpacingAndOrigin( lower_bound, upper_bound );
}


#ifdef AXOM_MINT_USE_SIDRE

//------------------------------------------------------------------------------
UniformMesh::UniformMesh( sidre::Group* group, const std::string& topo ) :
  StructuredMesh( group, topo )
{
  SLIC_ERROR_IF( m_type != STRUCTURED_UNIFORM_MESH,
                 "supplied Sidre group does not correspond to a UniformMesh!" );

  blueprint::getUniformMesh( m_ndims, m_origin, m_h, getCoordsetGroup() );
}

//------------------------------------------------------------------------------
UniformMesh::UniformMesh( int dimension,
                          const double* lower_bound,
                          const double* upper_bound,
                          const IndexType* extent,
                          sidre::Group* group,
                          const std::string& topo,
                          const std::string& coordset ) :
  StructuredMesh( STRUCTURED_UNIFORM_MESH, dimension, extent, group, topo,
                  coordset )
{
  SLIC_ERROR_IF( extent == nullptr, "supplied extent buffer is null" );
  SLIC_ERROR_IF( lower_bound == nullptr, "supplied null for lower_bound" );
  SLIC_ERROR_IF( upper_bound == nullptr, "supplied null for upper_bound" );

  // STEP 0: initialize mesh
  setSpacingAndOrigin( lower_bound, upper_bound );

  // STEP 1: populate sidre
  blueprint::setUniformMesh( m_ndims, m_origin, m_h, getCoordsetGroup() );
}

//------------------------------------------------------------------------------
UniformMesh::UniformMesh( const double* lower_bound,
                          const double* upper_bound,
                          sidre::Group* group,
                          const std::string& topo,
                          const std::string& coordset,
                          IndexType Ni,
                          IndexType Nj,
                          IndexType Nk ) :
  StructuredMesh( STRUCTURED_UNIFORM_MESH, Ni, Nj, Nk, group, topo, coordset )
{
  SLIC_ERROR_IF( lower_bound == nullptr, "supplied null for lower_bound" );
  SLIC_ERROR_IF( upper_bound == nullptr, "supplied null for upper_bound" );

  // STEP 0: initialize mesh
  setSpacingAndOrigin( lower_bound, upper_bound );

  // STEP 1: populate sidre
  blueprint::setUniformMesh( m_ndims, m_origin, m_h, getCoordsetGroup() );
}

#endif

//------------------------------------------------------------------------------
void UniformMesh::getNode( IndexType nodeID, double* node ) const
{
  SLIC_ASSERT( 0 <= nodeID && nodeID < getNumberOfNodes() );
  SLIC_ASSERT( node != nullptr );

  IndexType i = -1;
  IndexType j = -1;
  IndexType k = -1;
  switch ( m_ndims )
  {
  case 1:
    node[ 0 ] = evaluateCoordinate( nodeID, I_DIRECTION );
    break;
  case 2:
    getNodeGridIndex( nodeID, i, j );
    node[ 0 ] = evaluateCoordinate( i, I_DIRECTION );
    node[ 1 ] = evaluateCoordinate( j, J_DIRECTION );
    break;
  default:
    SLIC_ASSERT( m_ndims == 3 );
    getNodeGridIndex( nodeID, i, j, k );
    node[ 0 ] = evaluateCoordinate( i, I_DIRECTION );
    node[ 1 ] = evaluateCoordinate( j, J_DIRECTION );
    node[ 2 ] = evaluateCoordinate( k, K_DIRECTION );
  } // END switch()

}

void UniformMesh::setSpacingAndOrigin( const double* lo, const double* hi )
{
  SLIC_ASSERT( lo != nullptr );
  SLIC_ASSERT( hi != nullptr );

  for ( int dim = 0 ; dim < m_ndims ; ++dim )
  {
    m_origin[ dim ] = lo[ dim ];
    double dx   = hi[ dim ] - lo[ dim ];
    SLIC_ERROR_IF( utilities::isNearlyEqual( dx, 0.0 ) || dx < 0.0,
                   "supplied invalid bounds!" );
    m_h[ dim ] = dx / getCellDimension( dim );
  }
}

} /* namespace mint */
} /* namespace axom */
