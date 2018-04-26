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

#include "mint/Mesh.hpp"

// axom includes
#include "axom/Types.hpp"

#include "mint/FieldData.hpp"
#include "mint/ParticleMesh.hpp"
#include "mint/UnstructuredMesh.hpp"

#ifdef MINT_USE_SIDRE
#include "sidre/sidre.hpp"
#endif

#include "mint/blueprint.hpp" // for blueprint functions

// C/C++ includes
#include <cstring> // for strcmp()

namespace axom
{
namespace mint
{


//------------------------------------------------------------------------------
Mesh::Mesh( int ndims, int type ) :
  m_ndims( ndims ),
  m_type( type ),
  m_block_idx( -1 ),
  m_part_idx( -1 ),
  m_explicit_coords( false ),
  m_explicit_connectivity( false ),
  m_has_mixed_topology( false ),
#ifdef MINT_USE_SIDRE
  m_group( AXOM_NULLPTR ),
  m_topology()
#endif
{
  SLIC_ERROR_IF( !validMeshType(), "invalid mesh type=" << m_type );
  SLIC_ERROR_IF( !validDimension(), "invalid mesh dimension=" << m_ndims );
  allocateFieldData();
}

#ifdef MINT_USE_SIDRE

//------------------------------------------------------------------------------
Mesh::Mesh( sidre::Group* group, const std::string& topo ) :
  m_ndims( -1 ),
  m_type( UNDEFINED_MESH ),
  m_block_idx( -1 ),
  m_part_idx( -1 ),
  m_explicit_coords( false ),
  m_explicit_connectivity( false ),
  m_has_mixed_topology( false ),
  m_group( group ),
  m_topology( topo )
{
  SLIC_ERROR_IF( m_group==AXOM_NULLPTR, "NULL sidre group" );
  SLIC_ERROR_IF( ! blueprint::validRootGroup( m_group ),
                "root group does not conform to blueprint" );

  if ( m_group->hasChildGroup("state") )
  {
    sidre::Group* state_group = m_group->getGroup( "state" );
    if ( state_group->hasChildView( "block_id" ) )
    {
      m_block_idx = state_group->getView( "block_id" )->getScalar();
    }

    if ( state_group->hasChildView( "partition_id" ) )
    {
      m_part_idx = state_group->getView( "partition_id" )->getScalar();
    }
  }

  blueprint::getMeshTypeAndDimension( m_type, m_ndims, m_group, m_topology );
  m_topology = getTopologyGroup( )->getName( );

  SLIC_ERROR_IF( !validMeshType(), "invalid mesh type=" << m_type );
  SLIC_ERROR_IF( !validDimension(), "invalid mesh dimension=" << m_ndims );

  allocateFieldData();
}

Mesh::Mesh( sidre::Group* group ) : Mesh( group, "" )
{}

Mesh::Mesh( int ndims, int type, sidre::Group* group, const std::string& topo,
            const std::string& coordset ) :
  m_ndims( ndims ),
  m_type( type ),
  m_block_idx( -1 ),
  m_part_idx( -1 ),
  m_explicit_coords( false ),
  m_explicit_connectivity( false ),
  m_has_mixed_topology( false ),
  m_group( group ),
  m_topology( topo )
{
  SLIC_ERROR_IF( !validMeshType(), "invalid mesh type=" << m_type );
  SLIC_ERROR_IF( !validDimension(), "invalid mesh dimension=" << m_ndims );
  SLIC_ERROR_IF( m_group==AXOM_NULLPTR, "NULL sidre group" );
  SLIC_ERROR_IF( m_group->getNumGroups() != 0, "group is not empty!" );
  SLIC_ERROR_IF( m_group->getNumViews() != 0, "group is not empty!" );

  // provide default names if not specified
  m_topology = ( topo.empty() )? "t1" : topo;
  std::string coordset_name = ( coordset.empty() )? "c1" : coordset;

  sidre::Group* state_group = m_group->createGroup( "state" );
  state_group->createView( "block_id" )->setScalar( m_block_idx );
  state_group->createView( "partition_id" )->setScalar( m_part_idx );

  // create the coordset group
  m_group->createGroup( "coordsets" )->createGroup( coordset_name );

  // create the topology group for this mesh and link to the coordset
  sidre::Group* t =
      m_group->createGroup( "topologies" )->createGroup( m_topology);
  t->createView( "coordset" )->setString( coordset_name );
  t->createView( "type" )->setString( "undefined" );

  // create the fields group  for this mesh
  m_group->createGroup( "fields" );

  allocateFieldData();
}

//------------------------------------------------------------------------------
Mesh::Mesh( int ndims, int type, sidre::Group* group ) :
  Mesh( ndims, type, group, "", "" )
{}

//------------------------------------------------------------------------------
sidre::Group* Mesh::getCoordsetGroup( )
{
  const sidre::Group* g =
      blueprint::getCoordsetGroup( m_group, getTopologyGroup() );
  return (  const_cast< sidre::Group* >( g ) );
}

//------------------------------------------------------------------------------
sidre::Group* Mesh::getTopologyGroup( )
{
  const sidre::Group* g = blueprint::getTopologyGroup( m_group, m_topology );
  return ( const_cast< sidre::Group* >( g ) );
}

#endif

//------------------------------------------------------------------------------
Mesh::~Mesh()
{
  deallocateFieldData();
}

//------------------------------------------------------------------------------
void Mesh::setBlockId( int ID )
{
  m_block_idx = ID;
#ifdef MINT_USE_SIDRE
  if ( hasSidreGroup() )
  {
    sidre::Group* state_group = m_group->getGroup( "state" );
    SLIC_ASSERT( state_group != AXOM_NULLPTR );
    sidre::View* block_view = state_group->getView( "block_id" );
    SLIC_ASSERT( block_view != AXOM_NULLPTR );
    block_view->setScalar( m_block_idx );
  }
#endif
}

//------------------------------------------------------------------------------
void Mesh::setPartitionId( int ID )
{
  m_part_idx = ID;
#ifdef MINT_USE_SIDRE
  if ( hasSidreGroup() )
  {
    sidre::Group* state_group = m_group->getGroup( "state" );
    SLIC_ASSERT( state_group != AXOM_NULLPTR );
    sidre::View* partition_view = state_group->getView( "partition_id" );
    SLIC_ASSERT( partition_view != AXOM_NULLPTR );
    partition_view->setScalar( m_part_idx );
  }
#endif
}

//------------------------------------------------------------------------------
void Mesh::allocateFieldData( )
{
#ifdef MINT_USE_SIDRE
  if ( hasSidreGroup() )
  {
    sidre::Group* fields_group = m_group->getGroup( "fields");
    SLIC_ASSERT( fields_group != AXOM_NULLPTR );

    for ( int assoc=0; assoc < NUM_FIELD_ASSOCIATIONS; ++assoc )
    {
      m_mesh_fields[ assoc ] = new FieldData( assoc, fields_group, m_topology );
    }

    return;
  }
#endif

  for ( int assoc=0; assoc < NUM_FIELD_ASSOCIATIONS; ++assoc )
  {
    m_mesh_fields[ assoc ] = new FieldData( assoc );
  }

}

//------------------------------------------------------------------------------
void Mesh::deallocateFieldData( )
{
  for ( int i=0; i < NUM_FIELD_ASSOCIATIONS; ++i )
  {
    SLIC_ASSERT( m_mesh_fields[ i ] != AXOM_NULLPTR );

    delete m_mesh_fields[ i ];
    m_mesh_fields[ i ] = AXOM_NULLPTR;
  }
}

//------------------------------------------------------------------------------
Mesh* Mesh::getMesh( sidre::Group* group, const std::string& topo )
{
  SLIC_ERROR_IF( group==AXOM_NULLPTR, "supplied group is null" );

  int mesh_type = UNDEFINED_MESH;
  int dimension = -1;
  blueprint::getMeshTypeAndDimension( mesh_type, dimension, group, topo );

  Mesh* m = AXOM_NULLPTR;
  switch ( mesh_type )
  {
  case STRUCTURED_CURVILINEAR_MESH:
    SLIC_ERROR( "!!! NOT IMPLEMENTED YET !!!" );
    // TODO: implement this
    break;
  case STRUCTURED_RECTILINEAR_MESH:
    SLIC_ERROR( "!!! NOT IMPLEMENTED YET !!!" );
    // TODO: implement this
    break;
  case STRUCTURED_UNIFORM_MESH:
    SLIC_ERROR( "!!! NOT IMPLEMENTED YET !!!" );
    // TODO: implement this
    break;
  case UNSTRUCTURED_MESH:
    SLIC_ERROR( "!!! NOT IMPLEMENTED YET !!!" );
    // TODO: implement this
    break;
  case PARTICLE_MESH:
    m = new ParticleMesh( const_cast< sidre::Group* >( group ), topo );
    break;
  default:
    SLIC_ERROR( "undefined mesh_type [" << mesh_type << "]\n" );
  } // END switch

  SLIC_ASSERT( m != AXOM_NULLPTR );
  return ( m );
}

//------------------------------------------------------------------------------
Mesh* Mesh::getMesh( sidre::Group* group )
{
  return ( Mesh::getMesh( group, "" ) );
}

} /* namespace mint */
} /* namespace axom */
