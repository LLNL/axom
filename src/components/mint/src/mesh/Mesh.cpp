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

#ifdef MINT_USE_SIDRE
#include "sidre/sidre.hpp"
#endif

// C/C++ includes
#include <cstring> // for strcmp()

namespace axom
{
namespace mint
{


//------------------------------------------------------------------------------
Mesh::Mesh( int ndims, int type, int blockId, int partId ) :
  m_ndims( ndims ),
  m_type( type ),
  m_block_idx( blockId ),
  m_part_idx( partId ),
  m_num_cells( 0 ),
  m_num_faces( 0 ),
  m_num_edges( 0 ),
  m_num_nodes( 0 ),
  m_explicit_coords( false ),
  m_explicit_connectivity( false ),
  m_has_mixed_topology( false )
#ifdef MINT_USE_SIDRE
  , m_group( AXOM_NULLPTR ),
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
  m_block_idx( 0 ),
  m_part_idx( 0 ),
  m_num_cells( 0 ),
  m_num_faces( 0 ),
  m_num_edges( 0 ),
  m_num_nodes( 0 ),
  m_explicit_coords( false ),
  m_explicit_connectivity( false ),
  m_has_mixed_topology( false ),
  m_group( group ),
  m_topology( topo )
{
  SLIC_ERROR_IF( m_group==AXOM_NULLPTR, "NULL sidre group" );

  if ( m_group->hasChildView("state") )
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
  
  const bool hasCoordsets  = m_group->hasChildGroup( "coordsets" );
  const bool hasTopologies = m_group->hasChildGroup( "topologies" );
  const bool hasFields     = m_group->hasChildGroup( "fields" );

  SLIC_ERROR_IF( !hasCoordsets, "sidre::Group " <<  m_group->getPathName() << 
                 " is missing coordsets group!" );
  SLIC_ERROR_IF( !hasTopologies, "sidre::Group " <<  m_group->getPathName() << 
                 " is missing topologies group!" );
  SLIC_ERROR_IF( !hasFields, "sidre::Group " <<  m_group->getPathName() << 
                 " is missing fields group!" );

  sidre::Group* topologies_group = m_group->getGroup( "topologies" );
  SLIC_ASSERT( topologies_group != AXOM_NULLPTR );

  sidre::Group* coordsets_group = m_group->getGroup( "coordsets" );
  SLIC_ASSERT( coordsets_group != AXOM_NULLPTR );

  if ( m_topology == "" )
  {
    SLIC_ERROR_IF( topologies_group->getNumGroups() != 1, 
              "No topology was specified so the topologies group must have " <<
              "a single child group not " << topologies_group->getNumGroups() );
    m_topology = topologies_group->getGroup(0)->getName();
  }
  SLIC_ERROR_IF( topologies_group->hasChildGroup( m_topology ), 
                 "No topology named " << m_topology << " found in " << 
                  topologies_group->getPathName() );
  sidre::Group* topology = topologies_group->getGroup( m_topology );
  SLIC_ASSERT( topology != AXOM_NULLPTR );

  SLIC_ERROR_IF( topology->hasChildView( "coordset" ), 
                 topology->getPathName() << " has no associated coordset.");

  sidre::View* coordset_view = topology->getView( "coordset" );
  SLIC_ASSERT( coordset_view != AXOM_NULLPTR );
  SLIC_ERROR_IF( coordset_view->isString(), 
                 "topology coordset view needs to hold a string." );

  const char* coordset_name = coordset_view->getString();
  SLIC_ERROR_IF( coordsets_group->hasChildGroup( coordset_name ), 
                 "No coordset named " << coordset_name << " found in " << 
                  coordsets_group->getPathName() );

  sidre::Group* coords = coordsets_group->getGroup( coordset_name );
  SLIC_ASSERT( coords != AXOM_NULLPTR );

  const char* coord_type = coords->getView( "type" )->getString();
  SLIC_ASSERT( coord_type != AXOM_NULLPTR );

  const char* topo_type = topology->getView( "type" )->getString();
  SLIC_ASSERT( topo_type != AXOM_NULLPTR );

  if ( strcmp( topo_type, "uniform" )==0 )
  {
    m_type  = UNIFORM_MESH;
    m_ndims = coords->getGroup("origin")->getNumViews();

  } // END if UNIFORM MESH
  else if ( strcmp( topo_type, "rectilinear" )== 0 )
  {
    m_type  = RECTILINEAR_MESH;
    m_ndims = coords->getGroup( "values" )->getNumViews();

  } // END if RECTILINEAR_MESH
  else if ( strcmp( topo_type, "structured" )==0 )
  {
    m_type  = STRUCTURED_MESH;
    m_ndims = coords->getGroup( "values" )->getNumViews();

  } // END if STRUCTURED_MESH
  else if ( strcmp( topo_type, "particle" )==0 )
  {
    // NOTE: currently the blueprint doesn't provide a topology type
    // that indicates particles
    m_type  = PARTICLE_MESH;
    m_ndims = coords->getGroup( "values" )->getNumViews();

  } // END if PARTICLE_MESH
  else if ( strcmp( topo_type, "unstructured" )==0 )
  {
    // check if this is a particle mesh stored as an unstructured mesh
    const char* shape = topology->getView("elements/shape")->getString();
    if ( strcmp( shape, "point" ) == 0 )
    {
      m_type = PARTICLE_MESH;
    }
    else 
    {
      m_type = UNSTRUCTURED_MESH;
    }

  } // END if UNSTRUCTURED_MESH
  else
  {
    m_type = UNDEFINED_MESH;
    SLIC_ERROR( "invalid blueprint mesh coord_type=[" << coord_type << "] " <<
                "topo_type=[" << topo_type << "] " );
  }

  SLIC_ERROR_IF( !validMeshType(), "invalid mesh type=" << m_type );
  SLIC_ERROR_IF( !validDimension(), "invalid mesh dimension=" << m_ndims );

  allocateFieldData();
}

//------------------------------------------------------------------------------
Mesh::Mesh( sidre::Group* group, const std::string& topo, const std::string& coordset,
            int ndims, int type,
            int blockId, int partId ) :
  m_ndims( ndims ),
  m_type( type ),
  m_block_idx( blockId ),
  m_part_idx( partId ),
  m_num_cells( 0 ),
  m_num_faces( 0),
  m_num_edges( 0 ),
  m_num_nodes( 0 ),
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

  sidre::Group* state_group = m_group->createGroup( "state_group" );
  state_group->createView( "block_id" )->setScalar( m_block_idx );
  state_group->createView( "partition_id" )->setScalar( m_part_idx );
  m_group->createGroup( "coordsets" )->createGroup( coordset );
  m_group->createGroup( "topologies" )->createGroup( topo );
  m_group->createGroup( "fields" );

  allocateFieldData();
}

#endif

//------------------------------------------------------------------------------
Mesh::~Mesh()
{
  deallocateFieldData();
}

//------------------------------------------------------------------------------
void Mesh::getMeshNode( IndexType nodeIdx, double* node ) const
{
  // TODO: implement this
}

//------------------------------------------------------------------------------
void Mesh::getMeshCell( IndexType cellIdx, IndexType* cell ) const
{
  // TODO: implement this
}

//------------------------------------------------------------------------------
CellType Mesh::getMeshCellType( IndexType cellIdx ) const
{
  // TODO: implement this
  return -1;
}

//------------------------------------------------------------------------------
void Mesh::allocateFieldData( )
{
#ifdef MINT_USE_SIDRE

  if ( hasSidreGroup() )
  {
    sidre::Group* fields_group = ( m_group->hasChildGroup("fields") ?
            m_group->getGroup( "fields") : m_group->createGroup( "fields") );
    SLIC_ASSERT( fields_group != AXOM_NULLPTR );

    for ( int i=0; i < NUM_FIELD_ASSOCIATIONS; ++i )
    {
      m_mesh_fields[ i ] = new FieldData( i, fields_group, m_topology );
    }
  }
  else
  {
    for ( int i=0; i < NUM_FIELD_ASSOCIATIONS; ++i )
    {
       m_mesh_fields[ i ] = new FieldData( i );
    }
  }

#else

  for ( int i=0; i < NUM_FIELD_ASSOCIATIONS; ++i )
  {
    m_mesh_fields[ i ] = new FieldData( i );
  }

#endif

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

} /* namespace mint */
} /* namespace axom */
