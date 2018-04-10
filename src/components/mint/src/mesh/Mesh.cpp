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
  m_coordinates( AXOM_NULLPTR ),
  m_explicit_coords( false ),
  m_explicit_connectivity( false ),
  m_has_mixed_topology( false ),
#ifdef MINT_USE_SIDRE
  m_group( AXOM_NULLPTR ),
  m_fields_group( AXOM_NULLPTR ),
  m_coordsets_group( AXOM_NULLPTR ),
  m_topologies_group( AXOM_NULLPTR )
#endif
{
  SLIC_ERROR_IF( !validMeshType(), "invalid mesh type=" << m_type );
  SLIC_ERROR_IF( !validDimension(), "invalid mesh dimension=" << m_ndims );
}

#ifdef MINT_USE_SIDRE
//------------------------------------------------------------------------------
Mesh::Mesh( sidre::Group* group ) :
  m_ndims( -1 ),
  m_type( mint::UNDEFINED_MESH ),
  m_block_idx( 0 ),
  m_part_idx( 0 ),
  m_num_cells( 0 ),
  m_num_faces( 0 ),
  m_num_edges( 0 ),
  m_num_nodes( 0 ),
  m_coordinates( AXOM_NULLPTR ),
  m_explicit_coords( false ),
  m_explicit_connectivity( false ),
  m_has_mixed_topology( false ),
  m_group( group ),
  m_fields_group( AXOM_NULLPTR ),
  m_coordsets_group( AXOM_NULLPTR ),
  m_topologies_group( AXOM_NULLPTR )
{
  SLIC_ERROR_IF( m_group==AXOM_NULLPTR, "NULL sidre group" );

  if ( m_group->hasChildView( "block_id" ) )
  {
    m_block_idx = m_group->getView( "block_id" )->getScalar();
  }

  if ( m_group->hasChildView( "partition_id" ) )
  {
    m_part_idx = m_group->getView( "partition_id" )->getScalar();
  }

  const bool hasCoordsets  = m_group->hasChildGroup( "coordsets" );
  const bool hasTopologies = m_group->hasChildGroup( "topologies" );
  const bool hasFields     = m_group->hasChildGroup( "fields" );

  SLIC_ERROR_IF( ! hasCoordsets, "missing coordsets group!" );
  SLIC_ERROR_IF( ! hasTopologies, "missing topologies group!" );
  SLIC_ERROR_IF( ! hasFields, "missing fields group!" );

  m_coordsets_group  = m_group->getGroup( "coordsets" );
  m_topologies_group = m_group->getGroup( "topologies" );
  m_fields_group     = m_group->getGroup( "fields" );

  detectMeshTypeAndDimension( );
}

//------------------------------------------------------------------------------
Mesh::Mesh( sidre::Group* group,
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
  m_coordinates( AXOM_NULLPTR ),
  m_explicit_coords( false ),
  m_explicit_connectivity( false ),
  m_has_mixed_topology( false ),
  m_group( group ),
  m_fields_group( AXOM_NULLPTR ),
  m_coordsets_group( AXOM_NULLPTR ),
  m_topologies_group( AXOM_NULLPTR )
{
  SLIC_ERROR_IF( !validMeshType(), "invalid mesh type=" << m_type );
  SLIC_ERROR_IF( !validDimension(), "invalid mesh dimension=" << m_ndims );
  SLIC_ERROR_IF( m_group==AXOM_NULLPTR, "NULL sidre group" );
  SLIC_ERROR_IF( m_group->getNumGroups() != 0, "group is not empty!" );
  SLIC_ERROR_IF( m_group->getNumViews() != 0, "group is not empty!" );

  m_group->createView( "block_id" )->setScalar( m_block_idx );
  m_group->createView( "partition_id" )->setScalar( m_part_idx );
  m_coordsets_group  = m_group->createGroup( "coordsets/c1" )->getParent();
  m_topologies_group = m_group->createGroup( "topologies/t1" )->getParent();
  m_fields_group     = m_group->createGroup( "fields" );
}

//------------------------------------------------------------------------------
void Mesh::detectMeshTypeAndDimension( )
{
  SLIC_ASSERT( m_coordsets_group != AXOM_NULLPTR );
  SLIC_ASSERT( m_topologies_group != AXOM_NULLPTR );

  // NOTE: the code assumes one coordset and group
  SLIC_ERROR_IF( m_coordsets_group->getNumGroups() > 1,
                 "code currently assumes one coordset per Mesh object!" );
  SLIC_ERROR_IF( m_topologies_group->getNumGroups() > 1,
                 "code currently assumes one topology per Mesh object!" );

  sidre::Group* coords = m_coordsets_group->getGroup( 0 );
  SLIC_ASSERT( coords != AXOM_NULLPTR );

  sidre::Group* topo = m_topologies_group->getGroup( 0 );
  SLIC_ASSERT( topo != AXOM_NULLPTR );

  const char* coord_type = coords->getView( "type" )->getString();
  SLIC_ASSERT( coord_type != AXOM_NULLPTR );

  const char* topo_type = topo->getView( "type" )->getString();
  SLIC_ASSERT( topo_type != AXOM_NULLPTR );

  if ( strcmp( topo_type, "uniform" )==0 )
  {
    m_type  = mint::UNIFORM_MESH;
    m_ndims = coords->getGroup("origin")->getNumViews();

  } // END if UNIFORM MESH
  else if ( strcmp( topo_type, "rectilinear" )== 0 )
  {
    m_type  = mint::RECTILINEAR_MESH;
    m_ndims = coords->getGroup( "values" )->getNumViews();

  } // END if RECTILINEAR_MESH
  else if ( strcmp( topo_type, "structured" )==0 )
  {
    m_type  = mint::STRUCTURED_MESH;
    m_ndims = coords->getGroup( "values" )->getNumViews();

  } // END if STRUCTURED_MESH
  else if ( strcmp( topo_type, "particle" )==0 )
  {
    // NOTE: currently the blue-print doesn't provide a topology type
    // that indicates particles
    m_type  = mint::PARTICLE_MESH;
    m_ndims = coords->getGroup( "values" )->getNumViews();

  } // END if PARTICLE_MESH
  else if ( strcmp( topo_type, "unstructured" )==0 )
  {
    // check if this is a particle mesh stored as an unstructured mesh
    const char* shape = topo->getView("elements/shape")->getString();
    m_type = ( strcmp(shape,"point")==0 ) ?
                  mint::PARTICLE_MESH : mint::UNSTRUCTURED_MESH;
    m_ndims = coords->getGroup( "values" )->getNumViews();

  } // END if UNSTRUCTURED_MESH
  else
  {
    m_type = mint::UNDEFINED_MESH;
    SLIC_ERROR( "invalid blueprint mesh coord_type=[" << coord_type << "] " <<
                "topo_type=[" << topo_type << "] " );
  }

  SLIC_ERROR_IF( !validMeshType(), "invalid mesh type=" << m_type );
  SLIC_ERROR_IF( !validDimension(), "invalid mesh dimension=" << m_ndims );
}

#endif

//------------------------------------------------------------------------------
Mesh::~Mesh()
{

  if ( m_coordinates != AXOM_NULLPTR )
  {
    delete m_coordinates;
    m_coordinates = AXOM_NULLPTR;
  }

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
int Mesh::getMeshCellType( IndexType cellIdx ) const
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
    SLIC_ASSERT( fields_group->getParent()==m_group );

    for ( int i=0; i < NUM_FIELD_ASSOCIATIONS; ++i )
    {
      m_mesh_fields[ i ] = new mint::FieldData( i, fields_group );
    }
  }
  else
  {
    for ( int i=0; i < NUM_FIELD_ASSOCIATIONS; ++i )
    {
       m_mesh_fields[ i ] = new mint::FieldData( i );
    }
  }

#else

  for ( int i=0; i < NUM_FIELD_ASSOCIATIONS; ++i )
  {
    m_mesh_fields[ i ] = new mint::FieldData( i );
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
