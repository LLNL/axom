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

#include "mint/blueprint.hpp" // for blueprint functions

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
  SLIC_ERROR_IF( ! blueprint::validRootGroup( m_group ),
                "root group does not conform to blueprint" );

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

  blueprint::getMeshTypeAndDimension( m_type, m_ndims, m_group, m_topology );

  SLIC_ERROR_IF( !validMeshType(), "invalid mesh type=" << m_type );
  SLIC_ERROR_IF( !validDimension(), "invalid mesh dimension=" << m_ndims );

  allocateFieldData();
}

//------------------------------------------------------------------------------
Mesh::Mesh( sidre::Group* group ) : Mesh( group, "" )
{

}

//------------------------------------------------------------------------------
Mesh::Mesh( int ndims, int type, int blockId, int partId,
             sidre::Group* group,
             const std::string& topo,
             const std::string& coordset ) :
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

  // provide default names if not specified
  m_topology = ( topo.empty() )? "t1" : topo;
  std::string coordset_name = ( coordset.empty() )? "c1" : coordset;

  sidre::Group* state_group = m_group->createGroup( "state_group" );
  state_group->createView( "block_id" )->setScalar( m_block_idx );
  state_group->createView( "partition_id" )->setScalar( m_part_idx );
  m_group->createGroup( "coordsets" )->createGroup( coordset_name );
  m_group->createGroup( "topologies" )->createGroup( m_topology);
  m_group->createGroup( "fields" );

  allocateFieldData();
}

//------------------------------------------------------------------------------
Mesh::Mesh( int ndims, int type,
            int blockId, int partId, sidre::Group* group) :
                Mesh( ndims, type, blockId, partId, group, "", "" )
{

}

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
double* Mesh::getCoordinateArray( int dim )
{
  SLIC_ERROR_IF( !hasExplicitCoordinates(),
        "mesh of type [" << m_type << "] does not have explicit coordinates" );
  SLIC_ERROR_IF( ( (dim >= 0) && ( dim < getDimension() ) ),
    "requested coordinate array dim=[" << dim << "] on a mesh of dimension [" <<
  getDimension() << "]" );

  switch ( m_type )
  {
  case PARTICLE_MESH:
    SLIC_ERROR( "!!! NOT IMPLEMENTED YET !!!" );
    // TODO: implement this
    break;
  case STRUCTURED_MESH:
    SLIC_ERROR( "!!! NOT IMPLEMENTED YET !!!" );
    // TODO: implement this
    break;
  case RECTILINEAR_MESH:
    SLIC_ERROR( "!!! NOT IMPLEMENTED YET !!!" );
    // TODO: implement this
    break;
  case UNSTRUCTURED_MESH:
    SLIC_ERROR( "!!! NOT IMPLEMENTED YET !!!" );
    // TODO: implement this
    break;
  default:
    SLIC_ERROR( "Undefined MeshType!" );
  } // END switch

  return AXOM_NULLPTR;
}

//------------------------------------------------------------------------------
inline const double* Mesh::getCoordinateArray( int dim ) const
{
  SLIC_ERROR_IF( !hasExplicitCoordinates(),
       "mesh of type [" << m_type << "] does not have explicit coordinates" );
  SLIC_ERROR_IF( ( (dim >= 0) && ( dim < getDimension() ) ),
   "requested coordinate array dim=[" << dim << "] on a mesh of dimension [" <<
   getDimension() << "]" );

  switch ( m_type )
  {
  case PARTICLE_MESH:
    SLIC_ERROR( "!!! NOT IMPLEMENTED YET !!!" );
    // TODO: implement this
    break;
  case STRUCTURED_MESH:
    SLIC_ERROR( "!!! NOT IMPLEMENTED YET !!!" );
    // TODO: implement this
    break;
  case RECTILINEAR_MESH:
    SLIC_ERROR( "!!! NOT IMPLEMENTED YET !!!" );
    // TODO: implement this
    break;
  case UNSTRUCTURED_MESH:
    SLIC_ERROR( "!!! NOT IMPLEMENTED YET !!!" );
    // TODO: implement this
    break;
  default:
    SLIC_ERROR( "Undefined MeshType!" );
  } // END switch

  return AXOM_NULLPTR;
}

//------------------------------------------------------------------------------
void Mesh::getMeshNode( IndexType nodeIdx, double* node ) const
{
  // sanity checks
  SLIC_ASSERT( node != AXOM_NULLPTR );
  SLIC_ASSERT( nodeIdx >= 0 && nodeIdx < m_num_nodes );

  switch ( m_type )
  {
  case PARTICLE_MESH:
    SLIC_ERROR( "!!! NOT IMPLEMENTED YET !!!" );
    // TODO: implement this
    break;
  case UNIFORM_MESH:
    SLIC_ERROR( "!!! NOT IMPLEMENTED YET !!!" );
    // TODO: implement this
    break;
  case STRUCTURED_MESH:
    SLIC_ERROR( "!!! NOT IMPLEMENTED YET !!!" );
    // TODO: implement this
    break;
  case RECTILINEAR_MESH:
    SLIC_ERROR( "!!! NOT IMPLEMENTED YET !!!" );
    // TODO: implement this
    break;
  case UNSTRUCTURED_MESH:
    SLIC_ERROR( "!!! NOT IMPLEMENTED YET !!!" );
    // TODO: implement this
    break;
  default:
    SLIC_ERROR( "Undefined MeshType!" );
  } // END switch

}

//------------------------------------------------------------------------------
void Mesh::getMeshCell( IndexType cellIdx, IndexType* cell ) const
{
  // sanity checks
  SLIC_ASSERT( cell != AXOM_NULLPTR );
  SLIC_ASSERT( cellIdx >= 0 && cellIdx < m_num_cells );

  switch ( m_type )
  {
  case PARTICLE_MESH:
    cell[ 0 ] = cellIdx;
    break;
  case UNIFORM_MESH:
    SLIC_ERROR( "!!! NOT IMPLEMENTED YET !!!" );
    // TODO: implement this
    break;
  case STRUCTURED_MESH:
    SLIC_ERROR( "!!! NOT IMPLEMENTED YET !!!" );
    // TODO: implement this
    break;
  case RECTILINEAR_MESH:
    SLIC_ERROR( "!!! NOT IMPLEMENTED YET !!!" );
    // TODO: implement this
    break;
  case UNSTRUCTURED_MESH:
    SLIC_ERROR( "!!! NOT IMPLEMENTED YET !!!" );
    // TODO: implement this
    break;
  default:
    SLIC_ERROR( "Undefined MeshType!" );
  } // END switch

}

//------------------------------------------------------------------------------
CellType Mesh::getMeshCellType( IndexType cellIdx ) const
{
  CellType type = UNDEFINED_CELL;

  switch ( m_type )
  {
  case PARTICLE_MESH:
    type = VERTEX;
    break;
  case UNIFORM_MESH:       /* intentional fall-through */
  case STRUCTURED_MESH:    /* intentional fall-through */
  case RECTILINEAR_MESH:
    type = (getDimension()==3)? HEX : ( (getDimension()==2)? QUAD:SEGMENT );
    break;
  case UNSTRUCTURED_MESH:
    {
      // TODO: implement this
      SLIC_ERROR( "!!! NOT IMPLEMENTED YET !!!" );
    }
    break;
  default:
    SLIC_ERROR( "Undefined MeshType!" );
  }

  return type;
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

//------------------------------------------------------------------------------
Mesh* Mesh::getMesh( const sidre::Group* group, const std::string& topo )
{
  SLIC_ERROR_IF( group==AXOM_NULLPTR, "supplied group is null" );

  int mesh_type = UNDEFINED_MESH;
  int dimension = -1;
  blueprint::getMeshTypeAndDimension( mesh_type, dimension, group, topo );

  Mesh* m = AXOM_NULLPTR;
  switch ( mesh_type )
  {
  case STRUCTURED_MESH:
    SLIC_ERROR( "!!! NOT IMPLEMENTED YET !!!" );
    // TODO: implement this
    break;
  case RECTILINEAR_MESH:
    SLIC_ERROR( "!!! NOT IMPLEMENTED YET !!!" );
    // TODO: implement this
    break;
  case UNIFORM_MESH:
    SLIC_ERROR( "!!! NOT IMPLEMENTED YET !!!" );
    // TODO: implement this
    break;
  case UNSTRUCTURED_MESH:
    SLIC_ERROR( "!!! NOT IMPLEMENTED YET !!!" );
    // TODO: implement this
    break;
  case PARTICLE_MESH:
    SLIC_ERROR( "!!! NOT IMPLEMENTED YET !!!" );
    // TODO: implement this
    break;
  default:
    SLIC_ERROR( "undefined mesh_type [" << mesh_type << "]\n" );
  } // END switch

  SLIC_ASSERT( m != AXOM_NULLPTR );
  return ( m );
}

//------------------------------------------------------------------------------
Mesh* Mesh::getMesh( const sidre::Group* group )
{
  return ( Mesh::getMesh( group, "" ) );
}

} /* namespace mint */
} /* namespace axom */
