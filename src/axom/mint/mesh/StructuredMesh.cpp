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
#include "axom/mint/mesh/StructuredMesh.hpp"
#include "axom/mint/mesh/MeshTypes.hpp"
#include "axom/mint/mesh/blueprint.hpp"

#include <cstring>             /* for memcpy() */
#include <limits>

namespace axom
{
namespace mint
{



//------------------------------------------------------------------------------
// HELPER METHODS
//------------------------------------------------------------------------------
namespace
{

bool validStructuredMeshType( int type )
{
  return ( (type==STRUCTURED_CURVILINEAR_MESH) ||
           (type==STRUCTURED_RECTILINEAR_MESH) ||
           (type==STRUCTURED_UNIFORM_MESH)
           );
}

inline int dim ( const IndexType& AXOM_NOT_USED( Ni ),
                 const IndexType& Nj,
                 const IndexType& Nk  )
{
  return ( (Nk >= 1) ? 3 : ( (Nj >= 1) ? 2 : 1 ) );
}

} /* end anonymous namespace */

//------------------------------------------------------------------------------
// IMPLEMENTATION
//------------------------------------------------------------------------------

void StructuredMesh::setGlobalNodeExtent( int ndims, const int64* extent )
{
  SLIC_ASSERT( 0 < ndims && ndims <= 3 );
  SLIC_ASSERT( extent != AXOM_NULLPTR );

  std::memset( m_global_node_extent, 0, sizeof( m_global_node_extent ) );

  int64 numNodes = 1;
  for ( int dim = 0; dim < ndims; ++dim )
  {
    const int64 min = extent[ 2 * dim ];
    const int64 max = extent[ 2 * dim + 1 ];
    SLIC_ASSERT( min <= max );
    m_global_node_extent[ 2 * dim ] = min;
    m_global_node_extent[ 2 * dim + 1 ] = max;
    numNodes *= max - min;
  }

  SLIC_ASSERT( numNodes == getNumberOfNodes() );

#ifdef MINT_USE_SIDRE
  if ( hasSidreGroup() )
  {
    blueprint::setGlobalNodeExtent( getCoordsetGroup(), m_global_node_extent );
  }
#endif
}

StructuredMesh::StructuredMesh( int meshType, int dimension,
                                const IndexType* node_ext) :
  Mesh( dimension, meshType )
{
  SLIC_ERROR_IF( !validStructuredMeshType( m_type ),
                 "invalid structured mesh type!" );
  SLIC_ERROR_IF( node_ext == AXOM_NULLPTR, "invalid extent." );

  memcpy( m_node_extent, node_ext, m_ndims * sizeof( IndexType ) );

  structuredInit();
}

//------------------------------------------------------------------------------
StructuredMesh::StructuredMesh( int meshType, IndexType Ni, IndexType Nj, 
                                IndexType Nk ) :
  Mesh( dim( Ni, Nj, Nk ), meshType )
{
  SLIC_ERROR_IF( !validStructuredMeshType( m_type ),
                 "invalid structured mesh type!" );

  SLIC_ERROR_IF( Ni <= 1, "Ni must be greater or equal to 2" );
  m_node_extent[0] = Ni;
  if ( m_ndims > 1 )
  {
    SLIC_ERROR_IF( Nj <= 1, "Nj must be greater or equal to 2" );
    m_node_extent[1] = Nj;
  }
  if ( m_ndims > 2 )
  {
    SLIC_ERROR_IF( Nk <= 1, "Nk must be greater or equal to 2" );
    m_node_extent[2] = Nk;
  }

  structuredInit();
}

#ifdef MINT_USE_SIDRE

//------------------------------------------------------------------------------
StructuredMesh::StructuredMesh( sidre::Group* group, const std::string& topo ) :
  Mesh( group, topo )
{
  SLIC_ERROR_IF( !validStructuredMeshType( m_type ),
                 "invalid structured mesh type!" );

  blueprint::getStructuredMesh( m_ndims, m_node_extent, m_global_node_extent,
                                getCoordsetGroup() );
  structuredInit();
}

//------------------------------------------------------------------------------
StructuredMesh::StructuredMesh( int meshType, int dimension, 
                                const IndexType* node_ext, sidre::Group* group,
                                const std::string& topo,
                                const std::string& coordset ) :
  Mesh( dimension, meshType, group, topo, coordset )
{
  SLIC_ERROR_IF( !validStructuredMeshType( m_type ),
                 "invalid structured mesh type!" );
  SLIC_ERROR_IF( node_ext == AXOM_NULLPTR, "invalid extent." );

  memcpy( m_node_extent, node_ext, m_ndims * sizeof( IndexType ) );

  std::string topo_type;
  if ( meshType == STRUCTURED_UNIFORM_MESH )
  {
    topo_type = "uniform";
  }
  else if ( meshType == STRUCTURED_RECTILINEAR_MESH )
  {
    topo_type = "rectilinear";
  }
  else
  {
    topo_type = "structured";
  }

  blueprint::initializeTopologyGroup( m_group, m_topology, m_coordset,
                                      topo_type );
  SLIC_ERROR_IF( !blueprint::isValidTopologyGroup( getTopologyGroup() ),
                 "invalid topology group!" );

  blueprint::setStructuredMesh( m_ndims, m_node_extent, m_global_node_extent,
                                getCoordsetGroup() );
  structuredInit();
}

//------------------------------------------------------------------------------
StructuredMesh::StructuredMesh( int meshType, IndexType Ni, IndexType Nj,
                                IndexType Nk, sidre::Group* group,
                                const std::string& topo,
                                const std::string& coordset ) :
  Mesh( dim( Ni, Nj, Nk ), meshType, group, topo, coordset )
{
  SLIC_ERROR_IF( !validStructuredMeshType( m_type ),
                 "invalid structured mesh type!" );

  SLIC_ERROR_IF( Ni <= 1, "Ni must be greater or equal to 2" );
  m_node_extent[0] = Ni;
  if ( m_ndims > 1 )
  {
    SLIC_ERROR_IF( Nj <= 1, "Nj must be greater or equal to 2" );
    m_node_extent[1] = Nj;
  }
  if ( m_ndims > 2 )
  {
    SLIC_ERROR_IF( Nk <= 1, "Nk must be greater or equal to 2" );
    m_node_extent[2] = Nk;
  }

  std::string topo_type;
  if ( meshType == STRUCTURED_UNIFORM_MESH )
  {
    topo_type = "uniform";
  }
  else if ( meshType == STRUCTURED_RECTILINEAR_MESH )
  {
    topo_type = "rectilinear";
  }
  else
  {
    topo_type = "structured";
  }
  
  blueprint::initializeTopologyGroup( m_group, m_topology, m_coordset,
                                      topo_type );
  SLIC_ERROR_IF( !blueprint::isValidTopologyGroup( getTopologyGroup() ),
                 "invalid topology group!" );

  blueprint::setStructuredMesh( m_ndims, m_node_extent, m_global_node_extent,
                                getCoordsetGroup() );
  structuredInit();
}

#endif

//------------------------------------------------------------------------------
void StructuredMesh::structuredInit()
{
  for ( int dim = 0; dim < m_ndims; ++dim )
  {
    SLIC_ERROR_IF( m_node_extent[ dim ] < 2, "invalid extent" );
  }

  /* Initialize the node meta data. */
  m_node_jp = ( m_ndims > 1 ) ? getNodeExtent( 0 ) : 
                                          std::numeric_limits<IndexType>::max();
  m_node_kp = ( m_ndims > 2 ) ? m_node_jp * getNodeExtent( 1 ) : 
                                          std::numeric_limits<IndexType>::max();

  /* Initialize the cell meta data */
  for ( int dim = 0; dim < m_ndims; ++dim )
  {
    m_cell_extent[ dim ] = getNodeExtent( dim ) - 1;
  }

  m_cell_jp = ( m_ndims > 1 ) ? getCellExtent( 0 ) : 
                                          std::numeric_limits<IndexType>::max();
  m_cell_kp = ( m_ndims > 2 ) ? m_cell_jp * getCellExtent( 1 ) :
                                          std::numeric_limits<IndexType>::max();

  /* Build the cell to node offsets. */
  m_cell_node_offsets[ 0 ] = 0;
  m_cell_node_offsets[ 1 ] = 1;
  m_cell_node_offsets[ 2 ] = 1 + nodeJp();
  m_cell_node_offsets[ 3 ] = nodeJp();

  m_cell_node_offsets[ 4 ] = nodeKp();
  m_cell_node_offsets[ 5 ] = 1 + nodeKp();
  m_cell_node_offsets[ 6 ] = 1 + nodeJp() + nodeKp();
  m_cell_node_offsets[ 7 ] = nodeJp() + nodeKp();

  /* Initialize the face meta data */
  if ( m_ndims == 2 )
  {
    m_total_faces[ 0 ] = getNodeExtent( 0 ) * getCellExtent( 1 );
    m_total_faces[ 1 ] = getCellExtent( 0 ) * getNodeExtent( 1 );
  }
  else if ( m_ndims == 3 )
  {
    m_total_faces[ 0 ] = getNodeExtent( 0 ) * getCellExtent( 1 ) * getCellExtent( 2 );
    m_total_faces[ 1 ] = getCellExtent( 0 ) * getNodeExtent( 1 ) * getCellExtent( 2 );
    m_total_faces[ 2 ] = getCellExtent( 0 ) * getCellExtent( 1 ) * getNodeExtent( 2 );
  }

  /* Initialize the edge meta data */
  if ( m_ndims < 3 )
  {
    m_num_edges = 0;
  }
  else
  {
    m_num_edges = getCellExtent( 0 ) * getNodeExtent( 1 ) * getNodeExtent( 2 )
                + getCellExtent( 0 ) * getNodeExtent( 2 ) * getNodeExtent( 1 )
                + getCellExtent( 1 ) * getNodeExtent( 2 ) * getNodeExtent( 0 );
  }

  /* Initialize the fields */
  m_mesh_fields[ NODE_CENTERED ]->setResizeRatio( 0.0 );
  m_mesh_fields[ CELL_CENTERED ]->setResizeRatio( 0.0 );
  m_mesh_fields[ FACE_CENTERED ]->setResizeRatio( 0.0 );
  m_mesh_fields[ EDGE_CENTERED ]->setResizeRatio( 0.0 );

  m_mesh_fields[ NODE_CENTERED ]->resize( getNumberOfNodes() );
  m_mesh_fields[ CELL_CENTERED ]->resize( getNumberOfCells() );
  m_mesh_fields[ FACE_CENTERED ]->resize( getNumberOfFaces() );
  m_mesh_fields[ EDGE_CENTERED ]->resize( getNumberOfEdges() );
}

}   /* end namespace mint */
}   /* end namespace axom */
