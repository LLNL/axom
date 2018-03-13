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
#include "mint/FieldAssociation.hpp"

#ifdef MINT_USE_SIDRE
#include "sidre/sidre.hpp"
#endif

// C/C++ includes
#include <cstddef>

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
  m_cell_data( CELL_CENTERED ),
  m_face_data( FACE_CENTERED ),
  m_edge_data( EDGE_CENTERED ),
  m_node_data( NODE_CENTERED ),
#ifdef MINT_USE_SIDRE
  m_group( AXOM_NULLPTR ),
#endif
  m_num_cells( AXOM_NULLPTR ),
  m_cell_capacity( AXOM_NULLPTR ),
  m_cell_resize_ratio( AXOM_NULLPTR ),
  m_num_faces( AXOM_NULLPTR ),
  m_face_capacity( AXOM_NULLPTR ),
  m_face_resize_ratio( AXOM_NULLPTR ),
  m_num_edges( AXOM_NULLPTR ),
  m_edge_capacity( AXOM_NULLPTR ),
  m_edge_resize_ratio( AXOM_NULLPTR ),
  m_num_nodes( AXOM_NULLPTR ),
  m_node_capacity( AXOM_NULLPTR ),
  m_node_resize_ratio( AXOM_NULLPTR )
{
  SLIC_ERROR_IF( m_ndims < 0 || m_ndims > 3, "invalid dimension" );

  m_num_cells = new IndexType(0);
  m_cell_capacity = new IndexType(0);
  m_cell_resize_ratio = new double(0.0);

  m_num_faces = new IndexType(0);
  m_face_capacity = new IndexType(0);
  m_face_resize_ratio = new double(0.0);

  m_num_edges = new IndexType(0);
  m_edge_capacity = new IndexType(0);
  m_edge_resize_ratio = new double(0.0);

  m_num_nodes = new IndexType(0);
  m_node_capacity = new IndexType(0);
  m_node_resize_ratio = new double(0.0);
}

#ifdef MINT_USE_SIDRE
//------------------------------------------------------------------------------
Mesh::Mesh( sidre::Group* group ) :
  m_ndims( 0 ),
  m_type( 0 ),
  m_block_idx( 0 ),
  m_part_idx( 0 ),
  m_cell_data( CELL_CENTERED ),
  m_face_data( FACE_CENTERED ),
  m_edge_data( EDGE_CENTERED ),
  m_node_data( NODE_CENTERED ),
  m_group( group ),
  m_num_cells( AXOM_NULLPTR ),
  m_cell_capacity( AXOM_NULLPTR ),
  m_cell_resize_ratio( AXOM_NULLPTR ),
  m_num_faces( AXOM_NULLPTR ),
  m_face_capacity( AXOM_NULLPTR ),
  m_face_resize_ratio( AXOM_NULLPTR ),
  m_num_edges( AXOM_NULLPTR ),
  m_edge_capacity( AXOM_NULLPTR ),
  m_edge_resize_ratio( AXOM_NULLPTR ),
  m_num_nodes( AXOM_NULLPTR ),
  m_node_capacity( AXOM_NULLPTR ),
  m_node_resize_ratio( AXOM_NULLPTR )
{
  SLIC_ERROR_IF( m_ndims < 0 || m_ndims > 3, "");
  SLIC_ERROR_IF( m_group == AXOM_NULLPTR, "");
  SLIC_ERROR_IF( m_group->getNumGroups() == 0, "");
  SLIC_ERROR_IF( m_group->getNumViews() == 0, "");
  SLIC_ERROR_IF( !m_group->hasChildView( "ndims" ), "");
  SLIC_ERROR_IF( !m_group->hasChildView( "type" ), "");
  SLIC_ERROR_IF( !m_group->hasChildView( "block_idx" ), "");
  SLIC_ERROR_IF( !m_group->hasChildView( "part_idx" ), "");
  SLIC_ERROR_IF( !m_group->hasChildView( "num_cells" ), "");
  SLIC_ERROR_IF( !m_group->hasChildView( "cell_resize_ratio" ), "");
  SLIC_ERROR_IF( !m_group->hasChildView( "num_faces" ), "");
  SLIC_ERROR_IF( !m_group->hasChildView( "face_resize_ratio" ), "");
  SLIC_ERROR_IF( !m_group->hasChildView( "num_edges" ), "");
  SLIC_ERROR_IF( !m_group->hasChildView( "edge_resize_ratio" ), "");
  SLIC_ERROR_IF( !m_group->hasChildView( "num_nodes" ), "");
  SLIC_ERROR_IF( !m_group->hasChildView( "node_resize_ratio" ), "");

  sidre::View* view = m_group->getView( "ndims" );
  SLIC_ERROR_IF( !view->isScalar(), "" );
  m_ndims = view->getData();

  view = m_group->getView( "type" );
  SLIC_ERROR_IF( !view->isScalar(), "" );
  m_type = view->getData();

  view = m_group->getView( "block_idx" );
  SLIC_ERROR_IF( !view->isScalar(), "" );
  m_block_idx = view->getData();

  view = m_group->getView( "part_idx" );
  SLIC_ERROR_IF( !view->isScalar(), "" );
  m_part_idx = view->getData();

  view = m_group->getView( "num_cells" );
  SLIC_ERROR_IF( !view->isScalar(), "" );
  m_num_cells = static_cast< IndexType* >( view->getVoidPtr() );

  view = m_group->getView( "cell_resize_ratio" );
  SLIC_ERROR_IF( !view->isScalar(), "" );
  m_cell_resize_ratio = static_cast< double* >( view->getVoidPtr() );

  view = m_group->getView( "num_faces" );
  SLIC_ERROR_IF( !view->isScalar(), "" );
  m_num_faces = static_cast< IndexType* >( view->getVoidPtr() );

  view = m_group->getView( "face_resize_ratio" );
  SLIC_ERROR_IF( !view->isScalar(), "" );
  m_face_resize_ratio = static_cast< double* >( view->getVoidPtr() );

  view = m_group->getView( "num_edges" );
  SLIC_ERROR_IF( !view->isScalar(), "" );
  m_num_edges = static_cast< IndexType* >( view->getVoidPtr() );

  view = m_group->getView( "edge_resize_ratio" );
  SLIC_ERROR_IF( !view->isScalar(), "" );
  m_edge_resize_ratio = static_cast< double* >( view->getVoidPtr() );

  view = m_group->getView( "num_nodes" );
  SLIC_ERROR_IF( !view->isScalar(), "" );
  m_num_nodes = static_cast< IndexType* >( view->getVoidPtr() );

  view = m_group->getView( "node_resize_ratio" );
  SLIC_ERROR_IF( !view->isScalar(), "" );
  m_node_resize_ratio = static_cast< double* >( view->getVoidPtr() );
}

//------------------------------------------------------------------------------
Mesh::Mesh( sidre::Group* group, int ndims, int type, int blockId,
            int partId ) :
  m_ndims( ndims ),
  m_type( type ),
  m_block_idx( blockId ),
  m_part_idx( partId ),
  m_cell_data( CELL_CENTERED ),
  m_face_data( FACE_CENTERED ),
  m_edge_data( EDGE_CENTERED ),
  m_node_data( NODE_CENTERED ),
  m_group( group ),
  m_num_cells( AXOM_NULLPTR ),
  m_cell_capacity( AXOM_NULLPTR ),
  m_cell_resize_ratio( AXOM_NULLPTR ),
  m_num_faces( AXOM_NULLPTR ),
  m_face_capacity( AXOM_NULLPTR ),
  m_face_resize_ratio( AXOM_NULLPTR ),
  m_num_edges( AXOM_NULLPTR ),
  m_edge_capacity( AXOM_NULLPTR ),
  m_edge_resize_ratio( AXOM_NULLPTR ),
  m_num_nodes( AXOM_NULLPTR ),
  m_node_capacity( AXOM_NULLPTR ),
  m_node_resize_ratio( AXOM_NULLPTR )
{
  SLIC_ERROR_IF( m_ndims < 0 || m_ndims > 3, "" );
  SLIC_ERROR_IF( m_group == AXOM_NULLPTR, "" );
  SLIC_ERROR_IF( m_group->getNumGroups() != 0, "" );
  SLIC_ERROR_IF( m_group->getNumViews() != 0, "" );

  m_group->createView( "ndims" )->setScalar( m_ndims );
  m_group->createView( "type" )->setScalar( m_type );
  m_group->createView( "block_idx" )->setScalar( m_block_idx );
  m_group->createView( "part_idx" )->setScalar( m_part_idx );

  IndexType zero = 0;
  double zero_f = 0.0;
  m_num_cells = static_cast< IndexType* >(
    m_group->createView( "num_cells" )
    ->setScalar( zero )->getVoidPtr() );
  m_cell_resize_ratio = static_cast< double* >(
    m_group->createView( "cell_resize_ratio")
    ->setScalar( zero_f )->getVoidPtr() );
  m_num_faces = static_cast< IndexType* >(
    m_group->createView( "num_faces" )
    ->setScalar( zero )->getVoidPtr() );
  m_face_resize_ratio = static_cast< double* >(
    m_group->createView( "face_resize_ratio")
    ->setScalar( zero_f )->getVoidPtr() );
  m_num_edges = static_cast< IndexType* >(
    m_group->createView( "num_edges" )
    ->setScalar( zero )->getVoidPtr() );
  m_edge_resize_ratio = static_cast< double* >(
    m_group->createView( "edge_resize_ratio")
    ->setScalar( zero_f )->getVoidPtr() );
  m_num_nodes = static_cast< IndexType* >(
    m_group->createView( "num_nodes" )
    ->setScalar( zero )->getVoidPtr() );
  m_node_resize_ratio = static_cast< double* >(
    m_group->createView( "node_resize_ratio")
    ->setScalar( zero_f )->getVoidPtr() );
}
#endif

//------------------------------------------------------------------------------
Mesh::~Mesh()
{
#ifdef MINT_USE_SIDRE
  if ( m_group != AXOM_NULLPTR )
  {
    m_group = AXOM_NULLPTR;
    m_num_cells = AXOM_NULLPTR;
    m_cell_capacity = AXOM_NULLPTR;
    m_cell_resize_ratio = AXOM_NULLPTR;
    m_num_faces = AXOM_NULLPTR;
    m_face_capacity = AXOM_NULLPTR;
    m_face_resize_ratio = AXOM_NULLPTR;
    m_num_edges = AXOM_NULLPTR;
    m_edge_capacity = AXOM_NULLPTR;
    m_edge_resize_ratio = AXOM_NULLPTR;
    m_num_nodes = AXOM_NULLPTR;
    m_node_capacity = AXOM_NULLPTR;
    m_node_resize_ratio = AXOM_NULLPTR;
    return;
  }
#endif

  delete m_num_cells;
  delete m_cell_capacity;
  delete m_cell_resize_ratio;

  delete m_num_faces;
  delete m_face_capacity;
  delete m_face_resize_ratio;

  delete m_num_edges;
  delete m_edge_capacity;
  delete m_edge_resize_ratio;

  delete m_num_nodes;
  delete m_node_capacity;
  delete m_node_resize_ratio;
}



} /* namespace mint */
} /* namespace axom */
