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

} /* end anonymous namespace */

//------------------------------------------------------------------------------
// IMPLEMENTATION
//------------------------------------------------------------------------------
StructuredMesh::StructuredMesh( int meshType, int dimension, const int64* ext) :
  Mesh( dimension, meshType ),
  m_extent( new Extent(dimension, ext) )
{
  SLIC_ERROR_IF( !validStructuredMeshType( m_type ),
                 "invalid structured mesh type!" );
  initializeFields();
}

//------------------------------------------------------------------------------
StructuredMesh::StructuredMesh( int meshType, int dimension ) :
  Mesh( dimension, meshType ),
  m_extent( nullptr )
{
  SLIC_ERROR_IF( !validStructuredMeshType( m_type ),
                 "invalid structured mesh type!" );
}

//------------------------------------------------------------------------------
StructuredMesh::~StructuredMesh( )
{
  delete m_extent;
  m_extent = nullptr;
}

//------------------------------------------------------------------------------
void StructuredMesh::initializeFields()
{
  SLIC_ERROR_IF( m_extent==nullptr, "null extent for mesh!" );

  m_mesh_fields[ NODE_CENTERED ]->setResizeRatio( 1.0 );
  m_mesh_fields[ CELL_CENTERED ]->setResizeRatio( 1.0 );
  m_mesh_fields[ FACE_CENTERED ]->setResizeRatio( 1.0 );
  m_mesh_fields[ EDGE_CENTERED ]->setResizeRatio( 1.0 );

  m_mesh_fields[ NODE_CENTERED ]->resize( getNumberOfNodes() );
  m_mesh_fields[ CELL_CENTERED ]->resize( getNumberOfCells() );
}

#ifdef MINT_USE_SIDRE

//------------------------------------------------------------------------------
StructuredMesh::StructuredMesh( sidre::Group* group, const std::string& topo ) :
  Mesh( group, topo ),
  m_extent( nullptr )
{}

//------------------------------------------------------------------------------
StructuredMesh::StructuredMesh( int meshType, int dimension,
                                sidre::Group* group,
                                const std::string& topo,
                                const std::string& coordset ) :
  Mesh( dimension, meshType, group, topo, coordset ),
  m_extent( nullptr )
{}

#endif

}   /* end namespace mint */
}   /* end namespace axom */
