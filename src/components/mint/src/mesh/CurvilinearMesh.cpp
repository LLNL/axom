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

#include "CurvilinearMesh.hpp"
#include "MeshType.hpp"

// axom includes
#include "slic/slic.hpp"

// C/C++ includes
#include <cstddef> // for definition of AXOM_NULLPTR

namespace axom
{
namespace mint
{

CurvilinearMesh::CurvilinearMesh() :
  StructuredMesh( MINT_UNDEFINED_MESH, -1, AXOM_NULLPTR ),
  m_coordinates( AXOM_NULLPTR )
{}

//------------------------------------------------------------------------------
CurvilinearMesh::CurvilinearMesh( int ndims, int ext[6] ) :
  StructuredMesh( MINT_STRUCTURED_CURVILINEAR_MESH, ndims, ext ),
  m_coordinates( new MeshCoordinates( ndims, m_extent->getNumNodes() ) )
{}

//------------------------------------------------------------------------------
CurvilinearMesh::CurvilinearMesh( int ndims, int ext[6], int blockId, 
                                  int partId ):
  StructuredMesh( MINT_STRUCTURED_CURVILINEAR_MESH, ndims, ext, blockId,
                  partId),
  m_coordinates( new MeshCoordinates( ndims, m_extent->getNumNodes() ) )
{}

//------------------------------------------------------------------------------
CurvilinearMesh::CurvilinearMesh( int ndims, int ext[6], int nodeCapacity ):
  StructuredMesh( MINT_STRUCTURED_CURVILINEAR_MESH, ndims, ext ),
  m_coordinates( new MeshCoordinates( ndims, nodeCapacity ) )
{
  SLIC_ERROR_IF( nodeCapacity < m_extent->getNumNodes(), "Node capacity " <<
                    "isn't large enough to hold all the nodes in the extent." );
}

//------------------------------------------------------------------------------
CurvilinearMesh::CurvilinearMesh( int ndims, int ext[6], int nodeCapacity,
                                  int blockId, int partId ):
  StructuredMesh( MINT_STRUCTURED_CURVILINEAR_MESH, ndims, ext, blockId,
                  partId),
  m_coordinates( new MeshCoordinates( ndims, nodeCapacity ) )
{
  SLIC_ERROR_IF( nodeCapacity < m_extent->getNumNodes(), "Node capacity " <<
                    "isn't large enough to hold all the nodes in the extent." );
}

//------------------------------------------------------------------------------
CurvilinearMesh::~CurvilinearMesh()
{
  delete m_coordinates;
  m_coordinates = AXOM_NULLPTR;
}

} /* namespace mint */
} /* namespace axom */
