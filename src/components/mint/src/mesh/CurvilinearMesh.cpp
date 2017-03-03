/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */

#include "CurvilinearMesh.hpp"
#include "MeshType.hpp"

// axom includes
#include "slic/slic.hpp"

// C/C++ includes
#include <cstddef> // for definition of AXOM_NULLPTR

namespace axom {
namespace mint {

CurvilinearMesh::CurvilinearMesh():
  StructuredMesh( MINT_UNDEFINED_MESH, -1, AXOM_NULLPTR ),
  m_coordinates( AXOM_NULLPTR )
{}

//------------------------------------------------------------------------------
CurvilinearMesh::CurvilinearMesh( int ndims, int ext[6] ):
  StructuredMesh( MINT_STRUCTURED_CURVILINEAR_MESH, ndims, ext ),
  m_coordinates( new MeshCoordinates(ndims,m_extent->getNumNodes()) )
{}

//------------------------------------------------------------------------------
CurvilinearMesh::CurvilinearMesh( int ndims, int ext[6],
                                  int blockId, int partId ):
  StructuredMesh( MINT_STRUCTURED_CURVILINEAR_MESH, ndims, ext, blockId,
                  partId),
  m_coordinates( new MeshCoordinates(ndims,m_extent->getNumNodes()) )
{}

//------------------------------------------------------------------------------
CurvilinearMesh::~CurvilinearMesh()
{
  delete m_coordinates;
  m_coordinates = AXOM_NULLPTR;
}

} /* namespace mint */
} /* namespace axom */
