/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */

#include "mint/RectilinearMesh.hpp"

#include "common/CommonTypes.hpp"

namespace axom {
namespace mint {

RectilinearMesh::RectilinearMesh() :
        StructuredMesh(MINT_UNDEFINED_MESH,-1,ATK_NULLPTR),
        m_coordinates( ATK_NULLPTR )
{

}

//------------------------------------------------------------------------------
RectilinearMesh::RectilinearMesh( int dimension, int ext[6] ) :
       StructuredMesh( MINT_STRUCTURED_RECTILINEAR_MESH, dimension, ext )
{
  int ndims[3];
  this->getDimensions( ndims );
  m_coordinates = new MeshCoordinates( dimension, ndims );
}

//------------------------------------------------------------------------------
RectilinearMesh::RectilinearMesh( int dimension, int ext[6],
                                  int blockId, int partId ) :
 StructuredMesh( MINT_STRUCTURED_RECTILINEAR_MESH, dimension, ext, blockId, partId )
{
  int ndims[3];
  this->getDimensions( ndims );
  m_coordinates = new MeshCoordinates( dimension, ndims );
}

//------------------------------------------------------------------------------
RectilinearMesh::~RectilinearMesh()
{
  delete m_coordinates;
  m_coordinates = ATK_NULLPTR;
}

} /* namespace mint */
} /* namespace axom */
