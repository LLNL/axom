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

#include "mint/RectilinearMesh.hpp"

#include "axom/Types.hpp"

namespace axom
{
namespace mint
{

RectilinearMesh::RectilinearMesh() :
  StructuredMesh(MINT_UNDEFINED_MESH,-1,AXOM_NULLPTR),
  m_coordinates( AXOM_NULLPTR )
{}

//------------------------------------------------------------------------------
RectilinearMesh::RectilinearMesh( int dimension, int ext[6] ) :
  StructuredMesh( MINT_STRUCTURED_RECTILINEAR_MESH, dimension, ext )
{
  int ext_size[3];
  this->getExtentSize( ext_size );
  m_coordinates = new MeshCoordinates( dimension, ext_size );
}

//------------------------------------------------------------------------------
RectilinearMesh::RectilinearMesh( int dimension, int ext[6],
                                  int blockId, int partId ) :
  StructuredMesh( MINT_STRUCTURED_RECTILINEAR_MESH, dimension, ext, blockId,
                  partId )
{
  int ext_size[3];
  this->getExtentSize( ext_size );
  m_coordinates = new MeshCoordinates( dimension, ext_size );
}

//------------------------------------------------------------------------------
RectilinearMesh::~RectilinearMesh()
{
  delete m_coordinates;
  m_coordinates = AXOM_NULLPTR;
}

} /* namespace mint */
} /* namespace axom */
