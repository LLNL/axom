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

#include "mint/StructuredMesh.hpp"
#include "mint/MeshType.hpp"
#include "slic/slic.hpp"

namespace axom
{
namespace mint
{

StructuredMesh::StructuredMesh() :
  Mesh(-1,MINT_UNDEFINED_MESH,-1,-1),
  m_extent( AXOM_NULLPTR )
{
// TODO Auto-generated constructor stub

}

//------------------------------------------------------------------------------
StructuredMesh::StructuredMesh( int meshType, int ndims, const int ext[ 6 ]) :
  Mesh( ndims, meshType, 0, 0 ),
  m_extent( new Extent< int >( ndims, ext ) )
{}

//------------------------------------------------------------------------------
StructuredMesh::StructuredMesh( int meshType, int ndims, const int ext[ 6 ],
                                int blockId, int partId) :
  Mesh( ndims, meshType, blockId, partId ),
  m_extent( new Extent< int >( ndims, ext ) )

{}

//------------------------------------------------------------------------------
StructuredMesh::~StructuredMesh()
{
  delete m_extent;
  m_extent = AXOM_NULLPTR;
}

} /* namespace mint */
} /* namespace axom */
