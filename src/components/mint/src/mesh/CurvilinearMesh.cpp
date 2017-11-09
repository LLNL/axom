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

//------------------------------------------------------------------------------
CurvilinearMesh::CurvilinearMesh( int ndims, const globalIndex ext[6] ):
  StructuredMesh( MINT_STRUCTURED_CURVILINEAR_MESH, ndims, ext ),
  m_coordinates( ndims, m_extent.getNumNodes(), 0.0 )
{ m_coordinates.setSize( m_extent.getNumNodes() ); }

//------------------------------------------------------------------------------
CurvilinearMesh::CurvilinearMesh( int ndims, const globalIndex ext[6], 
                                  int blockId, int partId ):
  StructuredMesh( MINT_STRUCTURED_CURVILINEAR_MESH, ndims, ext, blockId,
                  partId),
  m_coordinates( ndims, m_extent.getNumNodes(), 0.0  )
{ m_coordinates.setSize( m_extent.getNumNodes() ); }

//------------------------------------------------------------------------------
CurvilinearMesh::CurvilinearMesh():
  StructuredMesh( MINT_UNDEFINED_MESH, -1, AXOM_NULLPTR, -1, -1 ),
  m_coordinates( 0, 0, 0 )
{}


} /* namespace mint */
} /* namespace axom */
