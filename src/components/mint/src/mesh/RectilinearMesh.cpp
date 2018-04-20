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
#include "mint/config.hpp"         /* For IndexType */
#include "mint/MeshTypes.hpp"          /* For MINT_STRUCTURED_RECTILINEAR_MESH */

namespace axom
{
namespace mint
{

//------------------------------------------------------------------------------
RectilinearMesh::RectilinearMesh( int dimension, int64 ext[6] ) :
  StructuredMesh( mint::RECTILINEAR_MESH, dimension, ext )
{
  IndexType ext_size[3];
  this->getExtentSize( ext_size );

  for ( int dim = 0 ; dim < m_ndims ; ++dim )
  {
    m_coordinates[ dim ] =
      new Array< double >( ext_size[ dim ], 1, ext_size[ dim ] );
    m_coordinates[ dim ]->setResizeRatio( 0.0 );
  }

  for (int dim = m_ndims ; dim < 3 ; ++dim )
  {
    m_coordinates[ dim ] = AXOM_NULLPTR;
  }

}

//------------------------------------------------------------------------------
RectilinearMesh::~RectilinearMesh()
{
  for ( int dim = 0 ; dim < 3 ; ++dim )
  {
    if ( m_coordinates[ dim ] != AXOM_NULLPTR )
    {
      delete m_coordinates[ dim ];
      m_coordinates[ dim ] = AXOM_NULLPTR;
    }
  }
}

} /* namespace mint */
} /* namespace axom */
