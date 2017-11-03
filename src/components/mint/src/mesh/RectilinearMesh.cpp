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
#include "mint/DataTypes.hpp"         /* For localIndex */
#include "mint/MeshType.hpp"          /* For MINT_STRUCTURED_RECTILINEAR_MESH */

#include "axom/Types.hpp"             /* For AXOM_NULLPTR */

namespace axom {
namespace mint {

//------------------------------------------------------------------------------
RectilinearMesh::RectilinearMesh( int dimension, globalIndex ext[6] ):
  StructuredMesh( MINT_STRUCTURED_RECTILINEAR_MESH, dimension, ext )
{
  localIndex ext_size[3];
  this->getExtentSize( ext_size );

  for ( int dim = 0; dim < m_ndims; ++dim ) {
    m_coordinates[ dim ].setSize( ext_size[ dim ] );
    m_coordinates[ dim ].setResizeRatio( 0.0 );
  }

  for (int dim = m_ndims; dim < 3; ++dim ) {
    m_coordinates[ dim ].setCapacity( 0 );
    m_coordinates[ dim ].setResizeRatio( 0.0 );
  }
}

//------------------------------------------------------------------------------
RectilinearMesh::RectilinearMesh( int dimension, globalIndex ext[6],
                                  int blockId, int partId ):
  StructuredMesh( MINT_STRUCTURED_RECTILINEAR_MESH, dimension, ext, blockId,
                  partId ) 
{
  localIndex ext_size[3];
  this->getExtentSize( ext_size );

  for ( int dim = 0; dim < m_ndims; ++dim ) {
    m_coordinates[ dim ].setSize( ext_size[ dim ] );
    m_coordinates[ dim ].setResizeRatio( 0.0 );
  }

  for (int dim = m_ndims; dim < 3; ++dim ) {
    m_coordinates[ dim ].setCapacity( 0 );
    m_coordinates[ dim ].setResizeRatio( 0.0 );
  }
}

RectilinearMesh::RectilinearMesh():
  StructuredMesh( MINT_UNDEFINED_MESH, -1, AXOM_NULLPTR )
{}

} /* namespace mint */
} /* namespace axom */
