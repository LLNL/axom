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

namespace axom {
namespace mint {

//------------------------------------------------------------------------------
RectilinearMesh::RectilinearMesh( int dimension, int ext[6] ) :
  StructuredMesh( MINT_STRUCTURED_RECTILINEAR_MESH, dimension, ext )
{
  int ext_size[3];
  this->getExtentSize( ext_size );

  for ( int dim = 0; dim < 3; ++dim ) {
    m_coordinates[ dim ].setCapacity(0);
    m_coordinates[ dim ].setResizeRatio(0.0);
  }

  for ( int dim = 0; dim < dimension; ++dim ) {
    m_coordinates[ dim ].setSize( ext_size[ dim ] );
  }
}

//------------------------------------------------------------------------------
RectilinearMesh::RectilinearMesh( int dimension, int ext[6],
                                  int blockId, int partId ) :
  StructuredMesh( MINT_STRUCTURED_RECTILINEAR_MESH, dimension, ext, blockId,
                  partId ) 
{
  int ext_size[3];
  this->getExtentSize( ext_size );

  for ( int dim = 0; dim < 3; ++dim ) {
    m_coordinates[ dim ].setCapacity(0);
    m_coordinates[ dim ].setResizeRatio(0.0);
  }

  for ( int dim = 0; dim < dimension; ++dim ) {
    m_coordinates[ dim ].setSize( ext_size[ dim ] );
  }
}

//------------------------------------------------------------------------------
const double* RectilinearMesh::getCoordinateArray( int dim ) const {
  SLIC_ASSERT( dim >= 0 && dim < this->getDimension() );
  return m_coordinates[ dim ].getData();
}

} /* namespace mint */
} /* namespace axom */
