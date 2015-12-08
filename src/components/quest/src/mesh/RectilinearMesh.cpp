/*
 * $Id$
 */

/*!
 *******************************************************************************
 * \file RectilinearMesh.cxx
 *
 * \date Sep 26, 2015
 * \author George Zagaris (zagaris2@llnl.gov)
 *******************************************************************************
 */

#include "quest/RectilinearMesh.hpp"

#include "common/CommonTypes.hpp"

namespace meshtk {

RectilinearMesh::RectilinearMesh() :
        StructuredMesh(UNDEFINED_MESH,-1,ATK_NULLPTR),
        m_coordinates( ATK_NULLPTR )
{

}

//------------------------------------------------------------------------------
RectilinearMesh::RectilinearMesh( int dimension, int ext[6] ) :
       StructuredMesh( STRUCTURED_RECTILINEAR_MESH, dimension, ext )
{
  int ndims[3];
  this->getDimensions( ndims );
  m_coordinates = new MeshCoordinates( dimension, ndims );
}

//------------------------------------------------------------------------------
RectilinearMesh::RectilinearMesh( int dimension, int ext[6],
                                  int blockId, int partId ) :
 StructuredMesh( STRUCTURED_RECTILINEAR_MESH, dimension, ext, blockId, partId )
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

} /* namespace meshtk */
