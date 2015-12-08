/*
 * $Id$
 */

/*!
 *******************************************************************************
 * \file CurvilinearMesh.cxx
 *
 * \date Sep 20, 2015
 * \author George Zagaris (zagaris2@llnl.gov)
 *******************************************************************************
 */

#include "CurvilinearMesh.hpp"

#include <cassert> // for assert()
#include <cstddef> // for definition of ATK_NULLPTR

namespace meshtk {

CurvilinearMesh::CurvilinearMesh() :
        StructuredMesh( UNDEFINED_MESH, -1, ATK_NULLPTR ),
        m_coordinates( ATK_NULLPTR )
{

}

//------------------------------------------------------------------------------
CurvilinearMesh::CurvilinearMesh( int ndims, int ext[6] ) :
        StructuredMesh( STRUCTURED_CURVILINEAR_MESH, ndims, ext ),
        m_coordinates( new MeshCoordinates(ndims,m_extent->getNumNodes()) )
{

}

//------------------------------------------------------------------------------
CurvilinearMesh::CurvilinearMesh( int ndims, int ext[6],
                                  int blockId, int partId ) :
     StructuredMesh( STRUCTURED_CURVILINEAR_MESH, ndims, ext, blockId, partId),
     m_coordinates( new MeshCoordinates(ndims,m_extent->getNumNodes()) )
{

}

//------------------------------------------------------------------------------
CurvilinearMesh::~CurvilinearMesh()
{
  delete m_coordinates;
  m_coordinates = ATK_NULLPTR;
}

} /* namespace meshtk */
