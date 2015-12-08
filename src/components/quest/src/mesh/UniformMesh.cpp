/*
 * $Id$
 */

/*!
 *******************************************************************************
 * \file UniformMesh.cxx
 *
 * \date Sep 26, 2015
 * \author George Zagaris (zagaris2@llnl.gov)
 *******************************************************************************
 */

#include "quest/UniformMesh.hpp"
#include "quest/MeshType.hpp"

#include "common/CommonTypes.hpp"

#include <algorithm> // for std::fill()

namespace meshtk
{

UniformMesh::UniformMesh() :
    StructuredMesh(UNDEFINED_MESH,-1,ATK_NULLPTR)
{
  std::fill(m_origin,m_origin+3,0.0);
  std::fill(m_h,m_h+3,1.0);
}

//------------------------------------------------------------------------------
UniformMesh::UniformMesh( int dimension,
                          const double origin[3],
                          const double h[3],
                          const int ext[6] ) :
     StructuredMesh( STRUCTURED_UNIFORM_MESH,dimension,const_cast<int*>(ext) )
{
  std::fill( m_origin, m_origin+3, 0.0 );
  std::fill( m_h, m_h+3 ,1.0 );

  memcpy( m_origin, origin, dimension*sizeof(double) );
  memcpy( m_h, h, dimension*sizeof(double) );
}

//------------------------------------------------------------------------------
UniformMesh::UniformMesh( int dimension,
                          const double origin[3],
                          const double h[3],
                          const int ext[6],
                          int blockId,
                          int partitionId ) :
    StructuredMesh( STRUCTURED_UNIFORM_MESH,dimension,const_cast<int*>(ext),
                    blockId, partitionId )
{
  std::fill( m_origin, m_origin+3, 0.0 );
  std::fill( m_h, m_h+3 ,1.0 );

  memcpy( m_origin, origin, dimension*sizeof(double) );
  memcpy( m_h, h, dimension*sizeof(double) );
}

//------------------------------------------------------------------------------
UniformMesh::~UniformMesh()
{
  // TODO Auto-generated destructor stub
}


} /* namespace meshtk */
