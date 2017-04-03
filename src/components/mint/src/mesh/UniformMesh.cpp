/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */

#include "mint/UniformMesh.hpp"
#include "mint/MeshType.hpp"

#include "axom/Types.hpp"

#include <algorithm> // for std::fill()

namespace axom {
namespace mint {

UniformMesh::UniformMesh():
  StructuredMesh( MINT_UNDEFINED_MESH,-1,AXOM_NULLPTR)
{
  std:: fill( m_origin, m_origin+3, 0.0);
  std:: fill( m_h,      m_h+3,      1.0);
}

//------------------------------------------------------------------------------
UniformMesh::UniformMesh( int dimension,
                          const double origin[3],
                          const double h[3],
                          const int ext[6] ):
  StructuredMesh( MINT_STRUCTURED_UNIFORM_MESH,dimension,
                  const_cast< int* >(ext) )
{
  std:: fill( m_origin, m_origin+3, 0.0 );
  std:: fill( m_h,      m_h+3 ,     1.0 );

  memcpy( m_origin, origin, dimension*sizeof(double) );
  memcpy( m_h,      h,      dimension*sizeof(double) );
}

//------------------------------------------------------------------------------
UniformMesh::UniformMesh( int dimension,
                          const double origin[3],
                          const double h[3],
                          const int ext[6],
                          int blockId,
                          int partitionId ):
  StructuredMesh( MINT_STRUCTURED_UNIFORM_MESH,dimension,
                  const_cast< int* >(ext),
                  blockId, partitionId )
{
  std:: fill( m_origin, m_origin+3, 0.0 );
  std:: fill( m_h,      m_h+3 ,     1.0 );

  memcpy( m_origin, origin, dimension*sizeof(double) );
  memcpy( m_h,      h,      dimension*sizeof(double) );
}

//------------------------------------------------------------------------------
UniformMesh::~UniformMesh()
{
  // TODO Auto-generated destructor stub
}

} /* namespace mint */
} /* namespace axom */
