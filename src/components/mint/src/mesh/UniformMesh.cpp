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

#include "mint/UniformMesh.hpp"

#include "axom_utils/Utilities.hpp"       /* for utilities::abs */
#include "mint/MeshTypes.hpp"              /* for MINT_STRUCTURED_UNIFORM_MESH */
#include "mint/config.hpp"             /* for IndexType */

#include <algorithm>                      /* for std::fill() */

namespace axom
{
namespace mint
{

//------------------------------------------------------------------------------
UniformMesh::UniformMesh( int dimension, const double origin[3],
                          const double h[3],
                          const int64 ext[6] ) :
  StructuredMesh( mint::STRUCTURED_UNIFORM_MESH, dimension, ext )
{
  std::fill( m_origin, m_origin + 3, 0.0 );
  std::fill( m_h,      m_h + 3,      1.0 );

  memcpy( m_origin, origin, dimension * sizeof( double ) );
  memcpy( m_h,      h,      dimension * sizeof( double ) );
}

//------------------------------------------------------------------------------
UniformMesh::UniformMesh( int dimension, const int64 ext[6],
                          const double lower_bound[3],
                          const double upper_bound[3] ) :
  StructuredMesh( mint::STRUCTURED_UNIFORM_MESH, dimension, ext )

{
  std::fill( m_origin, m_origin + 3, 0.0 );
  std::fill( m_h,      m_h + 3,      1.0 );

  double h[ dimension ];
  double dim_length;
  for ( int dim = 0 ; dim < dimension ; ++dim )
  {
    dim_length = utilities::abs( lower_bound[ dim ] - upper_bound[ dim ] );
    h[ dim ] = dim_length / ( m_extent->size( dim ) - 1.0 );
  }

  memcpy( m_origin, lower_bound,  dimension * sizeof( double ) );
  memcpy( m_h,      h,            dimension * sizeof( double ) );
}

} /* namespace mint */
} /* namespace axom */
