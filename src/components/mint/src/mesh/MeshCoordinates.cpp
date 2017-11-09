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

#include "mint/MeshCoordinates.hpp"     /* for MeshCoordinates */
#include "mint/DataTypes.hpp"           /* for localIndex */
#include "mint/Vector.hpp"              /* for Vector */

namespace axom
{
namespace mint
{


//------------------------------------------------------------------------------
MeshCoordinates::MeshCoordinates( int dimension, localIndex capacity,
                                  double resize_ratio ) :
  m_ndims( dimension )
{
  for ( int dim = 0 ; dim < m_ndims ; ++dim )
  {
    m_coordinates[ dim ].setCapacity( capacity );
    m_coordinates[ dim ].setResizeRatio( resize_ratio );
  }

  for (int dim = m_ndims ; dim < 3 ; ++dim )
  {
    m_coordinates[ dim ].setCapacity( 0 );
    m_coordinates[ dim ].setResizeRatio( 0.0 );
  }
}

//------------------------------------------------------------------------------
MeshCoordinates::MeshCoordinates() :
  m_ndims(0)
{
  for (int dim = 0 ; dim < 3 ; ++dim )
  {
    m_coordinates[ dim ].setCapacity( 0 );
    m_coordinates[ dim ].setResizeRatio( 0.0 );
  }
}

}   /* end namespace mint */
}   /* end namespace axom */
