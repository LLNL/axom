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

#include "mint/Extent.hpp"     /* for Extent */
#include "mint/config.hpp"     /* for int64 */

#include <cstring>             /* for memcpy() */

namespace axom
{
namespace mint
{

//------------------------------------------------------------------------------
Extent::Extent( int ndims, const int64* ext ) :
    m_ndims( ndims ),
    m_numnodes( 1 ),
    m_numcells( 1 )
{
  SLIC_ERROR_IF( ext==AXOM_NULLPTR, "provided extent is null!" );
  SLIC_ERROR_IF( !(m_ndims >= 1 && m_ndims <= 3),
                 "provided dimension is invalid!"  );

  // zero out extent information
  memset(m_sizes, 0, 3 * sizeof(IndexType) );
  memset(m_extent, 0, 6 * sizeof( int64 ) );

  // copy in the user-supplied extent
  memcpy(m_extent, ext, 2 * ndims * sizeof( int64 )  );

  for ( int i=0; i < m_ndims; ++i )
  {
    const int64 ilo = ext[ i*2   ];
    const int64 ihi = ext[ i*2+1 ];
    SLIC_ERROR_IF( (ihi < ilo), "invalid extent!" );

    m_sizes[ i ]  = ihi - ilo + 1;
    m_numnodes   *= m_sizes[ i ];
    m_numcells   *= ( ( m_sizes[ i ] > 1 )? m_sizes[ i ] - 1 : 1 );
  }

  // compute strides
  m_jp = ( m_ndims > 1) ? m_sizes[ 0 ] : 0;
  m_kp = ( m_ndims > 2) ? m_jp * m_sizes[ 1 ] : 0;

  // build the cell offsets
  buildCellOffsets();
}

}   /* end namespace mint */
}   /* end namespace axom */
