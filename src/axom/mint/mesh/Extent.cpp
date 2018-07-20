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

#include "axom/mint/mesh/Extent.hpp"     /* for Extent */
#include "axom/mint/config.hpp"     /* for int64 */

#include <cstring>             /* for memcpy() */

namespace axom
{
namespace mint
{

//------------------------------------------------------------------------------
Extent::Extent( int ndims, const int64* ext ) :
  m_ndims( ndims )
{
  SLIC_ERROR_IF( ext == AXOM_NULLPTR, "provided extent is null!" );
  SLIC_ERROR_IF( !(m_ndims >= 1 && m_ndims <= 3), 
                 "provided dimension is invalid!"  );

  // copy in the user-supplied extent
  std::memcpy( m_extent, ext, 2 * ndims * sizeof( int64 )  );

  for ( int i = 0; i < m_ndims; ++i )
  {
    const int64 ilo = ext[ i * 2 ];
    const int64 ihi = ext[ i * 2 + 1 ];
    SLIC_ERROR_IF( ihi <= ilo, "invalid extent!" );

    m_sizes[ i ]  = ihi - ilo + 1;
  }

  // compute strides
  m_node_jp = ( m_ndims > 1 ) ? m_sizes[ 0 ] : 0;
  m_node_kp = ( m_ndims > 2 ) ? m_node_jp * m_sizes[ 1 ] : 0;

  m_cell_jp = ( m_ndims > 1 ) ? m_sizes[ 0 ] - 1 : 0;
  m_cell_kp = ( m_ndims > 2 ) ? m_cell_jp * ( m_sizes[ 2 ] - 1 ) : 0;

  // calculate the number of nodes, cells, faces and edges.
  calculateNumberOfNodes();
  calculateNumberOfCells();
  calculateNumberOfFaces();
  calculateNumberOfEdges();

  // build the cell offsets
  buildCellOffsets();
}

//------------------------------------------------------------------------------
inline void Extent::buildCellOffsets()
{
  m_cell_offsets[ 0 ] = 0;
  m_cell_offsets[ 1 ] = 1;
  m_cell_offsets[ 2 ] = 1 + m_node_jp;
  m_cell_offsets[ 3 ] = m_node_jp;

  m_cell_offsets[ 4 ] = m_node_kp;
  m_cell_offsets[ 5 ] = 1 + m_node_kp;
  m_cell_offsets[ 6 ] = 1 + m_node_jp + m_node_kp;
  m_cell_offsets[ 7 ] = m_node_jp + m_node_kp;
}

//------------------------------------------------------------------------------
void Extent::calculateNumberOfNodes()
{
  m_numnodes = 1;

  for ( int dim = 0; dim < m_ndims; ++dim )
  {
    m_numnodes *= m_sizes[ dim ];
  }
}

//------------------------------------------------------------------------------
void Extent::calculateNumberOfCells()
{
  m_numcells = 1;

  for ( int dim = 0; dim < m_ndims; ++dim )
  {
    m_numcells *= m_sizes[ dim ] - 1;
  }
}

//------------------------------------------------------------------------------
void Extent::calculateNumberOfFaces()
{
  if ( m_ndims == 1 )
  {
    m_numfaces = 0;
  }
  else if ( m_ndims == 2 )
  {
    m_numfaces = (m_sizes[ 0 ] - 1) * m_sizes[ 1 ]
               + (m_sizes[ 1 ] - 1) * m_sizes[ 0 ];
  }
  else
  {
    m_numfaces = (m_sizes[ 0 ] - 1) * (m_sizes[ 1 ] - 1) * m_sizes[ 2 ]
               + (m_sizes[ 0 ] - 1) * (m_sizes[ 2 ] - 1) * m_sizes[ 1 ]
               + (m_sizes[ 1 ] - 1) * (m_sizes[ 2 ] - 1) * m_sizes[ 0 ];

  }
}

//------------------------------------------------------------------------------
void Extent::calculateNumberOfEdges()
{
  if ( m_ndims < 3 )
  {
    m_numedges = 0;
  }
  else
  {
    m_numedges = (m_sizes[ 0 ] - 1) * m_sizes[ 1 ] * m_sizes[ 2 ]
               + (m_sizes[ 1 ] - 1) * m_sizes[ 0 ] * m_sizes[ 2 ]
               + (m_sizes[ 2 ] - 1) * m_sizes[ 0 ] * m_sizes[ 1 ];
  }
}


}   /* end namespace mint */
}   /* end namespace axom */
