#include "mint/Extent.hpp"        /* for Extent */
#include "mint/DataTypes.hpp"     /* for globalIndex */

#include <cstring>                /* for memcpy() */

namespace axom
{
namespace mint
{


//------------------------------------------------------------------------------
Extent::Extent( int ndims, const globalIndex * ext ) :
  m_ndims( ndims )
{
  SLIC_ASSERT( ndims >= 1 && ndims <= 3 );

  // zero out all extents
  memset(m_extent, 0, 6 * sizeof( globalIndex ) );

  // copy in the user-supplied extent
  memcpy(m_extent, ext, 2 * ndims * sizeof( globalIndex )  );

  // compute strides
  m_jp = 0;
  m_kp = 0;
  if ( ndims > 1 )
  {
    m_jp = size( 0 );
  }

  if ( ndims > 2 )
  {
    m_kp = m_jp * size( 1 );
  }

  buildCellOffsets();
}

//------------------------------------------------------------------------------
Extent::Extent() :
  m_ndims(0),
  m_jp(0),
  m_kp(0)
{
  memset(m_extent, 0, 6 * sizeof( globalIndex ) );
  buildCellOffsets();
}

}   /* end namespace mint */
}   /* end namespace axom */
