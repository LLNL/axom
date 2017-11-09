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

#include "mint/StructuredMesh.hpp"
#include "mint/Mesh.hpp"                    /* For Mesh */
#include "mint/DataTypes.hpp"               /* For localIndex */
#include "axom/Types.hpp"                   /* For AXOM_NULLPTR */
#include "mint/Extent.hpp"                  /* For Extent */
#include "mint/CellType.hpp"                /* For MINT_HEX, MINT_QUAD, ect. */

namespace axom
{
namespace mint
{

//------------------------------------------------------------------------------
StructuredMesh::StructuredMesh( int meshType, int ndims,
                                const globalIndex ext[ 6 ]) :
  Mesh( ndims, meshType, 0, 0 ),
  m_extent( ndims, ext )
{}

//------------------------------------------------------------------------------
StructuredMesh::StructuredMesh( int meshType, int ndims,
                                const globalIndex ext[ 6 ], int blockId,
                                int partId) :
  Mesh( ndims, meshType, blockId, partId ),
  m_extent( ndims, ext )
{}

//------------------------------------------------------------------------------
int StructuredMesh::getMeshCellType( localIndex AXOM_NOT_USED( cellIdx ) ) const
{
  const int dim = this->getDimension();
  if ( dim == 3 )
  {
    return MINT_HEX;
  }
  else if ( dim == 2 )
  {
    return MINT_QUAD;
  }
  return MINT_SEGMENT;
}

}   /* end namespace mint */
}   /* end namespace axom */
