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

#include "CurvilinearMesh.hpp"
#include "MeshTypes.hpp"

// axom includes
#include "slic/slic.hpp"

// C/C++ includes
#include <cstddef> // for definition of AXOM_NULLPTR

namespace axom
{
namespace mint
{

//------------------------------------------------------------------------------
CurvilinearMesh::CurvilinearMesh( int ndims, const int64 ext[6] ) :
  StructuredMesh( STRUCTURED_CURVILINEAR_MESH, ndims, ext ),
  m_coordinates( ndims, m_extent->getNumNodes() )
{
//  TODO: ???
//  m_coordinates.setSize( m_extent.getNumNodes() );
}

} /* namespace mint */
} /* namespace axom */
