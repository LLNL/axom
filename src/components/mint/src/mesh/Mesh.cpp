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

#include "mint/Mesh.hpp"

// axom includes
#include "axom/Types.hpp"
#include "mint/FieldData.hpp"

// C/C++ includes
#include <cstddef>

namespace axom
{
namespace mint
{

Mesh::Mesh() :
  m_ndims(2),
  m_type(MINT_UNDEFINED_MESH),
  m_block_idx(0),
  m_part_idx(0),
  m_cell_data(),
  m_face_data(),
  m_edge_data(),
  m_node_data()
{}

//------------------------------------------------------------------------------
Mesh::Mesh( int ndims, int type, int blockId, int partId ) :
  m_ndims( ndims ),
  m_type( type ),
  m_block_idx( blockId ),
  m_part_idx( partId ),
  m_cell_data(),
  m_face_data(),
  m_edge_data(),
  m_node_data()
{}

//------------------------------------------------------------------------------
Mesh::~Mesh()
{}



} /* namespace mint */
} /* namespace axom */
