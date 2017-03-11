/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */

#include "mint/Mesh.hpp"

// axom includes
#include "axom/CommonTypes.hpp"
#include "mint/FieldData.hpp"

// C/C++ includes
#include <cstddef>

namespace axom {
namespace mint {

Mesh::Mesh():
  m_ndims(2),
  m_block_idx(0),
  m_part_idx(0),
  m_type(MINT_UNDEFINED_MESH),
  m_node_data( new FieldData() ),
  m_cell_data( new FieldData() ),
  m_face_data( new FieldData() )
{}

//------------------------------------------------------------------------------
Mesh::Mesh(int ndims, int type, int blockId, int partId ):
  m_ndims( ndims ),
  m_block_idx( blockId ),
  m_part_idx( partId ),
  m_type( type ),
  m_node_data( new FieldData() ),
  m_cell_data( new FieldData() ),
  m_face_data( new FieldData() )
{}

//------------------------------------------------------------------------------
Mesh::~Mesh()
{
  m_node_data->clear();
  delete m_node_data;
  m_node_data = AXOM_NULLPTR;

  m_cell_data->clear();
  delete m_cell_data;
  m_cell_data = AXOM_NULLPTR;

  m_face_data->clear();
  delete m_face_data;
  m_face_data = AXOM_NULLPTR;
}

} /* namespace mint */
} /* namespace axom */
