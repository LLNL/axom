// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/quest/readers/ProEReader.hpp"

// Axom includes
#include "axom/core/utilities/Utilities.hpp"
#include "axom/mint/mesh/CellTypes.hpp"
#include "axom/slic/interface/slic.hpp"

// C/C++ includes
#include <fstream>
#include <map>

//------------------------------------------------------------------------------
//      ProEReader Implementation
//------------------------------------------------------------------------------
namespace axom
{
namespace quest
{
ProEReader::ProEReader() : m_fileName(""), m_num_nodes(0), m_num_tets(0) { }

//------------------------------------------------------------------------------
ProEReader::~ProEReader() { this->clear(); }

//------------------------------------------------------------------------------
void ProEReader::clear()
{
  m_num_nodes = 0;
  m_num_tets = 0;
  m_num_unique_nodes = 0;
  m_nodes.clear();
}

//------------------------------------------------------------------------------
int ProEReader::read()
{
  constexpr int NUM_NODES_PER_TET = 4;
  constexpr int NUM_COMPS_PER_NODE = 3;

  std::string junk;
  int id;
  int tet_nodes[NUM_NODES_PER_TET];

  struct coordinate
  {
    double comp[NUM_COMPS_PER_NODE];
  } cur_coord;

  std::ifstream ifs(m_fileName.c_str());

  if(!ifs.is_open())
  {
    SLIC_WARNING("Cannot open the provided Pro/E file [" << m_fileName << "]");
    return (-1);
  }

  // Remove any comments from header
  while(ifs.peek() == '#')
  {
    std::getline(ifs, junk);
  }

  // Initialize number of nodes and tetrahedra.
  // 4 nodes per tet, 3 components per node.
  ifs >> m_num_unique_nodes >> m_num_tets;
  m_num_nodes = m_num_tets * NUM_NODES_PER_TET;
  m_nodes.reserve(m_num_nodes * NUM_COMPS_PER_NODE);

  // Read nodes
  std::map<int, coordinate> nodes;
  for(int i = 0; i < m_num_unique_nodes; i++)
  {
    ifs >> id >> cur_coord.comp[0] >> cur_coord.comp[1] >> cur_coord.comp[2];
    nodes[id] = {cur_coord.comp[0], cur_coord.comp[1], cur_coord.comp[2]};
  }

  // Initialize nodes
  for(int i = 0; i < m_num_tets; i++)
  {
    ifs >> id >> tet_nodes[0] >> tet_nodes[1] >> tet_nodes[2] >> tet_nodes[3];

    for(int j = 0; j < NUM_NODES_PER_TET; j++)
    {
      for(int k = 0; k < NUM_COMPS_PER_NODE; k++)
      {
        m_nodes.push_back(nodes[tet_nodes[j]].comp[k]);
      }
    }
  }

  ifs.close();
  return (0);
}

//------------------------------------------------------------------------------
void ProEReader::getMesh(axom::mint::UnstructuredMesh<mint::SINGLE_SHAPE>* mesh)
{
  /* Sanity checks */
  SLIC_ERROR_IF(mesh == nullptr, "supplied mesh is null!");
  SLIC_ERROR_IF(static_cast<axom::IndexType>(m_nodes.size()) != m_num_nodes * 3,
                "nodes vector size doesn't match expected size!");
  SLIC_ERROR_IF(mesh->getDimension() != 3, "Pro/E reader expects a 3D mesh!");
  SLIC_ERROR_IF(mesh->getCellType() != mint::TET,
                "Pro/E reader expects a tetrahedra mesh!");

  // pre-allocate space to store the mesh
  if(!mesh->isExternal())
  {
    mesh->resize(m_num_nodes, m_num_tets);
  }

  SLIC_ERROR_IF(mesh->getNumberOfNodes() != m_num_nodes,
                "mesh number of nodes does not match the number of nodes in "
                "the Pro/E file!");
  SLIC_ERROR_IF(mesh->getNumberOfCells() != m_num_tets,
                "mesh number of cells does not match number of tetrahedra in "
                "the Pro/E file!");

  double* x = mesh->getCoordinateArray(mint::X_COORDINATE);
  double* y = mesh->getCoordinateArray(mint::Y_COORDINATE);
  double* z = mesh->getCoordinateArray(mint::Z_COORDINATE);

  // Load the vertices into the mesh
  for(axom::IndexType i = 0; i < m_num_nodes; ++i)
  {
    const axom::IndexType offset = i * 3;
    x[i] = m_nodes[offset];
    y[i] = m_nodes[offset + 1];
    z[i] = m_nodes[offset + 2];
  }

  // Load the tetrahedra.  Note that the indices are implicitly defined.
  axom::IndexType* conn = mesh->getCellNodesArray();
  for(axom::IndexType i = 0; i < m_num_tets; ++i)
  {
    const axom::IndexType offset = i * 4;
    conn[offset] = offset;
    conn[offset + 1] = offset + 1;
    conn[offset + 2] = offset + 2;
    conn[offset + 3] = offset + 3;
  }
}

}  // end namespace quest
}  // end namespace axom
