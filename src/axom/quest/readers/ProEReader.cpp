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
  m_nodes.clear();
}

//------------------------------------------------------------------------------
int ProEReader::read()
{
  std::string id;
  std::string tet_nodes[4];

  struct coordinate
  {
    double comp[3];
  } cur_coord;

  std::ifstream ifs(m_fileName.c_str());

  if(!ifs.is_open())
  {
    SLIC_WARNING("Cannot open the provided Pro/E file [" << m_fileName << "]");
    return (-1);
  }

  // Initialize number of nodes and tetrahedra
  ifs >> m_num_nodes >> m_num_tets;
  m_nodes.reserve(m_num_tets * 4 * 3);

  // Read nodes
  std::map<std::string, coordinate> nodes;
  for(int i = 0; i < m_num_nodes; i++)
  {
    ifs >> id >> cur_coord.comp[0] >> cur_coord.comp[1] >> cur_coord.comp[2];
    nodes[id] = {cur_coord.comp[0], cur_coord.comp[1], cur_coord.comp[2]};
  }

  // Initialize nodes
  for(int i = 0; i < m_num_tets; i++)
  {
    ifs >> id >> tet_nodes[0] >> tet_nodes[1] >> tet_nodes[2] >> tet_nodes[3];
    int tet_offset = i * 12;

    for(int j = 0; j < 4; j++)
    {
      int vert_offset = j * 3;
      for(int k = 0; k < 3; k++)
      {
        m_nodes[tet_offset + vert_offset + k] = nodes[tet_nodes[j]].comp[k];
      }
    }
  }

  ifs.close();
  return (0);
}

//------------------------------------------------------------------------------
void ProEReader::getMesh(axom::mint::UnstructuredMesh<mint::SINGLE_SHAPE>* mesh)
{ }

}  // end namespace quest
}  // end namespace axom
