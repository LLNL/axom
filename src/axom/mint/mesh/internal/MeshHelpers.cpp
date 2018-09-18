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

#include "axom/mint/mesh/internal/MeshHelpers.hpp"
#include "axom/mint/mesh/Mesh.hpp"

#include <algorithm>
#include <map>
#include <vector>

namespace axom
{
namespace mint
{
namespace internal
{

//------------------------------------------------------------------------------
// IMPLEMENTATION OF FREE METHODS
//------------------------------------------------------------------------------
std::string join_ints_into_string(int count,
                                  m::IndexType * values,
                                  char sep)
{
  std::stringstream joined;

  for (int i = 0; i < count; ++i)
  {
    if (i > 0) { joined << sep; }
    joined << values[i];
  }

  return joined.str();
}


//------------------------------------------------------------------------------
bool initFaces(m::Mesh * m,
               int & facecount,
               m::IndexType *& f2c,
               m::IndexType *& c2f,
               m::IndexType *& c2foffsets)
{
  bool success = true;

  typedef std::pair< m::CellType, std::vector<m::IndexType> > FaceTypeAndCells;
  typedef std::map< std::string, FaceTypeAndCells > FaceBuilderType;
  typedef std::map< m::IndexType, std::string > IDtoKeyType;
  FaceBuilderType workface;

  // Iterate over each cell.
  const int cellcount = m->getNumberOfCells();
  for (int c = 0; c < cellcount; ++c)
  {
    // Step 1. For every cell, get the nodes.
    m::IndexType nodes[m::MAX_NUM_NODES];
    m->getCellNodeIDs(c, nodes);
    const m::CellType celltype = m->getCellType(c);
    const m::CellInfo thisCell = m::getCellInfo(celltype);

    int base = 0;
    m::IndexType face_nodes[m::MAX_ONE_FACE_NODES];
    for (int f = 0; f < thisCell.num_faces; ++f)
    {
      // Step 2. The cell nodes will be in "VTK Order," as specified by
      // https://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf.
      // This routine selects the face nodes as listed in the registered
      // cell info to make sure that face normals point outward.
      // 
      // For every face, sort the node IDs and join them together in a
      // string with '.' delimiters.  This is the key for an associative
      // array whose value is a list of cell IDs.
      int num_face_nodes = thisCell.face_nodecount[f];
      const m::IndexType * face_node_offsets = thisCell.face_nodes + base;
      for (int fn = 0; fn < num_face_nodes; ++fn)
      {
        face_nodes[fn] = nodes[face_node_offsets[fn]];
      }
      base += num_face_nodes;

      std::sort(face_nodes, face_nodes + num_face_nodes);
      std::string face_hash = 
        join_ints_into_string(num_face_nodes, face_nodes, '.');

      FaceTypeAndCells cells_of_face;
      if (workface.count(face_hash) > 0)
      {
        cells_of_face = workface[face_hash];
      }
      else
      {
        cells_of_face.first = thisCell.face_types[f];
      }
      cells_of_face.second.push_back(c);
      workface[face_hash] = cells_of_face;
    }
  }

  // Step 3. For each face, record its incident cells.
  // Here we use ConnectivityArray::reserve() and then append() each
  // face.  This means the kth inserted face gets ID k.
  int faceID = 0;
  IDtoKeyType keys;
  f2c = new m::IndexType[2 * workface.size()];
  for (FaceBuilderType::value_type v : workface)
  {
    int faceNodeCount = v.second.second.size();
    m::IndexType * faceNodes = v.second.second.data();

    if (faceNodeCount < 1 || faceNodeCount > 2)
    {
      success = false;
    }
    else
    {
      f2c[2*faceID] = faceNodes[0];
      f2c[2*faceID + 1] = (faceNodeCount > 1? faceNodes[1]: -1);
    }

    keys[faceID] = v.first;

    faceID += 1;
  }

  // Record how many faces we have in this mesh.
  facecount = faceID;

  // Step 4. Now that we have face IDs, record cell-to-face relation.
  typedef std::map< m::IndexType, std::vector<m::IndexType> >
    CellFaceBuilderType;

  CellFaceBuilderType cell_to_face;
  int cellFaceCount = 0;
  for (int f = 0; f < facecount; ++f)
  {
    int faceCellCount = workface[keys[f]].second.size();
    std::vector<m::IndexType> cell_faces;
    for (int c = 0; c < faceCellCount; ++c)
    {
      if (cell_to_face.count(c) > 0)
      {
        cell_faces = cell_to_face[c];
      }
      cell_faces.push_back(f);
      cellFaceCount += 1;
      cell_to_face[c] = cell_faces;
    }
  }

  // Step 4b. Put cell-to-face relation into output arrays.
  c2f = new m::IndexType[cellFaceCount];
  c2foffsets = new m::IndexType[cell_to_face.size() + 1];
  cellFaceCount = 0;

  for (int cellID = 0; cellID < cellcount; ++cellID)
  {
    // For the current cell ID, set its type and value offset
    //  cellTypes[cellID] = m_cell_connectivity->getIDType(cellID);
    c2foffsets[cellID] = cellFaceCount;

    if (cell_to_face.count(cellID) > 0)
    {
      // If we have faces for the current cell (which we certainly
      // *SHOULD* have), copy them in.
      std::vector<m::IndexType> & theCells = cell_to_face[cellID];
      m::IndexType * faceIDs = theCells.data();
      int thisCellFaceCount = theCells.size();

      // Copy in values
      for (int f = 0; f < thisCellFaceCount; ++f)
      {
        c2f[cellFaceCount + f] = faceIDs[f];
      }

      // Maintain offset
      cellFaceCount += thisCellFaceCount;
    }        
  }
  c2foffsets[cellcount] = cellFaceCount;

  return success;
}

} /* namespace internal */
} /* namespace mint */
} /* namespace axom */
