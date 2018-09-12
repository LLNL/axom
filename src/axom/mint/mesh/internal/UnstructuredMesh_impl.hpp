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

#ifndef MINT_UNSTRUCTUREDMESH_IMPL_HPP_
#define MINT_UNSTRUCTUREDMESH_IMPL_HPP_

#include "axom/mint/mesh/UnstructuredMesh.hpp"
#include "axom/mint/mesh/MeshTypes.hpp"
#include "axom/mint/mesh/blueprint.hpp"

#include <unordered_map>
#include <sstream>

namespace
{
  static inline
  std::string join_ints_into_string(int count,
                                    axom::mint::IndexType * values,
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
} // end anonymous namespace

namespace axom
{
namespace mint
{

//------------------------------------------------------------------------------
// IMPLEMENTATION
//------------------------------------------------------------------------------
  template < Topology TOPO >
  void UnstructuredMesh<TOPO>::initializeFaceConnectivity()
  {
    typedef std::pair< CellType, std::vector<int> > FaceTypeAndCells;
    typedef std::unordered_map< std::string, FaceTypeAndCells > FaceBuilderType;
    FaceBuilderType workface;

    // Step 0. Allocate arrays to hold face association
    m_cell_to_face = new FaceConnectivity( );
    m_face_to_cell = new FaceConnectivity( );

    // Iterate over each cell.
    int cellcount = getNumberOfCells();
    for (int c = 0; c < cellcount; ++c)
    {
      // Step 1. For every cell, get the nodes.
      const IndexType * nodes = getCellNodeIDs(c);
      const CellType celltype = getCellType(c);
      CellInfo thisCell = getCellInfo(celltype);

      int base = 0;
      IndexType face_nodes[MAX_ONE_FACE_NODES];
      for (int f = 0; f < thisCell.num_faces; ++f)
      {
        // Step 2. The face nodes will be in "VTK Order," as specified by
        // https://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf.
        // For every face, sort the node IDs and join them together in a
        // string with '.' delimiters.  This is the key for an associative
        // array whose value is a list of cell IDs.
        int num_face_nodes = thisCell.face_nodecount[f];
        const int * face_node_offsets = thisCell.face_nodes + base;
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
    m_face_to_cell->reserve(workface.size());
    for (FaceBuilderType::value_type v : workface)
    {
      CellType faceType = v.second.first;
      int faceNodeCount = v.second.second.size();
      int * faceNodes = v.second.second.data();
      m_face_to_cell->append(faceNodes, faceNodeCount, faceType);
    }

    // Step 4. Now that we have face IDs, tie faces to cells.
    typedef std::map< IndexType, std::vector<IndexType> > CellFaceBuilderType;

    CellFaceBuilderType cell_to_face;
    int faceCount = m_face_to_cell->getNumberOfIDs();
    int cellFaceCount = 0;
    for (int f = 0; f < faceCount; ++f)
    {
      IndexType faceCellCount = m_face_to_cell->getNumberOfValuesForID(f);
      std::vector<IndexType> cell_faces;
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

    // In contrast to Step 3, here we use ConnectivityArray::resize()
    // and then set() each cell with faces.  This is because the cell IDs
    // in this container have to match IDs already established in
    // m_cell_connectivity.
    SLIC_ASSERT( (int)m_cell_connectivity->getNumberOfIDs() ==
                 (int)cell_to_face.size() );
    m_cell_to_face->resize(cell_to_face.size(), cellFaceCount);
    cellFaceCount = 0;

    // Because ConnectivityArray::set() won't set the ID type or the number
    // of values for each ID, we have to work with the raw arrays.
    CellType * cellTypes = m_cell_to_face->getTypePtr();
    IndexType * cellFaceIDs = m_cell_to_face->getValuePtr();
    IndexType * cellFaceOffsets = m_cell_to_face->getOffsetPtr();

    int cellCount = getNumberOfCells();
    for (int cellID = 0; cellID < cellCount; ++cellID)
    {
      // For the current cell ID, set its type and value offset
      cellTypes[cellID] = m_cell_connectivity->getIDType(cellID);
      cellFaceOffsets[cellID] = cellFaceCount;

      if (cell_to_face.count(cellID) > 0)
      {
        // If we have faces for the current cell (which we certainly
        // *SHOULD* have), copy them in.
        std::vector<IndexType> & theCells = cell_to_face[cellID];
        IndexType * faceIDs = theCells.data();
        int faceCount = theCells.size();

        // Copy in values
        for (int f = 0; f < faceCount; ++f)
        {
          cellFaceIDs[cellFaceCount + f] = faceIDs[f];
        }

        // Maintain offset
        cellFaceCount += faceCount;
      }        
    }
    cellFaceOffsets[cell_to_face.size()] = cellFaceCount;
  }

  
}   /* end namespace mint */
}   /* end namespace axom */

#endif /* MINT_UNSTRUCTUREDMESH_IMPL_HPP_ */
