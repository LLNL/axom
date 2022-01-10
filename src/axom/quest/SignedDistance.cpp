// Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/quest/SignedDistance.hpp"

#include "axom/mint/mesh/UnstructuredMesh.hpp"

namespace axom
{
namespace quest
{
namespace detail
{
AXOM_HOST_DEVICE mint::CellType UcdMeshData::getCellType(IndexType cellId) const
{
  if(shape_type == mint::SINGLE_SHAPE)
  {
    SLIC_ASSERT(single_cell_type != mint::UNDEFINED_CELL);
    return single_cell_type;
  }
  else
  {
    SLIC_ASSERT(cell_types != nullptr);
    return cell_types[cellId];
  }
}

AXOM_HOST_DEVICE const IndexType* UcdMeshData::getCellNodeIDs(IndexType cellId,
                                                              int& nnodes) const
{
  SLIC_ASSERT(cells_to_nodes != nullptr);
  int cellBegin;
  if(shape_type == mint::SINGLE_SHAPE)
  {
    SLIC_ASSERT(nodes_per_cell != -1);
    cellBegin = cellId * nodes_per_cell;
    nnodes = nodes_per_cell;
  }
  else
  {
    SLIC_ASSERT(cell_node_offsets != nullptr);
    cellBegin = cell_node_offsets[cellId];
    nnodes = cell_node_offsets[cellId + 1] - cell_node_offsets[cellId];
  }
  return cells_to_nodes + cellBegin;
}

bool SD_GetUcdMeshData(const mint::Mesh* surfaceMesh, UcdMeshData& outSurfData)
{
  using SingleShapeMesh = mint::UnstructuredMesh<mint::SINGLE_SHAPE>;
  using MixedShapeMesh = mint::UnstructuredMesh<mint::MIXED_SHAPE>;

  if(typeid(*surfaceMesh) == typeid(const SingleShapeMesh))
  {
    auto pSingleSurfMesh = static_cast<const SingleShapeMesh*>(surfaceMesh);
    outSurfData.shape_type = mint::SINGLE_SHAPE;
    outSurfData.single_cell_type = pSingleSurfMesh->getCellType();
    outSurfData.cells_to_nodes = pSingleSurfMesh->getCellNodesArray();
    outSurfData.nodes_per_cell = pSingleSurfMesh->getNumberOfCellNodes();
    outSurfData.cell_node_offsets = nullptr;
  }
  else if(typeid(*surfaceMesh) == typeid(const MixedShapeMesh))
  {
    auto pMixSurfMesh = static_cast<const MixedShapeMesh*>(surfaceMesh);
    outSurfData.shape_type = mint::MIXED_SHAPE;
    outSurfData.cell_types = pMixSurfMesh->getCellTypesArray();
    outSurfData.cells_to_nodes = pMixSurfMesh->getCellNodesArray();
    outSurfData.nodes_per_cell = -1;
    outSurfData.cell_node_offsets = pMixSurfMesh->getCellNodesOffsetsArray();
  }
  else
  {
    SLIC_ASSERT_MSG(false, "Mesh is not an unstructured surface mesh.");
    return false;
  }
  return true;
}

}  // end namespace detail
}  // end namespace quest
}  // end namespace axom
