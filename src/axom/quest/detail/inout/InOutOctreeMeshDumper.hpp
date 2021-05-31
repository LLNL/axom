// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file InOutOctreeMeshDumper.hpp
 *
 * \brief Defines helper class to write meshes for InOutOctree instances
 */

#ifndef INOUT_OCTREE_MESHDUMPER__HXX_
#define INOUT_OCTREE_MESHDUMPER__HXX_

#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/slam.hpp"
#include "axom/primal.hpp"
#include "axom/spin.hpp"

#include "BlockData.hpp"
#include "MeshWrapper.hpp"

#include "fmt/fmt.hpp"

namespace axom
{
namespace quest
{
// Predeclare InOutOctree class
template <int DIM>
class InOutOctree;

namespace detail
{
template <int DIM>
class InOutOctreeMeshDumper
{
public:
  using InOutOctreeType = InOutOctree<DIM>;
  using CellIndexSet = typename InOutOctreeType::CellIndexSet;

  using OctreeBaseType = typename InOutOctreeType::OctreeBaseType;
  using OctreeLevels = typename OctreeBaseType::OctreeLevels;
  using BlockIndex = typename OctreeBaseType::BlockIndex;

  using SpacePt = typename InOutOctreeType::SpacePt;
  using GridPt = typename InOutOctreeType::GridPt;
  using VertexIndex = typename InOutOctreeType::VertexIndex;
  using CellIndex = typename InOutOctreeType::CellIndex;
  using CellVertIndices = typename MeshWrapper<DIM>::CellVertIndices;
  using GeometricBoundingBox = typename InOutOctreeType::GeometricBoundingBox;
  using SpaceCell = typename InOutOctreeType::SpaceCell;

  using LeafVertMap = slam::Map<slam::Set<VertexIndex>, VertexIndex>;
  using LeafIntMap = slam::Map<slam::Set<axom::IndexType>, axom::IndexType>;
  using LeafGridPtMap = slam::Map<slam::Set<axom::IndexType>, GridPt>;

  using DebugMesh = mint::UnstructuredMesh<mint::MIXED_SHAPE>;

  using ColorsMap = std::map<InOutBlockData::LeafColor, int>;

  using GridPtHash = spin::PointHash<typename GridPt::CoordType>;
  using GridIntMap = std::unordered_map<GridPt, int, GridPtHash>;

public:
  InOutOctreeMeshDumper(const InOutOctreeType& octree)
    : m_octree(octree)
    , m_generationState(m_octree.m_generationState)
  {
    // Create a small lookup table to map block colors to ints
    m_colorsMap[InOutBlockData::White] = -1;
    m_colorsMap[InOutBlockData::Gray] = 0;
    m_colorsMap[InOutBlockData::Black] = 1;
  }

  /**
   *  Generates a hexahedral VTK mesh with all neighboring blocks where one is
   * inside and the other is outside
   *  \note By construction, there should be no such pairs in a valid
   * InOutOctree mesh.
   */
  void dumpDifferentColoredNeighborsMeshVTK(const std::string& name) const
  {
    if(m_generationState < InOutOctreeType::INOUTOCTREE_LEAVES_COLORED)
    {
      SLIC_INFO("Need to generate octree colors before visualizing them.");
      return;
    }

    using LevelGridIntMap = slam::Map<slam::Set<>, GridIntMap>;
    LevelGridIntMap diffBlocks(&(m_octree.m_levels));

    int totalBlocks = 0;

    // Iterate through the octree leaves
    // looking for neighbor blocks with different labelings
    for(int lev = m_octree.maxLeafLevel() - 1; lev >= 0; --lev)
    {
      diffBlocks[lev] = GridIntMap(0, GridPtHash());

      const auto& levelLeafMap = m_octree.getOctreeLevel(lev);
      auto itEnd = levelLeafMap.end();
      for(auto it = levelLeafMap.begin(); it != itEnd; ++it)
      {
        const BlockIndex block(it.pt(), lev);
        const InOutBlockData& data = *it;
        if(data.isLeaf() && data.color() != InOutBlockData::Gray)
        {
          for(int i = 0; i < block.numFaceNeighbors(); ++i)
          {
            BlockIndex neighborBlk =
              m_octree.coveringLeafBlock(block.faceNeighbor(i));
            if(neighborBlk != BlockIndex::invalid_index())
            {
              const InOutBlockData& neighborData = m_octree[neighborBlk];
              switch(neighborData.color())
              {
              case InOutBlockData::Black:  // intentional fallthrough
              case InOutBlockData::White:
                if(data.color() != neighborData.color())
                {
                  diffBlocks[lev][block.pt()] = 1;
                  diffBlocks[neighborBlk.level()][neighborBlk.pt()] = 1;
                }
                break;
              case InOutBlockData::Gray:
              case InOutBlockData::Undetermined:
                break;
              }
            }
          }
        }
      }
      totalBlocks += diffBlocks[lev].size();
    }

    // Add all such blocks to a vector
    std::vector<BlockIndex> blocks;
    blocks.reserve(totalBlocks);

    for(int lev = m_octree.maxLeafLevel() - 1; lev >= 0; --lev)
    {
      for(auto&& blk : diffBlocks[lev])
      {
        blocks.emplace_back(BlockIndex(blk.first, lev));
      }
    }

    // Generate a VTK mesh with these blocks
    dumpOctreeMeshBlocks(name, blocks, false);
  }

  void dumpLocalOctreeMeshesForVertex(const std::string& name,
                                      VertexIndex vIdx) const
  {
    std::string vertStr = fmt::format("{}vertex_{}", name, vIdx);

    BlockIndex vertexBlock = m_octree.m_vertexToBlockMap[vIdx];

    // Dump a mesh for the vertex's containing block
    std::string blockStr = fmt::format("{}_block_{}_{}_{}",
                                       vertStr,
                                       vertexBlock.pt()[0],
                                       vertexBlock.pt()[1],
                                       vertexBlock.pt()[2]);
    std::vector<BlockIndex> blocks;
    blocks.push_back(vertexBlock);
    dumpOctreeMeshBlocks(blockStr, blocks, true);

    // Dump a mesh for the incident triangles
    std::string triStr = fmt::format("{}_triangles", vertStr);
    std::vector<CellIndex> tris;

    CellIndexSet triSet = m_octree.leafCells(vertexBlock, m_octree[vertexBlock]);
    for(int i = 0; i < triSet.size(); ++i)
    {
      CellIndex tIdx = triSet[i];
      CellVertIndices tv = m_octree.m_meshWrapper.cellVertexIndices(tIdx);
      for(int j = 0; j < tv.size(); ++j)
      {
        if(tv[j] == vIdx) tris.push_back(tIdx);
      }
    }
    dumpTriangleMesh(triStr, tris, true);
  }

  void dumpLocalOctreeMeshesForTriangle(const std::string& name,
                                        CellIndex tIdx) const
  {
    std::string triStr = fmt::format("{}triangle_{}", name, tIdx);

    // Dump a triangle mesh with the single triangle
    std::vector<CellIndex> tris;
    tris.push_back(tIdx);
    dumpTriangleMesh(triStr, tris, true);

    // Dump a hex mesh of all blocks that index this triangle
    std::vector<BlockIndex> blocks;
    for(int lev = 0; lev < m_octree.m_levels.size(); ++lev)
    {
      const auto& levelLeafMap = m_octree.getOctreeLevel(lev);
      auto itEnd = levelLeafMap.end();
      for(auto it = levelLeafMap.begin(); it != itEnd; ++it)
      {
        if(it->isLeaf() && it->hasData())
        {
          BlockIndex leafblk(it.pt(), lev);
          CellIndexSet triSet = m_octree.leafCells(leafblk, *it);

          bool found = false;
          for(int i = 0; !found && i < triSet.size(); ++i)
          {
            if(triSet[i] == tIdx)
            {
              blocks.push_back(leafblk);
              found = true;
            }
          }
        }
      }
    }
    dumpOctreeMeshBlocks(fmt::format("{}_blocks", triStr), blocks, true);
  }

  void dumpLocalOctreeMeshesForBlock(const std::string& name,
                                     const BlockIndex& block) const
  {
    // Dump a mesh with the single block
    std::string blockStr = fmt::format("{}block_{}_{}_{}",
                                       name,
                                       block.pt()[0],
                                       block.pt()[1],
                                       block.pt()[2]);

    std::vector<BlockIndex> blocks;
    blocks.push_back(block);
    dumpOctreeMeshBlocks(blockStr, blocks, true);

    // Dump a mesh with the indexed triangles for this block
    const InOutBlockData& blkData = m_octree[block];
    if(blkData.isLeaf() && blkData.hasData())
    {
      std::vector<CellIndex> tris;
      CellIndexSet triSet = m_octree.leafCells(block, blkData);
      for(int i = 0; i < triSet.size(); ++i)
      {
        tris.push_back(triSet[i]);
      }

      dumpTriangleMesh(fmt::format("{}_triangles", blockStr), tris, true);
    }
  }

  /** Generates a hexahedral VTK mesh with all octree blocks   */
  void dumpOctreeMeshVTK(const std::string& name) const
  {
    std::vector<BlockIndex> blocks;

    // Create an stl vector of all leaf blocks
    for(int lev = 0; lev < m_octree.m_levels.size(); ++lev)
    {
      const auto& levelLeafMap = m_octree.getOctreeLevel(lev);
      auto itEnd = levelLeafMap.end();
      for(auto it = levelLeafMap.begin(); it != itEnd; ++it)
      {
        if(it->isLeaf()) blocks.push_back(BlockIndex(it.pt(), lev));
      }
    }
    SLIC_INFO("Dump vtk:: Octree has " << blocks.size() << " leaves.");

    dumpOctreeMeshBlocks(name, blocks, false);
  }

  /** Generates a VTK mesh with all triangles in the mesh */
  void dumpTriMeshVTK(const std::string& name) const
  {
    const int numElts = m_octree.m_meshWrapper.numMeshCells();

    std::vector<CellIndex> tris;
    tris.reserve(numElts);

    for(int i = 0; i < numElts; ++i)
    {
      tris.push_back(i);
    }
    SLIC_INFO("Dump vtk:: Mesh has " << numElts << " triangles.");

    dumpTriangleMesh(name, tris, false);
  }

private:
  void dumpOctreeMeshBlocks(const std::string& name,
                            const std::vector<BlockIndex>& blocks,
                            bool shouldLogBlocks = false) const
  {
    if(blocks.empty()) return;

    // Dump an octree mesh containing all blocks
    std::string fName = fmt::format("{}.vtk", name);

    DebugMesh* debugMesh = new DebugMesh(3, 8 * blocks.size(), blocks.size());
    const bool hasCells =
      (m_generationState >= InOutOctreeType::INOUTOCTREE_ELEMENTS_INSERTED);
    const bool hasColors =
      (m_generationState >= InOutOctreeType::INOUTOCTREE_LEAVES_COLORED);

    // Allocate Slam Maps for the field data
    slam::PositionSet<> leafSet(blocks.size());

    LeafVertMap leafVertID(&leafSet);
    LeafVertMap leafVertID_unique(&leafSet);
    LeafIntMap leafTriCount(&leafSet);
    LeafIntMap leafColors(&leafSet);
    LeafIntMap leafLevel(&leafSet);
    LeafGridPtMap leafPoint(&leafSet);

    // Iterate through blocks -- and set the field data
    int leafCount = 0;
    for(auto it = blocks.begin(); it < blocks.end(); ++it)
    {
      const BlockIndex& block = *it;

      // Add the hex to the mesh
      addOctreeBlock(debugMesh, block, shouldLogBlocks);

      const InOutBlockData& leafData = m_octree[block];

      int vIdx = leafData.hasData() ? m_octree.leafVertex(block, leafData)
                                    : MeshWrapper<DIM>::NO_VERTEX;

      leafVertID[leafCount] = vIdx;
      leafLevel[leafCount] = block.level();
      leafPoint[leafCount] = block.pt();

      if(hasCells)
      {
        leafVertID_unique[leafCount] = m_octree.blockIndexesVertex(vIdx, block)
          ? vIdx
          : MeshWrapper<DIM>::NO_VERTEX;

        leafTriCount[leafCount] =
          leafData.hasData() ? m_octree.leafCells(block, leafData).size() : 0;
      }

      if(hasColors)
      {
        leafColors[leafCount] = m_colorsMap.at(leafData.color());
      }

      leafCount++;
    }

    // Add the fields to the mint mesh
    axom::IndexType* vertID = addIntField(debugMesh, "vertID");
    axom::IndexType* lLevel = addIntField(debugMesh, "level");

    axom::IndexType* blockCoord[3];
    blockCoord[0] = addIntField(debugMesh, "block_x");
    blockCoord[1] = addIntField(debugMesh, "block_y");
    blockCoord[2] = addIntField(debugMesh, "block_z");

    for(int i = 0; i < leafSet.size(); ++i)
    {
      vertID[i] = leafVertID[i];
      lLevel[i] = leafLevel[i];

      blockCoord[0][i] = leafPoint[i][0];
      blockCoord[1][i] = leafPoint[i][1];
      blockCoord[2][i] = leafPoint[i][2];
    }

    if(hasCells)
    {
      axom::IndexType* uniqVertID = addIntField(debugMesh, "uniqVertID");
      axom::IndexType* triCount = addIntField(debugMesh, "triCount");

      for(int i = 0; i < leafSet.size(); ++i)
      {
        uniqVertID[i] = leafVertID_unique[i];
        triCount[i] = leafTriCount[i];
      }
    }

    if(hasColors)
    {
      axom::IndexType* colors = addIntField(debugMesh, "colors");
      for(int i = 0; i < leafSet.size(); ++i) colors[i] = leafColors[i];
    }

    mint::write_vtk(debugMesh, fName);

    delete debugMesh;
    debugMesh = nullptr;
  }

  void dumpTriangleMesh(const std::string& name,
                        const std::vector<CellIndex>& tris,
                        bool shouldLogTris = false) const
  {
    std::string fName = fmt::format("{}.vtk", name);

    DebugMesh* debugMesh = new DebugMesh(3, 3 * tris.size(), tris.size());

    for(auto it = tris.begin(); it < tris.end(); ++it)
    {
      CellIndex tIdx = *it;
      addTriangle(debugMesh, tIdx, shouldLogTris);
    }

    // Add fields to the triangle mesh
    int numTris = tris.size();

    // Index of each triangle within the mesh
    axom::IndexType* triIdx = addIntField(debugMesh, "triangle_index");

    // Indices of the three boundary vertices of this triangle
    axom::IndexType* vertIdx[3];
    vertIdx[0] = addIntField(debugMesh, "vertex_index_0");
    vertIdx[1] = addIntField(debugMesh, "vertex_index_1");
    vertIdx[2] = addIntField(debugMesh, "vertex_index_2");

    for(int i = 0; i < numTris; ++i)
    {
      CellIndex tIdx = tris[i];
      triIdx[i] = tIdx;

      CellVertIndices tv = m_octree.m_meshWrapper.cellVertexIndices(tIdx);
      vertIdx[0][i] = tv[0];
      vertIdx[1][i] = tv[1];
      vertIdx[2][i] = tv[2];
    }

    // other possible fields on triangles
    // -- number of blocks that index this triangle (blockCount)?
    // -- vertex field for number of triangles incident in the vertex (vtCount)?

    mint::write_vtk(debugMesh, fName);

    delete debugMesh;
    debugMesh = nullptr;
  }

private:
  axom::IndexType* addIntField(DebugMesh* mesh, const std::string& name) const
  {
    axom::IndexType* fld =
      mesh->createField<axom::IndexType>(name, mint::CELL_CENTERED);
    SLIC_ASSERT(fld != nullptr);
    return fld;
  }

  double* addRealField(DebugMesh* mesh, const std::string& name) const
  {
    double* fld = mesh->createField<double>(name, mint::CELL_CENTERED);
    SLIC_ASSERT(fld != nullptr);
    return fld;
  }

  void addTriangle(DebugMesh* mesh, const CellIndex& tIdx, bool shouldLogTris) const
  {
    SpaceCell triPos = m_octree.m_meshWrapper.cellPositions(tIdx);

    axom::IndexType vStart = mesh->getNumberOfNodes();
    mesh->appendNode(triPos[0][0], triPos[0][1], triPos[0][2]);
    mesh->appendNode(triPos[1][0], triPos[1][1], triPos[1][2]);
    mesh->appendNode(triPos[2][0], triPos[2][1], triPos[2][2]);

    axom::IndexType data[3];
    for(int i = 0; i < 3; ++i)
    {
      data[i] = vStart + i;
    }

    mesh->appendCell(data, mint::TRIANGLE);

    // Log the triangle info as primal code to simplify adding a test for this case
    if(shouldLogTris)
    {
      SLIC_INFO("// Triangle " << tIdx << "\n\t"
                               << "TriangleType tri("
                               << "PointType::make_point" << triPos[0] << ","
                               << "PointType::make_point" << triPos[1] << ","
                               << "PointType::make_point" << triPos[2] << ");");
    }
  }

  void addOctreeBlock(DebugMesh* mesh,
                      const BlockIndex& block,
                      bool shouldLogBlocks) const
  {
    GeometricBoundingBox blockBB = m_octree.blockBoundingBox(block);

    axom::IndexType vStart = mesh->getNumberOfNodes();

    const SpacePt& bMin = blockBB.getMin();
    const SpacePt& bMax = blockBB.getMax();

    mesh->appendNode(bMin[0], bMin[1], bMin[2]);
    mesh->appendNode(bMax[0], bMin[1], bMin[2]);
    mesh->appendNode(bMax[0], bMax[1], bMin[2]);
    mesh->appendNode(bMin[0], bMax[1], bMin[2]);

    mesh->appendNode(bMin[0], bMin[1], bMax[2]);
    mesh->appendNode(bMax[0], bMin[1], bMax[2]);
    mesh->appendNode(bMax[0], bMax[1], bMax[2]);
    mesh->appendNode(bMin[0], bMax[1], bMax[2]);

    axom::IndexType data[8];
    for(int i = 0; i < 8; ++i) data[i] = vStart + i;

    mesh->appendCell(data, mint::HEX);

    // Log bounding box info to simplify adding a test for this case
    if(shouldLogBlocks)
    {
      static int counter = 0;
      SLIC_INFO("// Block index " << block);
      SLIC_INFO("BoundingBoxType box"
                << ++counter << "(PointType::make_point" << bMin << ","
                << "PointType::make_point" << bMax << ");");
    }
  }

private:
  const InOutOctreeType& m_octree;
  typename InOutOctreeType::GenerationState m_generationState;

  ColorsMap m_colorsMap;
};

}  // namespace detail
}  // namespace quest
}  // namespace axom

#endif  // INOUT_OCTREE_MESHDUMPER__HXX_
