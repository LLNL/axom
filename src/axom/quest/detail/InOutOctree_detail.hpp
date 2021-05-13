// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file InOutOctree_detail.hpp
 *
 * \brief Defines helper classes for the InOutOctree.
 */

#ifndef INOUT_OCTREE_DETAIL__HXX_
#define INOUT_OCTREE_DETAIL__HXX_

#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/slam.hpp"
#include "axom/primal.hpp"
#include "axom/mint.hpp"
#include "axom/spin.hpp"

#include "inout/BlockData.hpp"
#include "inout/MeshWrapper.hpp"

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

template <int DIM>
class InOutOctreeValidator
{
public:
  using InOutOctreeType = InOutOctree<DIM>;
  using CellIndexSet = typename InOutOctreeType::CellIndexSet;

  using OctreeBaseType = typename InOutOctreeType::OctreeBaseType;
  using OctreeLevels = typename OctreeBaseType::OctreeLevels;
  using BlockIndex = typename OctreeBaseType::BlockIndex;

  using SpacePt = typename InOutOctreeType::SpacePt;
  using VertexIndex = typename InOutOctreeType::VertexIndex;
  using CellIndex = typename InOutOctreeType::CellIndex;
  using CellVertIndices = typename MeshWrapper<DIM>::CellVertIndices;
  using GeometricBoundingBox = typename InOutOctreeType::GeometricBoundingBox;

public:
  InOutOctreeValidator(const InOutOctreeType& octree)
    : m_octree(octree)
    , m_generationState(m_octree.m_generationState)
  { }

  void checkAllLeavesColored() const
  {
    SLIC_DEBUG(
      "--Checking that all leaves have a color -- black, white and gray");

    for(int lev = 0; lev < m_octree.m_levels.size(); ++lev)
    {
      checkAllLeavesColoredAtLevel(lev);
    }
  }

  void checkAllLeavesColoredAtLevel(int level) const
  {
    const auto& levelLeafMap = m_octree.getOctreeLevel(level);
    auto itEnd = levelLeafMap.end();
    for(auto it = levelLeafMap.begin(); it != itEnd; ++it)
    {
      if(!it->isLeaf()) continue;

      SLIC_ASSERT_MSG(
        it->isColored(),
        fmt::format(
          "Error after coloring level {} leaf block {} was not colored.",
          level,
          BlockIndex(it.pt(), level)));
    }
  }

  void checkEachVertexIsIndexed() const
  {
    SLIC_DEBUG("--Checking that each vertex is in a leaf block of the tree.");

    const VertexIndex numVertices = m_octree.m_meshWrapper.numMeshVertices();
    for(VertexIndex i = 0; i < numVertices; ++i)
    {
      const SpacePt& pos = m_octree.m_meshWrapper.vertexPosition(i);
      BlockIndex vertBlock = m_octree.findLeafBlock(pos);
      const InOutBlockData& leafData = m_octree[vertBlock];

      VertexIndex vertInBlock = m_octree.leafVertex(vertBlock, leafData);
      AXOM_DEBUG_VAR(vertInBlock);

      // Check that we can find the leaf block indexing each vertex from its position
      SLIC_ASSERT_MSG(
        leafData.hasData() && vertInBlock == i,
        fmt::format(
          " Vertex {} at position {} was not indexed in block {} with bounding "
          "box {} (point {} contained in block)."
          " \n\n *** \n Leaf data: {} \n ***",
          i,
          pos,
          vertBlock,
          m_octree.blockBoundingBox(vertBlock),
          (m_octree.blockBoundingBox(vertBlock).contains(pos) ? "is" : "is NOT"),
          leafData));

      // Check that our cached value of the vertex's block is accurate
      SLIC_ASSERT_MSG(
        vertBlock == m_octree.m_vertexToBlockMap[i],
        fmt::format(
          "Cached block for vertex {} differs from its found block. "
          "\n\t -- cached block {} -- is leaf? {}"
          "\n\t -- actual block {} -- is leaf? {}"
          "\n\t -- vertex set size: {}",
          i,
          m_octree.m_vertexToBlockMap[i],
          (m_octree[m_octree.m_vertexToBlockMap[i]].isLeaf() ? "yes" : "no"),
          vertBlock,
          (m_octree[vertBlock].isLeaf() ? "yes" : "no"),
          numVertices));
    }
  }

  void checkTrianglesReferencedInBoundaryVertexBlocks() const
  {
    SLIC_DEBUG(
      "--Checking that each triangle is referenced by the leaf blocks "
      "containing its vertices.");

    const axom::IndexType numTriangles = m_octree.m_meshWrapper.numMeshCells();
    for(axom::IndexType tIdx = 0; tIdx < numTriangles; ++tIdx)
    {
      CellVertIndices tvRel = m_octree.m_meshWrapper.cellVertexIndices(tIdx);
      for(int j = 0; j < tvRel.size(); ++j)
      {
        VertexIndex vIdx = tvRel[j];
        BlockIndex vertBlock = m_octree.m_vertexToBlockMap[vIdx];
        const InOutBlockData& leafData = m_octree[vertBlock];

        // Check that this triangle is referenced here.
        bool foundTriangle = false;
        CellIndexSet leafTris = m_octree.leafCells(vertBlock, leafData);
        for(int k = 0; !foundTriangle && k < leafTris.size(); ++k)
        {
          if(leafTris[k] == tIdx) foundTriangle = true;
        }

        SLIC_ASSERT_MSG(
          foundTriangle,
          fmt::format("Did not find triangle with index {} and vertices [{}] "
                      "in block {} containing vertex {}"
                      " \n\n *** "
                      "\n Leaf data: {} "
                      "\n\t Triangles in block? [{}] "
                      "\n ***",
                      tIdx,
                      fmt::join(tvRel, ", "),
                      vertBlock,
                      vIdx,
                      leafData,
                      fmt::join(leafTris.begin(), leafTris.end(), ", ")));
      }
    }
  }

  void checkBlockIndexingConsistency() const
  {
    // Check that internal blocks have no triangle / vertex
    //       and leaf blocks satisfy the conditions above
    SLIC_DEBUG(
      "--Checking that internal blocks have no data, and that leaves satisfy "
      "all PM conditions");

    const double bb_scale_factor = m_octree.m_boundingBoxScaleFactor;

    for(int lev = 0; lev < m_octree.m_levels.size(); ++lev)
    {
      const auto& levelLeafMap = m_octree.getOctreeLevel(lev);
      auto itEnd = levelLeafMap.end();
      for(auto it = levelLeafMap.begin(); it != itEnd; ++it)
      {
        const BlockIndex block(it.pt(), lev);
        const InOutBlockData& data = *it;

        if(!data.isLeaf())
        {
          SLIC_ASSERT(!data.hasData());
        }
        else  // leaf block
        {
          if(data.hasData())
          {
            VertexIndex vIdx = m_octree.leafVertex(block, data);
            CellIndexSet triSet = m_octree.leafCells(block, data);
            for(int i = 0; i < triSet.size(); ++i)
            {
              CellIndex tIdx = triSet[i];

              // Check that vIdx is one of this triangle's vertices
              CellVertIndices tvRel =
                m_octree.m_meshWrapper.cellVertexIndices(tIdx);

              SLIC_ASSERT_MSG(
                m_octree.m_meshWrapper.incidentInVertex(tvRel, vIdx),
                fmt::format("All triangles in a gray block must be incident "
                            " in a common vertex, but triangle {} with "
                            "vertices [{}] in block {}"
                            " is not incident in vertex {}",
                            tIdx,
                            fmt::join(tvRel, ", "),
                            block,
                            vIdx));

              // Check that this triangle intersects the bounding box of the
              // block
              GeometricBoundingBox blockBB = m_octree.blockBoundingBox(block);
              blockBB.expand(bb_scale_factor);
              SLIC_ASSERT_MSG(
                m_octree.blockIndexesElementVertex(tIdx, block) ||
                  intersect(m_octree.m_meshWrapper.cellPositions(tIdx), blockBB),
                fmt::format("Triangle {} was indexed in block {}"
                            " but it does not intersect the block."
                            "\n\tBlock bounding box: {}"
                            "\n\tTriangle positions: {}"
                            "\n\tTriangle vertex indices: [{}]"
                            "\n\tLeaf vertex is: {}"
                            "\n\tLeaf triangles: {} ({})",
                            tIdx,
                            block,
                            blockBB,
                            m_octree.m_meshWrapper.cellPositions(tIdx),
                            fmt::join(tvRel, ", "),
                            vIdx,
                            fmt::join(triSet, ", "),
                            triSet.size()));
            }
          }
        }
      }
    }
  }

  void checkNeighboringBlockColors() const
  {
    SLIC_DEBUG("--Checking that inside blocks do not neighbor outside blocks");
    for(int lev = m_octree.maxLeafLevel() - 1; lev >= 0; --lev)
    {
      const auto& levelLeafMap = m_octree.getOctreeLevel(lev);
      auto itEnd = levelLeafMap.end();
      for(auto it = levelLeafMap.begin(); it != itEnd; ++it)
      {
        const BlockIndex block(it.pt(), lev);
        const InOutBlockData& data = *it;

        if(data.isLeaf() && data.color() != InOutBlockData::Gray)
        {
          // Leaf does not yet have a color... try to find its color from
          // same-level face neighbors
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
                SLIC_CHECK_MSG(data.color() == neighborData.color(),
                               fmt::format("Problem at block {}  with data {} "
                                           "--- neighbor {} with data {}."
                                           " Neighboring leaves that are not "
                                           "gray must have the same color.",
                                           block,
                                           data,
                                           neighborBlk,
                                           neighborData));
                break;
              case InOutBlockData::Gray:  // intentional fallthrough
              case InOutBlockData::Undetermined:
                break;
              }
            }
          }
        }
      }
    }
  }

  void checkValid() const
  {
    // We are assumed to be valid before we insert the vertices
    if(m_generationState < InOutOctreeType::INOUTOCTREE_VERTICES_INSERTED)
      return;

    // Iterate through the tree
    // Internal blocks should not have associated vertex data
    // Leaf block consistency depends on 'color'
    //      Black or White blocks should not have vertex data
    //      Gray blocks should have a vertex reference; it may or may not be
    // located within the block
    //          The sum of vertices located within a block should equal the
    // number of mesh vertices.
    //      Gray blocks should have one or more triangles.
    //          All triangles should be incident in a common vertex -- which
    // equals the indexed vertex, if it exists.

    if(m_generationState > InOutOctreeType::INOUTOCTREE_VERTICES_INSERTED)
    {
      checkEachVertexIsIndexed();
    }

    if(m_generationState >= InOutOctreeType::INOUTOCTREE_ELEMENTS_INSERTED)
    {
      checkTrianglesReferencedInBoundaryVertexBlocks();
      checkBlockIndexingConsistency();
    }

    if(m_generationState >= InOutOctreeType::INOUTOCTREE_LEAVES_COLORED)
    {
      checkAllLeavesColored();
      checkNeighboringBlockColors();
    }
  }

private:
  const InOutOctreeType& m_octree;
  typename InOutOctreeType::GenerationState m_generationState;
};

template <int DIM>
class InOutOctreeStats
{
public:
  using InOutOctreeType = InOutOctree<DIM>;
  using CellIndexSet = typename InOutOctreeType::CellIndexSet;

  using OctreeBaseType = typename InOutOctreeType::OctreeBaseType;
  using OctreeLevels = typename OctreeBaseType::OctreeLevels;
  using BlockIndex = typename OctreeBaseType::BlockIndex;

  using LeafCountMap = slam::Map<slam::Set<>, int>;
  using TriCountMap = slam::Map<slam::Set<>, int>;
  using CardinalityVTMap = slam::Map<slam::Set<>, int>;

  using LogHistogram = std::map<int, int>;
  using MinMaxRange = primal::BoundingBox<double, 1>;
  using LengthType = MinMaxRange::PointType;
  using LogRangeMap = std::map<int, MinMaxRange>;

  /** A simple struct to track totals within the octree levels */
  struct Totals
  {
    /** Default constructor to set everything to 0 */
    Totals()
      : blocks(0)
      , leaves(0)
      , leavesWithVert(0)
      , triangleRefCount(0)
      , whiteBlocks(0)
      , blackBlocks(0)
      , grayBlocks(0)
    { }

    int blocks;
    int leaves;
    int leavesWithVert;
    int triangleRefCount;
    int whiteBlocks;
    int blackBlocks;
    int grayBlocks;
  };

public:
  InOutOctreeStats(const InOutOctreeType& octree)
    : m_octree(octree)
    , m_generationState(m_octree.m_generationState)
    , m_levelBlocks(&m_octree.m_levels)
    , m_levelLeaves(&m_octree.m_levels)
    , m_levelLeavesWithVert(&m_octree.m_levels)
    , m_levelTriangleRefCount(&m_octree.m_levels)
    , m_levelWhiteBlockCount(&m_octree.m_levels)
    , m_levelBlackBlockCount(&m_octree.m_levels)
    , m_levelGrayBlockCount(&m_octree.m_levels)
  {
    if(m_generationState >= InOutOctreeType::INOUTOCTREE_ELEMENTS_INSERTED)
    {
      m_triCount = TriCountMap(&m_octree.m_meshWrapper.elementSet());
    }

    // Iterate through blocks -- count the numbers of internal and leaf blocks
    for(int lev = 0; lev < m_octree.m_levels.size(); ++lev)
    {
      const auto& levelLeafMap = m_octree.getOctreeLevel(lev);

      m_levelBlocks[lev] = levelLeafMap.numBlocks();
      m_levelLeaves[lev] = levelLeafMap.numLeafBlocks();
      m_levelLeavesWithVert[lev] = 0;
      m_levelTriangleRefCount[lev] = 0;
      m_levelWhiteBlockCount[lev] = 0;
      m_levelBlackBlockCount[lev] = 0;
      m_levelGrayBlockCount[lev] = 0;

      auto itEnd = levelLeafMap.end();
      for(auto it = levelLeafMap.begin(); it != itEnd; ++it)
      {
        const InOutBlockData& blockData = *it;
        BlockIndex block(it.pt(), lev);

        if(blockData.isLeaf())
        {
          if(blockData.hasData())
          {
            ++m_levelLeavesWithVert[lev];

            if(m_generationState >= InOutOctreeType::INOUTOCTREE_ELEMENTS_INSERTED)
            {
              m_levelTriangleRefCount[lev] +=
                m_octree.leafCells(block, blockData).size();

              BlockIndex blk(it.pt(), lev);
              CellIndexSet tris = m_octree.leafCells(blk, blockData);
              for(int i = 0; i < tris.size(); ++i)
              {
                ++m_triCount[tris[i]];
              }
            }
          }

          if(m_generationState >= InOutOctreeType::INOUTOCTREE_LEAVES_COLORED)
          {
            switch(blockData.color())
            {
            case InOutBlockData::Black:
              ++m_levelBlackBlockCount[lev];
              break;
            case InOutBlockData::White:
              ++m_levelWhiteBlockCount[lev];
              break;
            case InOutBlockData::Gray:
              ++m_levelGrayBlockCount[lev];
              break;
            case InOutBlockData::Undetermined:
              break;
            }
          }
        }
      }

      m_totals.blocks += m_levelBlocks[lev];
      m_totals.leaves += m_levelLeaves[lev];
      m_totals.leavesWithVert += m_levelLeavesWithVert[lev];
      m_totals.triangleRefCount += m_levelTriangleRefCount[lev];
      m_totals.whiteBlocks += m_levelWhiteBlockCount[lev];
      m_totals.blackBlocks += m_levelBlackBlockCount[lev];
      m_totals.grayBlocks += m_levelGrayBlockCount[lev];
    }
  }

  /** Generates a string summarizing information about the leaves and blocks of
     the octree */
  std::string blockDataStats() const
  {
    std::stringstream sstr;

    for(int lev = 0; lev < m_octree.m_levels.size(); ++lev)
    {
      if(m_levelBlocks[lev] > 0)
      {
        sstr << fmt::format(
          "\t Level {} has {} blocks -- {} internal; {} leaves ({}% w/ vert);",
          lev,
          m_levelBlocks[lev],
          m_levelBlocks[lev] - m_levelLeaves[lev],
          m_levelLeaves[lev],
          integerPercentage(m_levelLeavesWithVert[lev], m_levelLeaves[lev]));

        if(m_generationState >= InOutOctreeType::INOUTOCTREE_LEAVES_COLORED)
        {
          sstr << fmt::format(
            " Leaf counts: {} Black, {} White, {} Gray w/ {} triangle refs.",
            m_levelBlackBlockCount[lev],
            m_levelWhiteBlockCount[lev],
            m_levelGrayBlockCount[lev],
            m_levelTriangleRefCount[lev]);
        }
        //sstr <<"Hash load factor: "
        //     << this->m_leavesLevelMap[ lev ].load_factor()
        //     << " -- max lf: " << this->m_leavesLevelMap[ lev
        // ].max_load_factor();
        sstr << "\n";
      }
    }

    return sstr.str();
  }

  /** Generates a string summarizing information about the mesh elements indexed
     by the octree */
  std::string meshDataStats() const
  {
    std::stringstream sstr;

    double meshNumTriangles = m_octree.m_meshWrapper.numMeshCells();

    sstr << fmt::format(
      "  Mesh has {} vertices."
      "\n  Octree has {} blocks; {} internal; {} leaves ({}% w/ vert); ",
      meshNumTriangles,
      m_totals.blocks,
      m_totals.blocks - m_totals.leaves,
      m_totals.leaves,
      integerPercentage(m_totals.leavesWithVert, m_totals.leaves));

    if(m_generationState >= InOutOctreeType::INOUTOCTREE_ELEMENTS_INSERTED)
    {
      sstr << fmt::format(
        " \n\t There were {} triangle references "
        " (avg. {} refs per triangle).",
        m_totals.triangleRefCount,
        (m_totals.triangleRefCount / meshNumTriangles));
    }

    return sstr.str();
  }

  std::string triangleCountHistogram() const
  {
    std::stringstream sstr;

    // Generate and output a histogram of the bucket counts on a lg-scale
    LogHistogram triCountHist;  // Create histogram of edge lengths (log scale)
    LogRangeMap triCountRange;

    int numElems = m_octree.m_meshWrapper.numMeshCells();

    for(int i = 0; i < numElems; ++i)
    {
      LengthType count(m_triCount[i]);
      int expBase2;
      std::frexp(m_triCount[i], &expBase2);
      triCountHist[expBase2]++;
      triCountRange[expBase2].addPoint(count);
    }

    std::stringstream triCountStr;
    triCountStr << "\tTriangle index count "
                << "(lg-arithmic bins for number of references per triangle):";
    for(auto it = triCountHist.begin(); it != triCountHist.end(); ++it)
    {
      triCountStr << fmt::format("\n\t exp: {}\t count: {}\tRange: {}",
                                 it->first,
                                 it->second,
                                 triCountRange[it->first]);
    }

    return triCountStr.str();
  }

  std::string vertexCardinalityHistogram() const
  {
    std::stringstream sstr;

    using CellVertIndices = typename MeshWrapper<DIM>::CellVertIndices;

    // Generate and output histogram of VT relation
    CardinalityVTMap cardVT(&m_octree.m_meshWrapper.vertexSet());

    int numElems = m_octree.m_meshWrapper.numMeshCells();
    for(int i = 0; i < numElems; ++i)
    {
      CellVertIndices tvRel = m_octree.m_meshWrapper.cellVertexIndices(i);
      cardVT[tvRel[0]]++;
      cardVT[tvRel[1]]++;
      cardVT[tvRel[2]]++;
    }

    using LinHistogram = std::map<int, int>;
    LinHistogram vtCardHist;
    int numVerts = m_octree.m_meshWrapper.numMeshVertices();
    for(int i = 0; i < numVerts; ++i)
    {
      LengthType count(cardVT[i]);
      vtCardHist[cardVT[i]]++;
    }

    sstr << "\tCardinality VT relation histogram (linear): ";
    for(auto it = vtCardHist.begin(); it != vtCardHist.end(); ++it)
    {
      sstr << fmt::format("\n\t exp: {}\t count: {}", it->first, it->second);
    }

    return sstr.str();
  }

  std::string summaryStats() const
  {
    std::stringstream octreeStatsStr;

    octreeStatsStr << fmt::format(
      "*** {} octree summary *** \n",
      (m_generationState >= InOutOctreeType::INOUTOCTREE_ELEMENTS_INSERTED
         ? "PM"
         : "PR"));

    octreeStatsStr << blockDataStats() << "\n" << meshDataStats();

    if(m_generationState >= InOutOctreeType::INOUTOCTREE_LEAVES_COLORED)
    {
      octreeStatsStr << fmt::format("\n\tColors: {} Black, {} White, {} Gray",
                                    m_totals.blackBlocks,
                                    m_totals.whiteBlocks,
                                    m_totals.grayBlocks);
    }

    if(m_generationState >= InOutOctreeType::INOUTOCTREE_ELEMENTS_INSERTED)
    {
      octreeStatsStr << "\n"
                     << triangleCountHistogram() << "\n"
                     << vertexCardinalityHistogram();
    }

    return octreeStatsStr.str();
  }

private:
  int integerPercentage(double val, double size) const
  {
    return (size > 0) ? static_cast<int>((100. * val) / size) : 0;
  }

private:
  const InOutOctreeType& m_octree;
  typename InOutOctreeType::GenerationState m_generationState;

  LeafCountMap m_levelBlocks;
  LeafCountMap m_levelLeaves;
  LeafCountMap m_levelLeavesWithVert;
  LeafCountMap m_levelTriangleRefCount;

  LeafCountMap m_levelWhiteBlockCount;
  LeafCountMap m_levelBlackBlockCount;
  LeafCountMap m_levelGrayBlockCount;

  TriCountMap m_triCount;

  Totals m_totals;
};

}  // end namespace detail
}  // namespace quest
}  // namespace axom

#endif  // INOUT_OCTREE_DETAIL__HXX_
