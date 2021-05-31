// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file InOutOctreeValidator.hpp
 *
 * \brief Defines helper class to validate an InOutOctree instance
 */

#ifndef INOUT_OCTREE_VALIDATOR__HXX_
#define INOUT_OCTREE_VALIDATOR__HXX_

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

}  // namespace detail
}  // namespace quest
}  // namespace axom

#endif  // INOUT_OCTREE_VALIDATOR__HXX_
