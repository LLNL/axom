// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MIR_VIEWS_SHAPES_HPP_
#define AXOM_MIR_VIEWS_SHAPES_HPP_

#include "axom/core/AxomArrayView.hpp"

#include <conduit/conduit.hpp>
#include <iostream>

namespace axom
{
namespace mir
{
namespace views
{

/*
      3
      *
     /|\         face 0: 0,2,1
    / | \        face 1: 0,1,3
   /  |  \       face 2: 1,2,3
  /   |   \      face 3: 2,0,3
0*----|----* 2
  \   |   /      edge 0: 0,2
   \  |  /       edge 1: 2,1
    \ | /        edge 2: 1,0
     \|/         edge 3: 1,3
      *          edge 4: 3,0
       1         edge 5: 2,3

 */
template <typename IndexT>
struct TetShape
{
  using IndexType = IndexT;

  AXOM_HOST_DEVICE constexpr bool is_polyhedral() { return false; }

  AXOM_HOST_DEVICE constexpr static IndexType numberOfNodes() { return 4; }
  AXOM_HOST_DEVICE constexpr static IndexType numberOfNodesInFace(int /*faceIndex*/) { return 3; }
  AXOM_HOST_DEVICE constexpr static IndexType numberOfFaces() { return 4; }

  AXOM_HOST_DEVICE constexpr static IndexType numberOfEdges() { return 6; }
  AXOM_HOST_DEVICE constexpr static IndexType zoneOffset(int zoneIndex) { return numberOfNodes() * zoneIndex; }

  AXOM_HOST_DEVICE constexpr static IndexType faces[][3] = {{0,2,1}, {0,1,3}, {1,2,3}, {2,0,3}};
  AXOM_HOST_DEVICE constexpr static IndexType edges[][2] = {{0,2}, {2,1}, {1,0}, {1,3}, {3,0}, {2,3}};

  AXOM_HOST_DEVICE TetShape(const axom::ArrayView<IndexType> &ids) : m_ids(ids), m_faceIds()
  {
  }

  AXOM_HOST_DEVICE axom::ArrayView<IndexType> getIds() const { return m_ids; }

  AXOM_HOST_DEVICE axom::ArrayView<IndexType> getFace(int faceIndex) const
  {
    m_faceIds[0] = m_ids[faces[faceIndex][0]];
    m_faceIds[1] = m_ids[faces[faceIndex][1]];
    m_faceIds[2] = m_ids[faces[faceIndex][2]];
    return axom::ArrayView<IndexType>(m_faceIds.m_data, 3);
  }
private:
  axom::ArrayView<IndexType> m_ids;
  mutable axom::StackArray<IndexType, 3> m_faceIds;
};

/*
      3*------------* 2   
      /|           /|     
     / |          / |
    /  |         /  |
  7*------------*6  |
   |   |        |   |
   |   |        |   |
   |  0*--------|---* 1
   |  /         |  /
   | /          | /
   |/           |/
   *------------*
   4            5

 */
template <typename IndexT>
struct HexShape
{
  using IndexType = IndexT;

  AXOM_HOST_DEVICE constexpr bool is_polyhedral() { return false; }

  AXOM_HOST_DEVICE constexpr static IndexType numberOfNodes() { return 8; }
  AXOM_HOST_DEVICE constexpr static IndexType numberOfNodesInFace(int /*faceIndex*/) { return 4; }

  AXOM_HOST_DEVICE constexpr static IndexType numberOfFaces() { return 6; }
  AXOM_HOST_DEVICE constexpr static IndexType numberOfEdges() { return 12; }
  AXOM_HOST_DEVICE constexpr static IndexType zoneOffset(int zoneIndex) { return numberOfNodes() * zoneIndex; }

  AXOM_HOST_DEVICE constexpr static IndexType faces[][4] = {
    {0, 3, 2, 1}, {0, 1, 5, 4}, {1, 2, 6, 5}, {2, 3, 7, 6}, {3, 0, 4, 7}, {4, 5, 6, 7}};
  AXOM_HOST_DEVICE constexpr static IndexType edges[][2] = {{0,3}, {3,2}, {2,1}, {1,0}, {1,5}, {5,4}, {4,0}, {2,6}, {6,5}, {3,7}, {7,6}, {7,4}};


  AXOM_HOST_DEVICE HexShape(const axom::ArrayView<IndexType> &ids) : m_ids(ids), m_faceIds()
  {
  }

  AXOM_HOST_DEVICE axom::ArrayView<IndexType> getIds() const { return m_ids; }

  AXOM_HOST_DEVICE axom::ArrayView<IndexType> getFace(int faceIndex) const
  {
    m_faceIds[0] = m_ids[faces[faceIndex][0]];
    m_faceIds[1] = m_ids[faces[faceIndex][1]];
    m_faceIds[2] = m_ids[faces[faceIndex][2]];
    m_faceIds[3] = m_ids[faces[faceIndex][3]];
    return axom::ArrayView<IndexType>(m_faceIds.m_data, 4);
  }
private:
  axom::ArrayView<IndexType> m_ids;
  mutable axom::StackArray<IndexType, 4> m_faceIds;
};

} // end namespace views
} // end namespace mir
} // end namespace axom

#endif
