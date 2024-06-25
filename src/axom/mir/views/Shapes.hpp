// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MIR_VIEWS_SHAPES_HPP_
#define AXOM_MIR_VIEWS_SHAPES_HPP_

#include "axom/core/ArrayView.hpp"

#include <conduit/conduit.hpp>
#include <iostream>

/**
 NOTES: I use the shapes for 2 things:
        1. As type traits to help the topology views traverse connectivity
        2. As a shape type that contains information for a given zone so we can do a few useful things with it.

 Q: I was going to ask whether a shape's getFace() should return a shape itself such as a tri/quad shape.
    That won't work for zones like pyr, wdg that have faces with different shapes.

 Q: Is it worth templating the shape's index type?
 */

namespace axom
{
namespace mir
{
namespace views
{

// TODO: PointTraits
// TODO: LineTraits

/*
  3*-----------* 2
   |           |
   |           |
   |           |
   |           |
   |           |
  0*-----------* 1

 */
template <typename IndexT>
struct TriTraits
{
  using IndexType = IndexT;

  AXOM_HOST_DEVICE constexpr static int  id() { return 1 << 3; }
  AXOM_HOST_DEVICE constexpr static bool is_polyhedral() { return false; }
  AXOM_HOST_DEVICE constexpr static bool is_variable_size() { return false; }

  AXOM_HOST_DEVICE constexpr static IndexType dimension() { return 2; }

  AXOM_HOST_DEVICE constexpr static IndexType numberOfNodes() { return 3; }
  AXOM_HOST_DEVICE constexpr static IndexType numberOfNodesInFace(int /*faceIndex*/) { return 3; }
  AXOM_HOST_DEVICE constexpr static IndexType maxNodesInFace() { return 3; }

  AXOM_HOST_DEVICE constexpr static IndexType numberOfFaces() { return 1; }
  AXOM_HOST_DEVICE constexpr static IndexType numberOfEdges() { return 3; }
  AXOM_HOST_DEVICE constexpr static IndexType zoneOffset(int zoneIndex) { return numberOfNodes() * zoneIndex; }

  AXOM_HOST_DEVICE constexpr static IndexType faces[][4] = {{0, 1, 2, 3}};
  AXOM_HOST_DEVICE constexpr static IndexType edges[][2] = {{0,1}, {1,2}, {2,3}, {3,0}};
};

/*
  3*-----------* 2
   |           |
   |           |
   |           |
   |           |
   |           |
  0*-----------* 1

 */
template <typename IndexT>
struct QuadTraits
{
  using IndexType = IndexT;

  AXOM_HOST_DEVICE constexpr static int  id() { return 1 << 4; }
  AXOM_HOST_DEVICE constexpr static bool is_polyhedral() { return false; }
  AXOM_HOST_DEVICE constexpr static bool is_variable_size() { return false; }

  AXOM_HOST_DEVICE constexpr static IndexType dimension() { return 2; }

  AXOM_HOST_DEVICE constexpr static IndexType numberOfNodes() { return 4; }
  AXOM_HOST_DEVICE constexpr static IndexType numberOfNodesInFace(int /*faceIndex*/) { return 4; }
  AXOM_HOST_DEVICE constexpr static IndexType maxNodesInFace() { return 4; }

  AXOM_HOST_DEVICE constexpr static IndexType numberOfFaces() { return 1; }
  AXOM_HOST_DEVICE constexpr static IndexType numberOfEdges() { return 4; }
  AXOM_HOST_DEVICE constexpr static IndexType zoneOffset(int zoneIndex) { return numberOfNodes() * zoneIndex; }

  AXOM_HOST_DEVICE constexpr static IndexType faces[][4] = {{0, 1, 2, 3}};
  AXOM_HOST_DEVICE constexpr static IndexType edges[][2] = {{0,1}, {1,2}, {2,3}, {3,0}};
};

// TODO: PolygonTraits
/*
  
n-1 *-... * 2
    |     |
    |     |
  0 *-----* 1

 */
template <typename IndexT>
struct PolygonShape
{
  using IndexType = IndexT;

  AXOM_HOST_DEVICE constexpr static int  id() { return 1 << 5; }
  AXOM_HOST_DEVICE constexpr static bool is_polyhedral() { return false; }
  AXOM_HOST_DEVICE constexpr static bool is_variable_size() { return true; }

  AXOM_HOST_DEVICE constexpr static IndexType dimension() { return 2; }

  /**
   * \brief Construct a shape.
   */
  AXOM_HOST_DEVICE PolygonShape(const axom::ArrayView<IndexType> &ids) : m_ids(ids)
  {
    assert(m_ids.size() == numberOfNodes());
  }

  /**
   * \brief Get the ids that make up this shape.
   *
   * \return A view containing the ids that make up this shape.
   */
  AXOM_HOST_DEVICE axom::ArrayView<IndexType> getIds() const { return m_ids; }

  /**
   * \brief Get the ids for the requested face.
   *
   * \param faceIndex The index of the desired face.
   *
   * \return An array view (wrapping m_faceIds) that contains the ids for the face.
   */
  AXOM_HOST_DEVICE axom::ArrayView<IndexType> getFace(int faceIndex) const
  {
    return m_ids;
  }
private:
  axom::ArrayView<IndexType> m_ids;
};


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
struct TetTraits
{
  using IndexType = IndexT;

  AXOM_HOST_DEVICE constexpr static int  id() { return 1 << 6; }
  AXOM_HOST_DEVICE constexpr static bool is_polyhedral() { return false; }
  AXOM_HOST_DEVICE constexpr static bool is_variable_size() { return false; }

  AXOM_HOST_DEVICE constexpr static IndexType dimension() { return 3; }

  AXOM_HOST_DEVICE constexpr static IndexType numberOfNodes() { return 4; }
  AXOM_HOST_DEVICE constexpr static IndexType numberOfNodesInFace(int /*faceIndex*/) { return 3; }
  AXOM_HOST_DEVICE constexpr static IndexType maxNodesInFace() { return 3; }

  AXOM_HOST_DEVICE constexpr static IndexType numberOfFaces() { return 4; }
  AXOM_HOST_DEVICE constexpr static IndexType numberOfEdges() { return 6; }
  AXOM_HOST_DEVICE constexpr static IndexType zoneOffset(int zoneIndex) { return numberOfNodes() * zoneIndex; }

  AXOM_HOST_DEVICE constexpr static IndexType faces[][3] = {{0,2,1}, {0,1,3}, {1,2,3}, {2,0,3}};
  AXOM_HOST_DEVICE constexpr static IndexType edges[][2] = {{0,2}, {2,1}, {1,0}, {1,3}, {3,0}, {2,3}};
};

/*

3*-----------* 2  face 0: 3,2,1,0
 |\         /|    face 1: 0,1,4
 | \       / |    face 2: 1,2,4
 |  \     /  |    face 3: 2,3,4
 |   \   /   |    face 4: 3,0,4
 |    \ /    |
 |     * 4   |    edge 0: 0,1
 |    / \    |    edge 1: 1,2
 |   /   \   |    edge 2: 2,3
 |  /     \  |    edge 3: 3,0
 | /       \ |    edge 4: 0,4
 |/         \|    edge 5: 1,4
0*-----------* 1  edge 6: 2,4
                  edge 7: 3,4

 */
template <typename IndexT>
struct PyramidTraits
{
  using IndexType = IndexT;

  AXOM_HOST_DEVICE constexpr int  id() { return 1 << 7; }
  AXOM_HOST_DEVICE constexpr bool is_polyhedral() { return false; }
  AXOM_HOST_DEVICE constexpr static bool is_variable_size() { return false; }

  AXOM_HOST_DEVICE constexpr static IndexType dimension() { return 3; }

  AXOM_HOST_DEVICE constexpr static IndexType numberOfNodes() { return 5; }
  AXOM_HOST_DEVICE constexpr static IndexType numberOfNodesInFace(int faceIndex) { return faceIndex == 0 ? 4 : 3; }
  AXOM_HOST_DEVICE constexpr static IndexType maxNodesInFace() { return 4; }

  AXOM_HOST_DEVICE constexpr static IndexType numberOfFaces() { return 5; }
  AXOM_HOST_DEVICE constexpr static IndexType numberOfEdges() { return 8; }
  AXOM_HOST_DEVICE constexpr static IndexType zoneOffset(int zoneIndex) { return numberOfNodes() * zoneIndex; }

  AXOM_HOST_DEVICE constexpr static int faces[][4] = {{3,2,1,0}, {0,1,4,-1}, {1,2,4,-1}, {2,3,4,-1}, {3,0,4,-1}};
  AXOM_HOST_DEVICE constexpr static IndexType edges[][2] = {{0,1}, {1,2}, {2,3}, {3,0}, {0,4}, {1,4}, {2,4}, {3,4}};
};

/*

3*---------* 5  face 0: 0,2,1
 |\       /|    face 1: 3,4,5
 | \     / |    face 2: 0,1,4,3
 |  \   /  |    face 3: 1,2,5,4
 |   \ /   |    face 4: 2,0,3,5
 |    *4   |
 |    |    |    edge 0: 0,1
0*----|----* 2  edge 1: 1,2
  \   |   /     edge 2: 2,0
   \  |  /      edge 3: 3,4
    \ | /       edge 4: 4,5
     \|/        edge 5: 5,3
      * 1       edge 6: 0,3
                edge 7: 1,4
                edge 8: 2,3

 */
template <typename IndexT>
struct WedgeTraits
{
  using IndexType = IndexT;

  AXOM_HOST_DEVICE constexpr static int  id() { return 1 << 8; }

  AXOM_HOST_DEVICE constexpr static bool is_polyhedral() { return false; }
  AXOM_HOST_DEVICE constexpr static bool is_variable_size() { return false; }

  AXOM_HOST_DEVICE constexpr static IndexType dimension() { return 3; }

  AXOM_HOST_DEVICE constexpr static IndexType numberOfNodes() { return 6; }
  AXOM_HOST_DEVICE constexpr static IndexType numberOfNodesInFace(int faceIndex) { return (faceIndex == 0 || faceIndex == 1) ? 3 : 4; }
  AXOM_HOST_DEVICE constexpr static IndexType maxNodesInFace() { return 4; }

  AXOM_HOST_DEVICE constexpr static IndexType numberOfFaces() { return 5; }
  AXOM_HOST_DEVICE constexpr static IndexType numberOfEdges() { return 9; }
  AXOM_HOST_DEVICE constexpr static IndexType zoneOffset(int zoneIndex) { return numberOfNodes() * zoneIndex; }

  AXOM_HOST_DEVICE constexpr static int faces[][4] = {{0,2,1,-1}, {3,4,5,-1}, {0,1,4,3}, {1,2,5,4}, {2,0,3,5}};
  AXOM_HOST_DEVICE constexpr static IndexType edges[][2] = {{0,1}, {1,2}, {2,0}, {3,4}, {4,5}, {5,3}, {0,3}, {1,4}, {2,3}};
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
struct HexTraits
{
  using IndexType = IndexT;

  AXOM_HOST_DEVICE constexpr static int  id() { return 1 << 9; }

  AXOM_HOST_DEVICE constexpr static bool is_polyhedral() { return false; }
  AXOM_HOST_DEVICE constexpr static bool is_variable_size() { return false; }

  AXOM_HOST_DEVICE constexpr static IndexType dimension() { return 3; }

  AXOM_HOST_DEVICE constexpr static IndexType numberOfNodes() { return 8; }
  AXOM_HOST_DEVICE constexpr static IndexType numberOfNodesInFace(int /*faceIndex*/) { return 4; }
  AXOM_HOST_DEVICE constexpr static IndexType maxNodesInFace() { return 4; }

  AXOM_HOST_DEVICE constexpr static IndexType numberOfFaces() { return 6; }
  AXOM_HOST_DEVICE constexpr static IndexType numberOfEdges() { return 12; }
  AXOM_HOST_DEVICE constexpr static IndexType zoneOffset(int zoneIndex) { return numberOfNodes() * zoneIndex; }

  AXOM_HOST_DEVICE constexpr static IndexType faces[][4] = {
    {0, 3, 2, 1}, {0, 1, 5, 4}, {1, 2, 6, 5}, {2, 3, 7, 6}, {3, 0, 4, 7}, {4, 5, 6, 7}};
  AXOM_HOST_DEVICE constexpr static IndexType edges[][2] = {{0,3}, {3,2}, {2,1}, {1,0}, {1,5}, {5,4}, {4,0}, {2,6}, {6,5}, {3,7}, {7,6}, {7,4}};
};

/**
 * \brief This class extends the ShapeTraits with object state so it can represent a zone.
 */
template <typename IndexT, typename ShapeTraits>
struct Shape : public ShapeTraits<IndexT>
{
  /**
   * \brief Construct a shape.
   */
  AXOM_HOST_DEVICE Shape(const axom::ArrayView<IndexType> &ids) : m_ids(ids), m_faceIds()
  {
    assert(m_ids.size() == numberOfNodes());
  }

  /**
   * \brief Get the ids that make up this shape.
   *
   * \return A view containing the ids that make up this shape.
   */
  AXOM_HOST_DEVICE axom::ArrayView<IndexType> getIds() const { return m_ids; }

  /**
   * \brief Get the ids for the requested face.
   *
   * \param faceIndex The index of the desired face.
   *
   * \return An array view (wrapping m_faceIds) that contains the ids for the face.
   */
  AXOM_HOST_DEVICE std::enable_if<dimension() == 2, axom::ArrayView<IndexType>>
  getFace(int faceIndex) const
  {
    return m_ids;
  }

  AXOM_HOST_DEVICE std::enable_if<dimension() == 3, axom::ArrayView<IndexType>>
  getFace(int faceIndex) const
  {
    const auto nnodes = numberOfNodesInFace(faceIndex);
    for(IndexType i = 0; i < nnodes; i++)
      m_faceIds[i] = m_ids[faces[faceIndex][i]];
    return axom::ArrayView<IndexType>(m_faceIds.m_data, nnodes);
  }
private:
  axom::ArrayView<IndexType> m_ids;
  mutable axom::StackArray<IndexType, maxNodesInFace()> m_faceIds;
};

// Make some concrete shape classes based on the shape traits.
template <typename IndexT>
using TriShape = Shape<IndexT, TriTraits<IndexT>>;

template <typename IndexT>
using QuadShape = Shape<IndexT, QuadTraits<IndexT>>;

template <typename IndexT>
using TetShape = Shape<IndexT, TetTraits<IndexT>>;

template <typename IndexT>
using PyramidShape = Shape<IndexT, PyramidTraits<IndexT>>;

template <typename IndexT>
using WedgeShape = Shape<IndexT, WedgeTraits<IndexT>>;

template <typename IndexT>
using HexShape = Shape<IndexT, HexTraits<IndexT>>;


} // end namespace views
} // end namespace mir
} // end namespace axom

#endif
