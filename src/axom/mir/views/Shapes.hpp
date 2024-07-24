// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MIR_VIEWS_SHAPES_HPP_
#define AXOM_MIR_VIEWS_SHAPES_HPP_

#include "axom/core/ArrayView.hpp"

#include <conduit/conduit.hpp>
#include <conduit/conduit_blueprint_mesh_utils.hpp>

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

/*

  0*-----------* 1

 */
template <typename IndexT>
struct LineTraits
{
  using IndexType = IndexT;

  AXOM_HOST_DEVICE constexpr static int  id() { return 1 << 2; }
  AXOM_HOST_DEVICE constexpr static bool is_polyhedral() { return false; }
  AXOM_HOST_DEVICE constexpr static bool is_variable_size() { return false; }

  AXOM_HOST_DEVICE constexpr static IndexType dimension() { return 1; }

  AXOM_HOST_DEVICE constexpr static IndexType numberOfNodes() { return 2; }
  AXOM_HOST_DEVICE constexpr static IndexType numberOfNodesInFace(int /*faceIndex*/) { return 2; }
  AXOM_HOST_DEVICE constexpr static IndexType maxNodesInFace() { return 2; }

  AXOM_HOST_DEVICE constexpr static IndexType numberOfFaces() { return 1; }
  AXOM_HOST_DEVICE constexpr static IndexType numberOfEdges() { return 1; }
  AXOM_HOST_DEVICE constexpr static IndexType zoneOffset(int zoneIndex) { return numberOfNodes() * zoneIndex; }

  AXOM_HOST_DEVICE constexpr static IndexType faces[][2] = {{0, 1}};
  AXOM_HOST_DEVICE constexpr static axom::StackArray<IndexType, 2> getEdge(int edgeIndex)
  {
    return axom::StackArray<IndexType, 2>{0,1};
  }

  AXOM_HOST_DEVICE constexpr static const char *name() { return "line"; }
};

/*
  2*
   |\
   | \
   |  \
   |   \
   |    \
  0*-----* 1

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

  AXOM_HOST_DEVICE constexpr static IndexType faces[][3] = {{0, 1, 2}};
  AXOM_HOST_DEVICE constexpr static axom::StackArray<IndexType, 2> getEdge(int edgeIndex)
  {
    const axom::StackArray<IndexType, 2> edges[] = {{0,1}, {1,2}, {2,0}};
    return edges[edgeIndex];
  }

  AXOM_HOST_DEVICE constexpr static const char *name() { return "tri"; }
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

  AXOM_HOST_DEVICE constexpr static axom::StackArray<IndexType, 2> getEdge(int edgeIndex)
  {
    const axom::StackArray<IndexType, 2> edges[] = {{0,1}, {1,2}, {2,3}, {3,0}};
    return edges[edgeIndex];
  }

  AXOM_HOST_DEVICE constexpr static const char *name() { return "quad"; }
};

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
  }

  AXOM_HOST_DEVICE constexpr static IndexType numberOfFaces() { return 1; }

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
  AXOM_HOST_DEVICE axom::ArrayView<IndexType> getFace(int /*faceIndex*/) const
  {
    return m_ids;
  }

  AXOM_HOST_DEVICE axom::StackArray<IndexType, 2> getEdge(int edgeIndex) const
  {
    const auto p0 = edgeIndex % m_ids.size();
    const auto p1 = (edgeIndex + 1) % m_ids.size();
    return axom::StackArray<IndexType, 2>{p0, p1};
  }

  AXOM_HOST_DEVICE constexpr static const char *name() { return "polygon"; }

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
  \   |   /      edge 0: 0,1
   \  |  /       edge 1: 1,2
    \ | /        edge 2: 2,0
     \|/         edge 3: 0,3
      *          edge 4: 1,3
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

  AXOM_HOST_DEVICE constexpr static axom::StackArray<IndexType, 2> getEdge(int edgeIndex)
  {
    const axom::StackArray<IndexType, 2> edges[] = {{0, 1}, {1, 2}, {2, 0},{0, 3}, {1, 3}, {2, 3}};
    return edges[edgeIndex];
  }

  AXOM_HOST_DEVICE constexpr static const char *name() { return "tet"; }
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

  AXOM_HOST_DEVICE constexpr static int id() { return 1 << 7; }
  AXOM_HOST_DEVICE constexpr static bool is_polyhedral() { return false; }
  AXOM_HOST_DEVICE constexpr static bool is_variable_size() { return false; }

  AXOM_HOST_DEVICE constexpr static IndexType dimension() { return 3; }

  AXOM_HOST_DEVICE constexpr static IndexType numberOfNodes() { return 5; }
  AXOM_HOST_DEVICE constexpr static IndexType numberOfNodesInFace(int faceIndex) { return faceIndex == 0 ? 4 : 3; }
  AXOM_HOST_DEVICE constexpr static IndexType maxNodesInFace() { return 4; }

  AXOM_HOST_DEVICE constexpr static IndexType numberOfFaces() { return 5; }
  AXOM_HOST_DEVICE constexpr static IndexType numberOfEdges() { return 8; }
  AXOM_HOST_DEVICE constexpr static IndexType zoneOffset(int zoneIndex) { return numberOfNodes() * zoneIndex; }

  AXOM_HOST_DEVICE constexpr static int faces[][4] = {{3,2,1,0}, {0,1,4,-1}, {1,2,4,-1}, {2,3,4,-1}, {3,0,4,-1}};

  AXOM_HOST_DEVICE constexpr static axom::StackArray<IndexType, 2> getEdge(int edgeIndex)
  {
    const axom::StackArray<IndexType, 2> edges[] = {{0,1}, {1,2}, {2,3}, {3,0}, {0,4}, {1,4}, {2,4}, {3,4}};
    return edges[edgeIndex];
  }

  AXOM_HOST_DEVICE constexpr static const char *name() { return "pyramid"; }
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

  AXOM_HOST_DEVICE constexpr static axom::StackArray<IndexType, 2> getEdge(int edgeIndex)
  {
    const axom::StackArray<IndexType, 2> edges[] = {{0,1}, {1,2}, {2,0}, {3,4}, {4,5}, {5,3}, {0,3}, {1,4}, {2,3}};
    return edges[edgeIndex];
  }

  AXOM_HOST_DEVICE constexpr static const char *name() { return "wedge"; }
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

  AXOM_HOST_DEVICE constexpr static axom::StackArray<IndexType, 2> getEdge(int edgeIndex)
  {
    const axom::StackArray<IndexType, 2> edges[] = {{0, 1}, {1, 2}, {3, 2}, {3, 0}, {4, 5}, {5, 6}, {6, 7}, {7, 4}, {0, 4}, {1, 5}, {3, 7}, {2, 6}};
    return edges[edgeIndex];
  }

  AXOM_HOST_DEVICE constexpr static const char *name() { return "hex"; }
};

/**
 * \brief This class extends the ShapeTraits with object state so it can represent a zone.
 */
template <typename ShapeTraits>
struct Shape : public ShapeTraits
{
  using IndexType = typename ShapeTraits::IndexType;

  /**
   * \brief Construct a shape.
   */
  AXOM_HOST_DEVICE Shape(const axom::ArrayView<IndexType> &ids) : m_ids(ids), m_faceIds()
  {
    SLIC_ASSERT(m_ids.size() == ShapeTraits::numberOfNodes());
  }

  /**
   * \brief Get a specific id that makes up this shape.
   *
   * \return The i'th id that makes up this shape.
   */
  AXOM_HOST_DEVICE IndexType getId(size_t index) const { return m_ids[index]; }

  /**
   * \brief Get the ids that make up this shape.
   *
   * \return A view containing the ids that make up this shape.
   */
  AXOM_HOST_DEVICE axom::ArrayView<IndexType> getIds() const { return m_ids; }

  /**
   * \brief Get the unique ids that make up this shape. For basic shapes, assume they are unique.
   *
   * \return A view containing the ids that make up this shape.
   */
  AXOM_HOST_DEVICE axom::ArrayView<IndexType> getUniqueIds() const { return m_ids; }

  /**
   * \brief Get the ids for the requested face.
   *
   * \param faceIndex The index of the desired face.
   *
   * \return An array view (wrapping m_faceIds) that contains the ids for the face.
   */
  AXOM_HOST_DEVICE
  axom::ArrayView<IndexType>
  getFace(int faceIndex) const
  {
    if constexpr(ShapeTraits::dimension() == 2)
      return m_ids;
    else
    {  
      const auto nnodes = ShapeTraits::numberOfNodesInFace(faceIndex);
      for(IndexType i = 0; i < nnodes; i++)
        m_faceIds[i] = m_ids[ShapeTraits::faces[faceIndex][i]];
      return axom::ArrayView<IndexType>(m_faceIds.m_data, nnodes);
    }
  }
private:
  axom::ArrayView<IndexType> m_ids;
  mutable axom::StackArray<IndexType, ShapeTraits::maxNodesInFace()> m_faceIds;
};

// Make some concrete shape classes based on the shape traits.
template <typename IndexT = axom::IndexType>
using LineShape = Shape<LineTraits<IndexT>>;

template <typename IndexT = axom::IndexType>
using TriShape = Shape<TriTraits<IndexT>>;

template <typename IndexT = axom::IndexType>
using QuadShape = Shape<QuadTraits<IndexT>>;

template <typename IndexT = axom::IndexType>
using TetShape = Shape<TetTraits<IndexT>>;

template <typename IndexT = axom::IndexType>
using PyramidShape = Shape<PyramidTraits<IndexT>>;

template <typename IndexT = axom::IndexType>
using WedgeShape = Shape<WedgeTraits<IndexT>>;

template <typename IndexT = axom::IndexType>
using HexShape = Shape<HexTraits<IndexT>>;


/**
 * \brief This is a shape that can act as any of the other shapes.
 *
 * \note This is a substitute for polymorphism so we can run on device.
 */
template <typename IndexT>
struct VariableShape
{
  using IndexType = IndexT;

  /**
   * \brief Constructor
   *
   * \param shapeId The shape id that describes the points.
   * \param ids The ids that describe the shape.
   */
  AXOM_HOST_DEVICE
  VariableShape(int shapeId, const axom::ArrayView<IndexType> &ids) : m_shapeId(shapeId), m_ids(ids)
  {
  }

  /**
   * \brief Returns the shape id of the actual shape represented by the variable shape.
   * \return The actual shape represented.
   */
  AXOM_HOST_DEVICE int id() const { return m_shapeId; }

  AXOM_HOST_DEVICE constexpr static bool is_polyhedral() { return false; }
  AXOM_HOST_DEVICE constexpr static bool is_variable_size() { return true; }

  AXOM_HOST_DEVICE IndexType dimension() const
  {
    int dim = 2;
    switch(m_shapeId)
    {
    case LineShape<IndexType>::id():    dim = LineShape<IndexType>::dimension(); break;
    case TriShape<IndexType>::id():     dim = TriShape<IndexType>::dimension(); break;
    case QuadShape<IndexType>::id():    dim = QuadShape<IndexType>::dimension(); break;
    case PolygonShape<IndexType>::id(): dim = PolygonShape<IndexType>::dimension(); break;
    case TetShape<IndexType>::id():     dim = TetShape<IndexType>::dimension(); break;
    case PyramidShape<IndexType>::id(): dim = PyramidShape<IndexType>::dimension(); break;
    case WedgeShape<IndexType>::id():   dim = WedgeShape<IndexType>::dimension(); break;
    case HexShape<IndexType>::id():     dim = HexShape<IndexType>::dimension(); break;
    }
    return dim;
  }

  AXOM_HOST_DEVICE IndexType numberOfNodes() const
  {
    IndexType nnodes = 0;
    switch(m_shapeId)
    {
    case LineShape<IndexType>::id():    nnodes = LineShape<IndexType>::numberOfNodes(); break;
    case TriShape<IndexType>::id():     nnodes = TriShape<IndexType>::numberOfNodes(); break;
    case QuadShape<IndexType>::id():    nnodes = QuadShape<IndexType>::numberOfNodes(); break;
    case PolygonShape<IndexType>::id(): nnodes = PolygonShape<IndexType>::numberOfNodes(); break;
    case TetShape<IndexType>::id():     nnodes = TetShape<IndexType>::numberOfNodes(); break;
    case PyramidShape<IndexType>::id(): nnodes = PyramidShape<IndexType>::numberOfNodes(); break;
    case WedgeShape<IndexType>::id():   nnodes = WedgeShape<IndexType>::numberOfNodes(); break;
    case HexShape<IndexType>::id():     nnodes = HexShape<IndexType>::numberOfNodes(); break;
    }
    return nnodes;
  }

  AXOM_HOST_DEVICE IndexType numberOfNodesInFace(int faceIndex) const
  {
    IndexType nnodes = 0;
    switch(m_shapeId)
    {
    case LineShape<IndexType>::id():    nnodes = LineShape<IndexType>::numberOfNodesInFace(faceIndex); break;
    case TriShape<IndexType>::id():     nnodes = TriShape<IndexType>::numberOfNodesInFace(faceIndex); break;
    case QuadShape<IndexType>::id():    nnodes = QuadShape<IndexType>::numberOfNodesInFace(faceIndex); break;
    case PolygonShape<IndexType>::id(): nnodes = PolygonShape<IndexType>::numberOfNodesInFace(faceIndex); break;
    case TetShape<IndexType>::id():     nnodes = TetShape<IndexType>::numberOfNodesInFace(faceIndex); break;
    case PyramidShape<IndexType>::id(): nnodes = PyramidShape<IndexType>::numberOfNodesInFace(faceIndex); break;
    case WedgeShape<IndexType>::id():   nnodes = WedgeShape<IndexType>::numberOfNodesInFace(faceIndex); break;
    case HexShape<IndexType>::id():     nnodes = HexShape<IndexType>::numberOfNodesInFace(faceIndex); break;
    }
    return nnodes;
  }

  AXOM_HOST_DEVICE IndexType maxNodesInFace() const
  {
    IndexType nnodes = 0;
    switch(m_shapeId)
    {
    case LineShape<IndexType>::id():    nnodes = LineShape<IndexType>::maxNodesInFace(); break;
    case TriShape<IndexType>::id():     nnodes = TriShape<IndexType>::maxNodesInFace(); break;
    case QuadShape<IndexType>::id():    nnodes = QuadShape<IndexType>::maxNodesInFace(); break;
    case PolygonShape<IndexType>::id(): nnodes = PolygonShape<IndexType>::maxNodesInFace(); break;
    case TetShape<IndexType>::id():     nnodes = TetShape<IndexType>::maxNodesInFace(); break;
    case PyramidShape<IndexType>::id(): nnodes = PyramidShape<IndexType>::maxNodesInFace(); break;
    case WedgeShape<IndexType>::id():   nnodes = WedgeShape<IndexType>::maxNodesInFace(); break;
    case HexShape<IndexType>::id():     nnodes = HexShape<IndexType>::maxNodesInFace(); break;
    }
    return nnodes;
  }

  AXOM_HOST_DEVICE IndexType numberOfFaces() const
  {
    IndexType nfaces = 0;
    switch(m_shapeId)
    {
    case LineShape<IndexType>::id():    nfaces = LineShape<IndexType>::numberOfFaces(); break;
    case TriShape<IndexType>::id():     nfaces = TriShape<IndexType>::numberOfFaces(); break;
    case QuadShape<IndexType>::id():    nfaces = QuadShape<IndexType>::numberOfFaces(); break;
    case PolygonShape<IndexType>::id(): nfaces = PolygonShape<IndexType>::numberOfFaces(); break;
    case TetShape<IndexType>::id():     nfaces = TetShape<IndexType>::numberOfFaces(); break;
    case PyramidShape<IndexType>::id(): nfaces = PyramidShape<IndexType>::numberOfFaces(); break;
    case WedgeShape<IndexType>::id():   nfaces = WedgeShape<IndexType>::numberOfFaces(); break;
    case HexShape<IndexType>::id():     nfaces = HexShape<IndexType>::numberOfFaces(); break;
    }
    return nfaces;
  }

  AXOM_HOST_DEVICE IndexType numberOfEdges() const
  {
    IndexType nedges = 0;
    switch(m_shapeId)
    {
    case LineShape<IndexType>::id():    nedges = LineShape<IndexType>::numberOfEdges(); break;
    case TriShape<IndexType>::id():     nedges = TriShape<IndexType>::numberOfEdges(); break;
    case QuadShape<IndexType>::id():    nedges = QuadShape<IndexType>::numberOfEdges(); break;
    case PolygonShape<IndexType>::id(): nedges = PolygonShape<IndexType>::numberOfEdges(); break;
    case TetShape<IndexType>::id():     nedges = TetShape<IndexType>::numberOfEdges(); break;
    case PyramidShape<IndexType>::id(): nedges = PyramidShape<IndexType>::numberOfEdges(); break;
    case WedgeShape<IndexType>::id():   nedges = WedgeShape<IndexType>::numberOfEdges(); break;
    case HexShape<IndexType>::id():     nedges = HexShape<IndexType>::numberOfEdges(); break;
    }
    return nedges;
  }

  AXOM_HOST_DEVICE axom::StackArray<IndexType, 2> getEdge(int edgeIndex)
  {
    axom::StackArray<IndexType, 2> edge;
    switch(m_shapeId)
    {
    case LineShape<IndexType>::id():    edge = LineShape<IndexType>::getEdge(edgeIndex); break;
    case TriShape<IndexType>::id():     edge = TriShape<IndexType>::getEdge(edgeIndex); break;
    case QuadShape<IndexType>::id():    edge = QuadShape<IndexType>::getEdge(edgeIndex); break;
    case PolygonShape<IndexType>::id(): edge = PolygonShape<IndexType>::getEdge(edgeIndex); break;
    case TetShape<IndexType>::id():     edge = TetShape<IndexType>::getEdge(edgeIndex); break;
    case PyramidShape<IndexType>::id(): edge = PyramidShape<IndexType>::getEdge(edgeIndex); break;
    case WedgeShape<IndexType>::id():   edge = WedgeShape<IndexType>::getEdge(edgeIndex); break;
    case HexShape<IndexType>::id():     edge = HexShape<IndexType>::getEdge(edgeIndex); break;
    }
    return edge;
  }

  /**
   * \brief Get a specific id that makes up this shape.
   *
   * \return The i'th id that makes up this shape.
   */
  AXOM_HOST_DEVICE IndexType getId(size_t index) const { return m_ids[index]; }

  /**
   * \brief Get the ids that make up this shape.
   *
   * \return A view containing the ids that make up this shape.
   */
  AXOM_HOST_DEVICE axom::ArrayView<IndexType> getIds() const { return m_ids; }

  AXOM_HOST_DEVICE constexpr static const char *name() { return "mixed"; }
private:
  int m_shapeId;
  axom::ArrayView<IndexType> m_ids;
};

/**
 * \brief Given a shape name (matches Blueprint shape name), return the Shape id() value.
 *
 * \param name The shape name.
 *
 * \return The shape id that matches the name, or 0 if there is no match.
 */
template <typename IndexT>
IndexT shapeNameToID(const std::string &name)
{
  IndexT id = 0;
  if(name == LineShape<IndexT>::name())
    id = LineShape<IndexT>::id();
  else if(name == TriShape<IndexT>::name())
    id = TriShape<IndexT>::id();
  else if(name == QuadShape<IndexT>::name())
    id = QuadShape<IndexT>::id();
  else if(name == PolygonShape<IndexT>::name())
    id = PolygonShape<IndexT>::id();
  else if(name == TetShape<IndexT>::name())
    id = TetShape<IndexT>::id();
  else if(name == PyramidShape<IndexT>::name())
    id = PyramidShape<IndexT>::id();
  else if(name == WedgeShape<IndexT>::name())
    id = WedgeShape<IndexT>::id();
  else if(name == HexShape<IndexT>::name())
    id = HexShape<IndexT>::id();
  return id;
}


} // end namespace views
} // end namespace mir
} // end namespace axom

#endif
